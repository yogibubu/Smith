#!/usr/bin/env python3
"""Reproduce the local-pseudosymmetry validation reported in the SMITH paper."""

from __future__ import annotations

import json
import shutil
import subprocess
import sys
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np

from _matrix_checkout import add_matrix_packages_to_path, matrix_root

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "local_salc_validation.json"
ATOMS = ("C", "H", "H", "H", "C", "F", "O", "H")
C1_METHYL_COORDINATES = np.array(
    [
        (0.0, 0.0, 0.0),
        (-0.629312, 0.629312, 0.629312),
        (-0.629312, -0.629312, 0.629312),
        (-0.629312, 0.0, -0.889982),
        (1.54, 0.0, 0.0),
        (2.14, 1.05, 0.11),
        (2.31, -0.87, 0.34),
        (2.95, -0.78, 0.95),
    ],
    dtype=float,
)
HIGH_COORDINATION_GROUPS = {
    "TRIGONAL_BIPYRAMIDAL": "D3h",
    "SQUARE_PYRAMIDAL": "C4v",
    "OCTAHEDRAL": "Oh",
    "TRIGONAL_PRISMATIC": "D3h",
    "PENTAGONAL_BIPYRAMIDAL": "D5h",
    "CAPPED_OCTAHEDRAL": "C3v",
    "SQUARE_ANTIPRISMATIC": "D4d",
    "DODECAHEDRAL_LIKE": "D2d",
    "TRICAPPED_TRIGONAL_PRISMATIC": "D3h",
    "CAPPED_SQUARE_ANTIPRISMATIC": "C4v",
}


def _git_revision(root: Path) -> str:
    result = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=root,
        check=True,
        capture_output=True,
        text=True,
    )
    return result.stdout.strip()


def _write_xyz(path: Path, title: str, atoms: tuple[str, ...], coordinates: np.ndarray) -> None:
    path.write_text(
        "\n".join(
            [str(len(atoms)), title]
            + [
                f"{atom} {xyz[0]:.12f} {xyz[1]:.12f} {xyz[2]:.12f}"
                for atom, xyz in zip(atoms, coordinates)
            ]
        )
        + "\n",
        encoding="utf-8",
    )


def _local_lines(definition, *, domain: str | None = None) -> tuple[str, ...]:
    lines = tuple(
        line
        for line in definition.reduction_diagnostics.skipped_dependent_details
        if "LOCAL_SALC" in line and (domain is None or f"DOMAIN={domain}" in line)
    )
    return lines


def _field_values(lines: tuple[str, ...], field: str) -> list[str]:
    marker = f"{field}="
    return sorted(
        {
            token.removeprefix(marker)
            for line in lines
            for token in line.split()
            if token.startswith(marker)
        }
    )


def _build(source: Path, target: Path):
    from matrix_chem import preprocess_to_enriched_xyz, write_validation_section
    from matrix_neo import write_gicforge_build_sections

    preprocess_to_enriched_xyz(source, target)
    write_validation_section(target)
    return write_gicforge_build_sections(target, local_salc=True, symmetrize=True)


def _system_record(name: str, definition) -> dict[str, object]:
    lines = _local_lines(definition)
    return {
        "system": name,
        "point_group": definition.point_group,
        "rank": definition.rank,
        "target_rank": definition.target_rank,
        "symmetry_method": definition.symmetry_diagnostics.method,
        "local_salc_records": len(lines),
        "local_groups": _field_values(lines, "GROUP"),
        "local_domains": _field_values(lines, "DOMAIN"),
        "local_irreps": _field_values(lines, "LOCAL_IRREP"),
        "kept": sum("STATUS=KEPT" in line for line in lines),
        "pruned": sum("STATUS=PRUNED" in line for line in lines),
    }


def main() -> None:
    add_matrix_packages_to_path()
    matrix = matrix_root()
    fixture_root = matrix / "tests" / "fixtures" / "test_molecules" / "molecules"
    with TemporaryDirectory(prefix="smith-local-salc-") as temporary:
        work = Path(temporary)
        c1_source = work / "c1_methyl.xyz"
        _write_xyz(c1_source, "globally C1, locally C3v methyl probe", ATOMS, C1_METHYL_COORDINATES)
        c1_definition = _build(c1_source, work / "c1_methyl.xyzin")

        angle = 0.731
        axis = np.array((1.0, -2.0, 0.5), dtype=float)
        axis /= np.linalg.norm(axis)
        cross = np.array(
            [
                (0.0, -axis[2], axis[1]),
                (axis[2], 0.0, -axis[0]),
                (-axis[1], axis[0], 0.0),
            ],
            dtype=float,
        )
        rotation = (
            np.eye(3) * np.cos(angle)
            + (1.0 - np.cos(angle)) * np.outer(axis, axis)
            + np.sin(angle) * cross
        )
        rotated_source = work / "c1_methyl_rotated.xyz"
        _write_xyz(
            rotated_source,
            "rigidly rotated C1 methyl probe",
            ATOMS,
            C1_METHYL_COORDINATES @ rotation.T,
        )
        rotated_definition = _build(rotated_source, work / "c1_methyl_rotated.xyzin")

        systems = [_system_record("C1 methyl probe", c1_definition)]
        definitions = {}
        for label, source_name in (
            ("benzene", "benzene.inp"),
            ("toluene", "toluene.inp"),
            ("azulene", "azulene.inp"),
        ):
            definition = _build(fixture_root / source_name, work / f"{label}.xyzin")
            definitions[label] = definition
            systems.append(_system_record(label, definition))

        c1_lines = _local_lines(c1_definition, domain="CENTER:1")
        rotated_lines = _local_lines(rotated_definition, domain="CENTER:1")
        benzene_ring = _local_lines(definitions["benzene"], domain="RING:1")
        toluene_ring = _local_lines(definitions["toluene"], domain="RING:1")
        azulene_ring = tuple(
            line
            for line in _local_lines(definitions["azulene"])
            if "DOMAIN=RING:" in line
        )

        from matrix_engines import gicforge_fortran_layout, run_legacy_gicforge
        from matrix_neo.runtime.gicforge_python import (
            _LOCAL_COORDINATION_TEMPLATES,
            build_gicforge_python_model,
        )

        fortran_layout = gicforge_fortran_layout(matrix)
        fortran_available = shutil.which("gfortran") is not None or fortran_layout.legacy_executable.is_file()
        fortran: dict[str, object] = {"available": fortran_available}
        if fortran_available:
            run = run_legacy_gicforge(
                work / "fortran_c1_methyl",
                atoms=ATOMS,
                coordinates_angstrom=C1_METHYL_COORDINATES,
                point_group="C1",
                keywords=("GNIC", "BMAT", "LOCSALC"),
                repo_root=matrix,
            )
            fortran.update(
                {
                    "gic_count": len(run.gic_labels),
                    "final_counts": list(run.final_counts),
                    "c3v_diagnostic_count": run.provout.count("GROUP=C3v"),
                    "a1_first_count": run.provout.count("A1_FIRST=YES"),
                }
            )

        high_coordination = []
        for coordination in range(5, 10):
            for template_index, template in enumerate(
                _LOCAL_COORDINATION_TEMPLATES[coordination]
            ):
                expected_group = HIGH_COORDINATION_GROUPS[template.name]
                directions = np.asarray(template.directions, dtype=float)
                directions /= np.linalg.norm(directions, axis=1)[:, None]
                variants = {}
                for distorted in (False, True):
                    working_directions = directions.copy()
                    if distorted:
                        working_directions[0] += np.array(
                            (0.025, -0.015, 0.010), dtype=float
                        )
                        working_directions[0] /= np.linalg.norm(working_directions[0])
                    coordinates = np.vstack(
                        (np.zeros((1, 3), dtype=float), 1.6 * working_directions)
                    )
                    atoms = ("Fe",) + ("H",) * coordination
                    model = build_gicforge_python_model(
                        atoms,
                        coordinates,
                        local_salc=True,
                    )
                    python_angles = "\n".join(
                        coordinate.diagnostic
                        for coordinate in model.coordinates
                        if "KIND=ANGLE" in coordinate.diagnostic
                    )
                    variant = {
                        "python_rank": len(model.coordinates),
                        "python_target_rank": model.target_rank,
                        "python_group_recognized": f"GROUP={expected_group}"
                        in python_angles,
                    }
                    if fortran_available:
                        high_run = run_legacy_gicforge(
                            work
                            / (
                                f"fortran_cn{coordination}_template{template_index}_"
                                f"distorted{int(distorted)}"
                            ),
                            atoms=atoms,
                            coordinates_angstrom=coordinates,
                            point_group="C1",
                            keywords=("GNIC", "BMAT", "LOCSALC"),
                            repo_root=matrix,
                        )
                        fortran_angles = "\n".join(
                            line
                            for line in high_run.provout.splitlines()
                            if "LOCAL_SALC DOMAIN=CENTER:" in line
                            and "KIND=ANGLE" in line
                        )
                        variant.update(
                            {
                                "fortran_rank": high_run.final_counts[-1],
                                "fortran_group_recognized": f"GROUP={expected_group}"
                                in fortran_angles,
                                "fortran_a1_first": "A1_FIRST=YES" in fortran_angles,
                            }
                        )
                    variants["distorted" if distorted else "ideal"] = variant
                high_coordination.append(
                    {
                        "coordination": coordination,
                        "template": template.name,
                        "expected_local_group": expected_group,
                        "target_rank": 3 * coordination - 3,
                        "variants": variants,
                    }
                )

        output = {
            "method": "SMITH/SONIC local-pseudosymmetry SALC validation",
            "matrix_revision": _git_revision(matrix),
            "python": sys.version.split()[0],
            "settings": {
                "local_salc": True,
                "global_symmetrization": True,
                "effective_atomic_number_tolerance": 5.0e-4,
                "radial_distance_tolerance_angstrom": 1.0e-3,
            },
            "checks": {
                "c1_methyl_rotation_invariant": c1_lines == rotated_lines,
                "c1_methyl_global_point_group": c1_definition.point_group,
                "c1_methyl_local_center_groups": _field_values(c1_lines, "GROUP"),
                "c1_methyl_local_center_irreps": _field_values(c1_lines, "LOCAL_IRREP"),
                "benzene_ring_groups": _field_values(benzene_ring, "GROUP"),
                "benzene_ring_operations": _field_values(benzene_ring, "OPERATIONS"),
                "toluene_ring_groups": _field_values(toluene_ring, "GROUP"),
                "toluene_ring_operations": _field_values(toluene_ring, "OPERATIONS"),
                "azulene_ring_domains": _field_values(azulene_ring, "DOMAIN"),
                "azulene_shared_edges": _field_values(azulene_ring, "SHARED_EDGES"),
            },
            "fortran_control": fortran,
            "high_coordination": high_coordination,
            "systems": systems,
            "notes": [
                "Local A1 labels refer to a center/ring/bond domain; molecular irreps remain controlled by the frozen point-group projector.",
                "The non-A1 labels identify a deterministic orthonormal complement and do not claim a complete local-irrep decomposition.",
                "Azulene ring-stretch ownership records the single shared fused edge only once.",
                "Coordination numbers 5--9 are tested for two idealized templates each, both ideally and after a small angular distortion, in Python and Fortran when available.",
            ],
        }
        DATA.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
