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


def _git_dirty(root: Path) -> bool:
    result = subprocess.run(
        ["git", "status", "--porcelain", "--untracked-files=no"],
        cwd=root,
        check=True,
        capture_output=True,
        text=True,
    )
    return bool(result.stdout.strip())


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


def _build(source: Path, target: Path, *, local_salc_settings=None):
    from matrix_chem import preprocess_to_enriched_xyz, write_validation_section
    from matrix_neo import write_gicforge_build_sections

    preprocess_to_enriched_xyz(source, target)
    write_validation_section(target)
    return write_gicforge_build_sections(
        target,
        local_salc=True,
        local_salc_settings=local_salc_settings,
        symmetrize=True,
    )


def _generic_coordination_geometry(coordination: int, radius: float = 1.9) -> np.ndarray:
    directions = []
    golden_angle = np.pi * (3.0 - np.sqrt(5.0))
    for index in range(coordination):
        z_coord = 1.0 - 2.0 * (index + 0.5) / coordination
        radial = np.sqrt(max(0.0, 1.0 - z_coord * z_coord))
        theta = index * golden_angle + 0.17
        directions.append((radial * np.cos(theta), radial * np.sin(theta), z_coord))
    return np.vstack((np.zeros((1, 3)), radius * np.asarray(directions, dtype=float)))


def _row_space_residual(left: np.ndarray, right: np.ndarray) -> float:
    def basis(matrix: np.ndarray) -> np.ndarray:
        _u, singular_values, vh = np.linalg.svd(matrix, full_matrices=False)
        tolerance = max(matrix.shape) * singular_values[0] * np.finfo(float).eps * 10.0
        return vh[: int(np.sum(singular_values > tolerance))]

    left_basis = basis(left)
    right_basis = basis(right)
    if left_basis.shape != right_basis.shape:
        return float("inf")
    projected = left_basis @ right_basis.T @ right_basis
    return float(np.linalg.norm(left_basis - projected) / np.sqrt(left_basis.shape[0]))


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
            LocalSALCSettings,
            _LOCAL_COORDINATION_TEMPLATES,
            _local_coordination_match,
            build_gicforge_python_model,
        )
        from matrix_neo.survibfit.pipeline import b_matrix_analytic

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
                        "python_template_frozen": "TEMPLATE_STATUS=FROZEN"
                        in python_angles,
                        "python_template_score": _field_values(
                            tuple(python_angles.splitlines()), "TEMPLATE_SCORE"
                        ),
                        "python_template_margin": _field_values(
                            tuple(python_angles.splitlines()), "TEMPLATE_MARGIN"
                        ),
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

        generic_high_coordination = []
        for coordination in (10, 12):
            coordinates = _generic_coordination_geometry(coordination)
            atoms = ("Fe",) + ("F",) * coordination
            model = build_gicforge_python_model(atoms, coordinates, local_salc=True)
            python_definition = model.to_definition()
            python_b = python_definition.u_matrix.T @ b_matrix_analytic(
                python_definition.primitives,
                coordinates,
            )
            python_diagnostics = "\n".join(
                coordinate.diagnostic for coordinate in model.coordinates
            )
            record = {
                "coordination": coordination,
                "target_rank": 3 * coordination - 3,
                "python_rank": len(model.coordinates),
                "python_generic_fallback": "TEMPLATE_STATUS=GENERIC"
                in python_diagnostics,
                "python_assignment_frozen": "ASSIGNMENT=FROZEN"
                in python_diagnostics,
            }
            if fortran_available:
                generic_run = run_legacy_gicforge(
                    work / f"fortran_generic_cn{coordination}",
                    atoms=atoms,
                    coordinates_angstrom=coordinates,
                    point_group="C1",
                    keywords=("GNIC", "BMAT", "LOCSALC"),
                    repo_root=matrix,
                )
                fortran_b = np.asarray(generic_run.b_matrix_rows, dtype=float)
                record.update(
                    {
                        "fortran_rank": generic_run.final_counts[-1],
                        "fortran_generic_fallback": "template GEN"
                        in generic_run.provout,
                        "python_fortran_row_space_residual": _row_space_residual(
                            python_b,
                            fortran_b,
                        ),
                    }
                )
            generic_high_coordination.append(record)

        first_template, second_template = _LOCAL_COORDINATION_TEMPLATES[6]
        left = np.asarray(first_template.directions, dtype=float)
        right = np.asarray(second_template.directions, dtype=float)
        template_tie_candidates = []
        for fraction in np.linspace(0.0, 1.0, 201):
            directions = (1.0 - fraction) * left + fraction * right
            directions /= np.linalg.norm(directions, axis=1)[:, None]
            coordinates = np.vstack((np.zeros((1, 3)), directions))
            match = _local_coordination_match(
                0,
                list(range(1, 7)),
                coords=coordinates,
                max_rms_cosine_error=1.0,
                min_score_margin=0.0,
            )
            template_tie_candidates.append((match.margin, fraction, coordinates))
        tie_margin, tie_fraction, tie_coordinates = min(
            template_tie_candidates,
            key=lambda item: item[0],
        )
        ambiguous_match = _local_coordination_match(
            0,
            list(range(1, 7)),
            coords=tie_coordinates,
            max_rms_cosine_error=1.0,
            min_score_margin=tie_margin + 1.0e-6,
        )

        chemistry_examples = []
        for label, source, settings, expected_template, expected_group in (
            (
                "distorted SF6",
                matrix / "examples/local_pseudosymmetry/sf6_distorted.xyz",
                LocalSALCSettings(
                    zeff_tolerance=0.1,
                    distance_tolerance_angstrom=0.05,
                ),
                "OCTAHEDRAL",
                "Oh",
            ),
            (
                "square-antiprismatic [ZrF8]4-",
                matrix / "examples/local_pseudosymmetry/zrf8_square_antiprismatic.xyz",
                LocalSALCSettings(),
                "SQUARE_ANTIPRISMATIC",
                "D4d",
            ),
        ):
            definition = _build(
                source,
                work / f"{source.stem}.xyzin",
                local_salc_settings=settings,
            )
            lines = _local_lines(definition, domain="CENTER:1")
            chemistry_examples.append(
                {
                    "system": label,
                    "rank": definition.rank,
                    "target_rank": definition.target_rank,
                    "expected_template": expected_template,
                    "expected_group": expected_group,
                    "templates": _field_values(lines, "TEMPLATE"),
                    "groups": _field_values(lines, "GROUP"),
                    "template_statuses": _field_values(lines, "TEMPLATE_STATUS"),
                    "template_scores": _field_values(lines, "TEMPLATE_SCORE"),
                    "template_margins": _field_values(lines, "TEMPLATE_MARGIN"),
                    "zeff_tolerance": settings.zeff_tolerance,
                    "distance_tolerance_angstrom": settings.distance_tolerance_angstrom,
                }
            )

        output = {
            "method": "SMITH/SONIC local-pseudosymmetry SALC validation",
            "matrix_revision": _git_revision(matrix),
            "matrix_worktree_dirty": _git_dirty(matrix),
            "python": sys.version.split()[0],
            "settings": {
                "local_salc": True,
                "global_symmetrization": True,
                "effective_atomic_number_tolerance": 5.0e-4,
                "radial_distance_tolerance_angstrom": 1.0e-3,
                "template_rms_threshold": 0.12,
                "template_min_margin": 0.02,
                "angle_class_tolerance": 0.02,
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
            "generic_high_coordination": generic_high_coordination,
            "template_ambiguity": {
                "coordination": 6,
                "interpolation_fraction": tie_fraction,
                "best_margin": tie_margin,
                "configured_margin": tie_margin + 1.0e-6,
                "status": ambiguous_match.status,
                "generic_fallback": ambiguous_match.template is None,
            },
            "chemical_examples": chemistry_examples,
            "systems": systems,
            "notes": [
                "Local A1 labels refer to a center/ring/bond domain; molecular irreps remain controlled by the frozen point-group projector.",
                "The non-A1 labels identify a deterministic orthonormal complement and do not claim a complete local-irrep decomposition.",
                "Azulene ring-stretch ownership records the single shared fused edge only once.",
                "Coordination numbers 5--9 are tested for two idealized templates each, both ideally and after a small angular distortion, in Python and Fortran when available.",
                "Coordination numbers 10 and 12 exercise the template-independent fallback and Python/Fortran Cartesian row-space agreement.",
                "A configurable near-tie demonstrates the AMBIGUOUS decision and generic fallback without changing the frozen SONIC identity.",
                "Distorted SF6 and square-antiprismatic [ZrF8]4- provide chemically named CN 6 and CN 8 examples.",
            ],
        }
        DATA.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
