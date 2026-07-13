#!/usr/bin/env python3
"""Prepare, run and summarize Gaussian ReadAllGIC optimization probes.

The committed JSON data are intentionally compact.  Raw Gaussian inputs, logs
and checkpoint files are written under ``calculations/readallgic_opt/``, which
is ignored by git.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
from pathlib import Path

from _matrix_checkout import add_matrix_packages_to_path, matrix_root


ROOT = Path(__file__).resolve().parents[1]
RUN_DIR = ROOT / "calculations" / "readallgic_opt"
DATA = ROOT / "data" / "readallgic_optimization_results.json"
DEFAULT_GAUSSIAN = Path("g16")
MOLECULES = (
    ("cyclohexane", "cyclohexane.inp"),
    ("norbornane", "norbornane.inp"),
    ("cubane", "cubane.inp"),
)
TOTAL_SYMMETRIC_MOLECULES = (("norbornane", "norbornane.inp"),)

SCF_RE = re.compile(r"SCF Done:\s+E\([^)]+\)\s+=\s+([-+0-9.Ee]+)")
STEP_RE = re.compile(r"Step number\s+(\d+)\s+out of a maximum")
RANK_RE = re.compile(r"NTRed=\s*(\d+).*?NRank=\s*(\d+).*?EigSml=([0-9.DdEe+-]+)")
NORMAL_TERMINATION = "Normal termination of Gaussian"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prepare", action="store_true", help="write Gaussian inputs")
    parser.add_argument("--run", action="store_true", help="run Gaussian after preparing inputs")
    parser.add_argument("--gaussian", type=Path, default=DEFAULT_GAUSSIAN)
    args = parser.parse_args()

    if args.prepare or args.run:
        prepare_inputs()
    if args.prepare and not args.run:
        return
    if args.run:
        run_gaussian(args.gaussian)

    DATA.parent.mkdir(parents=True, exist_ok=True)
    DATA.write_text(json.dumps(summarize_results(), indent=2) + "\n", encoding="utf-8")


def prepare_inputs() -> None:
    add_matrix_packages_to_path()

    from matrix_chem import preprocess_to_enriched_xyz, write_validation_section
    from matrix_gaussian import write_gicforge_gaussian_input
    from matrix_neo import gaussian_gic_lines_from_xyzin, write_gicforge_build_sections

    source_root = matrix_root() / "tests" / "fixtures" / "test_molecules" / "molecules"
    for name, source_name in MOLECULES:
        target_dir = RUN_DIR / name
        target_dir.mkdir(parents=True, exist_ok=True)
        xyzin = target_dir / f"{name}.xyzin"
        preprocess_to_enriched_xyz(source_root / source_name, xyzin)
        write_validation_section(xyzin)
        write_gicforge_build_sections(xyzin, symmetrize=False)

        unsupported = [
            line
            for line in gaussian_gic_lines_from_xyzin(xyzin)
            if "U(" in line or "OuPl" in line
        ]
        if unsupported:
            preview = "; ".join(unsupported[:3])
            raise RuntimeError(f"{name} has out-of-plane GIC lines: {preview}")

        write_gicforge_gaussian_input(
            xyzin,
            target_dir / f"{name}_sonic_readallgic.gjf",
            route="#p hf/sto-3g opt=(maxcycle=80) nosymm",
            title=f"{name} SONIC ReadAllGIC optimization probe",
            link0=(
                f"%chk={name}_sonic_readallgic.chk",
                "%nprocshared=4",
                "%mem=2GB",
            ),
        )
    for name, source_name in TOTAL_SYMMETRIC_MOLECULES:
        target_dir = RUN_DIR / f"{name}_total_symmetric"
        target_dir.mkdir(parents=True, exist_ok=True)
        xyzin = target_dir / f"{name}.xyzin"
        preprocess_to_enriched_xyz(source_root / source_name, xyzin)
        write_validation_section(xyzin)
        write_gicforge_build_sections(xyzin, symmetrize=True)
        write_gicforge_gaussian_input(
            xyzin,
            target_dir / f"{name}_sonic_a1_readallgic.gjf",
            route="#p hf/sto-3g opt=(readallgic,maxcycle=80)",
            title=f"{name} SONIC total-symmetric ReadAllGIC optimization probe",
            link0=(
                f"%chk={name}_sonic_a1_readallgic.chk",
                "%nprocshared=4",
                "%mem=2GB",
            ),
            total_symmetric_only=False,
            freeze_non_total=True,
        )


def run_gaussian(executable: Path) -> None:
    if not executable.exists() and executable.parent != Path("."):
        raise FileNotFoundError(executable)
    env = os.environ.copy()
    env.setdefault("GAUSS_SCRDIR", str(Path.home() / "Documents" / "GAUSSIAN_SCRATCH"))
    for name, _source_name in MOLECULES:
        target_dir = RUN_DIR / name
        input_path = target_dir / f"{name}_sonic_readallgic.gjf"
        log_path = target_dir / f"{name}_sonic_readallgic.log"
        chk_path = target_dir / f"{name}_sonic_readallgic.chk"
        for path in (log_path, chk_path, target_dir / "fort.7"):
            path.unlink(missing_ok=True)
        subprocess.run([str(executable), input_path.name], cwd=target_dir, env=env, check=True)
    for name, _source_name in TOTAL_SYMMETRIC_MOLECULES:
        target_dir = RUN_DIR / f"{name}_total_symmetric"
        input_path = target_dir / f"{name}_sonic_a1_readallgic.gjf"
        for path in (
            target_dir / f"{name}_sonic_a1_readallgic.log",
            target_dir / f"{name}_sonic_a1_readallgic.chk",
            target_dir / "fort.7",
        ):
            path.unlink(missing_ok=True)
        subprocess.run([str(executable), input_path.name], cwd=target_dir, env=env, check=True)


def summarize_results() -> dict[str, object]:
    rows = []
    for name, source_name in MOLECULES:
        xyzin = RUN_DIR / name / f"{name}.xyzin"
        log = RUN_DIR / name / f"{name}_sonic_readallgic.log"
        rows.append(summarize_one(name, source_name, xyzin, log))
    return {
        "method": "Gaussian 16 ReadAllGIC optimization with SONIC-generated coordinates",
        "electronic_structure": "HF/STO-3G",
        "route": "#p hf/sto-3g opt=(readallgic,maxcycle=80) nosymm output=pickett",
        "out_of_plane_gics": False,
        "systems": rows,
        "total_symmetric_probe": [
            summarize_one(
                name,
                source_name,
                RUN_DIR / f"{name}_total_symmetric" / f"{name}.xyzin",
                RUN_DIR / f"{name}_total_symmetric" / f"{name}_sonic_a1_readallgic.log",
                total_symmetric_only=True,
                route="#p hf/sto-3g opt=(readallgic,maxcycle=80) output=pickett",
                symmetry_enabled=True,
            )
            for name, source_name in TOTAL_SYMMETRIC_MOLECULES
        ],
    }


def summarize_one(
    name: str,
    source_name: str,
    xyzin: Path,
    log: Path,
    *,
    total_symmetric_only: bool = False,
    route: str | None = None,
    symmetry_enabled: bool = False,
) -> dict[str, object]:
    add_matrix_packages_to_path()

    from matrix_neo import gaussian_gic_lines_from_xyzin, read_gic_definition_from_xyzin
    from matrix_neo import total_symmetric_gic_names

    definition = read_gic_definition_from_xyzin(xyzin)
    total_names = total_symmetric_gic_names(definition)
    gic_lines = gaussian_gic_lines_from_xyzin(
        xyzin,
        freeze_non_total=total_symmetric_only,
    )
    text = log.read_text(encoding="utf-8", errors="replace")
    energies = [float(match.group(1).replace("D", "E")) for match in SCF_RE.finditer(text)]
    ranks = [
        (
            int(match.group(1)),
            int(match.group(2)),
            float(match.group(3).replace("D", "E")),
        )
        for match in RANK_RE.finditer(text)
    ]
    steps = [int(match.group(1)) for match in STEP_RE.finditer(text)]
    return {
        "name": name,
        "source": source_name,
        "point_group": definition.point_group,
        "sonic_rank": definition.rank,
        "target_rank": definition.target_rank,
        "gaussian_gic_count": len(gic_lines),
        "active_total_symmetric_gic_count": len(total_names) if total_symmetric_only else None,
        "contains_out_of_plane": any("U(" in line or "OuPl" in line for line in gic_lines),
        "total_symmetric_only": total_symmetric_only,
        "symmetry_enabled": symmetry_enabled,
        "route": route,
        "normal_termination": NORMAL_TERMINATION in text,
        "stationary_point": "Stationary point found" in text,
        "optimization_completed": "Optimization completed" in text,
        "steps": max(steps) if steps else None,
        "final_energy_hartree": energies[-1] if energies else None,
        "rank_pairs": [(rank[0], rank[1]) for rank in ranks],
        "minimum_eig_small": min((rank[2] for rank in ranks), default=None),
    }


if __name__ == "__main__":
    main()
