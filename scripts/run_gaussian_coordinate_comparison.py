#!/usr/bin/env python3
"""Prepare, run and summarize Gaussian optimization coordinate comparisons."""

from __future__ import annotations

import argparse
import json
import math
import os
import re
import subprocess
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from _matrix_checkout import add_matrix_packages_to_path, matrix_root


ROOT = Path(__file__).resolve().parents[1]
RUN_DIR = ROOT / "calculations" / "coordinate_comparison"
DATA = ROOT / "data" / "gaussian_coordinate_comparison.json"
FIGURE = ROOT / "figures" / "gaussian_coordinate_comparison.png"
DEFAULT_GAUSSIAN = "g16"
SYSTEMS = (
    ("water", "h2ocart.inp", None),
    ("hydrogen_peroxide", "h2o2zmat.inp", "h2o2.inp"),
    ("camphor", "camphor.inp", None),
    ("norbornane", "norbornane.inp", None),
)
MODES = ("default", "cartesian", "zmatrix", "sonic")
ROUTES = {
    "default": "#p hf/sto-3g opt=(redundant,maxcycle=80)",
    "cartesian": "#p hf/sto-3g opt=(cartesian,maxcycle=80)",
    "zmatrix": "#p hf/sto-3g opt=(z-matrix,maxcycle=80)",
    "sonic": "#p hf/sto-3g opt=(readallgic,calcfc,maxcycle=80)",
}
SCF_RE = re.compile(r"SCF Done:\s+E\([^)]+\)\s+=\s+([-+0-9.Ee]+)")
STEP_RE = re.compile(r"Step number\s+(\d+)\s+out of a maximum")
NORMAL_TERMINATION = "Normal termination of Gaussian"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prepare", action="store_true")
    parser.add_argument("--run", action="store_true")
    parser.add_argument("--gaussian", default=DEFAULT_GAUSSIAN)
    parser.add_argument("--plot", action="store_true")
    args = parser.parse_args()

    if args.prepare or args.run:
        prepare_inputs()
    if args.prepare and not args.run:
        return
    if args.run:
        run_gaussian(args.gaussian)
    if args.prepare or args.run:
        results = summarize_results()
        DATA.write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8")
    else:
        results = json.loads(DATA.read_text(encoding="utf-8"))
    plot_results(results)


def prepare_inputs() -> None:
    add_matrix_packages_to_path()
    from matrix_chem import preprocess_to_enriched_xyz, read_enriched_xyz, write_validation_section
    from matrix_gaussian import write_gicforge_gaussian_input
    from matrix_smith import write_gicforge_build_sections

    source_root = matrix_root() / "tests" / "fixtures" / "test_molecules" / "molecules"
    for name, source_name, zmat_source_name in SYSTEMS:
        target_dir = RUN_DIR / name
        target_dir.mkdir(parents=True, exist_ok=True)
        source = source_root / source_name
        xyzin = target_dir / f"{name}.xyzin"
        preprocess_to_enriched_xyz(source, xyzin)
        write_validation_section(xyzin)
        write_gicforge_build_sections(xyzin, symmetrize=False, improper_dihedrals=True)
        geometry = read_enriched_xyz(xyzin)
        write_cartesian_input(target_dir / f"{name}_default.gjf", geometry, ROUTES["default"])
        write_cartesian_input(target_dir / f"{name}_cartesian.gjf", geometry, ROUTES["cartesian"])
        if zmat_source_name is None:
            write_zmatrix_input(target_dir / f"{name}_zmatrix.gjf", geometry, ROUTES["zmatrix"])
        else:
            write_fixture_zmatrix_input(
                source_root / zmat_source_name,
                target_dir / f"{name}_zmatrix.gjf",
                ROUTES["zmatrix"],
                f"{name} Z-matrix coordinate comparison",
            )
        write_gicforge_gaussian_input(
            xyzin,
            target_dir / f"{name}_sonic.gjf",
            route=ROUTES["sonic"],
            title=f"{name} SONIC ReadAllGIC optimization",
            link0=(f"%chk={name}_sonic.chk", "%nprocshared=4", "%mem=2GB"),
            g16_compatibility=True,
        )


def write_cartesian_input(path: Path, geometry: object, route: str) -> None:
    lines = [
        f"%chk={path.stem}.chk",
        "%nprocshared=4",
        "%mem=2GB",
        route,
        "",
        f"{path.stem} coordinate comparison",
        "",
        f"{getattr(geometry, 'charge', None) or 0} {getattr(geometry, 'multiplicity', None) or 1}",
    ]
    for atom, (x, y, z) in zip(geometry.atoms, geometry.coordinates_angstrom):
        lines.append(f"{atom:2s} {x:15.8f} {y:15.8f} {z:15.8f}")
    lines.append("")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_fixture_zmatrix_input(source: Path, path: Path, route: str, title: str) -> None:
    source_lines = source.read_text(encoding="utf-8").splitlines()
    body_start = 0
    for idx, line in enumerate(source_lines):
        if line.strip().startswith("#"):
            body_start = idx + 1
            break
    body = [line.rstrip() for line in source_lines[body_start:]]
    for idx, line in enumerate(body):
        if re.match(r"^\s*-?\d+\s+\d+\s*$", line):
            body = body[idx:]
            break
    while body and not body[-1].strip():
        body.pop()
    lines = [
        f"%chk={path.stem}.chk",
        "%nprocshared=4",
        "%mem=2GB",
        route,
        "",
        title,
        "",
        *body,
        "",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_zmatrix_input(path: Path, geometry: object, route: str) -> None:
    coords = np.asarray(geometry.coordinates_angstrom, dtype=float)
    lines = [
        f"%chk={path.stem}.chk",
        "%nprocshared=4",
        "%mem=2GB",
        route,
        "",
        f"{path.stem} coordinate comparison",
        "",
        f"{getattr(geometry, 'charge', None) or 0} {getattr(geometry, 'multiplicity', None) or 1}",
    ]
    for idx, atom in enumerate(geometry.atoms):
        if idx == 0:
            lines.append(f"{atom}")
        elif idx == 1:
            lines.append(f"{atom} 1 R{idx}")
        elif idx == 2:
            lines.append(f"{atom} 2 R{idx} 1 A{idx}")
        else:
            lines.append(f"{atom} {idx} R{idx} {idx - 1} A{idx} {idx - 2} D{idx}")
    lines.extend(["", f"R1={_distance(coords, 1, 0):.8f}"])
    if len(geometry.atoms) > 2:
        lines.append(f"R2={_distance(coords, 2, 1):.8f}")
        lines.append(f"A2={_angle(coords, 2, 1, 0):.8f}")
    for idx in range(3, len(geometry.atoms)):
        lines.append(f"R{idx}={_distance(coords, idx, idx - 1):.8f}")
        lines.append(f"A{idx}={_angle(coords, idx, idx - 1, idx - 2):.8f}")
        lines.append(f"D{idx}={_dihedral(coords, idx, idx - 1, idx - 2, idx - 3):.8f}")
    lines.append("")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_gaussian(executable: str) -> None:
    env = os.environ.copy()
    env.setdefault("GAUSS_SCRDIR", str(Path.home() / "gaussian_scratch"))
    for name, _source, _zmat_source in SYSTEMS:
        target_dir = RUN_DIR / name
        for mode in MODES:
            input_path = target_dir / f"{name}_{mode}.gjf"
            for suffix in (".log", ".chk"):
                (target_dir / f"{name}_{mode}{suffix}").unlink(missing_ok=True)
            subprocess.run([executable, input_path.name], cwd=target_dir, env=env, check=False)


def summarize_results() -> dict[str, object]:
    add_matrix_packages_to_path()
    from matrix_smith import read_gic_definition_from_xyzin

    systems = []
    for name, source, zmat_source in SYSTEMS:
        xyzin = RUN_DIR / name / f"{name}.xyzin"
        definition = read_gic_definition_from_xyzin(xyzin)
        rows = []
        for mode in MODES:
            log = RUN_DIR / name / f"{name}_{mode}.log"
            rows.append(summarize_log(mode, log))
        systems.append(
            {
                "name": name,
                "source": source,
                "zmatrix_source": zmat_source,
                "point_group": definition.point_group,
                "sonic_rank": definition.rank,
                "target_rank": definition.target_rank,
                "modes": rows,
            }
        )
    return {
        "method": "Gaussian 16 coordinate-system optimization comparison",
        "electronic_structure": "Gaussian 16 HF/STO-3G",
        "systems": systems,
    }


def summarize_log(mode: str, path: Path) -> dict[str, object]:
    text = path.read_text(encoding="utf-8", errors="replace") if path.exists() else ""
    energies = [float(match.group(1).replace("D", "E")) for match in SCF_RE.finditer(text)]
    steps = [int(match.group(1)) for match in STEP_RE.finditer(text)]
    failed = not text or NORMAL_TERMINATION not in text
    return {
        "mode": mode,
        "normal_termination": not failed,
        "stationary_point": "Stationary point found" in text,
        "optimization_completed": "Optimization completed" in text,
        "steps": max(steps) if steps else None,
        "scf_cycles": len(energies),
        "final_energy_hartree": energies[-1] if energies else None,
    }


def plot_results(results: dict[str, object]) -> None:
    systems = results["systems"]
    labels = [str(system["name"]).replace("_", "\n") for system in systems]
    mode_labels = {
        "default": "Gaussian redundant",
        "cartesian": "Cartesian",
        "zmatrix": "Z-matrix",
        "sonic": "SONIC ReadAllGIC",
    }
    values = {}
    failures = {}
    for mode in MODES:
        mode_values = []
        mode_failures = []
        for system in systems:
            row = next(row for row in system["modes"] if row["mode"] == mode)
            failed = not row["normal_termination"]
            mode_values.append(np.nan if failed else row["steps"])
            mode_failures.append(failed)
        values[mode] = mode_values
        failures[mode] = mode_failures
    fig, ax = plt.subplots(figsize=(7.6, 4.25))
    x = np.arange(len(labels), dtype=float)
    width = 0.19
    colors = {
        "default": "#5b6670",
        "cartesian": "#b7534c",
        "zmatrix": "#d1a33b",
        "sonic": "#1f6f78",
    }
    for offset, mode in enumerate(MODES):
        xpos = x + (offset - 1.5) * width
        bars = ax.bar(
            xpos,
            values[mode],
            width,
            label=mode_labels[mode],
            color=colors[mode],
            edgecolor="white",
            linewidth=0.45,
        )
        for bar, value, failed in zip(bars, values[mode], failures[mode]):
            if failed:
                ax.scatter(
                    [bar.get_x() + bar.get_width() / 2.0],
                    [0.6],
                    marker="x",
                    s=38,
                    linewidths=1.2,
                    color=colors[mode],
                    zorder=4,
                )
            elif np.isfinite(value):
                ax.text(
                    bar.get_x() + bar.get_width() / 2.0,
                    float(value) + 0.35,
                    f"{int(value)}",
                    ha="center",
                    va="bottom",
                    fontsize=8,
                    color="#24272a",
                )
    ax.set_xticks(x, labels)
    ax.set_ylabel("optimization steps")
    ymax = max(
        float(value)
        for mode_values in values.values()
        for value in mode_values
        if np.isfinite(value)
    )
    ax.set_ylim(0, ymax + 5.0)
    ax.legend(frameon=False, ncol=4, loc="upper center", bbox_to_anchor=(0.5, 1.14))
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", color="#d7d9dc", linewidth=0.7)
    FIGURE.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(FIGURE, dpi=320)
    plt.close(fig)


def _distance(coords: np.ndarray, i: int, j: int) -> float:
    return float(np.linalg.norm(coords[i] - coords[j]))


def _angle(coords: np.ndarray, i: int, j: int, k: int) -> float:
    u = coords[i] - coords[j]
    v = coords[k] - coords[j]
    cosang = float(np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v)))
    return math.degrees(math.acos(max(-1.0, min(1.0, cosang))))


def _dihedral(coords: np.ndarray, i: int, j: int, k: int, m: int) -> float:
    p0, p1, p2, p3 = coords[[i, j, k, m]]
    b0 = -(p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2
    b1 /= np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return math.degrees(math.atan2(y, x))


if __name__ == "__main__":
    main()
