#!/usr/bin/env python3
"""Cyclohexane chair-to-boat export equivalence probe.

The script builds a fixed-frame puckering interpolation, writes both portable
Cartesian Gaussian inputs and SONIC ReadAllGIC inputs, optionally runs Gaussian,
and stores a compact JSON summary plus the manuscript figure.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from _matrix_checkout import add_matrix_packages_to_path, matrix_root


ROOT = Path(__file__).resolve().parents[1]
RUN_DIR = ROOT / "calculations" / "cyclohexane_puckering_equivalence"
DATA = ROOT / "data" / "cyclohexane_puckering_equivalence.json"
FIGURE = ROOT / "figures" / "cyclohexane_puckering_equivalence.png"
DEFAULT_GAUSSIAN = Path("g16")
ROUTE = "#p b3lyp/6-31+g* nosymm"
HARTREE_TO_KJMOL = 2625.499638
SCF_RE = re.compile(r"SCF Done:\s+E\([^)]+\)\s+=\s+([-+0-9.Ee]+)")


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

    data = summarize_results()
    DATA.parent.mkdir(parents=True, exist_ok=True)
    DATA.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    FIGURE.parent.mkdir(parents=True, exist_ok=True)
    plot(data)


def prepare_inputs() -> None:
    add_matrix_packages_to_path()

    from matrix_chem import preprocess_to_enriched_xyz, write_validation_section
    from matrix_gaussian import write_gicforge_gaussian_input
    from matrix_neo import gaussian_gic_lines_from_xyzin, write_gicforge_build_sections

    RUN_DIR.mkdir(parents=True, exist_ok=True)
    atoms, chair = read_cyclohexane_fixture()
    ring = carbon_ring_indices(atoms)
    hydrogens = attach_hydrogens(atoms, chair, ring)
    plane_center, normal, chair_z = ring_plane(chair, ring)
    if chair_z[0] < 0.0:
        normal *= -1.0
        chair_z *= -1.0
    chair_pattern = chair_z / np.sqrt(np.mean(chair_z**2))
    boat_pattern = np.asarray([1.0, -0.5, -0.5, 1.0, -0.5, -0.5], dtype=float)
    boat_pattern /= np.sqrt(np.mean(boat_pattern**2))
    amplitude = float(np.sqrt(np.mean(chair_z**2)))

    base = chair.copy()
    for local, atom_index in enumerate(ring):
        base[atom_index] = chair[atom_index] - chair_z[local] * normal
    for hydrogen, carbon in hydrogens.items():
        local = ring.index(carbon)
        base[hydrogen] = chair[hydrogen] - chair_z[local] * normal

    for index, lam in enumerate(np.linspace(0.0, 1.0, 5)):
        pattern = (1.0 - lam) * chair_pattern + lam * boat_pattern
        coords = base.copy()
        for local, atom_index in enumerate(ring):
            coords[atom_index] = base[atom_index] + amplitude * pattern[local] * normal
        for hydrogen, carbon in hydrogens.items():
            local = ring.index(carbon)
            coords[hydrogen] = base[hydrogen] + amplitude * pattern[local] * normal

        point_dir = RUN_DIR / f"p{index:02d}"
        point_dir.mkdir(parents=True, exist_ok=True)
        xyz = point_dir / f"cyclohexane_p{index:02d}.xyz"
        xyzin = point_dir / f"cyclohexane_p{index:02d}.xyzin"
        xyz.write_text(xyz_text(atoms, coords, f"cyclohexane chair-boat lambda={lam:.2f}"))
        preprocess_to_enriched_xyz(xyz, xyzin)
        write_validation_section(xyzin)
        definition = write_gicforge_build_sections(xyzin, symmetrize=False)
        gic_lines = gaussian_gic_lines_from_xyzin(xyzin)
        unsupported = [line for line in gic_lines if "U(" in line or "OuPl" in line]
        if unsupported:
            raise RuntimeError(f"Point {index} has out-of-plane GIC lines: {unsupported[:3]}")

        standard = point_dir / f"cyclohexane_p{index:02d}_cart.gjf"
        standard.write_text(
            gaussian_cartesian_input(
                atoms,
                coords,
                title=f"cyclohexane Cartesian puckering point {index}",
                chk=f"cyclohexane_p{index:02d}_cart.chk",
            ),
            encoding="utf-8",
        )
        write_gicforge_gaussian_input(
            xyzin,
            point_dir / f"cyclohexane_p{index:02d}_neo_gic.gjf",
            route=ROUTE,
            title=f"cyclohexane SONIC ReadAllGIC puckering point {index}",
            link0=(
                f"%chk=cyclohexane_p{index:02d}_neo_gic.chk",
                "%nprocshared=4",
                "%mem=3GB",
            ),
        )

        metadata = {
            "lambda": float(lam),
            "q3_chair_to_q2_boat": float(lam),
            "point_group": definition.point_group,
            "rank": definition.rank,
            "target_rank": definition.target_rank,
            "gic_count": len(gic_lines),
        }
        (point_dir / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")


def read_cyclohexane_fixture() -> tuple[tuple[str, ...], np.ndarray]:
    source = (
        matrix_root()
        / "tests"
        / "fixtures"
        / "test_molecules"
        / "molecules"
        / "cyclohexane.inp"
    )
    atoms: list[str] = []
    coords: list[tuple[float, float, float]] = []
    for line in source.read_text(encoding="utf-8").splitlines():
        parts = line.split()
        if len(parts) == 4 and parts[0] in {"C", "H", "6", "1"}:
            atoms.append("C" if parts[0] == "6" else "H" if parts[0] == "1" else parts[0])
            coords.append(tuple(float(value) for value in parts[1:]))
    return tuple(atoms), np.asarray(coords, dtype=float)


def carbon_ring_indices(atoms: tuple[str, ...]) -> tuple[int, ...]:
    return tuple(index for index, atom in enumerate(atoms) if atom == "C")


def attach_hydrogens(
    atoms: tuple[str, ...],
    coords: np.ndarray,
    ring: tuple[int, ...],
) -> dict[int, int]:
    attachments: dict[int, int] = {}
    for index, atom in enumerate(atoms):
        if atom != "H":
            continue
        nearest = min(ring, key=lambda carbon: float(np.linalg.norm(coords[index] - coords[carbon])))
        attachments[index] = nearest
    return attachments


def ring_plane(coords: np.ndarray, ring: tuple[int, ...]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    ring_coords = coords[list(ring)]
    center = ring_coords.mean(axis=0)
    _, _, vh = np.linalg.svd(ring_coords - center)
    normal = vh[-1]
    z = np.asarray([np.dot(coords[index] - center, normal) for index in ring], dtype=float)
    return center, normal, z


def xyz_text(atoms: tuple[str, ...], coords: np.ndarray, comment: str) -> str:
    rows = [str(len(atoms)), comment]
    rows.extend(
        f"{atom:2s} {xyz[0]:15.10f} {xyz[1]:15.10f} {xyz[2]:15.10f}"
        for atom, xyz in zip(atoms, coords)
    )
    return "\n".join(rows) + "\n"


def gaussian_cartesian_input(
    atoms: tuple[str, ...],
    coords: np.ndarray,
    *,
    title: str,
    chk: str,
) -> str:
    rows = [
        f"%chk={chk}",
        "%nprocshared=4",
        "%mem=3GB",
        ROUTE,
        "",
        title,
        "",
        "0 1",
    ]
    rows.extend(
        f"{atom:2s} {xyz[0]:15.10f} {xyz[1]:15.10f} {xyz[2]:15.10f}"
        for atom, xyz in zip(atoms, coords)
    )
    return "\n".join(rows) + "\n\n"


def run_gaussian(executable: Path) -> None:
    if not executable.exists():
        raise FileNotFoundError(executable)
    env = os.environ.copy()
    env.setdefault("GAUSS_SCRDIR", str(Path.home() / "Documents" / "GAUSSIAN_SCRATCH"))
    for point_dir in sorted(RUN_DIR.glob("p[0-9][0-9]")):
        for input_path in sorted(point_dir.glob("*.gjf")):
            log_path = input_path.with_suffix(".log")
            chk_path = input_path.with_suffix(".chk")
            log_path.unlink(missing_ok=True)
            chk_path.unlink(missing_ok=True)
            (point_dir / "fort.7").unlink(missing_ok=True)
            subprocess.run([str(executable), input_path.name], cwd=point_dir, env=env, check=True)


def summarize_results() -> dict[str, object]:
    points: list[dict[str, object]] = []
    for point_dir in sorted(RUN_DIR.glob("p[0-9][0-9]")):
        metadata = json.loads((point_dir / "metadata.json").read_text(encoding="utf-8"))
        cart_energy = parse_energy(point_dir / f"cyclohexane_{point_dir.name}_cart.log")
        gic_energy = parse_energy(point_dir / f"cyclohexane_{point_dir.name}_neo_gic.log")
        points.append(
            {
                **metadata,
                "cartesian_energy_hartree": cart_energy,
                "readallgic_energy_hartree": gic_energy,
                "delta_microhartree": (gic_energy - cart_energy) * 1.0e6,
            }
        )
    reference = min(float(point["cartesian_energy_hartree"]) for point in points)
    for point in points:
        point["relative_cartesian_kj_mol"] = (
            float(point["cartesian_energy_hartree"]) - reference
        ) * HARTREE_TO_KJMOL
        point["relative_readallgic_kj_mol"] = (
            float(point["readallgic_energy_hartree"]) - reference
        ) * HARTREE_TO_KJMOL
    return {
        "method": "B3LYP/6-31+G*",
        "route_cartesian": ROUTE,
        "route_readallgic": ROUTE + " geom=readallgic output=pickett",
        "path": "fixed-frame q3(chair) to q2(boat) cyclohexane puckering interpolation",
        "points": points,
        "max_abs_delta_microhartree": max(abs(float(p["delta_microhartree"])) for p in points),
    }


def parse_energy(path: Path) -> float:
    text = path.read_text(encoding="utf-8", errors="replace")
    matches = list(SCF_RE.finditer(text))
    if "Normal termination of Gaussian" not in text:
        raise RuntimeError(f"Gaussian did not terminate normally: {path}")
    if not matches:
        raise RuntimeError(f"No SCF energy found in {path}")
    return float(matches[-1].group(1).replace("D", "E"))


def plot(data: dict[str, object]) -> None:
    points = data["points"]
    lambdas = np.asarray([point["lambda"] for point in points], dtype=float)
    cart = np.asarray([point["relative_cartesian_kj_mol"] for point in points], dtype=float)
    gic = np.asarray([point["relative_readallgic_kj_mol"] for point in points], dtype=float)
    delta = np.asarray([point["delta_microhartree"] for point in points], dtype=float)

    fig, (top, bottom) = plt.subplots(
        2,
        1,
        figsize=(5.7, 4.1),
        sharex=True,
        gridspec_kw={"height_ratios": [3.0, 1.0]},
        constrained_layout=True,
    )
    top.plot(lambdas, cart, color="#1f6f78", lw=2.2, marker="o", label="Cartesian")
    top.plot(
        lambdas,
        gic,
        color="#f2a03d",
        lw=0.0,
        marker="s",
        markerfacecolor="none",
        markeredgewidth=1.7,
        label="SONIC ReadAllGIC",
    )
    top.set_ylabel(r"$\Delta E$ / kJ mol$^{-1}$")
    top.set_title("Cyclohexane chair-to-boat fixed puckering path")
    top.grid(True, color="0.88", linewidth=0.8)
    top.legend(frameon=True, fontsize=8)
    bottom.axhline(0.0, color="0.35", lw=0.8)
    bottom.plot(lambdas, delta, color="#b52b3a", lw=1.7, marker="o", ms=4)
    bottom.set_xlabel(r"path parameter $\lambda$: chair $q_3$ to boat $q_2$")
    bottom.set_ylabel(r"$\Delta$ / $\mu E_h$")
    bottom.grid(True, color="0.9", linewidth=0.8)
    fig.savefig(FIGURE, dpi=320)
    plt.close(fig)


if __name__ == "__main__":
    main()
