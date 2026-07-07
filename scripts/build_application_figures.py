#!/usr/bin/env python3
"""Build the application figures used by the SONIC manuscript.

The default mode redraws the figures from the compact JSON data stored in the
repository.  Use ``--from-qm`` to rebuild the JSON data from the local Gaussian
DV HF/STO-3G calculation files and the adjacent MATRIX checkout.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
import tempfile
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "application_results.json"
FIGURES = ROOT / "figures"
QM_DIR = ROOT / "calculations" / "quick_qm"
MATRIX_ROOT = Path("/Users/vincenzobarone/Documents/git/software/matrix")
BOHR_TO_ANGSTROM = 0.529177210903
HARTREE_TO_KJMOL = 2625.499638


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--from-qm",
        action="store_true",
        help="rebuild data from local Gaussian logs/FCHK and MATRIX SONIC",
    )
    args = parser.parse_args()

    if args.from_qm:
        data = build_data_from_qm()
        DATA.parent.mkdir(parents=True, exist_ok=True)
        DATA.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    else:
        data = json.loads(DATA.read_text(encoding="utf-8"))

    FIGURES.mkdir(parents=True, exist_ok=True)
    plot_hessian_heatmap(data)
    plot_puckering_scan(data)


def build_data_from_qm() -> dict[str, object]:
    add_matrix_packages_to_path()

    from matrix_chem import preprocess_to_enriched_xyz, write_validation_section
    from matrix_neo import build_gic_b_matrix_from_xyzin, write_gicforge_build_sections

    fchk = parse_fchk(QM_DIR / "cyclobutane_hfsto3g.fchk")
    atom_symbols = atomic_symbols(fchk["atomic_numbers"])
    coordinates_angstrom = fchk["coordinates_bohr"] * BOHR_TO_ANGSTROM
    hessian_cartesian = lower_triangle_to_matrix(
        fchk["force_constants"],
        3 * len(atom_symbols),
    )

    with tempfile.TemporaryDirectory(prefix="neo-figures-") as scratch:
        scratch_path = Path(scratch)
        xyz = scratch_path / "cyclobutane_hfsto3g.xyz"
        xyz.write_text(
            xyz_text(atom_symbols, coordinates_angstrom, "cyclobutane HF/STO-3G"),
            encoding="utf-8",
        )

        matrices: dict[str, np.ndarray] = {}
        definitions: dict[str, object] = {}
        for key, symmetrize in (("neo_contract", False), ("neo_symmetry", True)):
            xyzin = scratch_path / f"{key}.xyzin"
            preprocess_to_enriched_xyz(xyz, xyzin)
            write_validation_section(xyzin)
            definitions[key] = write_gicforge_build_sections(xyzin, symmetrize=symmetrize)
            matrices[key] = np.asarray(build_gic_b_matrix_from_xyzin(xyzin).rows, dtype=float)

    cartesian_vibrational = vibrational_cartesian_hessian(
        hessian_cartesian,
        coordinates_angstrom,
    )
    neo_contract = internal_hessian(hessian_cartesian, matrices["neo_contract"])
    neo_symmetry = internal_hessian(hessian_cartesian, matrices["neo_symmetry"])

    hessian_panels = {
        "cartesian": coupling_matrix(cartesian_vibrational),
        "neo_contract": coupling_matrix(neo_contract),
        "neo_symmetry": coupling_matrix(neo_symmetry),
    }
    hessian_metrics = {
        key: coupling_metrics(value)
        for key, value in hessian_panels.items()
    }

    scan = parse_puckering_scan(QM_DIR)
    min_energy = min(point["energy_hartree"] for point in scan)
    for point in scan:
        point["relative_kj_mol"] = (
            point["energy_hartree"] - min_energy
        ) * HARTREE_TO_KJMOL

    return {
        "method": "Gaussian HF/STO-3G",
        "system": "cyclobutane",
        "hessian": {
            "rank": 30,
            "point_group": getattr(definitions["neo_symmetry"], "point_group", "D2d"),
            "panels": {
                key: np.round(value, 6).tolist()
                for key, value in hessian_panels.items()
            },
            "metrics": hessian_metrics,
        },
        "puckering": {
            "coordinate": "alternating C-ring puckering amplitude / Angstrom",
            "points": scan,
        },
    }


def add_matrix_packages_to_path() -> None:
    package_root = MATRIX_ROOT / "packages"
    src_paths = [
        str(path / "src")
        for path in sorted(package_root.iterdir())
        if (path / "src").is_dir()
    ]
    for path in reversed(src_paths):
        if path not in sys.path:
            sys.path.insert(0, path)


def parse_fchk(path: Path) -> dict[str, np.ndarray | list[int]]:
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()

    def scalar_int(label: str) -> int:
        for line in lines:
            if line.startswith(label):
                return int(line.split()[-1])
        raise KeyError(label)

    def array(label: str, count: int, kind: type = float) -> np.ndarray:
        for index, line in enumerate(lines):
            if line.startswith(label):
                values: list[float | int] = []
                cursor = index + 1
                while len(values) < count:
                    tokens = lines[cursor].split()
                    if kind is int:
                        values.extend(int(token) for token in tokens)
                    else:
                        values.extend(float(token.replace("D", "E")) for token in tokens)
                    cursor += 1
                return np.asarray(values[:count], dtype=int if kind is int else float)
        raise KeyError(label)

    atom_count = scalar_int("Number of atoms")
    cartesian_count = 3 * atom_count
    force_constant_count = cartesian_count * (cartesian_count + 1) // 2
    return {
        "atomic_numbers": array("Atomic numbers", atom_count, int).astype(int),
        "coordinates_bohr": array("Current cartesian coordinates", cartesian_count).reshape(-1, 3),
        "force_constants": array("Cartesian Force Constants", force_constant_count),
    }


def atomic_symbols(atomic_numbers: np.ndarray) -> tuple[str, ...]:
    symbols = {1: "H", 6: "C"}
    return tuple(symbols[int(number)] for number in atomic_numbers)


def xyz_text(symbols: tuple[str, ...], coords: np.ndarray, comment: str) -> str:
    rows = [str(len(symbols)), comment]
    rows.extend(
        f"{symbol} {xyz[0]:.10f} {xyz[1]:.10f} {xyz[2]:.10f}"
        for symbol, xyz in zip(symbols, coords)
    )
    return "\n".join(rows) + "\n"


def lower_triangle_to_matrix(values: np.ndarray, size: int) -> np.ndarray:
    matrix = np.zeros((size, size), dtype=float)
    cursor = 0
    for row in range(size):
        for column in range(row + 1):
            matrix[row, column] = values[cursor]
            matrix[column, row] = values[cursor]
            cursor += 1
    return matrix


def vibrational_cartesian_hessian(hessian: np.ndarray, coords: np.ndarray) -> np.ndarray:
    center = coords.mean(axis=0)
    external_rows: list[np.ndarray] = []
    for axis in range(3):
        vector = np.zeros_like(coords)
        vector[:, axis] = 1.0
        external_rows.append(vector.ravel())
    for axis_vector in np.eye(3):
        external_rows.append(np.cross(np.tile(axis_vector, (len(coords), 1)), coords - center).ravel())

    external = np.asarray(external_rows)
    _, _, vt = np.linalg.svd(external, full_matrices=True)
    vibrational_basis = vt[6:].T
    return vibrational_basis.T @ hessian @ vibrational_basis


def internal_hessian(hessian: np.ndarray, b_matrix: np.ndarray) -> np.ndarray:
    inverse = b_matrix.T @ np.linalg.pinv(b_matrix @ b_matrix.T, rcond=1.0e-10)
    return inverse.T @ hessian @ inverse


def coupling_matrix(force_constants: np.ndarray) -> np.ndarray:
    diagonal_scale = np.sqrt(np.maximum(np.abs(np.diag(force_constants)), 1.0e-20))
    coupling = np.abs(force_constants) / diagonal_scale[:, None] / diagonal_scale[None, :]
    np.fill_diagonal(coupling, 0.0)
    return np.clip(coupling, 0.0, 1.0)


def coupling_metrics(coupling: np.ndarray) -> dict[str, float]:
    upper = coupling[np.triu_indices(coupling.shape[0], 1)]
    return {
        "mean_off_diagonal": round(float(upper.mean()), 5),
        "p95_off_diagonal": round(float(np.quantile(upper, 0.95)), 5),
        "max_off_diagonal": round(float(upper.max()), 5),
    }


def parse_puckering_scan(path: Path) -> list[dict[str, float]]:
    points: list[dict[str, float]] = []
    for log_path in sorted(path.glob("cyclobutane_hf_pucker_p*.log")):
        match = re.search(r"_p([0-9])p([0-9]+)\.log$", log_path.name)
        if match is None:
            continue
        q_value = float(f"{match.group(1)}.{match.group(2)}")
        if q_value > 0.12:
            continue
        energy = parse_scf_energy(log_path)
        points.append({"q_angstrom": q_value, "energy_hartree": energy})
    return sorted(points, key=lambda point: point["q_angstrom"])


def parse_scf_energy(path: Path) -> float:
    pattern = re.compile(r"SCF Done:\s+E\([^)]+\)\s+=\s+([-+0-9.Ee]+)")
    energy = None
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        match = pattern.search(line)
        if match is not None:
            energy = float(match.group(1).replace("D", "E"))
    if energy is None:
        raise ValueError(f"No SCF energy found in {path}")
    return energy


def plot_hessian_heatmap(data: dict[str, object]) -> None:
    hessian = data["hessian"]
    panels = hessian["panels"]
    metrics = hessian["metrics"]
    titles = [
        ("cartesian", "Cartesian vibrational"),
        ("neo_contract", "SONIC contract"),
        ("neo_symmetry", "Symmetry-adapted SONIC"),
    ]

    fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.35), constrained_layout=True)
    for ax, (key, title) in zip(axes, titles):
        matrix = np.asarray(panels[key], dtype=float)
        im = ax.imshow(matrix, cmap="magma_r", vmin=0.0, vmax=0.65, interpolation="nearest")
        ax.set_title(title, fontsize=10)
        ax.set_xticks([])
        ax.set_yticks([])
        metric = metrics[key]
        ax.text(
            0.03,
            0.96,
            f"mean={metric['mean_off_diagonal']:.3f}\np95={metric['p95_off_diagonal']:.3f}",
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=7.5,
            color="black",
            bbox={"boxstyle": "round,pad=0.2", "facecolor": "white", "edgecolor": "0.8", "alpha": 0.9},
        )
    colorbar = fig.colorbar(im, ax=axes, fraction=0.028, pad=0.02)
    colorbar.set_label(r"$|F_{ij}|/\sqrt{|F_{ii}F_{jj}|}$", fontsize=9)
    fig.suptitle("Cyclobutane HF/STO-3G Hessian: normalized off-diagonal coupling", fontsize=11)
    fig.savefig(FIGURES / "sonic_hessian_heatmap.png", dpi=320)
    plt.close(fig)


def plot_puckering_scan(data: dict[str, object]) -> None:
    points = data["puckering"]["points"]
    q = np.asarray([point["q_angstrom"] for point in points], dtype=float)
    energy = np.asarray([point["relative_kj_mol"] for point in points], dtype=float)

    fig, ax = plt.subplots(figsize=(5.4, 3.2), constrained_layout=True)
    ax.plot(q, energy, color="#1f6f78", lw=2.2)
    ax.scatter(q, energy, s=42, color="#f2a03d", edgecolor="#1f1f1f", linewidth=0.6, zorder=3)
    minimum = int(np.argmin(energy))
    ax.scatter([q[minimum]], [energy[minimum]], s=80, marker="D", color="#b52b3a", zorder=4)
    ax.axvline(q[minimum], color="#b52b3a", lw=1.0, ls="--", alpha=0.8)
    ax.set_xlabel(r"ring puckering amplitude $q$ / $\AA$", fontsize=10)
    ax.set_ylabel(r"$\Delta E$ / kJ mol$^{-1}$", fontsize=10)
    ax.set_title("Frozen-coordinate puckering branch", fontsize=10)
    ax.grid(True, color="0.88", linewidth=0.8)
    ax.set_xlim(-0.004, 0.124)
    ax.set_ylim(bottom=-0.4)
    ax.tick_params(labelsize=9)
    ax.text(
        0.98,
        0.94,
        "Gaussian\nHF/STO-3G",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=8,
        bbox={"boxstyle": "round,pad=0.25", "facecolor": "white", "edgecolor": "0.82", "alpha": 0.92},
    )
    fig.savefig(FIGURES / "sonic_puckering_scan.png", dpi=320)
    plt.close(fig)


if __name__ == "__main__":
    main()
