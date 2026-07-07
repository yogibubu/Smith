#!/usr/bin/env python3
"""Build a compact molecule panel for the SONIC manuscript."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
FIGURES = ROOT / "figures"

ATOM_COLORS = {
    "H": "#e8e8e8",
    "C": "#333333",
    "O": "#d94738",
    "Fe": "#d08a20",
}
ATOM_SIZES = {"H": 12, "C": 28, "O": 34, "Fe": 56}
COVALENT_RADII = {"H": 0.31, "C": 0.76, "O": 0.66, "Fe": 1.32}


def read_xyz_like(path: Path) -> tuple[list[str], np.ndarray]:
    lines = path.read_text(encoding="utf-8", errors="ignore").splitlines()
    n = int(lines[0].split()[0])
    symbols: list[str] = []
    coords: list[list[float]] = []
    for line in lines[2 : 2 + n]:
        parts = line.split()
        symbols.append(parts[0])
        coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
    return symbols, np.array(coords, dtype=float)


def read_bonds(path: Path) -> list[tuple[int, int]]:
    bonds: list[tuple[int, int]] = []
    in_bonds = False
    for raw in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        line = raw.strip()
        if line == "[BONDS]":
            in_bonds = True
            continue
        if in_bonds and (not line or line.startswith("[") or line.startswith("#")):
            break
        if in_bonds:
            parts = line.split()
            if len(parts) >= 2 and parts[0].isdigit() and parts[1].isdigit():
                bonds.append((int(parts[0]) - 1, int(parts[1]) - 1))
    return bonds


def infer_bonds(symbols: list[str], coords: np.ndarray) -> list[tuple[int, int]]:
    bonds: list[tuple[int, int]] = []
    for i in range(len(symbols)):
        for j in range(i + 1, len(symbols)):
            ri = COVALENT_RADII.get(symbols[i], 0.75)
            rj = COVALENT_RADII.get(symbols[j], 0.75)
            cutoff = 1.22 * (ri + rj) + 0.15
            if np.linalg.norm(coords[i] - coords[j]) <= cutoff:
                bonds.append((i, j))
    return bonds


def project(coords: np.ndarray, mode: str = "pca") -> np.ndarray:
    centered = coords - coords.mean(axis=0)
    if mode == "xz":
        xy = centered[:, [0, 2]]
    else:
        _, _, vt = np.linalg.svd(centered, full_matrices=False)
        xy = centered @ vt[:2].T
    span = np.ptp(xy, axis=0)
    scale = max(float(span.max()), 1e-6)
    return xy / scale


def draw_molecule(
    ax,
    symbols: list[str],
    coords: np.ndarray,
    bonds: list[tuple[int, int]],
    title: str,
    subtitle: str = "",
    projection: str = "pca",
) -> None:
    xy = project(coords, projection)
    for i, j in bonds:
        ax.plot(
            [xy[i, 0], xy[j, 0]],
            [xy[i, 1], xy[j, 1]],
            color="#6c6c6c",
            lw=1.15,
            solid_capstyle="round",
            zorder=1,
        )
    for symbol in sorted(set(symbols), key=lambda s: {"Fe": 0, "O": 1, "C": 2, "H": 3}.get(s, 9)):
        idx = [i for i, atom in enumerate(symbols) if atom == symbol]
        ax.scatter(
            xy[idx, 0],
            xy[idx, 1],
            s=ATOM_SIZES.get(symbol, 26),
            c=ATOM_COLORS.get(symbol, "#7a7a7a"),
            edgecolors="#1e1e1e" if symbol != "H" else "#aaaaaa",
            linewidths=0.35,
            zorder=2,
        )
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(title, fontsize=9.0, fontweight="bold", pad=2.5)
    if subtitle:
        ax.text(
            0.5,
            -0.03,
            subtitle,
            transform=ax.transAxes,
            ha="center",
            va="top",
            fontsize=7.2,
            color="#333333",
        )
    pad = 0.12
    ax.set_xlim(xy[:, 0].min() - pad, xy[:, 0].max() + pad)
    ax.set_ylim(xy[:, 1].min() - pad, xy[:, 1].max() + pad)


def main() -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    panels = [
        (
            "Norbornane",
            "bridged C$_{2v}$ test",
            ROOT / "calculations/gdv_oop_probe/norbornane.xyzin",
            None,
        ),
        (
            "Camphor",
            "rigid bridged ketone",
            ROOT / "calculations/g16_improper_fix/camphor.xyzin",
            None,
        ),
        (
            "Ferrocene",
            "special-center coordinates",
            ROOT / "calculations/special_coordinate_example/ferrocene.xyzin",
            None,
        ),
        (
            "Cyclohexane",
            "$\\lambda=0.00$ chair",
            ROOT / "calculations/cyclohexane_puckering_equivalence/p00/cyclohexane_p00.xyz",
            "xz",
        ),
        (
            "Cyclohexane",
            "$\\lambda=0.50$ path",
            ROOT / "calculations/cyclohexane_puckering_equivalence/p02/cyclohexane_p02.xyz",
            "xz",
        ),
        (
            "Cyclohexane",
            "$\\lambda=1.00$ boat",
            ROOT / "calculations/cyclohexane_puckering_equivalence/p04/cyclohexane_p04.xyz",
            "xz",
        ),
    ]
    fig, axes = plt.subplots(2, 3, figsize=(7.2, 5.15), dpi=320)
    for ax, (title, subtitle, path, projection) in zip(axes.flat, panels, strict=True):
        symbols, coords = read_xyz_like(path)
        bonds = infer_bonds(symbols, coords) if projection == "xz" else read_bonds(path)
        draw_molecule(ax, symbols, coords, bonds, title, subtitle, projection)

    fig.suptitle(
        "Representative molecular systems used in the validation probes",
        fontsize=10.5,
        fontweight="bold",
        y=0.985,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.955), h_pad=1.55, w_pad=0.1)
    fig.savefig(FIGURES / "validation_molecule_panel.png", bbox_inches="tight", pad_inches=0.03)


if __name__ == "__main__":
    main()
