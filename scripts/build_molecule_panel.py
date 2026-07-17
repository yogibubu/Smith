#!/usr/bin/env python3
"""Build a compact 3D molecule panel for the SONIC manuscript."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401


ROOT = Path(__file__).resolve().parents[1]
FIGURES = ROOT / "figures"

ATOM_COLORS = {
    "H": "#e8e8e8",
    "C": "#333333",
    "O": "#d94738",
    "Fe": "#d08a20",
    "Cl": "#3aa657",
    "Pd": "#8a86a8",
    "X": "#26a6a1",
}
ATOM_SIZES = {"H": 18, "C": 44, "O": 54, "Fe": 90, "Cl": 64, "Pd": 92, "X": 42}
COVALENT_RADII = {"H": 0.31, "C": 0.76, "O": 0.66, "Fe": 1.32, "Cl": 1.02, "Pd": 1.39}


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


def orient(coords: np.ndarray, mode: str = "pca") -> np.ndarray:
    centered = coords - coords.mean(axis=0)
    if mode == "cyclohexane":
        oriented = centered[:, [0, 2, 1]]
    else:
        _, _, vt = np.linalg.svd(centered, full_matrices=False)
        oriented = centered @ vt.T
    span = np.ptp(oriented, axis=0)
    scale = max(float(span.max()), 1e-6)
    return oriented / scale


def draw_molecule(
    ax,
    symbols: list[str],
    coords: np.ndarray,
    bonds: list[tuple[int, int]],
    title: str,
    subtitle: str = "",
    orientation: str = "pca",
    special_bonds: list[tuple[int, int]] | None = None,
) -> None:
    xyz = orient(coords, orientation)
    special = {tuple(sorted(pair)) for pair in (special_bonds or [])}
    for i, j in bonds:
        is_special = tuple(sorted((i, j))) in special
        ax.plot(
            [xyz[i, 0], xyz[j, 0]],
            [xyz[i, 1], xyz[j, 1]],
            [xyz[i, 2], xyz[j, 2]],
            color="#138f8a" if is_special else "#626262",
            lw=1.55 if is_special else 1.35,
            ls="--" if is_special else "-",
            solid_capstyle="round",
            zorder=1,
        )
    for symbol in sorted(
        set(symbols),
        key=lambda s: {"Pd": 0, "Fe": 1, "Cl": 2, "O": 3, "C": 4, "X": 5, "H": 6}.get(s, 9),
    ):
        idx = [i for i, atom in enumerate(symbols) if atom == symbol]
        ax.scatter(
            xyz[idx, 0],
            xyz[idx, 1],
            xyz[idx, 2],
            s=ATOM_SIZES.get(symbol, 26),
            c=ATOM_COLORS.get(symbol, "#7a7a7a"),
            edgecolors="#1e1e1e" if symbol not in {"H", "X"} else "#16837f" if symbol == "X" else "#aaaaaa",
            linewidths=0.35,
            depthshade=True,
            zorder=2,
        )
    ax.set_box_aspect((1, 1, 0.72))
    ax.axis("off")
    title_text = f"{title}\n{subtitle}" if subtitle else title
    ax.set_title(title_text, fontsize=8.8, fontweight="bold", pad=2.0, linespacing=1.35)
    pad = 0.12
    for setter, values in (
        (ax.set_xlim, xyz[:, 0]),
        (ax.set_ylim, xyz[:, 1]),
        (ax.set_zlim, xyz[:, 2]),
    ):
        setter(values.min() - pad, values.max() + pad)
    ax.view_init(elev=22, azim=-58)


def main() -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    panels = [
        (
            "Norbornane",
            "bridged C$_{2v}$ test",
            ROOT / "calculations/gdv_oop_probe/norbornane.xyzin",
            None,
            None,
        ),
        (
            "Camphor",
            "rigid bridged ketone",
            ROOT / "calculations/g16_improper_fix/camphor.xyzin",
            None,
            None,
        ),
        (
            "Ferrocene",
            "special-center coordinates",
            ROOT / "calculations/special_coordinate_example/ferrocene.xyzin",
            None,
            None,
        ),
        (
            "Formic acid--water",
            "interfragment coordinates",
            ROOT / "standalone/examples/formic_acid_water.xyzin",
            None,
            "weak_complex",
        ),
        (
            "$\\eta^3$-Allyl--PdCl",
            "ORACLE ligand center",
            ROOT / "standalone/examples/eta3_allyl_palladium.oracle.xyzin",
            None,
            "eta3_center",
        ),
        (
            "Cyclohexane",
            "$\\lambda=0.00$ chair",
            ROOT / "calculations/cyclohexane_puckering_equivalence/p00/cyclohexane_p00.xyz",
            "cyclohexane",
            None,
        ),
        (
            "Cyclohexane",
            "$\\lambda=0.50$ path",
            ROOT / "calculations/cyclohexane_puckering_equivalence/p02/cyclohexane_p02.xyz",
            "cyclohexane",
            None,
        ),
        (
            "Cyclohexane",
            "$\\lambda=1.00$ boat",
            ROOT / "calculations/cyclohexane_puckering_equivalence/p04/cyclohexane_p04.xyz",
            "cyclohexane",
            None,
        ),
    ]
    fig = plt.figure(figsize=(9.6, 5.55), dpi=320)
    axes = [fig.add_subplot(2, 4, index + 1, projection="3d") for index in range(8)]
    for ax, (title, subtitle, path, projection, special_kind) in zip(axes, panels, strict=True):
        symbols, coords = read_xyz_like(path)
        bonds = infer_bonds(symbols, coords) if projection == "cyclohexane" else read_bonds(path)
        special_bonds: list[tuple[int, int]] = []
        if special_kind == "weak_complex":
            # The two dashed contacts identify the intermolecular relation;
            # the protected contract itself contains relative translations and rotations.
            special_bonds = [(4, 5), (2, 7)]
            bonds.extend(special_bonds)
        elif special_kind == "eta3_center":
            # ORACLE supplies the centroid of the three allylic carbon atoms.
            center_index = len(symbols)
            symbols.append("X")
            coords = np.vstack((coords, coords[[2, 3, 4]].mean(axis=0)))
            special_bonds = [(0, center_index)]
            bonds.extend(special_bonds)
        draw_molecule(
            ax,
            symbols,
            coords,
            bonds,
            title,
            subtitle,
            projection,
            special_bonds=special_bonds,
        )

    fig.suptitle(
        "Representative molecular validation systems",
        fontsize=10.2,
        fontweight="bold",
        y=0.985,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.94), h_pad=1.85, w_pad=0.08)
    fig.savefig(FIGURES / "validation_molecule_panel.png", bbox_inches="tight", pad_inches=0.03)


if __name__ == "__main__":
    main()
