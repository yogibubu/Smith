#!/usr/bin/env python3
"""Build an artistic, text-free graphical TOC for the SMITH manuscript."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Arc, Circle, PathPatch, Polygon
from matplotlib.path import Path as MplPath
import numpy as np

from build_molecule_panel import ATOM_COLORS, infer_bonds, read_xyz_like


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "figures" / "graphical_toc.png"

NAVY = "#102d46"
SAPPHIRE = "#155d8c"
STEEL = "#95b6ca"
AMBER = "#e0a12b"
TEAL = "#23a89f"
VIOLET = "#8464bc"
BLUE = "#3f83c5"
WHITE = "#fcfbf8"


def _project(coords: np.ndarray) -> np.ndarray:
    centered = coords - coords.mean(axis=0)
    az, el = np.deg2rad(-20.0), np.deg2rad(28.0)
    rz = np.array([[np.cos(az), -np.sin(az), 0], [np.sin(az), np.cos(az), 0], [0, 0, 1]])
    rx = np.array([[1, 0, 0], [0, np.cos(el), -np.sin(el)], [0, np.sin(el), np.cos(el)]])
    xy = (centered @ (rz @ rx).T)[:, :2]
    return xy / max(float(np.ptp(xy, axis=0).max()), 1.0e-8)


def _molecule(ax, symbols: list[str], xy: np.ndarray, bonds: list[tuple[int, int]], *, center, scale, alpha=1.0):
    pts = xy * scale + np.asarray(center)
    for i, j in bonds:
        ax.plot(
            [pts[i, 0], pts[j, 0]],
            [pts[i, 1], pts[j, 1]],
            color="#647580",
            lw=2.0,
            alpha=0.78 * alpha,
            solid_capstyle="round",
            zorder=7,
        )
    order = sorted(range(len(symbols)), key=lambda i: 0 if symbols[i] != "H" else 1)
    for i in order:
        symbol = symbols[i]
        radius = 0.014 if symbol == "H" else 0.022 if symbol in {"C", "N"} else 0.027
        ax.add_patch(
            Circle(
                pts[i],
                radius,
                facecolor=ATOM_COLORS.get(symbol, "#788995"),
                edgecolor="#213746" if symbol != "H" else "#aebbc2",
                lw=0.75,
                alpha=alpha,
                zorder=9,
            )
        )
    return pts


def _bezier(ax, vertices, color, lw, alpha=1.0, zorder=5):
    codes = [MplPath.MOVETO] + [MplPath.CURVE4] * (len(vertices) - 1)
    ax.add_patch(
        PathPatch(
            MplPath(vertices, codes),
            fill=False,
            edgecolor=color,
            lw=lw,
            alpha=alpha,
            capstyle="round",
            zorder=zorder,
        )
    )


def _background(ax):
    nx, ny = 900, 480
    x = np.linspace(0, 1, nx)
    y = np.linspace(0, 1, ny)
    xx, yy = np.meshgrid(x, y)
    glow = np.exp(-(((xx - 0.50) / 0.20) ** 2 + ((yy - 0.53) / 0.48) ** 2))
    right = np.exp(-(((xx - 0.78) / 0.32) ** 2 + ((yy - 0.50) / 0.60) ** 2))
    base = np.ones((ny, nx, 3))
    warm = np.array([0.992, 0.988, 0.972])
    blue = np.array([0.79, 0.89, 0.95])
    base[:] = warm
    base = base * (1 - 0.20 * glow[..., None]) + blue * (0.20 * glow[..., None])
    base = base * (1 - 0.045 * right[..., None]) + np.array([0.91, 0.96, 0.98]) * (0.045 * right[..., None])
    ax.imshow(base, extent=(0, 1, 0, 1), origin="lower", aspect="auto", zorder=0)


def _forge(ax):
    # A faceted optical forge: part prism, part anvil, without literal machinery.
    back = Polygon(
        [(0.455, 0.27), (0.405, 0.38), (0.425, 0.72), (0.475, 0.84), (0.545, 0.78), (0.575, 0.50), (0.545, 0.28)],
        closed=True,
        facecolor=NAVY,
        edgecolor="#0b2438",
        lw=1.2,
        zorder=4,
    )
    ax.add_patch(back)
    facets = [
        ([(0.425, 0.72), (0.475, 0.84), (0.495, 0.52), (0.405, 0.38)], "#267caf", 0.76),
        ([(0.475, 0.84), (0.545, 0.78), (0.495, 0.52)], "#174f78", 0.88),
        ([(0.495, 0.52), (0.545, 0.78), (0.575, 0.50), (0.545, 0.28)], "#0f3d61", 0.93),
        ([(0.405, 0.38), (0.495, 0.52), (0.545, 0.28), (0.455, 0.27)], "#1a668e", 0.79),
    ]
    for vertices, color, alpha in facets:
        ax.add_patch(Polygon(vertices, closed=True, facecolor=color, edgecolor=STEEL, lw=0.55, alpha=alpha, zorder=5))
    ax.add_patch(Polygon([(0.39, 0.25), (0.56, 0.25), (0.61, 0.19), (0.43, 0.17), (0.36, 0.21)], closed=True,
                         facecolor="#153e59", edgecolor="#8eb9cf", lw=0.8, alpha=0.92, zorder=3))


def _matrix_lattice(ax):
    rng = np.random.default_rng(31)
    for row in range(7):
        for col in range(13):
            if rng.random() < 0.62:
                continue
            x = 0.395 + col * 0.017
            y = 0.075 + row * 0.015
            intensity = 0.35 + 0.55 * rng.random()
            ax.scatter([x], [y], s=5.0, marker="s", color=STEEL, alpha=intensity, linewidths=0, zorder=2)


def main() -> None:
    OUT.parent.mkdir(parents=True, exist_ok=True)
    symbols, coords = read_xyz_like(ROOT / "standalone/examples/saccharin.smith.xyz")
    bonds = infer_bonds(symbols, coords)
    xy = _project(coords)

    fig, ax = plt.subplots(figsize=(8.6, 4.65), dpi=360)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    _background(ax)

    # Many faint primitive traces approach the forge from the raw geometry.
    left_pts = _molecule(ax, symbols, xy, bonds, center=(0.205, 0.53), scale=0.31, alpha=0.92)
    rng = np.random.default_rng(17)
    for idx in rng.choice(len(bonds), size=min(12, len(bonds)), replace=False):
        i, j = bonds[int(idx)]
        mid = 0.5 * (left_pts[i] + left_pts[j])
        _bezier(
            ax,
            [tuple(left_pts[i]), tuple(mid + rng.normal(0, 0.035, 2)), tuple(mid + rng.normal(0, 0.035, 2)), tuple(left_pts[j])],
            "#97a4ab",
            0.85,
            alpha=0.34,
            zorder=3,
        )
    for y0 in np.linspace(0.35, 0.71, 7):
        _bezier(ax, [(0.31, y0), (0.36, y0 + 0.05), (0.39, 0.53), (0.435, 0.53)], "#7890a0", 0.8, alpha=0.30, zorder=2)

    _matrix_lattice(ax)
    _forge(ax)

    # A small number of ordered, family-colored coordinate ribbons emerge.
    for offset, color, width in ((0.10, AMBER, 3.4), (0.035, TEAL, 3.1), (-0.035, VIOLET, 3.0), (-0.10, BLUE, 2.8)):
        _bezier(
            ax,
            [(0.545, 0.52 + offset * 0.25), (0.62, 0.57 + offset), (0.69, 0.57 + offset), (0.735, 0.53 + offset * 0.35)],
            color,
            width,
            alpha=0.91,
            zorder=6,
        )

    right_pts = _molecule(ax, symbols, xy, bonds, center=(0.805, 0.53), scale=0.31, alpha=1.0)
    ring_atoms = [i for i, atom in enumerate(symbols) if atom in {"C", "N", "S"}]
    ring_center = right_pts[ring_atoms].mean(axis=0)
    ax.add_patch(Arc(ring_center, 0.31, 0.20, angle=8, theta1=8, theta2=174, color=VIOLET, lw=3.0, alpha=0.92, zorder=10))
    if bonds:
        i, j = bonds[0]
        ax.plot([right_pts[i, 0], right_pts[j, 0]], [right_pts[i, 1], right_pts[j, 1]], color=AMBER, lw=4.0, alpha=0.92, zorder=8)
    ax.add_patch(Arc(ring_center + np.array([0.11, -0.02]), 0.17, 0.17, theta1=35, theta2=305, color=TEAL, lw=2.7, alpha=0.90, zorder=10))

    fig.tight_layout(pad=0)
    fig.savefig(OUT, bbox_inches="tight", pad_inches=0.02, facecolor=WHITE)


if __name__ == "__main__":
    main()
