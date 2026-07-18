#!/usr/bin/env python3
"""Build the introductory SONIC atlas and the special-coordinate gallery."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Arc, FancyArrowPatch
import numpy as np

from build_molecule_panel import (
    ATOM_COLORS,
    ATOM_SIZES,
    infer_bonds,
    read_bonds,
    read_xyz_like,
)


ROOT = Path(__file__).resolve().parents[1]
FIGURES = ROOT / "figures"

NAVY = "#12304a"
AMBER = "#d89422"
TEAL = "#15958f"
VIOLET = "#7857a6"
BLUE = "#397ab7"
SOFT_GREY = "#6b7379"
WARM_WHITE = "#fcfbf8"


def _rotation(azimuth: float, elevation: float) -> np.ndarray:
    az = np.deg2rad(azimuth)
    el = np.deg2rad(elevation)
    rz = np.array(
        [[np.cos(az), -np.sin(az), 0.0], [np.sin(az), np.cos(az), 0.0], [0.0, 0.0, 1.0]]
    )
    rx = np.array(
        [[1.0, 0.0, 0.0], [0.0, np.cos(el), -np.sin(el)], [0.0, np.sin(el), np.cos(el)]]
    )
    return rz @ rx


def project(coords: np.ndarray, *, azimuth: float = -34.0, elevation: float = 24.0) -> np.ndarray:
    centered = coords - coords.mean(axis=0)
    rotated = centered @ _rotation(azimuth, elevation).T
    xy = rotated[:, :2]
    scale = max(float(np.ptp(xy, axis=0).max()), 1.0e-8)
    return xy / scale


def molecule(
    ax,
    symbols: list[str],
    coords: np.ndarray,
    bonds: list[tuple[int, int]],
    *,
    azimuth: float = -34.0,
    elevation: float = 24.0,
    fragment_split: int | None = None,
) -> np.ndarray:
    xy = project(coords, azimuth=azimuth, elevation=elevation)
    for i, j in bonds:
        color = "#7b858b"
        if fragment_split is not None and (i < fragment_split) != (j < fragment_split):
            color = TEAL
        ax.plot(
            [xy[i, 0], xy[j, 0]],
            [xy[i, 1], xy[j, 1]],
            color=color,
            lw=1.5,
            solid_capstyle="round",
            zorder=1,
        )
    for symbol in sorted(set(symbols), key=lambda item: ATOM_SIZES.get(item, 25), reverse=True):
        idx = [i for i, atom in enumerate(symbols) if atom == symbol]
        ax.scatter(
            xy[idx, 0],
            xy[idx, 1],
            s=np.asarray([ATOM_SIZES.get(symbol, 26)] * len(idx)) * 1.55,
            c=ATOM_COLORS.get(symbol, "#7a7a7a"),
            edgecolors="#26343e" if symbol != "H" else "#a7afb4",
            linewidths=0.45,
            zorder=3,
        )
    # Keep the full subplot box for headings; expand the data limits instead of
    # vertically collapsing a wide molecular panel.
    ax.set_aspect("equal", adjustable="datalim")
    ax.axis("off")
    ax.set_xlim(float(xy[:, 0].min()) - 0.24, float(xy[:, 0].max()) + 0.24)
    ax.set_ylim(float(xy[:, 1].min()) - 0.22, float(xy[:, 1].max()) + 0.22)
    return xy


def title(ax, heading: str, subtitle: str) -> None:
    ax.text(
        0.02,
        0.98,
        heading,
        transform=ax.transAxes,
        ha="left",
        va="top",
        color=NAVY,
        fontsize=10.0,
        fontweight="bold",
    )
    ax.text(
        0.02,
        0.885,
        subtitle,
        transform=ax.transAxes,
        ha="left",
        va="top",
        color=SOFT_GREY,
        fontsize=8.0,
    )


def double_arrow(ax, start: np.ndarray, end: np.ndarray, color: str, *, shrink: float = 0.08) -> None:
    vec = end - start
    length = max(float(np.linalg.norm(vec)), 1.0e-8)
    unit = vec / length
    a = start + shrink * unit
    b = end - shrink * unit
    ax.add_patch(
        FancyArrowPatch(
            a,
            b,
            arrowstyle="<->",
            mutation_scale=11,
            linewidth=1.8,
            color=color,
            zorder=5,
        )
    )


def center_marker(ax, xy: np.ndarray, color: str = TEAL) -> None:
    ax.scatter([xy[0]], [xy[1]], marker="D", s=34, c=color, edgecolors="white", linewidths=0.7, zorder=6)


def build_atlas() -> None:
    fig, axes = plt.subplots(1, 4, figsize=(10.0, 2.95), dpi=320, facecolor=WARM_WHITE)
    for ax in axes:
        ax.set_facecolor(WARM_WHITE)

    # Stretch and bend: water.
    path = ROOT / "calculations/coordinate_comparison/water/water.xyzin"
    symbols, coords = read_xyz_like(path)
    bonds = read_bonds(path) or infer_bonds(symbols, coords)
    xy = molecule(axes[0], symbols, coords, bonds, azimuth=0.0, elevation=90.0)
    double_arrow(axes[0], xy[0], xy[1], AMBER)
    v0 = xy[0] - xy[1]
    v2 = xy[2] - xy[1]
    theta0 = float(np.degrees(np.arctan2(v0[1], v0[0])))
    theta2 = float(np.degrees(np.arctan2(v2[1], v2[0])))
    if theta2 < theta0:
        theta2 += 360.0
    axes[0].add_patch(
        Arc(xy=xy[1], width=0.48, height=0.48, theta1=theta0, theta2=theta2, color=TEAL, lw=2.0, zorder=5)
    )
    title(axes[0], "Local primitives", "stretch and valence bend")

    # Curvilinear ring out-of-plane source.
    path = ROOT / "calculations/quick_qm/cyclobutane_hfsto3g.xyzin"
    symbols, coords = read_xyz_like(path)
    bonds = read_bonds(path) or infer_bonds(symbols, coords)
    xy = molecule(axes[1], symbols, coords, bonds, azimuth=-18.0, elevation=25.0)
    carbon = [i for i, atom in enumerate(symbols) if atom == "C"][:4]
    for order, idx in enumerate(carbon):
        dy = 0.15 if order % 2 == 0 else -0.15
        axes[1].annotate(
            "",
            xy=xy[idx] + np.array([0.0, dy]),
            xytext=xy[idx],
            arrowprops={"arrowstyle": "-|>", "color": VIOLET, "lw": 1.7},
            zorder=6,
        )
    title(axes[1], "Curvilinear ring source", "$U$ out-of-plane and $RPck$")

    # Fused-ring butterfly and puckering phase.
    path = ROOT / "standalone/examples/saccharin.smith.xyz"
    symbols, coords = read_xyz_like(path)
    bonds = infer_bonds(symbols, coords)
    xy = molecule(axes[2], symbols, coords, bonds, azimuth=-12.0, elevation=28.0)
    center = xy[[i for i, atom in enumerate(symbols) if atom in {"C", "N", "S"}]].mean(axis=0)
    axes[2].add_patch(
        Arc(xy=center, width=0.82, height=0.46, angle=8, theta1=8, theta2=172, color=VIOLET, lw=2.2, zorder=5)
    )
    axes[2].add_patch(
        FancyArrowPatch(
            center + np.array([-0.35, 0.02]),
            center + np.array([0.35, -0.02]),
            connectionstyle="arc3,rad=-0.32",
            arrowstyle="-|>",
            mutation_scale=10,
            color=AMBER,
            lw=1.7,
            zorder=5,
        )
    )
    title(axes[2], "Fused-ring coordinates", "puckering phase and butterfly")

    # Fragment pose.
    path = ROOT / "standalone/examples/formic_acid_water.smith.xyz"
    symbols, coords = read_xyz_like(path)
    bonds = infer_bonds(symbols, coords)
    xy = molecule(axes[3], symbols, coords, bonds, azimuth=0.0, elevation=0.0, fragment_split=5)
    c1 = xy[:5].mean(axis=0)
    c2 = xy[5:].mean(axis=0)
    center_marker(axes[3], c1, BLUE)
    center_marker(axes[3], c2, BLUE)
    double_arrow(axes[3], c1, c2, BLUE, shrink=0.03)
    axes[3].add_patch(
        Arc(xy=c2, width=0.37, height=0.37, theta1=30, theta2=300, color=VIOLET, lw=1.7, zorder=5)
    )
    axes[3].plot([xy[4, 0], xy[5, 0]], [xy[4, 1], xy[5, 1]], ls="--", lw=1.5, color=TEAL, zorder=2)
    title(axes[3], "Disconnected fragments", "relative translation and rotation")

    fig.suptitle(
        "A visual atlas of SONIC source families",
        x=0.5,
        y=0.99,
        fontsize=12.2,
        fontweight="bold",
        color=NAVY,
    )
    fig.tight_layout(rect=(0.01, 0.01, 0.99, 0.92), w_pad=0.35)
    fig.savefig(FIGURES / "sonic_coordinate_atlas.png", bbox_inches="tight", pad_inches=0.05, facecolor=fig.get_facecolor())


def _special_panel(ax, title_text: str, subtitle: str, path: Path, *, az: float, el: float) -> tuple[list[str], np.ndarray, np.ndarray]:
    symbols, coords = read_xyz_like(path)
    bonds = read_bonds(path) or infer_bonds(symbols, coords)
    xy = molecule(ax, symbols, coords, bonds, azimuth=az, elevation=el)
    title(ax, title_text, subtitle)
    return symbols, coords, xy


def build_special_gallery() -> None:
    grid = plt.subplots(2, 3, figsize=(9.4, 5.25), dpi=320, facecolor=WARM_WHITE)
    fig, axes_grid = grid
    axes = axes_grid.ravel()
    for ax in axes:
        ax.set_facecolor(WARM_WHITE)

    # Ferrocene: two ring centers and center--metal rows.
    symbols, coords, xy = _special_panel(
        axes[0],
        "Ferrocene",
        "ring centers / Fe distances",
        ROOT / "calculations/special_coordinate_example/ferrocene.xyzin",
        az=-32.0,
        el=22.0,
    )
    top = xy[[0, 2, 3, 4, 5]].mean(axis=0)
    bottom = xy[[6, 7, 8, 9, 10]].mean(axis=0)
    fe = xy[1]
    for c in (top, bottom):
        center_marker(axes[0], c)
        axes[0].plot([c[0], fe[0]], [c[1], fe[1]], ls="--", lw=1.8, color=TEAL, zorder=4)

    # Formic-acid--water: fragment centers and weak contacts.
    symbols, coords, xy = _special_panel(
        axes[1],
        "Formic acid--water",
        "fragment pose / pseudo-contacts",
        ROOT / "standalone/examples/formic_acid_water.smith.xyz",
        az=0.0,
        el=0.0,
    )
    c1, c2 = xy[:5].mean(axis=0), xy[5:].mean(axis=0)
    for c in (c1, c2):
        center_marker(axes[1], c, BLUE)
    axes[1].plot([c1[0], c2[0]], [c1[1], c2[1]], ls=":", lw=1.8, color=BLUE, zorder=4)
    for i, j in ((4, 5), (2, 7)):
        axes[1].plot([xy[i, 0], xy[j, 0]], [xy[i, 1], xy[j, 1]], ls="--", lw=1.6, color=TEAL, zorder=4)

    # Water dimer.
    symbols, coords, xy = _special_panel(
        axes[2],
        "Water dimer",
        "six relative rigid-body rows",
        ROOT / "standalone/examples/water_dimer.smith.xyz",
        az=0.0,
        el=0.0,
    )
    c1, c2 = xy[:3].mean(axis=0), xy[3:].mean(axis=0)
    center_marker(axes[2], c1, BLUE)
    center_marker(axes[2], c2, BLUE)
    double_arrow(axes[2], c1, c2, BLUE, shrink=0.025)
    axes[2].plot([xy[1, 0], xy[3, 0]], [xy[1, 1], xy[3, 1]], ls="--", lw=1.7, color=TEAL, zorder=4)

    # Benzene--water: ring center and fragment-center distance.
    symbols, coords, xy = _special_panel(
        axes[3],
        "Benzene--water",
        "ring center / fragment pose",
        ROOT / "standalone/examples/benzene_water.smith.xyz",
        az=-28.0,
        el=34.0,
    )
    ring = xy[:6].mean(axis=0)
    water = xy[12:].mean(axis=0)
    center_marker(axes[3], ring)
    center_marker(axes[3], water, BLUE)
    axes[3].plot([ring[0], water[0]], [ring[1], water[1]], ls="--", lw=1.8, color=TEAL, zorder=4)

    # eta3 ligand center.
    symbols, coords, xy = _special_panel(
        axes[4],
        r"$\eta^3$ allyl--PdCl",
        "supplied ligand center",
        ROOT / "standalone/examples/eta3_allyl_palladium.source.xyz",
        az=0.0,
        el=90.0,
    )
    ligand = xy[[2, 3, 4]].mean(axis=0)
    center_marker(axes[4], ligand)
    axes[4].plot([xy[0, 0], ligand[0]], [xy[0, 1], ligand[1]], ls="--", lw=2.0, color=TEAL, zorder=4)

    fig.suptitle(
        "Molecular structures carrying protected special coordinates",
        x=0.5,
        y=0.995,
        fontsize=12.0,
        fontweight="bold",
        color=NAVY,
    )
    # The sixth cell is a visual key, kept at the same scale as the structures.
    key = axes[5]
    key.axis("off")
    key.text(0.08, 0.82, "Protected objects", transform=key.transAxes, color=NAVY, fontsize=10.0, fontweight="bold")
    key.scatter([0.16], [0.60], transform=key.transAxes, marker="D", s=45, c=TEAL, edgecolors="white", linewidths=0.7)
    key.text(0.28, 0.60, "stored center", transform=key.transAxes, va="center", color=SOFT_GREY, fontsize=8.3)
    key.plot([0.10, 0.23], [0.41, 0.41], transform=key.transAxes, ls="--", lw=2.0, color=TEAL)
    key.text(0.28, 0.41, "center/contact row", transform=key.transAxes, va="center", color=SOFT_GREY, fontsize=8.3)
    key.add_patch(
        FancyArrowPatch(
            (0.10, 0.22),
            (0.23, 0.22),
            transform=key.transAxes,
            arrowstyle="<->",
            mutation_scale=11,
            linewidth=2.0,
            color=BLUE,
        )
    )
    key.text(0.28, 0.22, "relative fragment pose", transform=key.transAxes, va="center", color=SOFT_GREY, fontsize=8.3)

    fig.tight_layout(rect=(0.01, 0.015, 0.99, 0.91), h_pad=0.65, w_pad=0.35)
    fig.savefig(FIGURES / "special_coordinate_structures.png", bbox_inches="tight", pad_inches=0.05, facecolor=fig.get_facecolor())


def main() -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    build_atlas()
    build_special_gallery()


if __name__ == "__main__":
    main()
