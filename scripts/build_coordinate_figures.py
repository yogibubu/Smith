#!/usr/bin/env python3
"""Build the introductory SONIC atlas and the special-coordinate gallery."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Arc, FancyArrowPatch, FancyBboxPatch
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
    """Show the complete, inspectable record behind one representative SONIC."""
    fig = plt.figure(figsize=(10.0, 3.35), dpi=320, facecolor=WARM_WHITE)
    canvas = fig.add_axes((0, 0, 1, 1))
    canvas.set_xlim(0, 1)
    canvas.set_ylim(0, 1)
    canvas.axis("off")
    canvas.text(
        0.035,
        0.93,
        "Anatomy of an inspectable SONIC coordinate",
        ha="left",
        va="center",
        color=NAVY,
        fontsize=14.0,
        fontweight="bold",
    )
    canvas.text(
        0.036,
        0.865,
        "Every accepted row retains its chemical sources, coefficients, metadata, analytic derivative and Cartesian fingerprint.",
        ha="left",
        va="center",
        color=SOFT_GREY,
        fontsize=8.2,
    )

    lefts = (0.035, 0.278, 0.522, 0.765)
    widths = (0.205, 0.205, 0.205, 0.200)
    headings = (
        ("1", "TYPED SOURCES", TEAL),
        ("2", "FROZEN ROW", AMBER),
        ("3", "CONTRACT RECORD", BLUE),
        ("4", "CARTESIAN FINGERPRINT", VIOLET),
    )
    for left, width, (number, heading, accent) in zip(lefts, widths, headings, strict=True):
        canvas.add_patch(
            FancyBboxPatch(
                (left, 0.22),
                width,
                0.57,
                boxstyle="round,pad=0.009,rounding_size=0.018",
                facecolor="white",
                edgecolor="#d7dde1",
                linewidth=0.85,
            )
        )
        canvas.text(left + 0.018, 0.745, number, color="white", fontsize=7.4, fontweight="bold",
                    bbox={"boxstyle": "circle,pad=0.28", "facecolor": accent, "edgecolor": accent})
        canvas.text(left + 0.050, 0.746, heading, color=accent, fontsize=7.6, fontweight="bold", va="center")

    # The same cyclobutane fixture is used at both ends of the record.
    path = ROOT / "calculations/quick_qm/cyclobutane_hfsto3g.xyzin"
    symbols, coords = read_xyz_like(path)
    bonds = read_bonds(path) or infer_bonds(symbols, coords)
    source_ax = fig.add_axes((0.055, 0.315, 0.165, 0.34))
    xy = molecule(source_ax, symbols, coords, bonds, azimuth=-18.0, elevation=25.0)
    carbon = [i for i, atom in enumerate(symbols) if atom == "C"][:4]
    for order, idx in enumerate(carbon):
        dy = 0.14 if order % 2 == 0 else -0.14
        source_ax.annotate(
            "",
            xy=xy[idx] + np.array([0.0, dy]),
            xytext=xy[idx],
            arrowprops={"arrowstyle": "-|>", "color": TEAL, "lw": 1.45},
            zorder=6,
        )
    canvas.text(0.137, 0.285, r"$U_1\quad U_2\quad U_3\quad U_4$", ha="center", color=NAVY, fontsize=8.5)

    coeff_ax = fig.add_axes((0.305, 0.405, 0.150, 0.19))
    coeff_ax.bar(range(4), [0.5, -0.5, 0.5, -0.5], color=[AMBER, "#e8c785", AMBER, "#e8c785"], width=0.58)
    coeff_ax.axhline(0.0, color="#83909a", lw=0.65)
    coeff_ax.set_xticks(range(4), [r"$U_1$", r"$U_2$", r"$U_3$", r"$U_4$"], fontsize=6.7)
    coeff_ax.set_yticks([-0.5, 0.0, 0.5], [r"$-1/2$", "0", r"$1/2$"], fontsize=6.2)
    coeff_ax.spines[["top", "right", "left"]].set_visible(False)
    coeff_ax.tick_params(axis="y", length=0)
    canvas.text(0.380, 0.645, "RPck001", ha="center", color=NAVY, fontsize=11.0, fontweight="bold")
    canvas.text(0.380, 0.325, r"$q=\frac{1}{2}(U_1-U_2+U_3-U_4)$", ha="center", color=NAVY, fontsize=8.5)

    records = (
        ("family", "ring puckering"),
        ("symmetry", r"$\Gamma_k$ + phase"),
        ("unit", "radian"),
        ("protected", "true"),
        ("Wilson B", "analytic"),
        ("provenance", "frozen"),
    )
    y = 0.645
    for key, value in records:
        canvas.text(0.545, y, key, color=SOFT_GREY, fontsize=6.8, va="center")
        canvas.text(0.705, y, value, color=NAVY, fontsize=7.2, va="center", ha="right", fontweight="bold")
        canvas.plot([0.543, 0.707], [y - 0.027, y - 0.027], color="#edf0f2", lw=0.65)
        y -= 0.065

    motion_ax = fig.add_axes((0.785, 0.355, 0.160, 0.30))
    xy_motion = molecule(motion_ax, symbols, coords, bonds, azimuth=-18.0, elevation=25.0)
    for order, idx in enumerate(carbon):
        dy = 0.16 if order % 2 == 0 else -0.16
        motion_ax.annotate(
            "",
            xy=xy_motion[idx] + np.array([0.0, dy]),
            xytext=xy_motion[idx],
            arrowprops={"arrowstyle": "-|>", "color": VIOLET, "lw": 1.55},
            zorder=6,
        )
    canvas.text(0.865, 0.315, r"stored $-h\;/\;0\;/\;+h$ trajectory", ha="center", color=NAVY, fontsize=7.4)

    for x0, x1 in ((0.244, 0.272), (0.487, 0.516), (0.731, 0.759)):
        canvas.add_patch(FancyArrowPatch((x0, 0.505), (x1, 0.505), arrowstyle="-|>", mutation_scale=9,
                                         linewidth=1.1, color="#9aa6ae"))

    canvas.add_patch(
        FancyBboxPatch((0.035, 0.075), 0.930, 0.090, boxstyle="round,pad=0.006,rounding_size=0.012",
                       facecolor="#f0f6f7", edgecolor="#d6e7e7", linewidth=0.75)
    )
    audit_items = (
        "human-readable primitive expansion",
        "ordered coefficient row",
        "rank and symmetry diagnostics",
        "reproducible Cartesian motion",
    )
    for index, item in enumerate(audit_items):
        x = 0.060 + index * 0.231
        canvas.text(x, 0.120, "✓", color=TEAL, fontsize=9.0, fontweight="bold", va="center")
        canvas.text(x + 0.021, 0.120, item, color=NAVY, fontsize=6.9, va="center")

    for suffix in ("png", "pdf"):
        fig.savefig(FIGURES / f"sonic_coordinate_atlas.{suffix}", dpi=320, bbox_inches="tight",
                    pad_inches=0.04, facecolor=fig.get_facecolor())
    plt.close(fig)


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
