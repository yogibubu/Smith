#!/usr/bin/env python3
"""Build the graphical table-of-contents image for the SONIC manuscript."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "figures" / "graphical_toc.png"


def box(ax, xy, width, height, text, face, edge="#202020", size=12, weight="normal"):
    patch = FancyBboxPatch(
        xy,
        width,
        height,
        boxstyle="round,pad=0.028,rounding_size=0.035",
        linewidth=1.15,
        edgecolor=edge,
        facecolor=face,
    )
    ax.add_patch(patch)
    ax.text(
        xy[0] + width / 2,
        xy[1] + height / 2,
        text,
        ha="center",
        va="center",
        fontsize=size,
        fontweight=weight,
        color="#141414",
        linespacing=1.12,
    )


def arrow(ax, start, end, rad=0.0):
    ax.add_patch(
        FancyArrowPatch(
            start,
            end,
            arrowstyle="-|>",
            mutation_scale=16,
            linewidth=1.35,
            color="#333333",
            connectionstyle=f"arc3,rad={rad}",
            shrinkA=7,
            shrinkB=7,
        )
    )


def main() -> None:
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(8.6, 3.15), dpi=320)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    box(
        ax,
        (0.035, 0.50),
        0.245,
        0.30,
        "Molecular state\ngeometry | topology\nsymmetry | fragments",
        "#e9f2ff",
        size=8.8,
    )
    box(
        ax,
        (0.365, 0.47),
        0.27,
        0.36,
        "SONIC contract\nnon-redundant internals\nsparse B rows + provenance",
        "#f1eefb",
        size=9.0,
        weight="bold",
    )
    box(
        ax,
        (0.720, 0.50),
        0.245,
        0.30,
        "Portable use\nG16 GICs | GF/PED\nDVR | PES | refinement",
        "#edf7ee",
        size=8.8,
    )

    arrow(ax, (0.280, 0.650), (0.365, 0.650))
    arrow(ax, (0.635, 0.650), (0.720, 0.650))

    box(
        ax,
        (0.080, 0.12),
        0.24,
        0.22,
        "protected rows\nrings | fragments\npseudo-bonds",
        "#fff4df",
        edge="#6a5428",
        size=8.4,
    )
    box(
        ax,
        (0.380, 0.11),
        0.24,
        0.24,
        "block-local\nrank reduction\nsymmetry adaptation",
        "#f7e8ee",
        edge="#6d3f55",
        size=8.4,
    )
    box(
        ax,
        (0.680, 0.12),
        0.24,
        0.22,
        "validated export\nGaussian-readable GICs\nand Cartesian paths",
        "#e9f7f7",
        edge="#356667",
        size=8.4,
    )

    arrow(ax, (0.455, 0.47), (0.270, 0.34), rad=0.12)
    arrow(ax, (0.500, 0.47), (0.500, 0.35))
    arrow(ax, (0.545, 0.47), (0.730, 0.34), rad=-0.12)

    ax.text(
        0.5,
        0.94,
        "Automatic construction and validation of symmetry-oriented non-redundant internal coordinates",
        ha="center",
        va="center",
        fontsize=10.5,
        fontweight="bold",
        color="#111111",
    )
    fig.tight_layout(pad=0.05)
    fig.savefig(OUT, bbox_inches="tight", pad_inches=0.03)


if __name__ == "__main__":
    main()
