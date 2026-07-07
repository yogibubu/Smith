#!/usr/bin/env python3
"""Measure SMITH/SONIC construction and B-matrix scaling for the manuscript."""

from __future__ import annotations

import json
import math
import statistics
import sys
import time
import tracemalloc
from dataclasses import dataclass
from pathlib import Path
from tempfile import TemporaryDirectory

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
MATRIX_ROOT = Path("/Users/vincenzobarone/Documents/git/software/matrix")
DATA = ROOT / "data" / "construction_scaling_benchmark.json"
FIGURE = ROOT / "figures" / "construction_scaling_benchmark.png"
SOURCE_ROOT = MATRIX_ROOT / "tests" / "fixtures" / "test_molecules" / "molecules"

SYSTEMS = (
    ("PAH/fused rings", "benzene", "benzene.inp", True, "none"),
    ("PAH/fused rings", "naphthalene", "naphtalene.inp", True, "none"),
    ("PAH/fused rings", "anthracene", "anthracene.inp", True, "none"),
    ("PAH/fused rings", "phenanthrene", "phenantrene.inp", True, "none"),
    ("PAH/fused rings", "pyrene", "pyrene.inp", True, "none"),
    ("PAH/fused rings", "coronene", "coronene.inp", True, "none"),
    ("ring/cage topology", "cubane", "cubane.inp", True, "none"),
    ("ring/cage topology", "norbornane", "norbornane.inp", True, "none"),
    ("ring/cage topology", "cyclooctane", "cyclottane.inp", True, "none"),
    ("ring/cage topology", "camphor", "camphor.inp", False, "none"),
    ("ring/cage topology", "spiro", "spiro.inp", True, "none"),
    ("special centers", "ferrocene_d5h", "ferrocene.inp", True, "special-coordinates"),
    ("special centers", "ferrocene_d5d", "ferrocene_staggered.inp", True, "special-coordinates"),
)
REPEATS = 5


def add_matrix_packages_to_path() -> None:
    package_root = MATRIX_ROOT / "packages"
    for path in reversed(
        [str(item / "src") for item in sorted(package_root.iterdir()) if (item / "src").is_dir()]
    ):
        if path not in sys.path:
            sys.path.insert(0, path)


@dataclass(frozen=True)
class Trial:
    build_seconds: float
    b_seconds: float
    peak_mib: float
    atoms: int
    candidates: int
    primitives: int
    gics: int
    rank: int
    target_rank: int
    b_nnz: int
    b_density: float
    point_group: str


def run_trial(name: str, source_name: str, symmetrize: bool, fragment_mode: str) -> Trial:
    from matrix_chem import preprocess_to_enriched_xyz, read_enriched_xyz, write_validation_section
    from matrix_fragments import write_fragment_build_section, write_interaction_center_section
    from matrix_neo import build_gic_b_matrix, write_gicforge_build_sections

    with TemporaryDirectory(prefix=f"smith-bench-{name}-") as tmp:
        tmp_path = Path(tmp)
        xyzin = tmp_path / f"{name}.xyzin"
        preprocess_to_enriched_xyz(SOURCE_ROOT / source_name, xyzin)
        write_validation_section(xyzin)
        if fragment_mode == "special-coordinates":
            write_fragment_build_section(xyzin)
            write_interaction_center_section(xyzin)
        geometry = read_enriched_xyz(xyzin)
        tracemalloc.start()
        start = time.perf_counter()
        definition = write_gicforge_build_sections(
            xyzin,
            symmetrize=symmetrize,
            fragment_mode=fragment_mode,
        )
        build_seconds = time.perf_counter() - start
        start = time.perf_counter()
        b_matrix = build_gic_b_matrix(definition, coordinates_angstrom=geometry.coordinates_angstrom)
        b_seconds = time.perf_counter() - start
        _current, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        sparse = b_matrix.sparse_matrix()
        return Trial(
            build_seconds=build_seconds,
            b_seconds=b_seconds,
            peak_mib=peak / (1024.0 * 1024.0),
            atoms=len(geometry.atoms),
            candidates=definition.candidate_count,
            primitives=len(definition.primitives),
            gics=len(definition.gics),
            rank=definition.rank,
            target_rank=definition.target_rank,
            b_nnz=sparse.nnz,
            b_density=sparse.density,
            point_group=definition.point_group,
        )


def median(values: list[float]) -> float:
    return float(statistics.median(values))


def _fit_power_law(results: list[dict], x_key: str, y_key: str) -> dict[str, float | int | str]:
    points = [
        (float(item[x_key]), float(item[y_key]))
        for item in results
        if float(item[x_key]) > 0.0 and float(item[y_key]) > 0.0
    ]
    if len(points) < 3:
        return {"x": x_key, "y": y_key, "n": len(points), "slope": math.nan, "r2": math.nan}
    x = np.log10(np.asarray([point[0] for point in points], dtype=float))
    y = np.log10(np.asarray([point[1] for point in points], dtype=float))
    slope, intercept = np.polyfit(x, y, 1)
    fitted = slope * x + intercept
    ss_res = float(np.sum((y - fitted) ** 2))
    ss_tot = float(np.sum((y - float(np.mean(y))) ** 2))
    r2 = 1.0 if ss_tot == 0.0 else 1.0 - ss_res / ss_tot
    return {
        "x": x_key,
        "y": y_key,
        "n": len(points),
        "slope": float(slope),
        "intercept": float(intercept),
        "r2": float(r2),
    }


def _scaling_fits(results: list[dict]) -> dict[str, dict]:
    pah = [item for item in results if item["series"] == "PAH/fused rings"]
    regular = [item for item in results if item["series"] != "special centers"]
    return {
        "pah_build_vs_atoms": _fit_power_law(pah, "atoms", "build_time_median_s"),
        "pah_b_vs_b_nnz": _fit_power_law(pah, "b_nnz", "b_time_median_s"),
        "regular_build_vs_candidates": _fit_power_law(
            regular, "candidates", "build_time_median_s"
        ),
        "regular_memory_vs_primitives": _fit_power_law(
            regular, "primitive_rows_stored", "peak_memory_median_mib"
        ),
    }


def _plot_results(results: list[dict], fits: dict[str, dict]) -> None:
    FIGURE.parent.mkdir(parents=True, exist_ok=True)
    colors = {
        "PAH/fused rings": "#1f77b4",
        "ring/cage topology": "#2ca02c",
        "special centers": "#d62728",
    }
    markers = {
        "PAH/fused rings": "o",
        "ring/cage topology": "s",
        "special centers": "^",
    }
    plt.rcParams.update(
        {
            "font.size": 10,
            "axes.titlesize": 11,
            "axes.labelsize": 10,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "legend.fontsize": 9,
        }
    )
    fig, axes = plt.subplots(1, 3, figsize=(14.2, 4.25))
    panels = (
        ("candidates", "build_time_median_s", "Build time / s"),
        ("b_nnz", "b_time_median_s", r"$\mathbf{B}$ time / s"),
        ("primitive_rows_stored", "peak_memory_median_mib", "Peak memory / MiB"),
    )
    for axis, (x_key, y_key, ylabel) in zip(axes, panels):
        for series in colors:
            subset = [item for item in results if item["series"] == series]
            axis.scatter(
                [item[x_key] for item in subset],
                [item[y_key] for item in subset],
                label=series,
                color=colors[series],
                marker=markers[series],
                s=48,
                edgecolor="white",
                linewidth=0.5,
                zorder=3,
            )
        axis.set_xscale("log")
        axis.set_yscale("log")
        axis.grid(True, which="both", color="#e6e6e6", linewidth=0.6)
        axis.set_xlabel(
            {
                "candidates": "Primitive candidates",
                "b_nnz": r"Sparse $\mathbf{B}$ nonzeros",
                "primitive_rows_stored": "Stored primitive rows",
            }[x_key]
        )
        axis.set_ylabel(ylabel)
    axes[0].set_title(
        f"Build slope {fits['regular_build_vs_candidates']['slope']:.2f}"
    )
    axes[1].set_title(f"B slope {fits['pah_b_vs_b_nnz']['slope']:.2f}")
    axes[2].set_title(
        f"Memory slope {fits['regular_memory_vs_primitives']['slope']:.2f}"
    )
    axes[0].legend(frameon=False, loc="upper left")
    fig.tight_layout(w_pad=2.0)
    fig.savefig(FIGURE, dpi=300)


def main() -> None:
    add_matrix_packages_to_path()
    results = []
    for series, name, source_name, symmetrize, fragment_mode in SYSTEMS:
        trials = [run_trial(name, source_name, symmetrize, fragment_mode) for _ in range(REPEATS)]
        first = trials[0]
        results.append(
            {
                "series": series,
                "system": name,
                "source": source_name,
                "point_group": first.point_group,
                "symmetrize": symmetrize,
                "fragment_mode": fragment_mode,
                "atoms": first.atoms,
                "candidates": first.candidates,
                "primitive_rows_stored": first.primitives,
                "gics": first.gics,
                "rank": first.rank,
                "target_rank": first.target_rank,
                "b_nnz": first.b_nnz,
                "b_density": first.b_density,
                "build_time_median_s": median([trial.build_seconds for trial in trials]),
                "b_time_median_s": median([trial.b_seconds for trial in trials]),
                "peak_memory_median_mib": median([trial.peak_mib for trial in trials]),
                "build_time_trials_s": [trial.build_seconds for trial in trials],
                "b_time_trials_s": [trial.b_seconds for trial in trials],
                "peak_memory_trials_mib": [trial.peak_mib for trial in trials],
            }
        )
    fits = _scaling_fits(results)
    _plot_results(results, fits)
    output = {
        "method": "SMITH/SONIC construction scaling study",
        "repeats": REPEATS,
        "python": sys.version.split()[0],
        "notes": [
            "Median wall times from single-process Python runs on oracle.sns.it (2 x AMD EPYC 7543, 128 logical CPUs, 1.0 TiB RAM); values are implementation diagnostics, not hardware-independent constants.",
            "Build time includes coordinate construction, rank reduction, optional symmetry projection and writing the xyzin GIC sections.",
            "B time is one analytic Wilson B evaluation from an already frozen contract.",
            "Peak memory is tracemalloc peak during build plus B evaluation for a single trial.",
            "Power-law slopes are log-log least-squares fits and are meant to summarize the sampled implementation regime, not asymptotic proofs.",
        ],
        "figure": str(FIGURE.relative_to(ROOT)),
        "fits": fits,
        "systems": results,
    }
    DATA.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
