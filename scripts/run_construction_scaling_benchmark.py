#!/usr/bin/env python3
"""Measure SMITH/SONIC construction and B-matrix costs for manuscript Table."""

from __future__ import annotations

import json
import statistics
import sys
import time
import tracemalloc
from dataclasses import dataclass
from pathlib import Path
from tempfile import TemporaryDirectory

ROOT = Path(__file__).resolve().parents[1]
MATRIX_ROOT = Path("/Users/vincenzobarone/Documents/git/software/matrix")
DATA = ROOT / "data" / "construction_scaling_benchmark.json"
SOURCE_ROOT = MATRIX_ROOT / "tests" / "fixtures" / "test_molecules" / "molecules"

SYSTEMS = (
    ("benzene", "benzene.inp", True, "none"),
    ("norbornane", "norbornane.inp", True, "none"),
    ("camphor", "camphor.inp", False, "none"),
    ("ferrocene", "ferrocene.inp", True, "special-coordinates"),
    ("spiro", "spiro.inp", True, "none"),
    ("cyclooctane", "cyclottane.inp", True, "none"),
    ("coronene", "coronene.inp", True, "none"),
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


def main() -> None:
    add_matrix_packages_to_path()
    results = []
    for name, source_name, symmetrize, fragment_mode in SYSTEMS:
        trials = [run_trial(name, source_name, symmetrize, fragment_mode) for _ in range(REPEATS)]
        first = trials[0]
        results.append(
            {
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
    output = {
        "method": "SMITH/SONIC construction scaling microbenchmark",
        "repeats": REPEATS,
        "python": sys.version.split()[0],
        "notes": [
            "Median wall times from one local macOS run; values are implementation diagnostics, not hardware-independent constants.",
            "Build time includes coordinate construction, rank reduction, optional symmetry projection and writing the xyzin GIC sections.",
            "B time is one analytic Wilson B evaluation from an already frozen contract.",
            "Peak memory is tracemalloc peak during build plus B evaluation for a single trial.",
        ],
        "systems": results,
    }
    DATA.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
