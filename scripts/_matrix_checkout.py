"""Resolve an explicitly configured MATRIX source checkout."""

from __future__ import annotations

import os
from pathlib import Path
import sys


def matrix_root() -> Path:
    configured = os.environ.get("MATRIX_ROOT")
    if not configured:
        raise RuntimeError(
            "MATRIX_ROOT is not set. Define it as the path to your MATRIX checkout, "
            "for example: export MATRIX_ROOT=/path/to/MATRIX"
        )
    root = Path(configured).expanduser().resolve()
    if not (root / "packages").is_dir():
        raise RuntimeError(f"MATRIX_ROOT does not contain a packages directory: {root}")
    return root


def add_matrix_packages_to_path() -> Path:
    root = matrix_root()
    sources = [
        str(package / "src")
        for package in sorted((root / "packages").iterdir())
        if (package / "src").is_dir()
    ]
    for source in reversed(sources):
        if source not in sys.path:
            sys.path.insert(0, source)
    return root
