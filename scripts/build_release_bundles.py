#!/usr/bin/env python3
"""Build deterministic manuscript, reviewer, and standalone RC bundles."""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path
import re
import zipfile


ROOT = Path(__file__).resolve().parents[1]
VERSION = "0.1.0rc7"
ZIP_TIME = (2026, 7, 18, 12, 0, 0)


def _write_zip(target: Path, members: list[tuple[Path, str]], extras: dict[str, bytes] | None = None) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    payloads: list[tuple[str, bytes]] = []
    for source, archive_name in members:
        if not source.is_file():
            raise FileNotFoundError(source)
        payloads.append((archive_name, source.read_bytes()))
    payloads.extend((name, data) for name, data in (extras or {}).items())

    checksums = "".join(
        f"{hashlib.sha256(data).hexdigest()}  {name}\n"
        for name, data in sorted(payloads)
    ).encode("ascii")
    payloads.append(("SHA256SUMS", checksums))

    with zipfile.ZipFile(target, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as archive:
        for name, data in sorted(payloads):
            info = zipfile.ZipInfo(name, ZIP_TIME)
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = 0o100644 << 16
            archive.writestr(info, data)


def _figure_members() -> list[tuple[Path, str]]:
    text = (ROOT / "main.tex").read_text(encoding="utf-8")
    names = sorted(set(re.findall(r"figures/[A-Za-z0-9_./-]+\.(?:png|pdf|jpg)", text)))
    return [(ROOT / name, name) for name in names]


def build_arxiv_bundle() -> Path:
    style_root = ROOT / "arxiv" / "SMITH_arxiv_v1"
    members = [
        (ROOT / "main.tex", "main.tex"),
        (ROOT / "main.bbl", "main.bbl"),
        (ROOT / "acs-main.bib", "acs-main.bib"),
        (ROOT / "references.bib", "references.bib"),
        (style_root / "achemso.cls", "achemso.cls"),
        (style_root / "achemso-jctcce.cfg", "achemso-jctcce.cfg"),
        (style_root / "natmove.sty", "natmove.sty"),
        (
            ROOT / "output" / "pdf" / "SMITH_supporting_information.pdf",
            "anc/SMITH_supporting_information.pdf",
        ),
        *_figure_members(),
    ]
    target = ROOT / "arxiv" / "SMITH_manuscript_rc7.zip"
    _write_zip(target, members)
    return target


def build_reviewer_bundle() -> Path:
    members = [
        (ROOT / "output" / "pdf" / "SMITH_manuscript.pdf", "SMITH_manuscript.pdf"),
        (
            ROOT / "output" / "pdf" / "SMITH_supporting_information.pdf",
            "SMITH_supporting_information.pdf",
        ),
        (
            ROOT / "output" / "pdf" / "SMITH_Standalone_Manual.pdf",
            "SMITH_Standalone_Manual.pdf",
        ),
        (ROOT / "RELEASE_CANDIDATE.md", "RELEASE_CANDIDATE.md"),
    ]
    target = ROOT / "release" / "SMITH_review_bundle_rc7.zip"
    _write_zip(target, members)
    return target


def build_standalone_bundle(dist_dir: Path) -> Path:
    artifacts = sorted(dist_dir.glob(f"smith_sonic-{VERSION}*"))
    if not artifacts:
        raise FileNotFoundError(f"No smith_sonic-{VERSION} artifacts in {dist_dir}")
    members = [(path, f"dist/{path.name}") for path in artifacts]
    members.extend(
        [
            (ROOT / "standalone" / "README.md", "README.md"),
            (ROOT / "standalone" / "MANUAL.md", "MANUAL.md"),
            (
                ROOT / "output" / "pdf" / "SMITH_Standalone_Manual.pdf",
                "SMITH_Standalone_Manual.pdf",
            ),
            (ROOT / "RELEASE_CANDIDATE.md", "RELEASE_CANDIDATE.md"),
        ]
    )
    target = ROOT / "release" / "SMITH_standalone_rc7.zip"
    _write_zip(target, members)
    return target


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dist-dir", type=Path, required=True)
    args = parser.parse_args()
    targets = (
        build_arxiv_bundle(),
        build_reviewer_bundle(),
        build_standalone_bundle(args.dist_dir.resolve()),
    )
    for target in targets:
        print(target)


if __name__ == "__main__":
    main()
