from __future__ import annotations

import argparse
from importlib.resources import as_file, files
from pathlib import Path

from matrix_core import read_sectioned_lines, replace_section, section_content
from matrix_neo.definition import write_sonic_build_sections_from_cartesian
from matrix_neo.standalone import write_smith_build_sections_from_input

from . import __version__


REQUIRED_ORACLE_SECTIONS = ("VALIDATION", "TOPOLOGY", "SYNTHONS", "SYMMETRY")
PROVENANCE_SCHEMA = "matrix.smith.standalone.v1"
MATRIX_REVISION = "f943523dd5468d35c7ebdc5bfa9f7bb305afda7f"


def _has_oracle_state(path: Path) -> bool:
    lines = read_sectioned_lines(path)
    return all(section_content(lines, name) for name in REQUIRED_ORACLE_SECTIONS)


def _build(args: argparse.Namespace) -> int:
    source = Path(args.input)
    target = Path(args.output) if args.output is not None else source.with_suffix(".xyzin")
    has_oracle_state = _has_oracle_state(source)
    if args.require_oracle_state and not has_oracle_state:
        missing = [
            name
            for name in REQUIRED_ORACLE_SECTIONS
            if not section_content(read_sectioned_lines(source), name)
        ]
        raise SystemExit(
            "SMITH was asked to require a frozen ORACLE state; missing sections: "
            + ", ".join(missing)
        )

    if has_oracle_state:
        definition = write_sonic_build_sections_from_cartesian(source, target)
        profile = "ORACLE_STATE"
    else:
        definition = write_smith_build_sections_from_input(source, target)
        profile = "REDUCED_ORACLE"

    replace_section(
        target,
        "SMITH_PROVENANCE",
        [
            f"SCHEMA {PROVENANCE_SCHEMA}",
            f"SMITH_VERSION {__version__}",
            f"MATRIX_REVISION {MATRIX_REVISION}",
            f"PERCEPTION_PROFILE {profile}",
            "ORACLE_RELATION CONTINUOUS_PERCEPTION_DEVELOPED_FROM_PROXIMA",
        ],
    )
    print(
        f"Wrote {target} (profile={profile}, GICs={len(definition.gics)}, "
        f"rank={definition.rank})"
    )
    if profile == "REDUCED_ORACLE":
        print(
            "Note: the Cartesian input used the reduced bundled ORACLE perception profile. "
            "Use --require-oracle-state for production runs that must consume a separately "
            "validated ORACLE state."
        )
    return 0


def _inspect(args: argparse.Namespace) -> int:
    path = Path(args.xyzin)
    lines = read_sectioned_lines(path)
    provenance = section_content(lines, "SMITH_PROVENANCE")
    gics = section_content(lines, "GIC")
    print(f"path: {path}")
    print(f"oracle_state: {'complete' if _has_oracle_state(path) else 'incomplete'}")
    print(f"gic_section: {'present' if gics else 'missing'}")
    if provenance:
        print("provenance:")
        for line in provenance:
            print(f"  {line}")
    return 0


def _example(args: argparse.Namespace) -> int:
    resource = files("smith_sonic.examples").joinpath(f"{args.name}.smith.xyz")
    target = Path(args.output) if args.output is not None else Path(f"{args.name}.xyzin")
    with as_file(resource) as source:
        build_args = argparse.Namespace(
            input=source,
            output=target,
            require_oracle_state=False,
        )
        return _build(build_args)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="smith-sonic",
        description="Build a frozen SONIC coordinate contract with standalone SMITH.",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    subparsers = parser.add_subparsers(dest="command", required=True)

    build = subparsers.add_parser("build", help="Build a SONIC contract")
    build.add_argument("input", type=Path, help="SMITH extended XYZ or ORACLE-enriched xyzin")
    build.add_argument("output", type=Path, nargs="?", help="Output xyzin path")
    build.add_argument(
        "--require-oracle-state",
        action="store_true",
        help="Refuse plain XYZ and require frozen ORACLE validation/topology/synthon/symmetry sections",
    )
    build.set_defaults(func=_build)

    inspect = subparsers.add_parser("inspect", help="Inspect state and standalone provenance")
    inspect.add_argument("xyzin", type=Path)
    inspect.set_defaults(func=_inspect)

    example = subparsers.add_parser("example", help="Run an example shipped in the package")
    example.add_argument("name", choices=("water", "norbornane"))
    example.add_argument("output", type=Path, nargs="?", help="Output xyzin path")
    example.set_defaults(func=_example)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    return int(args.func(args))


if __name__ == "__main__":
    raise SystemExit(main())
