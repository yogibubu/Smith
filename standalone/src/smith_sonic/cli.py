from __future__ import annotations

import argparse
from importlib.resources import as_file, files
from pathlib import Path
import tempfile

from matrix_chem import preprocess_to_enriched_xyz, write_validation_section
from matrix_core import read_sectioned_lines, replace_section, section_content
from matrix_fragments import (
    build_fragment_definition_from_xyzin,
    write_fragment_build_section,
    write_interaction_center_section,
)
from matrix_neo import write_gic_report, write_gicforge_gaussian_input
from matrix_neo.definition import write_sonic_build_sections_from_cartesian
from matrix_neo.standalone import (
    _normalized_source_kind,
    _optional_bool,
    _optional_string,
    _pairs,
    _read_smith_input,
    _strings,
)

from . import __version__


REQUIRED_ORACLE_SECTIONS = (
    "VALIDATION",
    "TOPOLOGY",
    "SYNTHONS",
    "SYMMETRY",
    "PRIMITIVES",
)
PROVENANCE_SCHEMA = "matrix.smith.standalone.v1"
MATRIX_REVISION = "cf5fdcebb85a5035d5c0400d7cc2398a0580df66"
DEFAULT_G16_ROUTE = "#p hf/sto-3g opt=(readallgic,calcfc,maxcycle=80)"
EXAMPLES = {
    "water": ("water.smith.xyz", False),
    "norbornane": ("norbornane.smith.xyz", False),
    "formic-acid-water": ("formic_acid_water.smith.xyz", False),
    "eta3-allyl-palladium": ("eta3_allyl_palladium.oracle.xyzin", True),
}


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
        definition = write_sonic_build_sections_from_cartesian(
            source,
            target,
            improper_dihedrals=True,
        )
        profile = "ORACLE_STATE"
        fragment_profile = (
            "ORACLE_SUPPLIED"
            if section_content(read_sectioned_lines(source), "FRAGMENTS")
            else "NONE"
        )
    else:
        definition, fragment_count = _build_from_reduced_input(source, target)
        profile = "REDUCED_ORACLE"
        fragment_profile = "AUTO_CONNECTED_COMPONENTS" if fragment_count > 1 else "NONE"

    replace_section(
        target,
        "SMITH_PROVENANCE",
        [
            f"SCHEMA {PROVENANCE_SCHEMA}",
            f"SMITH_VERSION {__version__}",
            f"MATRIX_REVISION {MATRIX_REVISION}",
            f"PERCEPTION_PROFILE {profile}",
            f"FRAGMENT_PROFILE {fragment_profile}",
            "ORACLE_RELATION CONTINUOUS_PERCEPTION_DEVELOPED_FROM_PROXIMA",
        ],
    )
    report_path, gaussian_path = _write_sidecars(target)
    print(
        f"Wrote {target} (profile={profile}, GICs={len(definition.gics)}, "
        f"rank={definition.rank})"
    )
    print(f"Wrote coordinate report {report_path}")
    print(f"Wrote Gaussian 16 input {gaussian_path}")
    if profile == "REDUCED_ORACLE":
        print(
            "Note: the Cartesian input used the reduced bundled ORACLE perception profile. "
            "Use --require-oracle-state for production runs that must consume a separately "
            "validated ORACLE state."
        )
    return 0


def _build_from_reduced_input(source: Path, target: Path):
    xyz_lines, options = _read_smith_input(source)
    target.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="smith-sonic-") as scratch:
        geometry_source = Path(scratch) / "geometry.xyz"
        geometry_source.write_text("\n".join(xyz_lines) + "\n", encoding="utf-8")
        preprocess_to_enriched_xyz(
            geometry_source,
            target,
            source_kind=_normalized_source_kind(options.get("source_kind", "auto")),
        )
    write_validation_section(target)

    fragments = build_fragment_definition_from_xyzin(target)
    if len(fragments.fragments) > 1:
        write_fragment_build_section(target)
        write_interaction_center_section(target)

    g16_compatibility = bool(options.get("g16", True))
    improper_dihedrals = _optional_bool(options.get("improper_dihedrals"))
    if g16_compatibility:
        improper_dihedrals = True
    fragment_mode = _optional_string(options.get("fragment_mode"))
    if len(fragments.fragments) > 1 and fragment_mode is None:
        fragment_mode = "special-coordinates"
    definition = write_sonic_build_sections_from_cartesian(
        target,
        target,
        symmetrize=bool(options.get("symmetrize", False)),
        sycart=bool(options.get("sycart", False)),
        symmetry_group=_optional_string(options.get("symmetry_group")),
        improper_dihedrals=improper_dihedrals,
        fragment_mode=fragment_mode,
        xh_stretch_policy=_optional_string(options.get("xh_stretch_policy")),
        local_xh_bonds=_pairs(options.get("local_xh_bonds")),
        local_xh_classes=_strings(options.get("local_xh_classes")),
    )
    return definition, len(fragments.fragments)


def _write_sidecars(target: Path) -> tuple[Path, Path]:
    report_path = target.with_suffix(".smith.out")
    gaussian_path = target.with_suffix(".g16.gjf")
    write_gic_report(target, report_path)
    write_gicforge_gaussian_input(
        target,
        gaussian_path,
        route=DEFAULT_G16_ROUTE,
        title=f"{target.stem} SMITH/SONIC Gaussian 16 optimization",
        g16_compatibility=True,
    )
    return report_path, gaussian_path


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
    filename, require_oracle_state = EXAMPLES[args.name]
    resource = files("smith_sonic.examples").joinpath(filename)
    target = Path(args.output) if args.output is not None else Path(f"{args.name}.xyzin")
    with as_file(resource) as source:
        build_args = argparse.Namespace(
            input=source,
            output=target,
            require_oracle_state=require_oracle_state,
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
    example.add_argument("name", choices=tuple(EXAMPLES))
    example.add_argument("output", type=Path, nargs="?", help="Output xyzin path")
    example.set_defaults(func=_example)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    return int(args.func(args))


if __name__ == "__main__":
    raise SystemExit(main())
