# Standalone SMITH / SONIC manual

## 1. Scope

Standalone SMITH converts a Cartesian molecular description and its perceived
molecular state into a rank-complete SONIC internal-coordinate contract.  It is
the small, reproducible interface associated with the SMITH manuscript, not a
distribution of the complete MATRIX framework.

ORACLE and SMITH have distinct roles.  ORACLE performs continuous molecular
perception, including topology and cycles, symmetry and atom equivalence,
effective atomic numbers, synthons, and interaction centres.  It develops the
ideas introduced in PROXIMA and is in final testing; it will be released and
described separately when ready.  SMITH consumes a frozen ORACLE state and
constructs coordinate candidates, protected special coordinates, rank
reduction, symmetry adaptation, analytic Wilson rows, and the serialized SONIC
contract.

## 2. Requirements and installation

Use Python 3.11 or newer and Git.  For the private test repository, first accept
the GitHub invitation and authenticate Git on the command line:

```bash
gh auth login
gh auth setup-git
git clone https://github.com/yogibubu/Smith.git
cd Smith
git switch agent/oracle-boundary-standalone-smith
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install ./standalone
```

On Windows PowerShell, activate the environment with
`.venv\Scripts\Activate.ps1`.  The branch-switch command will no longer be
needed after the draft branch is merged.

## 3. Command summary

```text
smith-sonic --version
smith-sonic example NAME [OUTPUT]
smith-sonic build INPUT [OUTPUT] [--require-oracle-state]
smith-sonic inspect FILE
```

`example` runs an input installed with the package.  `build` accepts either a
plain extended XYZ input or an ORACLE-enriched `xyzin`.  `inspect` reports the
presence of the frozen state, GIC section, and provenance profile.

For a production state prepared and validated by ORACLE, require the boundary
explicitly:

```bash
smith-sonic build molecule.oracle.xyzin molecule.sonic.xyzin --require-oracle-state
```

This option refuses an input missing any of `VALIDATION`, `TOPOLOGY`,
`SYNTHONS`, or `SYMMETRY`.

## 4. Input profiles

### Frozen ORACLE state

This is the production interface.  SMITH preserves the molecular perception
sections and adds the SONIC construction.  The output provenance contains
`PERCEPTION_PROFILE ORACLE_STATE`.

### Reduced ORACLE convenience profile

A plain or extended XYZ is sufficient for small reproducibility examples.  The
packaged path performs a reduced perception pass and records
`PERCEPTION_PROFILE REDUCED_ORACLE`.  For disconnected components it also
constructs fragment definitions and six intermolecular translation/orientation
coordinates.  This profile is deliberately limited and is not a substitute for
the forthcoming ORACLE application.

Extended-XYZ directives follow the Cartesian block, for example:

```text
#SMITH
fragment_mode = special-coordinates
symmetrize = false
sycart = false
improper_dihedrals = false
```

The available advanced controls include `symmetry_group`, `xh_stretch_policy`,
`local_xh_bonds`, and `local_xh_classes`.  Omit them to use the pinned SMITH
defaults.

## 5. Packaged examples

Run the complete set with:

```bash
smith-sonic example water water.xyzin
smith-sonic example norbornane norbornane.xyzin
smith-sonic example formic-acid-water formic-acid-water.xyzin
smith-sonic example eta3-allyl-palladium eta3-allyl-palladium.xyzin
```

Expected results are:

| Example | Input profile | Target rank | Purpose |
|---|---:|---:|---|
| water | reduced | 3 | minimal nonlinear molecule |
| norbornane | reduced | 51 | bridged and cyclic topology |
| formic-acid-water | reduced + fragments | 18 | non-covalent two-fragment contract |
| eta3-allyl-palladium | frozen ORACLE state | 24 | protected metal-to-η3-centre coordinate |

The formic-acid–water output must contain three `FRAG_TRANSLATION` and three
`FRAG_ORIENTATION` primitives.  The η3 probe must retain an `ETA3_CENTER` over
the three allyl carbon atoms and create a protected `CENTER_ATOM_DISTANCE` from
Pd to that centre.

The η3 geometry is an idealized interface test, not an optimized or computed
chemical benchmark.  Its interaction centre is supplied explicitly in a
frozen ORACLE-state fixture.  It tests whether SMITH consumes the correct
contract; it does not claim that the reduced profile can perceive η3 bonding.

## 6. Reading an output

The most relevant sections are:

- `TOPOLOGY`, `SYNTHONS`, `SYMMETRY`, and `INTERACTION_CENTERS`: perceived
  molecular state consumed by SMITH;
- `GIC`: build options, target dimension, primitive candidates, protected
  special coordinates, and the selected generalized internal coordinates;
- `SMITH_PROVENANCE`: package version, pinned MATRIX revision, and perception
  and fragment profiles.

For a nonlinear system of `N` atoms without disconnected fragments, the usual
target is `3N - 6`.  Always use the serialized `TARGET_RANK` and reported final
rank when special fragment or centre coordinates are present.

## 7. Verification and problem reports

From the repository root, run:

```bash
python -m unittest discover -s standalone/tests -v
```

When reporting a problem, include the operating system and architecture,
Python version, exact command, complete terminal output, input file, and
generated `xyzin` if one was produced.  Do not include credentials or private
tokens.  During this test phase, open an issue in the private repository or
send the report directly to the project owner.

Gaussian export, optimization workflows, and the full ORACLE validation
surface remain outside this reduced package and are available only in the
corresponding MATRIX development environment.
