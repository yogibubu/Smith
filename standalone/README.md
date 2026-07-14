# Standalone SMITH / SONIC

This package builds a frozen SONIC internal-coordinate contract without
installing or exposing the full MATRIX command suite.  It pins the MATRIX
implementation used by the manuscript and installs only the core, chemical
perception, engine-interface, and SMITH coordinate packages needed by the
builder.

SMITH and ORACLE have separate scientific responsibilities:

- ORACLE performs continuous molecular perception.  It owns the molecular
  graph and cycle basis, point-group operations and atom permutations, atom
  equivalence, effective atomic number, and the charge, covalency,
  delocalization, strain, bond-order, and synthon descriptors, together with
  the redundant primitive/Wilson-B source.  ORACLE develops the ideas
  introduced in PROXIMA and is validated as an independent release candidate.
- SMITH consumes that frozen molecular and primitive state and constructs the SONIC
  coordinate families, protected rows, rank reduction, homogeneous symmetry
  adaptation, analytic Wilson rows, and serialized coordinate contract.

Two input profiles are available.  An ORACLE-enriched `xyzin` file is the
recommended production input.  A plain SMITH extended XYZ is accepted through
a reduced bundled ORACLE profile so that the paper examples can be reproduced
with one small input file.  This convenience profile is not presented as a
replacement or release of the full ORACLE application.

## Install from GitHub

Python 3.11 or newer and `git` are required.

```bash
export SMITH_ENV=/path/to/your/smith-venv
python -m venv "$SMITH_ENV"
source "$SMITH_ENV/bin/activate"
python -m pip install --upgrade pip
python -m pip install \
  "smith-sonic @ git+https://github.com/yogibubu/Smith.git@agent/oracle-boundary-standalone-smith#subdirectory=standalone"
```

The environment directory is selected by the installer.  The package contains
no absolute path to the developer's machine.
The branch name can be replaced by `main` after the release candidate is
merged.

For a local checkout:

```bash
python -m pip install ./standalone
```

## Run

The examples are installed inside the package, so they can be run even without
a source checkout:

```bash
smith-sonic example water water.xyzin
smith-sonic inspect water.xyzin
smith-sonic example norbornane norbornane.xyzin
smith-sonic example formic-acid-water formic-acid-water.xyzin
smith-sonic example eta3-allyl-palladium eta3-allyl-palladium.xyzin
```

The same input files are visible under `examples/` in the GitHub repository and
can be passed explicitly to `smith-sonic build`.

For a production state prepared and validated by ORACLE, make the boundary
explicit:

```bash
smith-sonic build molecule.oracle.xyzin molecule.sonic.xyzin --require-oracle-state
```

Every output receives a `#SMITH_PROVENANCE` section recording whether SMITH
consumed a complete `ORACLE_STATE` or invoked the `REDUCED_ORACLE` convenience
profile.  Gaussian export and the wider validation/optimizer commands remain
available in the full MATRIX distribution; the standalone package deliberately
limits its surface to SONIC contract construction and inspection.  Each build
also writes a `.smith.out` coordinate report and a `.g16.gjf` Gaussian 16 input.
The G16 profile is the default, and non-totally symmetric coordinates are
written as `Frozen`.

With `--require-oracle-state`, the required production boundary is
`VALIDATION`, `TOPOLOGY`, `SYNTHONS`, `SYMMETRY`, and `PRIMITIVES`. The last
section freezes the primitive ordering, reference values and Wilson-B
fingerprint that SMITH consumes before constructing SONIC.

The formic-acid–water example exercises all six intermolecular fragment
coordinates.  The eta3 allyl–palladium example consumes an explicit, frozen
ORACLE interaction centre; its idealized geometry is an interface test, not a
computed chemical benchmark.  See [`MANUAL.md`](MANUAL.md) for the complete
installation, input, example, output, and troubleshooting guide.
