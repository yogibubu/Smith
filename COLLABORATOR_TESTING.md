# Private release-candidate test for standalone SMITH

The purpose of this release-candidate test is to find installation, portability,
command-line, and documentation problems before freezing the manuscript and
standalone package.  No other application is part of this test.

## 1. Give the collaborator access

Keep `yogibubu/Smith` private during release-candidate testing.  Federico needs
read access only to this repository; all runtime sources are included.  In
GitHub open:

`Settings` -> `Collaborators` -> `Add people`

Invite Federico Lazzari
(`federico.lazzari@sns.it`).  If GitHub cannot resolve
the email address, ask him for his GitHub username.  He must accept the
invitation before cloning or installing the private repository.

## 2. Clone and install

The collaborator should install GitHub CLI, then run:

```bash
gh auth login
gh auth setup-git
export SMITH_CHECKOUT=/path/to/your/Smith
export SMITH_ENV=/path/to/your/smith-venv
git clone https://github.com/yogibubu/Smith.git "$SMITH_CHECKOUT"
cd "$SMITH_CHECKOUT"
git switch --detach v0.1.0rc7
python3 -m venv "$SMITH_ENV"
source "$SMITH_ENV/bin/activate"
python -m pip install --upgrade pip
python -m pip install "./standalone[test]"
```

The collaborator must replace both `/path/to/your/...` placeholders with
directories of their choice.  SMITH does not assume the owner's directory
layout.

On Windows PowerShell, activation is:

```powershell
<smith-env>\Scripts\Activate.ps1
```

The detached checkout is intentional: it guarantees that the test uses the
immutable release-candidate tag rather than a moving branch.

## 3. Run the packaged examples

```bash
smith-sonic --version
smith-sonic example water water.xyzin
smith-sonic inspect water.xyzin
smith-sonic example norbornane norbornane.xyzin
smith-sonic inspect norbornane.xyzin
smith-sonic example formic-acid-water formic-acid-water.xyzin
smith-sonic inspect formic-acid-water.xyzin
smith-sonic example water-dimer water-dimer.xyzin
smith-sonic inspect water-dimer.xyzin
smith-sonic example benzene-water benzene-water.xyzin
smith-sonic inspect benzene-water.xyzin
smith-sonic example eta3-allyl-palladium eta3-allyl-palladium.xyzin
smith-sonic inspect eta3-allyl-palladium.xyzin
```

Each example must also create a `.smith.out` report and `.g16.gjf` input beside
the requested `.xyzin` file.  The water report must contain both `state=ACTIVE`
and `state=FROZEN`, and the Gaussian input must contain `(Frozen)` rows.

Expected summaries are:

- water: 3 GICs, rank 3;
- norbornane: 51 GICs, rank 51;
- formic-acid--water: 18 GICs, rank 18, including three fragment translations
  and three fragment orientations;
- water dimer: 12 GICs, rank 12, including three fragment translations and
  three fragment orientations;
- benzene--water: 39 GICs, rank 39, including three fragment translations and
  three fragment orientations;
- eta3 allyl--palladium: 24 GICs, rank 24, including one protected
  centre--atom distance;
- the first five outputs: `PERCEPTION_PROFILE STANDALONE_MINIMAL`;
- the eta3 output: `PERCEPTION_PROFILE FROZEN_STATE`.

The source inputs are also visible under `standalone/examples/`.  Read
`standalone/MANUAL.md` before running the advanced examples.  The eta3 geometry
is an idealized interface probe and not a computed chemical benchmark.

## 4. Information to return

Before returning the report, run the complete standalone test from the
repository root:

```bash
python -m unittest discover -s standalone/tests -v
```

This confirms the four input profiles, both output formats, the non-covalent
examples and the supplied interaction-centre contract.

Please report:

- operating system, processor architecture, and Python version;
- whether authentication, cloning, and installation completed without manual
  intervention;
- the complete terminal output if a command fails;
- whether the standalone topology/primitive boundary and SMITH/SONIC
  coordinate construction are understandable;
- whether the commands and generated files are self-explanatory;
- the summaries printed for all six examples;
- whether the non-covalent output clearly represents six intermolecular
  degrees of freedom;
- whether the distinction between a supplied eta3 centre and the SMITH
  coordinate built from it is clear;
- whether the Gaussian 16 compatibility limitations and the authoritative
  native/human outputs are clear.

Treat the advanced cases as packaging and interface tests; do not add
application-level optimization or scan tests to the standalone acceptance gate.
