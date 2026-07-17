# Private release-candidate test for standalone SMITH

The purpose of this release-candidate test is to find installation, portability,
command-line, and documentation problems before freezing the manuscript and
standalone package.  The full ORACLE application is still in final testing and
is not part of this test.

## 1. Give the collaborator access

Keep `yogibubu/Smith` and `yogibubu/MATRIX` private.  Federico needs read
access to both repositories because the standalone wheel pins its MATRIX
components to an exact private commit.  In GitHub open:

`Settings` -> `Collaborators` -> `Add people`

Repeat the invitation in both repositories.  Invite Federico Lazzari
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
git switch --detach v0.1.0rc5
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
- eta3 allyl--palladium: 24 GICs, rank 24, including one protected
  centre--atom distance;
- the first three outputs: `PERCEPTION_PROFILE REDUCED_ORACLE`;
- the eta3 output: `PERCEPTION_PROFILE ORACLE_STATE`.

The source inputs are also visible under `standalone/examples/`.  Read
`standalone/MANUAL.md` before running the advanced examples.  The eta3 geometry
is an idealized interface probe and not a computed chemical benchmark.

## 4. Information to return

Before returning the report, run the complete consumer test from the repository
root:

```bash
SMITH_REQUIRE_DOWNSTREAM=1 python -m unittest discover -s standalone/tests -v
```

This confirms that the two advanced contracts are read directly by MORPHEUS
and converted into LINK SONIC coordinate models.

Please report:

- operating system, processor architecture, and Python version;
- whether authentication, cloning, and installation completed without manual
  intervention;
- the complete terminal output if a command fails;
- whether the distinction between reduced ORACLE perception and SMITH/SONIC
  coordinate construction is understandable;
- whether the commands and generated files are self-explanatory;
- the summaries printed for all four examples;
- whether the non-covalent output clearly represents six intermolecular
  degrees of freedom;
- whether the distinction between an ORACLE-supplied eta3 centre and the SMITH
  coordinate built from it is clear.

Do not broaden this test to the unreleased full ORACLE application.  Treat the
two advanced cases as packaging and interface tests; their inclusion in the
manuscript will be decided only after the collaborator's report.
