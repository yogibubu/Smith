# SMITH 0.1.0rc7 acceptance record

SMITH is ready for collaborator release-candidate testing when every item below
passes from a clean Python 3.11+ environment.

| Requirement | Artifact or verification |
|---|---|
| Installable standalone package | `standalone/pyproject.toml`; build wheel and install it in a clean virtual environment |
| Short complete manual | `standalone/MANUAL.md` and `output/pdf/SMITH_Standalone_Manual.pdf` |
| Reproducible examples | packaged water and norbornane examples plus the two advanced cases below |
| Non-covalent complex | `formic-acid-water`; target rank 18 and six fragment coordinates |
| η3 transition-metal complex | `eta3-allyl-palladium`; target rank 24 and protected Pd-to-η3-centre distance |
| Synthetic coordinate output | the `.smith.out` sidecar generated for every build |
| Gaussian 16 input | the default `.g16.gjf` sidecar; non-totally symmetric coordinates are `Frozen` |
| LINK and MORPHEUS consumption | `standalone/tests/test_downstream_consumers.py`, run separately with `MATRIX_ROOT` pointing to a compatible clean checkout |

## Clean verification

```bash
export SMITH_RC_ENV=/path/to/your/smith-rc-venv
export SMITH_WHEEL_DIR=/path/to/your/smith-wheel-output
python3 -m venv "$SMITH_RC_ENV"
source "$SMITH_RC_ENV/bin/activate"
python -m pip install --upgrade pip build
python -m pip install "./standalone[test]"
python -m pytest standalone/tests/test_examples.py
python -m build --wheel --outdir "$SMITH_WHEEL_DIR" standalone
```

The separate integration gate deliberately does not add MATRIX as a standalone
runtime dependency.  Point `MATRIX_ROOT` to a compatible clean checkout and
run:

```bash
MATRIX_ROOT=/absolute/path/to/MATRIX SMITH_REQUIRE_DOWNSTREAM=1 \
  python -m pytest standalone/tests/test_downstream_consumers.py
```

The pinned scientific implementation revision is
`bc6d62140aab5fbcef8dd6fa7bc6c468debba69a`. Change it only after repeating
the clean verification and all four packaged examples.

## Editorial gate

The self-contained code bundle, manuscript, manual and independent source
archive are regenerated from the matching SMITH tag `v0.1.0rc7`.  The embedded
implementation provenance above identifies the frozen MATRIX source snapshot;
no external MATRIX checkout is required at installation or run time.  New
features are not admitted into this candidate without repeating the complete
gate.
