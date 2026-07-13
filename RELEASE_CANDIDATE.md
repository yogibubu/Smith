# SMITH 0.1.0rc1 acceptance record

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
| LINK and MORPHEUS consumption | `standalone/tests/test_downstream_consumers.py`, run with the `test` extra |

## Clean verification

```bash
export SMITH_RC_ENV=/path/to/your/smith-rc-venv
export SMITH_WHEEL_DIR=/path/to/your/smith-wheel-output
python3 -m venv "$SMITH_RC_ENV"
source "$SMITH_RC_ENV/bin/activate"
python -m pip install --upgrade pip build
python -m pip install "./standalone[test]"
SMITH_REQUIRE_DOWNSTREAM=1 python -m unittest discover -s standalone/tests -v
python -m build --wheel --outdir "$SMITH_WHEEL_DIR" standalone
```

The pinned MATRIX revision is
`711a8283ba50a064ae36a438f8c6e5a24e1fe0e7`. Change it only after repeating
the clean verification and all four packaged examples.

## Editorial gate

After the collaborator confirms this candidate, the remaining submission work
is limited to a final manuscript reading, regeneration of the manual PDF and
arXiv source archive, and creation of the stable release tag.
