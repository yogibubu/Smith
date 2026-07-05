# NEO manuscript

Working LaTeX draft for the NEO coordinate-generation manuscript.

Build:

```bash
latexmk -pdf main.tex
```

The draft is intentionally source-close: the implementation-status table tracks
the current MATRIX/NEO code paths, including the distinction between the
production special-coordinate fragment model and the disabled pseudo-bond /
pseudo-cycle branch.
