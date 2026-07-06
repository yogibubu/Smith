# NEO manuscript

Working LaTeX draft for the NEO coordinate-generation manuscript.

Build:

```bash
latexmk -pdf main.tex
```

The draft is intentionally source-close: the implementation-status table tracks
the current NEO code paths, including the special-coordinate fragment model and
the active pseudo-bond / pseudo-cycle branch.

The GF/PED scaling add-on is tracked in
`data/gf_ped_scaling_probe.json`.  It records the cyclohexanol GDV
B3LYP/6-31G(d) Opt/Freq probe, the SQM-style family scale classes, and the
representative PED rows reported in the manuscript.
