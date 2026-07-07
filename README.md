# SMITH / SONIC manuscript

Working LaTeX draft for the SMITH coordinate-construction manuscript and the
SONIC coordinate family.

Build:

```bash
latexmk -pdf main.tex
```

The draft is intentionally source-close: the implementation-status table tracks
the current SMITH/SONIC code paths, including the special-coordinate fragment
model and the active pseudo-bond / pseudo-cycle branch.

The GF/PED scaling add-on is tracked in
`data/gf_ped_scaling_probe.json`.  It records the cyclohexanol Gaussian
B3LYP/6-31G(d) Opt/Freq probe, the SQM-style family scale classes, and the
representative PED rows reported in the manuscript.

The construction and B-matrix timing microbenchmark is tracked in
`data/construction_scaling_benchmark.json` and can be regenerated with:

```bash
PYTHONPATH=$(find /Users/vincenzobarone/Documents/git/software/matrix/packages -maxdepth 2 -type d -name src | paste -sd: -) \
  python scripts/run_construction_scaling_benchmark.py
```
