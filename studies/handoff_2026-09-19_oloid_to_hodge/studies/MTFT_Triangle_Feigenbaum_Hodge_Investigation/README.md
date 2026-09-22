# Reproduce the investigation

Start with `REPORT.md`. It contains the derivations, mathematical qualifications,
source references, and distinction between exact statements and numerical evidence.

From this directory, run:

```bash
python3 scaling/compute_scaling.py
python3 geometry/verify_geometry.py
python3 plot_results.py
```

The first two commands need only the Python standard library. The plot needs
Matplotlib. Generated results are written beside their programs; the plot is
`investigation.png` at the bundle root.

Contents:

- `scaling/`: 60/90-digit quadratic-map superstable centers through period 1024,
  residuals, and finite scaling ratios.
- `geometry/`: exact finite checks of discriminants, Hodge eigenspace dimensions,
  elliptic isogeny identities, and numerical AGM period calculations.
- `SOURCE_AUDIT.md`: the reviewed MTFT snapshot and bounded source-audit findings.
- `MANIFEST.json`: checksums of the distributed files, excluding the manifest itself.

The finite computations support the general proofs in the report. They are not
interval certificates, a proof of a new universality theorem, or a derivation of
physical constants. No MTFT installation is needed and no repository files were
modified by this investigation.
