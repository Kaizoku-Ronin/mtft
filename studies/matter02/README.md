# Matter-02 reproduction

This is the focused collinear-mode and bulk-coupling milestone, based on
MTFT 0.34.0 plus the unreleased Matter-01 additions. No release was published.

From the candidate repository:

```bash
python -m pip install -e .
python -m pip install pytest matplotlib
python studies/matter02/run_all.py
```

The runner prints verdicts and writes full logs and JSON under `results/`.
The archived run passes 148 cohomology, 34 Yukawa, and 38 photon arithmetic
gates; 133 focused Python tests; and 138 Node numerical/mock-DOM checks.
These counts describe finite checks, not a formal proof assistant validation
of sheaf cohomology. It does not run the complete historical release suite.

For just the standard-library exact geometry:

```bash
python studies/matter02/cohomology/cohomology_certificate.py
python studies/matter02/yukawa/yukawa_certificate.py
```

Public API: `python -m mtft.matter.collinear --output collinear_report.json`.
The photon audit needs the repository's existing SymPy dependency. The
scientific figure requires Matplotlib. Node.js is optional; its absence is
recorded as a skipped lab logic check. The prebuilt HTML needs no server or
external scripts. Browser layout was not tested here.

Read [the main report](../../docs/SM/MATTER02_REPORT.md), or open
`docs/SM/MATTER02_REPORT.md` from the repository root.
The main report names the chosen cusp, flat-twist assumptions, parent action
and opposite chirality conventions. The separate derivations and certificates
retain source hashes and exact rational witnesses.

The selected minor has determinant −1 in holomorphic coordinates and proves
rank at least three. It is a slice of an 18 × 44 × 10 interaction tensor;
no Higgs expectation value, selected physical families, kinetic Gram matrices
or measured masses have been derived. Physical photon masslessness remains
open; example axion matrices are explicitly INPUT probes.
