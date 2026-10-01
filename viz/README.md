# viz/ — MTFT visualization gallery

## Standard Model interaction atlas

Open [`sm_interaction_atlas.html`](sm_interaction_atlas.html) in a browser for
the prebuilt, offline viewer. It includes reference vertex diagrams, tensor
details, search, sector filters, source provenance and an all-unmapped MTFT
evidence ledger. No installation is needed to view this file.

Regenerate it after changing the reference or viewer:

```bash
python -m mtft.interactions render -o viz/sm_interaction_atlas.html --force
```

See the [atlas guide](../docs/SM/INTERACTION_ATLAS.md) for model scope and
evidence requirements. This is a tree-level interaction reference, not a
matter simulation.

## Geometry components

Five self-contained React components (hooks only, no external state, no
build config beyond React 18).  Render any of them in a React sandbox,
a claude.ai artifact, or your own app.

| Component | Shows |
|---|---|
| `mtft_enneper.jsx` | Interactive Weierstrass–Enneper minimal-surface explorer (Paper 32 — the modular Enneper surface) |
| `MTFT_HyperbolicTiling.jsx` | Poincaré-disk SL(2,ℤ) tessellation colored by the three X₀(143) newform orbits |
| `MTFT_MonsterFingerprint.jsx` | Deterministic avalanche "fingerprint" built from a 13-round mixer (a JS *lookalike* of MonsterHash, not the real algorithm) |
| `X0_143_BurningMandelbrot.jsx` | Burning-Ship/Mandelbrot hybrid modulated by the 200 Hecke traces (first 50 = the verified `TRACE_TOTALS_50`) |
| `MTFT_DrawnLoop.jsx` | The drawn-loop stage: sketch closed walks on the 56-triangle dual of X₀(143), read off the 29 exact intersection numbers with the Manin arcs, the cuspidal/Eisenstein split, and Zeno survival curves (exact rates 35/16, 10/3); interface data exported by `studies/du02_cycle_space_map.py` |

These live at the repo root rather than inside the package on purpose:
a `mtft/viz/` package directory would shadow the Python module
`mtft/viz.py`, and the wheel has no business shipping JSX.


## Matter-01 U(1) lab (unreleased)

[Open the standalone lab](matter01_u1_lab.html) and
[scientific overview](matter01_u1_overview.png). The lab visualizes a supplied-
parameter massive QED reference, its amplitude checks, and conditional
breathing/gauge normalization. It displays the MTFT vacuum gates still OPEN.
Its Feynman diagrams represent amplitude terms, not literal particle paths.

See the [Matter-01 report](../docs/SM/MATTER01_REPORT.md). Rebuild and verify
with `python studies/matter01/run_all.py` from the repository root.

## Matter-02 mode lab (unreleased)

[matter02_mode_lab.html](matter02_mode_lab.html) visualizes the specified cusp
mode counts, local fat jets and an exact residue minor. The selected minor
has rank three; the full physical mass matrix and photon remain open.
[Scientific overview](matter02_overview.png) and
[derivation](../docs/SM/MATTER02_REPORT.md) accompany it.

Reproduce: `python studies/matter02/run_all.py`. Node tests check logic with a
mock DOM; browser layout has not been verified.

## Matter-03 Higgs lab (unreleased)

[matter03_higgs_lab.html](matter03_higgs_lab.html) shows the single-arrow
obstruction, paired classical existence branch and unselected common scale.
It compares 279 split charged classes with 237 for the certified pair and
the universal lower bound 160 in the aligned two-arrow background. Separate
INPUT panels illustrate a finite electroweak quartic and positive metric
normalization; they do not solve the B3 metric or predict masses.

[Scientific overview](matter03_overview.png),
[report](../docs/SM/MATTER03_REPORT.md), and
[reproduction guide](../studies/matter03/README.md) accompany it. Reproduce with
`python studies/matter03/run_all.py`. Node/mock-DOM logic is checked; actual
browser layout has not been verified.
