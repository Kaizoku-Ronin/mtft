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

## Field lines and the force on the internal torus (unreleased, 2026-10-14)

[wilson_line_torus.html](wilson_line_torus.html) draws the one-loop force −G⁻¹∇V on
the Wilson line of PA-01's relative U(1), on 143a1's own torus of holonomies, for
the three even spin structures: the sink at the jet vacuum, the source at the
untwisting point and the two saddles that PG-01's first filing missed (PG-6,
CC-66). The torus is drawn in its kinetic metric, so the parallelogram has
143a1's shape. The page sums the series itself (|m|, |n| ≤ 7, within 5 × 10⁻⁶ of
`mtft.matter.vacuum_atlas.wilson_line_torus` at the half-periods); see section 5
of the [PG-01 report](../docs/SM/PG01_ENERGY_TWO_CLOCKS.md).

[field_line_islands.html](field_line_islands.html) is the companion picture from
fusion: the standard map as the local model of a field-line Poincaré map, with
islands at the rational rotational transforms and the Ford circles of the same
rationals, and a slider through Greene's K_c ≈ 0.9716. Magnetic field lines
preserve area and circulate; a force from a potential drains to its minimum —
the two pages show both kinds of field on a torus. Both pages are
self-contained (no network, light and dark themes); browser layout was checked
with headless Chromium.

## Sphereland (unreleased, 2026-10-24)

[sphereland/index.html](sphereland/index.html) is the visual companion to the
repository: twenty-five interactive instruments, from spheres and the laws of
motion through knots one dimension up and the modular torus to mass gaps. The
newest three are the modular flip between the confined and the fluid filling
of a thermal torus, the strips of analyticity behind the Yang–Mills and
Navier–Stokes problems, and the spectral gap of the modular curve. Every panel
is labelled by layer (structure, model, derived, hypothesis); MTFT readings are
folded into labelled boxes, and a reading-mode switch shows the essentials or
everything. The page is self-contained apart from web fonts, which fall back to
system fonts offline, and has light and dark themes; its layout was checked with
headless Chromium at desktop and phone widths.

Its numbers are computed in the page or certified by the scripts in
[sphereland/checks/](sphereland/checks/README.md), each of which ends with
`ALL CHECKS PASSED`. The page and its scripts sit in a subfolder on purpose:
`MANIFEST.in` ships only top-level `viz/*.html`, so neither is in the sdist or
the wheel.
