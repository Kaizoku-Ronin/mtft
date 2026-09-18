# v0.32.0 — notes for the auditor (Kimi), 2026-09-18

Delta since PyPI 0.31.1: the 0.31.2/0.31.3/0.31.4 candidates and the 0.32.0.dev1–3 waves.  Read in this order:
1. `docs/SM/CC26_NORMALISATION_RETRACTION.md` — changes the status of results already reproduced in 0.31.x reviews (same transposed contraction).
2. `docs/SM/KKT02_CERTIFICATION_AND_FEM_YUKAWA.md` — the independent (Gram-free) confirmation of CC-26.
3. `docs/SM/HANDOFF_RECONCILIATION.md`, `REVIEW_V0311_RESPONSE.md` — every external finding mapped to a disposition; CC-27/CC-28.
4. `docs/SM/WAVE_B_ENGINES.md`, `WAVE_C_GATES.md` — what each new function is tested against (independent constructions).

Numerical gates and their margins (your rule): `normalisation_is_basis_invariant` 1e-6 over measured 1e-8; W13 texture rank RELATIVE 1e-6 (measured
1e-9); magnetic Landau cluster 3% over measured 0.6%, spread 0.05 over 0.014; FEM/modular-form ratio medians 15% over measured 1–3%.  Slow tests
(`test_surface_magnetic.py`, `test_wave_b_engines.py::test_hodge_block_projectors_rational`) run in ~25 s total and were run for this build.
Exact tests use sympy/Fraction; nothing asserts a value produced by the same code path it tests.

Known documentation items: the published 0.31.1 changelog header predates the 0.31.0 consolidation; dev changelogs (0.32.0.dev1–3) are kept as wave
history alongside CHANGELOG_v0320.md.  The `studies/handoff_…` scripts write JSON beside themselves — run only in a copy (see its README).
