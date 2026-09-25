# MTFT v0.7.0 — Computation, the Jacobian Engine & Sonification

## 343 tests, all green (235 inherited + 108 new)

This release lands the staged v0.7.0 additions: the computation-theory
tier (five-primitive decomposition, Hecke-constrained Busy Beaver,
arithmetic Wick rotation), the Paper 30 Jacobian stiffness engine
re-based onto verified data, and the sonification module.  Everything
was integrated under the v0.6.1 audit-coalescence discipline: every new
module sources its X₀(143) data from `mtft.x0_143` / `mtft.constants`
instead of carrying local copies.

### New modules (5)

| Module | What it does |
|---|---|
| `arithmetic_machine` | Five-primitive decomposition of computation (Turing/λ/recursive), configuration-space geometry, computational stiffness μ_C, halting-surface topology, arithmetic entropy, computation-physics bridge. Ships with its 56-test suite |
| `busy_beaver` | Hecke-constrained Busy Beaver hierarchy BB_Fatou ≤ BB_g ≤ BB_Hecke ≤ BB on the genus-13 tape window, with CLI (`python -m mtft.busy_beaver`). Exact pins: BB_Hecke(1)=BB_Hecke(2)=0, BB_unc(1)=7, BB_unc(2)=13, BB_Hecke(3)=9 |
| `arithmetic_wick` | Laplace (mass-gap) vs Dirichlet (zeta) statistical ensembles on the shared holonomy weights w_n; the Mellin transform as arithmetic Wick rotation; critical-depth comparisons |
| `jacobian` | The Paper 30 3×3 Jacobian stiffness engine (port of `mtft_core.MTFT`), re-based onto the verified `ORBIT_TRACE_F1/F2/F3` tables |
| `music` | `VacuumSonifier`, `ModularScale` (supersingular / Farey / Hecke / Koide scales), `MonsterComposer` — deterministic, MonsterHash-driven composition; WAV output with stdlib `wave` |

New curated exports in `mtft/__init__.py` (Tier 9); version → 0.7.0.

### v0.7.1 flags settled early (independent audit computation)

- `L_VALUES["L'(f1, 1)"]` corrected 0.791 → **0.945696**: the audit's
  approximate-functional-equation computation (w = −1, converged by
  n = 100) agrees with PARI `ellL1(E,1) = 0.94570`; the old value was
  simply wrong.
- `EllipticCurve143a1.analytic_rank` corrected 0 → **1**: root number
  −1 forces odd analytic rank, L′(1) ≠ 0 forces exactly 1, and PARI's
  Mordell–Weil rank 1 (+ Kolyvagin) agrees.

### Audit-canon enforcement

- **Pre-audit `ORBIT_TRACES` removed** from `arithmetic_machine` and
  `busy_beaver`.  The staged copies disagreed with
  `ORBIT_TRACES_VERIFIED` at **all 9 primes** and were internally
  inconsistent with their own (correct) `HECKE_TRACES` at 7 of 9.  Both
  modules now re-export the verified table — single source of truth.
  (`busy_beaver.dominant_sector` and the sector-resolved statistics
  changed accordingly; the BB values themselves were never affected.)
- **`Y_C = 0.1812` → `CriticalDepths.y_conf = 0.18174`** in
  `arithmetic_wick` and `music` (the stale rounded value predated the
  audit).
- `music` dropped its unused, unsuffixed measured `DELTA_X`/`DELTA_Y`
  copies (the B9 provenance rule: measured values carry the
  `_MEASURED` suffix in `mtft.constants`).
- **`jacobian` data provenance**: per-orbit traces come from the
  verified n ≤ 50 tables instead of `mtft_core`'s reconstruction
  (exact at 6 primes + Hecke prime-power recurrences + a 40/60
  proportional split — wrong at 37 of the first 50 entries).  The n ≤ 50
  truncation is exact to ~1e−20 at y ≈ 0.18; a RuntimeWarning guards
  y < 0.08.  Effect on Paper 30's quoted numbers at y = 0.1812:

  |  | pre-audit (Paper 30) | verified data (v0.7.0) |
  |---|---|---|
  | λ₁, λ₂, λ₃ | 1.52628, 6.42705, 45.46824 | 1.531802, 6.591154, 45.466583 |
  | λ₃/λ₂ vs δ_x^meas = 7.2422 | 7.0745 (2.3%) | 6.8981 (4.8%) |
  | λ₂/λ₁ vs δ_F = 4.6692 | 4.2109 (9.8%) | 4.3029 (7.9%) |

  The eigenvalue hierarchy and couplings (e-μ ≈ 0.02, μ-τ ≈ 0.72) are
  unchanged; the λ₃/λ₂ ↔ δ_x agreement weakens, λ₂/λ₁ ↔ δ_F improves.
- New `tests/test_v070_canon.py`: cross-module identity locks (the
  re-exports must *be* the verified objects), trace-total consistency,
  canonical depths everywhere, no-phantom-roots regression.
- **`HECKE_TRACES` n = 1–200 independently re-verified** at integration
  time against the Eichler–Selberg trace form of the level-143 weight-2
  *newspace* (PARI/GP 2.15.4 `mftraceform([143,2],0)`): exact match, all
  200 entries.  (The full-cuspidal trace form differs by exactly twice
  the 11a1 eigenvalues — the oldspace — which confirms the decomposition
  independently.)

### Fixes

- `arithmetic_wick.laplace_ensemble`: the underflow branch constructed
  `LaplaceEnsemble(beta=y, …)` against a field named `y` — a guaranteed
  `TypeError` whenever Z underflowed (y ≳ 60).  Fixed + regression test.
- `tests/test_arithmetic_machine.py`: imports rewritten for the src
  layout (`from mtft.arithmetic_machine import …`); the old
  `sys.path` hack inserted the tests directory itself.
- `busy_beaver.bb_targeted` gained a `router_samples` knob (was a
  hardcoded 50,000/writer); dead `ROUTER_SAMPLES` code removed.
- ruff F401 pass over the new modules (20 unused imports removed).

### Packaging

- New optional extra **`[lhc]`** = `uproot>=5.0`, `awkward>=2.0` — the
  lazy dependencies of `lhcb_analysis`, previously undeclared.
- The staging README's proposed `[music]` extra was **not** added:
  `music` needs only numpy (already a core dependency) — WAV output
  uses the stdlib `wave` module.
- `CITATION.cff` updated to 0.7.0 (it had been stale at 0.5.0).
- Release plumbing: `publish.yml` now passes `skip-existing: true`
  (double-publish safety), and `release.ps1`'s local twine upload was
  removed — PyPI publishing is owned by the trusted-publishing workflow,
  which gates on the full test suite.
- Repo-root **`viz/`**: four self-contained React visualizations
  (Enneper surface, hyperbolic tiling, Monster fingerprint, Burning
  Mandelbrot).  Deliberately *outside* the package: a `mtft/viz/`
  package directory would shadow the existing `mtft/viz.py` module.
- **`scripts/pari/`**: the Paper 32/33-era PARI/GP provenance scripts
  with their verified output logs (banner lines carrying local paths
  trimmed).  Spot-checked against `x0_143.py`: field minpolys, Hecke
  charpoly roots, f₁ q-expansion, and cusp widths all match.

### Deliberately not included

- `tests/test_monsterhash_integrity.py` (the 5-battery suite) — the
  MonsterHash line now lives in its own repository; the suite is
  print-only (no asserts) and pins no generation vector.  Kept in
  staging.
- `mtft_core.py` as-is — superseded by `jacobian.py`; its constants
  block duplicated `mtft.constants` with pre-audit values.

### Flagged for v0.7.1

- `L_VALUES["L'(f1, 1)"] = 0.791` vs PARI `ellL1(E, 1) = 0.94570` —
  normalization to reconcile.
- `EllipticCurve143a1.analytic_rank = 0` sits inconsistently beside
  `root_number = −1` in `x0_143.py`.
