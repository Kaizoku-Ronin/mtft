# mtft v0.8.0 — The Legend, Tier 11, Visuals & Release-Pipeline Hardening

## Added — The Legend (`mtft.legend`)
A cartographic map of the API. Every entry carries four fields ordinary
docs cannot: **nature**, **AG primitive signature** (glyphs * / S ^ d for
ITERATE/DIVIDE/ASSEMBLE/EXTRACT/CURVE), **epistemic tag + exactness class**
(Df/Pp/Pr/Conj/Heur/Cert × EXACT/CERTIFIED/DIAGNOSTIC/PHENO/GIVEN), and
**upstream derivation links**. Because the pipeline has zero free
parameters, `trace NAME` walks any object's chain down to the integers.
```
python -m mtft.legend             # the map (tiers, glyphs, tags)
python -m mtft.legend trace alpha_inverse   # ...ends at the integers
python -m mtft.legend status Pr   # epistemic audit / filter
>>> mtft.what("dirichlet_curvature")
```

## Added — Tier 11: Certificates & Standards
- **`jc_counterexample`** — the degree-7 Jacobian Conjecture counterexample
  (Alpöge / Claude Fable 5, Jul 20 2026). A self-certifying, **dependency-free**
  machine certificate: its own exact multivariate-polynomial engine
  (`Fraction`s) re-derives det DF = -2, the depressed fiber cubic via the
  tautological identity, S₃ monodromy (disc = -4S²p₃), the escape wall, and
  the empty-fiber missed curve on every call. `mtft.jc_verify_all()`. JC is
  false for all n ≥ 3; the 7·6·4 = 168 Bézout coincidence is flagged AG-D5
  and dismissed. Heavy Gröbner cross-checks: `scripts/jc/`.
- **`estimator_standards`** — the A.7 log-log slope discipline from the L1/L2
  sessions: `binned_log_slope` (terminal-bin leverage guard) and
  `stride_resonance_check` (flags the γ₃ × 6-decade near-resonance, 54.994
  cycles, that manufactured phantom drift).

## Added — Visuals (`viz/`)
- `hero_stiffness.png` (GitHub-dark stiffness landscape, generated from the
  arithmetic) and `stiffness_navigator.html` (interactive 3D, plotly-CDN).
  Regenerate both with `python viz/make_hero.py`.

## Added — Research provenance (`scripts/`)
- `scripts/jc/` — the July 20 verification session verbatim (SymPy/Gröbner),
  Kimi-K3 re-verified, plus `MECHANISM.md`.
- `scripts/riemann/` — `conspiracy.py`, `crossover2.py`, and the
  `L1_L2_corrected_RH_equivalence_draft.md` (two-engine validated). The
  corrected equivalence: **RH ⟺ limsup |(κ^Λ − κ_Main)(2πy)^(−3/2)| < ∞**,
  superseding the four false κ ≥ 0 propositions (now Th 1).

## Changed
- **README.md** rewritten: escaped-bracket badge corruption fixed,
  auto-updating PyPI badge, hero visual, three-ensemble table, 15→16-tier
  module map, key-identities block, Legend and JC sections.
- **CITATION.cff** → 0.8.0, refreshed abstract.
- **`critical_ensemble.py`** docstrings: Cauchy radius documented as a
  *conditioning* guard (ξ′/ξ is analytic at s = -2; audit Addendum K).

## Fixed
- **The v0.7.2 silent-skip release failure.** Root cause: `skip-existing: true`
  + unbumped version strings → rebuilt 0.7.1 artifacts, PyPI skipped duplicates,
  run went green publishing nothing. `skip-existing` is now `false`; a workflow
  **version-consistency guard** (tag must match pyproject + `__init__` +
  CITATION) and a `release.ps1` preflight make recurrence impossible.
- `CITATION.cff` version drift (was 0.7.0 through two releases); a stray
  `(v0.7.2)` docstring header in `__init__.py`.

## Notes
- PyPI sequence is 0.7.1 → 0.8.0; v0.7.2 remains GitHub-only. All v0.7.2 code
  (critical ensemble, Tier 10) ships here alongside Tier 11 and the Legend.

## Audit fixes on release day (Kimi K3, independent pass)
- `viz/stiffness_navigator.html` shipped dead: the generator interpolated
  `list(np.round(...))` straight into JavaScript, so every `ys` entry read
  `np.float64(...)` — a `ReferenceError` before Plotly ever ran. Root-caused
  in `viz/make_hero.py` (`map(float, ...)`), both assets regenerated, page
  verified rendering in a live browser.
- README quick start referenced two APIs that do not exist
  (`mtft.NewformSpace`, `mtft.HosotaniMTFT(N=3).effective_potential()`);
  replaced with the real calls, and **every** code snippet in the README is
  now executed verbatim in the release gate.
- Legend examples corrected to the real engine surface (`mtft.X0(143)`,
  `mtft.HosotaniMTFT().find_vacuum()`); the fictional `NewformSpace` entry
  became `X0_143_engine`.
- Hadamard wording tightened to the certified bound: residual 4e-9 at s=3
  (per-s value), < 1e-5 certified across s ∈ [3,10] — the blanket "~1e-9"
  claim was the Addendum K fine print come home to roost.
- Module/tier counts corrected (35 modules, 16 tiers incl. Tier 11);
  test-count badge 386 → 401; CITATION date-released → ship date.
- `python -m mtft.legend` no longer emits a runpy RuntimeWarning (Legend
  exports are now lazy via PEP 562 `__getattr__`).
- `scripts/riemann/zeros.json` (first 12 zero ordinates, dps 30) now ships
  so the conspiracy/crossover provenance scripts run out of the folder.
