# MTFT v0.6.1 — Audit Coalescence Release

## 235 tests, all green (214 inherited + 21 new verified-data regression tests)

This release applies the independent-audit coalescence (July 2026): wrong
data corrected, superseded criteria replaced, and the verified X₀(143)
arithmetic spine made self-certifying.  No API removals — superseded
functions are retained with loud markers for historical continuity.

### Corrected data (bugs fixed)

| Bug | Fix |
|---|---|
| **B11 phantom polynomial** | `hecke_polynomial_f3_T2()` now returns the correct charpoly x⁶−10x⁴+2x³+24x²−7x−12 (6 real roots, Ramanujan-bounded). The old polynomial x⁶−x⁵−9x⁴+11x³+13x²−20x+8 was not a Hecke charpoly of this space; its complex pair 0.5732±0.3564i was an artifact |
| **Curve eigenvalues** | `CURVE_143A1.hecke_eigenvalues` corrected at p = 37 (−11), 41 (+10), 43 (−4), 47 (−4), 53 (+2), 59 (−1), 61 (−2); regression-locked by brute-force point counts |
| **B9 constants drift** | `DELTA_X_MEASURED = 7.2422`, `DELTA_Y_MEASURED = 3.3699` added; defined-vs-measured provenance documented |
| **B3 T_∞ factor 2** | `TORQUE_FULL = −ζ′(2)` documented as the *proved* Cesàro limit; `T_INF = TORQUE_FULL/2` as the physics-chain definition; locked by a Cesàro convergence test |
| **B2 ln\|M\| misprint** | `LN_MONSTER = 124.126423366` (AG's 124.01348 was a misprint); `GAUGE.alpha_inv_monster` = 137.0355 (3.5 ppm) |

### Superseded (retained with markers)

- **Tano Mass Formula** (`tano_mass_predictions`, `koide_angle_tano`, `A2_COMPLEX`):
  built on the phantom eigenvalue.  Functions still run; outputs carry
  `superseded: True`; `verify_complex_eigenvalue()` now verifies the
  *truth* (all roots real, phantom absent).
- **"RH ⟺ κ(y) ≥ 0"** (`rh_diagnostic`): false criterion (κ^Λ < 0
  unconditionally).  Replaced by the July 2026 Theorem 1 diagnostic:

### New: corrected RH equivalence (July 2026 Th 1) in `riemann.py`

- `skeleton_weights`, `skeleton_stiffness` — μ^Λ(y) with w_n = Σ_{dm=n} d²·mΛ(m)
- `delta_kappa_stable`, `normalized_oscillation` — Δκ in the stable normal
  form (Df 3), 𝒟(y) = Δκ·X^{−3/2}
- `envelope_slope`, `offline_quadruplet`, `corrected_rh_diagnostic` —
  envelope-slope test: ≈ 0 on true zeros (bounded), ≈ 1/2 − β₀ on
  synthetic off-line quadruplets (divergence).  Reproduces the draft's
  Appendix A.2/A.3 anchors.  Requires `mpmath` (new dependency).

### New: verified X₀(143) data in `x0_143.py`

`ORBIT_TRACES_VERIFIED` (p ≤ 23), full per-orbit traces to n = 50 with
trace-form totals (51/51 validated), coefficient-field polys +
discriminants + Galois groups (S₄, S₆), `ROOT_NUMBERS_LIST = (−1,+1,+1)`
re-derived from U_p eigenvalues, Frob₁₁ cycle type [2,4], and the
Rankin–Selberg coupling block: `rankin_selberg_Q(50)` recomputes
Q = +0.570292, Q_corr = +0.819186, strict sign rule 6/6 — converging to
the Session-4 NMAX = 1500 values (+0.587882, +0.813978).

### `falsify.py` honesty metadata (audit §3.1)

- `Prediction` gains `deviation_ppm`, `theory_tolerance_ppm`,
  `multiplicity`, `group`, and `theory_status`.
- `honest_report()`: σ-status vs theory-tolerance status side by side;
  correlated-duplicate groups (α_inv, α_s, charge, higgs_family,
  koide_family) with effective independent count 17/23.  The flagship
  α⁻¹ entry (#2) is judged against its pre-registered 1 ppm band
  (theory_status PASS) instead of the 12.1σ error-bar artifact.

### Tests

- `tests/test_x0_143_verified.py` (21 tests): point-count cross-check of
  the curve table, orbit-trace sum rules, q-expansion multiplicativity,
  prime-power recursion, charpoly Newton moments, all-roots-real +
  Ramanujan, exact Bareiss–Sylvester discriminants, unramifiedness at
  11/13, U_p rationality → root numbers, Rankin–Selberg reproduction,
  Cesàro identity.
- `tests/test_mtft.py`: the three phantom-asserting tests now encode the
  verified truth (supersession flags, real spectrum, historical
  continuity of the archived Koide angle).
- `examples/03_lepton_masses.py` narrative corrected (phantom retired).

### Note for v0.7.0 staging files (not yet in this repo)

`arithmetic_machine.py` and `busy_beaver.py` (Drive staging) hardcode a
wrong `ORBIT_TRACES` f₁ column (bug B1).  Fix by importing the verified
table: `from mtft.x0_143 import ORBIT_TRACES_VERIFIED`.
`mtft_core.py` is unaffected (computes a_p from the curve).
