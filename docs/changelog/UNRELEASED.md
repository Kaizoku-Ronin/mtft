# Unreleased — after v1.0.0

Everything after v0.34.0 shipped in v1.0.0 ([release notes](CHANGELOG_v100.md)). Since then:

## CC-71 (2026-10-29): the Tano metric of the diatomic gases — S°(X₂) over the mass of X₂

- `mtft.tano_metric` divided the standard entropy of H₂, N₂ and O₂, which is tabulated per mole of molecules, by the molar mass of one
  atom. That doubled their Tano values: H 129.64 → 64.82, N 13.68 → 6.84, O 12.82 → 6.41 J g⁻¹ K⁻¹. `Element` gains `atoms` (atoms per
  standard-state formula unit; default 1, and 2 for H, N, O), and `tano_metric` now divides by `atoms × molar_mass`. The values before the
  correction stay exposed as `tano_metric.TANO_METRIC_RETRACTED_CC71` (append-only). `seebeck_diagnostic` computes its own T_a and now
  uses the same normalisation; nothing changes for the metals it is meant for.
- Second route: the three gases counted from first principles reproduce the stored S° per mole of molecules: H₂ 130.682 against 130.68
  (ortho and para hydrogen; counting every J once instead gives 5.71 more, close to R ln 2), N₂ 191.564 against 191.61, O₂ 205.128 against
  205.14 (only odd rotational levels exist, because ¹⁶O has spin 0; R ln 3 for the electrons' triplet). So the division by 2M is the right one.
- Downstream, the default screen (`materials_screening()`, `python -m mtft screen`) moves with it: LaH₁₀ λ_eff 2.092 → 1.444 and T_c
  930.1 → 750.5 K, H₃S 743.1 → 598.7 K, and the O/Ba proxy for YBCO 125.6 → 109.8 K; MgB₂, Al and Pb are unchanged. The screen's other weakness, the form
  Θ e^(−1/λ) without the (1 + λ) renormalisation or μ*, is a model choice rather than part of this correction and stays as it is (DIAGNOSTIC).
- Found on 2026-10-28 in a review of Roger's Tano repository workbook, where the same doubling affects seven rows outside the package
  (H, N, O, F, Cl, Br, I). `tests/test_cc71_tano_diatomic.py` (11 tests).

## CC-72 (2026-11-15): the Burning Ship fold exponents are a caveat — Paper 4's off-axis table is not one cascade

- `mtft.constants.BETA_X = +0.29` and `BETA_Y = −0.21` (and `DELTA_X`, `DELTA_Y`, defined from them) were fitted to the five points of
  Paper 4's off-axis table, periods 48 to 768. The four gaps between those points are 2.42, 1.81, 2.08 and 1.34 × 10⁻⁴: they shrink 1.81
  times over three doublings, where a period-doubling cascade shrinks δ³ = 101.8 times, and the third gap is larger than the second. The
  five points are not consecutive doublings of one cascade, so the exponents are unconfirmed.
- No number changes and nothing is removed (append-only). `mtft.burning_ship.cc72_cascade_check()` states the comparison and lists the
  dependents (`AnisotropicScaling`, `satellite_mass_scale`, `MassRatios.m_tau_over_m_mu`); `mtft.burning_ship.real_axis_cascade()` is
  the control: on the real axis the Burning Ship's map is x² + c, and its superstable parameters do shrink by Feigenbaum's δ
  (δ₁₀ = 4.6691951, c∞ = −1.401155189092). `DELTA_X_MEASURED` and `DELTA_Y_MEASURED` (Paper 16) were not examined.
- What would settle it: an off-axis cascade search that confirms each doubling before any fit.
- Found on 2026-11-06 in the Modular Plane lab (record MP05, Room 11). The gaps are as transcribed there from Paper 4's table; the
  transcription is to be checked against the paper at audit. `tests/test_cc72_burning_ship_cascade.py` (4 tests).

## MA-01 (2026-11-15): `mtft.maass` — the Maass spectrum of Γ₀(143) by Hejhal's method, the cusps' scattering matrix, the Ford domain

A new subpackage, ported from the Modular Plane lab's records MP14 to MP20. Import it on demand (`import mtft.maass`); nothing is added
to the top-level namespace. Report: `docs/SM/MA01_MAASS_SPECTRUM.md`. It makes no statement about physics.

- `mtft.maass.kbessel`: K_{iR}(x) by a trapezoid rule (R ≤ 3), 2.1 × 10⁻¹⁴ against mpmath on the gate's sample; mpmath above R = 3.
- `mtft.maass.domain`: the reduction of a point into the Ford domain of Γ₀*(N) (Γ₀(N) with its Atkin–Lehner matrices, N squarefree),
  with the Atkin–Lehner characters and, for odd N, the sign of the η-twist χ = ν_η¹². Gated on Eisenstein series (no Maass form
  involved) and on words in S and T.
- `mtft.maass.hejhal`: the solver. `refine` finds the R at which a coefficient agrees on two horocycles; `sigma_min` is the pole-free
  detector. Gates: LMFDB's form 13.1 at level 13 (0.97081541), three forms at level 11, three at level 143, SL(2, ℤ)'s first two
  eigenvalues with mpmath's Bessel function (marked slow), and controls that a wrong sign, parity or R does not pass.
- `mtft.maass.spectrum`: the frozen table. 41 cusp forms below R = 2.0 (even 4, 4, 4, 7; odd 6, 6, 5, 5 by Atkin–Lehner sector) and 42
  forms of the η-twist below R = 1.8, each with 400 coefficients and SHA-1 provenance. The first eigenvalue is λ₁ = ¼ + 0.3738518157²
  = 0.3897652 (odd, signs (−,−)). Every eigenvalue was found again by the package's solver when the table was built, within 5 × 10⁻¹² of
  the lab's. `evaluate` gives a form's value anywhere in the upper half-plane; `rederive` repeats the solver for one form.
- `mtft.maass.scattering`: Φ(s) = φ(s) M₁₁(s) ⊗ M₁₃(s) and the Eisenstein series of the four cusps (mpmath). Φ(s)Φ(1 − s) = I,
  unitarity on Re s = ½, residue 1/(56π) at s = 1, Φ(½) = −I: every Eisenstein series vanishes at the bottom of the continuum.
- `mtft.maass.ford`: the Ford domain of Γ₀*(N) (17 tops, 24 vertices and 16 distinct heights at N = 143; lowest vertex at height
  0.009851) with a panel quadrature: the area 14π to 2 × 10⁻¹¹, norms, overlaps, triple products, participation ratios (the lab's seven
  values to 10⁻³). Products whose character is not trivial are zero on the curve by symmetry and are returned as zero without integrating.
- Ledger (`mtft.matter.ledger`), wave MA-01 in piece IV: MA-1 upgrades the diagnostic figure of VE-5 (whose text stands) to λ₁ =
  0.3897652, with SPEC-01's finite elements as the independent route (0.3925, and the next eight eigenvalues within 1.5 %, the double one
  included); MA-2, the tower, with the caveat below; MA-3, the edge of the continuum.
- **Caveat.** The lists' completeness (below R = 1.98 untwisted, 1.8 twisted) rests on the lab's trace formulas, whose geometric side
  needs class numbers from PARI. Those formulas are not in the package: completeness is a recorded result, not a gate.
- Provenance: `studies/maass_2026nov/` (the build script and the lab's scripts and number files). Tests: `test_maass_solver.py` (10),
  `test_maass_spectrum.py` (9), `test_maass_scattering.py` (5), `test_maass_ford.py` (5), `test_ma01_maass_rung.py` (3).

