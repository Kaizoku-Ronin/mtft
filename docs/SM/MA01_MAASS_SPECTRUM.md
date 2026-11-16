# MA-01 — the curve's own tower: the Maass spectrum of Γ₀(143) by Hejhal's method (2026-11-15)

Modules `mtft.maass.spectrum`, `mtft.maass.hejhal`, `mtft.maass.scattering`, `mtft.maass.ford`, `mtft.maass.trace`; tests
`tests/test_ma01_maass_rung.py`, `tests/test_maass_solver.py`, `tests/test_maass_spectrum.py`, `tests/test_maass_scattering.py`,
`tests/test_maass_ford.py`, `tests/test_maass_trace.py`. Follows VE-5,
which named the curve's first Kaluza–Klein rung from a finite-element number (SPEC-01, DIAGNOSTIC) and asked that GC-01's thresholds be
checked against it. This wave supplies the rung and the tower above it by a second, independent method. It makes no statement about physics.

## 1. The first rung (MA-1, NUMERICAL)

The first Maass cusp form of Γ₀(143) has R = 0.3738518157, λ₁ = ¼ + R² = 0.3897652. It is odd under the mirror x ↦ −x and has
Atkin–Lehner signs (−, −). Three things support the number:

- Hejhal's method on two horocycles (y = 0.009680 and 0.008874, 632 coefficients, 788 sample points): the two solutions agree to
  2 × 10⁻¹⁵ and the Hecke relations, which are never imposed, hold to 10⁻¹⁴.
- The package's own solver, started within 10⁻⁵ of the lab's value, lands within 5 × 10⁻¹³ of it; a wrong sign or parity does not pass.
- SPEC-01's finite elements (`mtft.surface.spectral.maass_candidates(143, h=0.22, nx=7, k=14)`, 28,447 nodes) give 0.3925, 0.7 % above.

In VE-5's units the rung is √λ₁/R_C = 0.6243/R_C, the same 0.079–0.087 M_KK at the precision quoted there. What changes is its class.
VE-5's own text is not edited.

## 2. The tower (MA-2, NUMERICAL; its completeness in section 4)

41 cusp forms lie below R = 2.0 (λ < 4.25): even 4, 4, 4, 7 and odd 6, 6, 5, 5 in the sectors (+,+), (+,−), (−,+), (−,−). The η-twist
χ = ν_η¹², which is −1 round all four cusps, has 42 forms below R = 1.8 and no continuum. Each form is frozen with 400 coefficients and
SHA-1 provenance, and every eigenvalue was found again by the package's solver when the table was built (largest change 5 × 10⁻¹²;
horocycles agreeing to 1.2 × 10⁻¹¹; Hecke relations to 2.3 × 10⁻¹¹).

The two routes, side by side, for the first nine eigenvalues:

| λ, Hejhal | λ, finite elements | difference |
|---|---|---|
| 0.3897652 | 0.3925 | +0.70 % |
| 0.4480075 | 0.4516 | +0.81 % |
| 0.5586374 | 0.5633 | +0.83 % |
| 0.9642333 | 0.9720 | +0.81 % |
| 0.9987261 | 1.0078 | +0.91 % |
| 1.0018827 | 1.0113 | +0.94 % |
| 1.0586600 | 1.0692 | +0.99 % |
| 1.1924826 | 1.2027 | +0.86 % |
| 1.1924826 | 1.2049 | +1.05 % |

The last two rows are one double eigenvalue: the level-13 form of LMFDB (13.1, R = 0.97081541), which is old at level 143 and appears
with both signs at 11. The finite elements see the pair; the solver, run at level 13, gives the same R to 10⁻⁹.

Whether the lists are complete is the business of the trace formula: section 4.

## 3. The edge of the continuum (MA-3, EXACT; established mathematics)

The scattering matrix of the four cusps is Φ(s) = φ(s) M₁₁(s) ⊗ M₁₃(s), φ(s) = ξ(2s − 1)/ξ(2s),
M_p(s) = (p^{2s} − 1)⁻¹ [[p − 1, p^s − p^{1−s}], [p^s − p^{1−s}, p − 1]]. Φ(s)Φ(1 − s) = I, Φ is unitary on Re s = ½, every entry has
residue 1/(56π) at s = 1, and Φ(½) = −I. So each of the four Eisenstein series vanishes identically at s = ½: the continuum begins at ¼
and its edge carries no state.

## 4. Are the lists complete? The trace formula (MA-4, NUMERICAL)

`mtft.maass.trace` carries the Selberg trace formula of X₀(143) sector by sector: the odd forms (no continuum), the even forms (with the
continuum of section 3 and, in the trivial sector, the constant function), and the η-twist's total. The geometric side is arithmetic:
the area, classes of closed geodesics and of glide reflections counted through binary quadratic forms, the fixed points of the
Atkin–Lehner involutions (0, 4 and 20 for W₁₁, W₁₃, W₁₄₃), their mirrors (of total length 4 log(12 + √143) for W₁₄₃: twice the Pell
geodesic), the cusps and the mirror's cusp ends. The weights are frozen for lengths up to 17 (4,912 traces for the untwisted coset),
recomputed with PARI and equal to the lab's tables.

For a test function h ≥ 0 every form adds to the spectral side, so geometry − list is the sum of h over the forms a list lacks. A missing
form at R would add h(R). With h(r) = 2 exp(−r²/2σ²):

| Lists | Residuals | No form missing below |
|---|---|---|
| Untwisted, each of the eight sectors and parities (σ = 0.43, 0.45, 0.5) | all positive, 3 × 10⁻⁶ to 1.1 × 10⁻³ | R = 1.959 (the weakest sector; the strongest reaches 2.21) |
| η-twisted, all sectors together | positive, 1.8 × 10⁻³ at σ = 0.43 | R = 1.607 |

No residual comes near 2, which an eigenvalue below ¼ would add, so the formula sees no exceptional eigenvalue. Controls: PSL(2, ℤ)'s whole
formula (seven pieces of size 0.06 to 3.7) cancels to 10⁻¹⁰; level 11's three forms are reproduced to 2 × 10⁻⁹; a form removed from a list
is located to 0.002; a spurious or misplaced form breaks positivity. The lab's records MP17 to MP19 used a narrower h (σ = 0.35, lengths
to 23.5) and reached R = 1.98 untwisted and 1.8 twisted; those stay recorded results.

**What this rests on.** The class numbers and regulators behind the weights come from PARI, which is conditional on GRH in general. 150
weights are re-derived in the tests from Dirichlet's finite sum, which is not; and an error in a class number large enough to matter here
would break the cancellation at level 1, which shares its discriminants with the untwisted coset of level 143.

## 5. What is not here

- No statement about masses, forces or thresholds: GC-01's check against the tower is still to be done (VE-5).
- No trace formula for one twisted sector at a time: the twisted bound is for the eight sectors together.
- No lattice comparison: the glueball table and its forecasts, with their misses, are filed apart.

Provenance: the Modular Plane lab, records MP14 to MP20 (12 to 14 November 2026); the lab's number files are under
`studies/maass_2026nov/lab/`.
