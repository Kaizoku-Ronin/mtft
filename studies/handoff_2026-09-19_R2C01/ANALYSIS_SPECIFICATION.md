# R2C-01 — six-dimensional chirality assignment audit

Date: 2026-09-18. Baseline: the attached MTFT 0.32.0 source distribution.

This specification fixes the finite experiment before executing its enumeration.
The expected elementary-scalar obstruction is already suggested by AXG-03; this
is a closure of the unassigned-chirality question, not a blind discovery study.

## Fixed inputs

- Global gauge group: U(3)c × U(2)L × U(1)a × U(1)b × U(1)d, direct product.
- Stack order: c, L, a, b, d. Ranks: 3, 2, 1, 1, 1.
- Flux degrees: 0, -3, 3, 3, 0.
- One complex six-dimensional Weyl field in each oriented bifundamental i<j;
  ten independent chirality signs, with no added chiral fields.
- Internal spin bundle S0; family divisor D consists of three distinct
  W13-fixed CM points with both signs of u represented. Stack line bundles:
  (O, O(-D), O(D), O(D), O). No new flat twists are inserted.
- Desired signed indices for cL, ca, cb, Ld, ad, bd:
  (3, -3, -3, -3, 3, 3), respectively.
- Ordinary local nonderivative Yukawa interactions with elementary 6D scalars
  of charges eL-ea and eL-eb. Both Dirac and transpose/charge-conjugated
  Lorentz-scalar contractions are checked.

## Computations

1. Enumerate all 2^10 assignments, retaining the complete sector ledger.
2. Separately count preservation of the six family indices, preservation of
   every nonzero M1 sector index, and allowance of all four scalar bilinears.
3. Construct each fermion's complete I8 from the Chern character and A-hat
   genus. Check every sector against a separate sum over Cartan weights.
4. Verify the flux pushforward against an independently constructed 4D
   anomaly ledger. Do not multiply 6D anomaly coefficients by the number of
   internal zero modes.
5. Identify irreducible gravitational and color-cubic terms, without assuming
   they can be cancelled by ordinary 2-form Green–Schwarz factorization.
6. As a clearly changed control, solve the flux equations for scalar-compatible
   family chiralities. Record changed Higgs bundles, index-zero modes and
   anomaly witnesses. Do not relabel such a control as M1.
7. Check that an odd Clifford insertion removes the scalar representation
   obstruction. This is an operator-level control, not a completed vector or
   gauge-Higgs parent.

## Boundaries

This audit does not enumerate additional representations, boundary/orbifold
projections, higher-dimensional inflow, new twists, higher-derivative operators,
or new gravity/tensor spectra. No supersymmetry, vacuum, GeV scale, UV completion
or physical particle identification is supplied by a sign assignment.
