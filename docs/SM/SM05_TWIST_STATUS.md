# SM-05 — the quadratic-character twist: what is established and what failed (2026-09-13)

## Established (EXACT)
- The three cuspidal 2-torsion classes are realised by eta quotients: 2t1 = div eta^(-11,1,11,-1),
  2t2 = div eta^(-13,13,1,-1), 2t3 = div eta^(-24,14,12,-2).  Only t3 has all-even exponents, so
  sqrt(v3) = eta(tau)^-12 eta(11tau)^7 eta(13tau)^6 eta(143tau)^-1 is a weight-0 eta quotient satisfying
  Newman's conditions, with multiplier character (11^6 13^5 / d) = chi_13(d).  Hence t3 = (−65,35,30,0) is
  the Dirichlet twist, and chi_13 is the only even quadratic character mod 143 (chi_{-11}, chi_{-143} are odd:
  dim S_4 = 0 for them).  The twists t1, t2 carry eta^{1/2}-type multipliers (non-Dirichlet).
- Twisted sections are chi_13-character forms divided by the same denominators: dim S_2(143,chi_13) = 12,
  S_4 = 40, S_6 = 68, and the condition counts 40 − 37 = 3, 68 − 50 = 18 give the twisted family and Higgs
  spaces with all conditions independent (smallest singular values 7.5e-4, 4.6e-4).  In particular the
  chi_13-twisted down/lepton sector is Brill–Noether pure and h^0(S0 (x) t3) = 0.
- Character forms are evaluated at CM points by Gamma_0(143)-only reduction with the chi_13(d) factor tracked
  on the composite integer matrix; the chain-rule derivative agrees with finite differences to 1e-8.

## Failed (recorded, not shipped)
The twisted Yukawa relation G_i G^t_j = f_K sum_k Y G'^t_k has expansion residual 8e-2 (untwisted control
3e-12).  Diagnosis: the quotients (G_i G^t_j)/f_K are weight-6 chi_13-forms (4e-6) but are NOT in the
constructed twisted 18-space, so at least one imposed condition is mis-stated for character forms — the cusp
orders at cusps 0 and 1/11, where chi_13 interacts with the cusp widths (expansions there are in Q(zeta_13),
and the first basis form vanishes identically through q_w^8 at cusp 0), are the suspects.  The correct
condition is ord_omega >= (6, 6, 0, 0) expressed in the right local parameter for the twisted forms; until it
is, the twisted tensor is not trusted.  The product test is doing exactly its job.

## Consequence for the physics
The mechanism that could lift the up/down/lepton degeneracy — a chi_13 twist separating L_b from L_a — is
mathematically available and pure; whether it produces Georgi–Jarlskog-type structure is still open, pending
the cusp bookkeeping fix.  No number from this section is a prediction.
