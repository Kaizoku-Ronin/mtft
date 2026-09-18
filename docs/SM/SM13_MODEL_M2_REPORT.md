# SM-13 — model M2: the mixed up sector reproduces the up-quark hierarchy (2026-09-14; DIAGNOSTIC)

## Assignment
Stack twists (L, c, a, b, d) = (0, 2, 2, 4, 2) in Z/6 (class 2 = chi_133, 4 = chi_100, 0 trivial):
   Q = (c, L-bar): cubic          u^c = (a, c-bar): untwisted        d^c = (b, c-bar): cubic
   L = (d, L-bar): cubic          e^c = (b, d-bar): cubic            nu^c = (a, d-bar): untwisted
   H_u in the chi_133 Higgs class (Q u^c = class 2), H_d in the chi_100 class (Q d^c = class 4); all triangles neutral;
   the shared-H_d constraint t_Q + t_dc = t_L + t_ec holds (2 + 2 = 2 + 2).
Spaces: chi_133 and chi_100 Higgs spaces from S_6(143, chi) with their own AL data (50 of 68 conditions, dim 18);
tensors pass the product test at 2.9e-8 (up) and 1.3e-8 (down/lepton); a character-violating pair fails it at > 1e-2
(control kept in the test).  HYM normalisation as in SM-08/09/10 (h = 0.2; low-height exclusion; Higgs Grams
regularised to 14 well-conditioned directions).

## Predictions over random Higgs directions (5th / 50th / 95th percentile) against observation
   up:     m_c/m_t   1.7e-3 / 7.1e-3 / 3.1e-2      obs 3e-3        m_u/m_t   2.2e-6 / 9.9e-6 / 5.0e-5     obs 6e-6
           44% of Higgs directions put BOTH up ratios within a factor 3 of observation.
   down:   m_s/m_b   1.8e-2 / 7.2e-2 / 3.4e-1      obs 2e-2        m_d/m_b   6.4e-4 / 2.8e-3 / 1.4e-2     obs 1e-3
   lepton: identical to down (same tensor, same H_d):                obs (6e-2, 3e-4): m_mu right, m_e heavy by 10.
   nu Dirac: (cubic L) x (untwisted nu^c) = the mixed class, hierarchy like the up sector (~1e-5) — not computed.

## What is exact and what is not
Exact: the three hierarchy classes of a fermion sector are the products of its two family classes — strong x strong
(m1/m3 ~ 1e-8, SM-09), mild x strong (~1e-5, this report and SM-12), mild x mild (~3e-3, SM-12) — with the family
classes set by the twist characters (certified spectra, SM-11); the character-neutrality algebra of the assignment.
Diagnostic: the ratio distributions (single mesh; Higgs direction a flat modulus).  Falsified within the class: any
assignment with equal down and lepton twists (identity), and the mixed lepton pairs (off by 40–700 either way).
The observed lepton first generation (3e-4) lies between the mild x mild and mild x strong classes.

## Tools (v0.30.9)
`rrspace.load_cubic_higgs`, `cubic_family_space_100`, `cubic_higgs_space(a)`, `mixed_yukawa(kind_A, kind_B, higgs_a)`;
frozen `x0143_cubic_higgs.npz` (1.3 MB).
