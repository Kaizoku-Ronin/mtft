# SM-10 — normalised down and lepton sectors of M1 (2026-09-14; DIAGNOSTIC)

## Method
Twisted sections (chi_13 for d^c and H_d, cubic 133 for L, sextic 56 for e^c) evaluated at the 467,520 quadrature
points by a vectorised Gamma_0(143)-only height-maximising reduction with the composite matrix tracked exactly
modulo 13 (validated against the per-point reducer to 1e-15).  The 5.5% of points in the cusp regions of 1/11
and 1/13, which Gamma_0 alone cannot raise (q-series non-convergent there), carry 0.6% of the family norms and
0.3% of the Higgs norms (measured with the correctly evaluated untwisted sections) and are excluded from every
sector consistently.  HYM metrics as in SM-08/09 (same divisor, flat torsion factor).  Shared H_d direction for
down and lepton.  Twisted Higgs Gram nearly singular (eigenvalues down to 7e-13): four numerically null
directions removed; results unchanged by the regularisation.

## HYM norm spectra (family Gram eigenvalues / max) — the structural result
   Q (untwisted)          (7.8e-5, 6.9e-4, 1)
   d^c (chi_13)           (1.9e-4, 1.1e-3, 1)
   e^c (sextic)           (2.6e-5, 3.0e-4, 1)
   L (cubic)              (8.2e-3, 1.1e-1, 1)
The strength of the family hierarchy is set by the twist character: untwisted, quadratic and sextic sectors
spread over 1e4; the cubic sector over 1e2.

## Mass ratios over random shared Higgs directions (5th / 50th / 95th percentiles)
   m_s/m_b     4e-5 / 2e-4 / 3e-3         m_d/m_b     1e-8 / 6e-8 / 4e-7
   m_mu/m_tau  4e-3 / 1.4e-2 / 7e-2       m_e/m_tau   4e-6 / 1.6e-5 / 9e-5
   (m_mu/m_tau)/(m_s/m_b)  4 / 65 / 670    (m_e/m_tau)/(m_d/m_b)  24 / 290 / 2900
Observed (GUT-scale orientation): m_s/m_b ~ 2e-2, m_d/m_b ~ 1e-3, m_mu/m_tau ~ 6e-2, m_e/m_tau ~ 3e-4;
Georgi–Jarlskog ratios ~ 3 and ~ 1/3.

## Reading
- The up and down quark sectors come out with the same over-strong hierarchy: second family a factor 10–100 low,
  first family essentially massless (1e-8 of the third).  The cubic-twisted lepton sector is realistic in its
  spread: m_mu/m_tau brackets the observed value (95th percentile 0.07), m_e/m_tau a factor ~5–20 low.
- The cubic twist therefore over-corrects the down–lepton relation: instead of leptons heavier than down quarks
  by ~3 (Georgi–Jarlskog), M1 gives factors 65 and 290 at the median.  The exact statement that survives is
  qualitative and direction-independent: the twist character controls the hierarchy strength, and a cubic twist
  on the singlet stack splits leptons from quarks in the observed DIRECTION.
- Where M1 fails is specific: the quark first generation.  In this class that points at the quark sectors needing
  a milder twist (a cubic or sextic twist between the a/b and c stacks), which is a finite search with the
  tools now in hand.
Status: DIAGNOSTIC (single mesh h = 0.2; twisted Higgs Gram ill-conditioned; low-height exclusion ~1%).
Tools: `hym.reduce_gamma0_vec` (validated), `hym.twisted_section_values`.
