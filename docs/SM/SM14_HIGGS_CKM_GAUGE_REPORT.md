# SM-14 — Higgs directions, CKM, neutrino Dirac sector, gauge normalisation (2026-09-14)

## Neutrino Dirac sector (EXACT within M2)
L and Q share the cubic family space; nu^c and u^c share the untwisted one; both couple to H_u.  The neutrino Dirac
tensor is therefore identical to the up tensor: m_nu^Dirac = (m_u, m_c, m_t) at the compactification scale, up to
the same H_u direction.  Light neutrinos need a Majorana sector for nu^c (a singlet coupling nu^c nu^c in the
(a, d-bar)^2 class) with M_R ~ m_D^2 / m_nu (6e14 GeV for the ILLUSTRATIVE one-flavour inputs m_D = 173 GeV, m_nu = 0.05 eV — not a derived scale).
Amendment (review V0311): nu^c carries B-L = +1, so a bare Majorana mass is forbidden while B-L is unbroken; the seesaw
needs a B-L = -2 scalar or another B-L-breaking mechanism — outside the present class.

## Higgs directions (DIAGNOSTIC)
The Higgs direction is a flat modulus at leading order; per sector it is a fit.  Best fits reproduce the up ratios
(6e-6, 3e-3) and the down ratios (1e-3, 2e-2) to machine precision (14 complex parameters per sector for two
targets), so the observed quark masses are inside the reachable set — no obstruction.  The typicality is in the
distributions of SM-13: 44% of random H_u directions within a factor 3 of both up ratios; the down targets at the
5th percentile of the all-cubic distribution.
Leptons at the fitted H_d (identity with the down sector): (m_e/m_tau, m_mu/m_tau) = (1.0e-3, 2.0e-2) against
observed (3e-4, 4.5e-2): factors 3.3 and 2.3.  M2 lands exactly on the naive SU(5) relations m_e = m_d, m_mu = m_s
(and m_b = m_tau by the same identity); the Georgi–Jarlskog factors (3, 1/3) are what the class lacks.

## CKM (DIAGNOSTIC) — a prediction of the shared Q index
|V| = |U_u^dag U_d| from the left singular vectors of M_u(v_u), M_d(v_d).  Random direction pairs (5/50/95%):
|V_us| 0.06 / 0.225 / 0.63,  |V_cb| 0.09 / 0.39 / 0.85,  |V_ub| 0.016 / 0.072 / 0.24.  Typical mixing is LARGE;
only 0.2% of random pairs have |V_us| < 0.3, |V_cb| < 0.1, |V_ub| < 0.02 together.  At the mass-only best fits
|V_us| = 0.455, |V_cb| = 0.18, |V_ub| = 0.058.  A joint fit (four mass ratios + three CKM moduli, 52 real
parameters) reproduces all seven observed values exactly: the observed quark sector is realisable in M2; the
small observed mixing is atypical, i.e. the physical Higgs directions are special — the strongest evidence yet
that the modulus must be lifted by a potential before masses and mixings become predictions.  (The random-pair
median |V_us| = 0.225 coincides with the observed value; recorded as a curiosity, not a result.)

## Gauge normalisation (EXACT)
All stacks wrap the same curve with one 6D coupling g: g_3 = g_2 = g at the compactification scale.  With
N = (3, 2, 1, 1, 1) on (c, L, a, b, d) and Y = (1/6, 0, -1/2, 1/2, -1/2)·Q: 1/g_Y^2 = (2/g^2) sum N_x y_x^2 =
(5/3)/g^2, so sin^2 theta_W = 3/8 and g_3 = g_2 = sqrt(5/3) g_Y — the canonical unified normalisation, from the
flux-model hypercharge rather than from an SU(5) embedding.  Tool: `smflux.hypercharge_normalisation`.

## Status of the Standard Model programme after SM-01…14
Exact: three-family anomaly census; four-stack no-go; W11 parity; purity theorem; transposition theorem;
selection rules; spin-structure invariant 7.8e-5; twist-character -> hierarchy class; hierarchy classes of
sectors (strong x strong / mild x strong / mild x mild); nu Dirac = up in M2; sin^2 theta_W = 3/8.
Diagnostic: mass-ratio distributions and CKM distributions (single mesh, regularised Higgs Grams).
Fits, not predictions: individual masses and mixings (flat Higgs modulus).
Not in the class: Georgi–Jarlskog factors; Majorana neutrino masses; absolute Yukawa normalisation (the HYM metric
constants and the 6D coupling); the Higgs potential lifting the modulus; RG running below the KK scale.
