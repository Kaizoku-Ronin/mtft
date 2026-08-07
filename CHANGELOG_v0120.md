# v0.12.0 — the peel wave (2026-08-07)

New modules: peel, lchannels, marked_gap, coset_reps, gl2_peel,
hodge_polarization, ledger_peel. Tests: tests/test_peel_v0120.py (13, all fast).
Records: studies/peel_2026aug/ (11 scripts).

* Index correction encoded: F(s) = -zeta(s-2) zeta'(s-1); pole s=3, residue
  T_inf. Certified bulk expansion (gamma via psi(2); skipped orders; odd zetas
  via trivial-zero derivatives). Skeleton expansion with zeta zeros; constant
  term = -zeta'(-1)/zeta(-1) = 1 - 12 ln A (Glaisher's birth certificate).
  Normalized-residual boundedness <=> RH, per Addendum M architecture.
* Even-character selection rule + exact SU(p) channel splits (Gauss-sum
  coefficients; 34.8-35.3 digits). Channel zeros located for chi5 and the
  three nontrivial even characters mod 13; veil table (cubic mod 13 first
  zero t = 2.2731: x2.0e7 vs zeta). Per-channel GRH boundedness statement.
  Mode minimum = the vacuum choosing L-channel signs/phases.
* Coset bijection: SU(p) channels = split-torus characters = principal-series
  parameters of PSL(2,p); discrete series gauge-invisible; 168 = 1+11+13+143
  with 143 = dim(St11 x St13).
* Rung-4 two-temperature identity: m_inf = marked-level spacing at beta_c = 3
  (the T_inf pole); R_inf reproduced; spectrum ordering 2,3,5,4,7,11,9,8,13
  predicted for PR-5.1 audit; sigma* crossing; sigma=3 site-measure audit flag.
* GL(2) peel of 143a1: conductor certified by exact integers; a_p by point
  counting (Satake table of the Hodge session confirmed); eps = -1 (10.5 vs
  0.7 digits); L'(1) = 0.9457; RANK READ = 1.000000 +- 7e-7 at five depths —
  third independent route after periods and root number. Cuspidality: no
  X^-2 term. Gamma-factor fingerprint: log terms at every order.
* Hodge polarization data module (eta certificates; 11 lambda ratios,
  status computed-not-certified per source session; CC-01/CC-02 summaries).
* Falsifier records: null-model decoy test + Prox protocol
  (studies/peel_2026aug/null_model_test.py, two_kinds.py).

## Addendum — chat-sweep landings (2026-08-08, K3 staging)

du03 dispersion wave landed: studies/du03_dispersion.py (final, 673
lines) + its JSON ledger and run log; tests/test_dynamical_units.py gains
the 6-test du03 tier (gate total 465 + 6 + 13 = 484). Headlines: cusp-well
route CLOSED (harmonic annihilation 7e-16); minimal coupling via the
eta-odd Fiedler profile; first-order dispersion vanishes by a parity
selection rule (V eta-ODD, comm_rel 2.0); the SECOND-order dispersion is
the real one (eta-even, parity branches; f1 w2_mean = -0.00693959); the
box's first ruler — systole of X0(143) = 2 arccosh(2) = 2.633915793849633
EXACT (trace-4 element [[-15,-2],[143,19]], trace 3 excluded mod 13);
disciplined ratio census finds no anomaly (p = 0.696); PHENO chi_g =
lambda1/(2 pi f_21cm) = 3.054870e-11 s (no claims).

Hodge session record -> studies/hodge_2026aug/ (early dispersion,
grok_triage, commutant, metric + CC-02 correction; the metric-independent
obstruction: joint admissible interactions are scalars only; anchor
count stays 2).

Encryption session record -> studies/crypto_2026aug/ (Hecke-CDH
eprint 2025/1681 broken twice — eigenbasis attack at 2.6e-15, and an
exact rational solve needing only commutativity; integral Hecke order
glued at congruence primes {2,3,443}; coset correspondence graph
Ramanujan at all good p <= 43).

Still open (per INTEGRATION_v0120 chat sweep): mtft.crypto jacobian_order
module + 93-prime CSV; periods.gp eta/Omega variants.

Gate note (K3, disclosed): two peel-tier tests now pin mp.workdps(50) —
test_critical_ensemble leaves ambient dps at 20-30 and the 1e-30
assertions must not depend on ambient precision. Historical tests
untouched.
