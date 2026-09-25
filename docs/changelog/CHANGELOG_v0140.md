# v0.14.0 — the promotion wave (2026-08-08)

Four new modules (1147 lines): moments, curvature, hecke, eisenstein.
Tests: tests/test_promotion.py (29). Everything was first proven in a
gated study; module docstrings carry the epistemic class of every claim
and name the study that certified it.

* moments (Tano weight closed forms): <w>_beta = -zeta'(beta+1) EXACT;
  <w^2>_beta = T(beta) = zeta''(v)C0 - 2 zeta'(v)C1 + zeta(v)C2
  (v = beta+2) EXACT closed form; susceptibility chi_w and
  Cov(log n, w) = zeta''(beta+1); <w^3> = U by the triple-Euler engine
  (Pr, dps-robust to 1.7e-35); Amari-Chentsov cumulant tensor; cold
  constants at the Hagedorn edge (T(1) = 1.70276979154901697001...,
  U(1) = 4.42947284842615649140232922679...).
* curvature: Brioschi machinery; the HESSIAN CANCELLATION THEOREM (the
  fourth cumulant contributes nothing to the curvature of any
  exponential family; Gaussian convention lock K = -1/2); the
  sign-changing profile of the (beta, lambda) manifold — flat at the
  Hagedorn wall with closed-form slope A = 0.423657463797093...,
  summit beta* = 4.593591164956 / K* = 1.1956959819919385, flat
  temperature beta_0 = 8.8565170425, cold dive K ~ -c (6/5)^beta with
  c = 0.270126465305425; finite_atom_curvature with the rigidity locks
  K = 1/4 on {1,2,3} and {1,2,3,4}, atom 5 flips the sign.
* hecke (Manin/Merel engine for X0(143)): exact 29-dim modular-symbol
  quotient / 26-dim cuspidal homology from P^1(Z/143); Merel matrices
  valid at ALL primes including the bad ones (images off P^1 dropped);
  charpoly(T_2|H_1) = x^2 (x+2)^4 g4^2 h6^2 exact — the four particle
  blocks [2, 4, 8, 12] (143a1, the level-11 oldspace ghost, the f2
  quartic, the f3 sextic); star involution splitting (1,1) (2,2)
  (4,4) (6,6); harmonic presence densities vanishing identically on
  the skeleton's self-loop edge.
* eisenstein: congruence moduli C = gcd_p det((p+1)I - T_p|block) —
  143a1: 1 (no congruence), 11a1 ghost: 5^4 (modulus 5 = Mazur's
  numerator((11-1)/12), Sturm-certified at bound 28), f2: 7^2, f3:
  12^2; U_11 = U_13 = -1 on 143a1; U_11 = 1, U_13 charpoly
  (x^2-4x+13)^2 on the ghost; independent Weierstrass point-count
  route. Corpus reconciliation documented: FIELD_POLY_F2 and g4 are
  the same quartic field (disc 1957, Frobenius patterns agree at all
  p < 400); FIELD_POLY_F3 = h6(-x) exactly. A2_COMPLEX remains exposed
  but is retracted under CC-01 (cleanup item, see below).

Studies shipped: studies/promotion_2026aug/ — the two provenance
chains (graph_uncertainty -> ribbon_embedding -> hecke_particles ->
eisenstein_congruences; w2_susceptibility -> w3_cumulants ->
curvature_tano_manifold) plus the standalone triangular_layers
(honest negative) and m7_graph_channel (study-only, M6/M7 arc):
nine scripts with their JSON ledgers, 76 gates total, all green.

## Auditor notes (K3, E2 round before integration)

* Independent certification: T(beta) closed form re-derived
  algebraically (Moebius coprime decomposition, sharing no steps with
  the module) AND against a 4M-sieve moment route (gap 1e-17 at
  beta=4, decaying as N^{1-beta}); COLD constants to 30 digits; U(4)
  vs sieve 7.8e-18; curvature anchors confirmed by an independent
  finite-difference Brioschi route; charpoly factorization re-computed
  with sympy; elliptic-block a_p checked against the auditor's own
  point counting; the f2/f3 norm-moduli 7 and 12 corroborated E2 by
  the v0.13.0 Jacobian point-count table (gcd of |A_i(F_p)| over 93
  primes) — a route sharing nothing with Manin symbols.
* F1 (resolved, module correct): the curvature study's gate JSON
  records K_12346 = 1.3549364825; the correct limit is
  1.3549368866023 (module value and auditor scan agree; the gate
  record is stale).
* F2 (resolved in the second drop): eisenstein_congruences.py and
  triangular_layers.py arrived and landed; both re-ran ledgers
  identical to the uploaded records (7/7 gates each). The M7 study
  (m7_graph_channel.py + ledger + transcript) arrived in the same
  drop and landed as study-only: the auditor's independent exact
  route (sympy rational arithmetic on the certified model primitives,
  mpmath 80-digit cross-check) confirms every certified claim —
  [V,T_p] nonzero in all 676 entries for p = 2, 3, 5; the full
  16-entry transmission table (off-block norm 561.042783); M
  symmetric / V G-self-adjoint / rational; iota*-even commutation;
  per-dim quietness (old 57.3402 < ell 74.973 < q6 120.4706 <
  q4 120.5288); degree-potential null control V = 3I exactly with
  vanishing commutator. Gate P7b (graph-distance potential),
  uncertified at drop time, is now certified on two routes: author's
  engine and auditor's route agree to all printed digits (off-block
  norm 16.772862, 676/676 commutator entries, G-self-adjoint).
  Landed with one disclosed one-line import patch (mtftpkg ->
  mtft.hecke) and a regenerated 8-gate ledger; the pre-revision
  ledger is preserved verbatim in the audit report. The M6 tension
  stands sharpened per the study's own verdict: the obstruction is
  the REALITY of the arithmetic; M8 (twisted homology / Hodge J)
  is the registered successor.
* F3 (deferred to the author): x0_143.A2_COMPLEX / A2_COMPLEX_CONJ
  remain exported; INTEGRATION recommends setting them to None with a
  RETRACTED_CC01 pointer (append-only discipline) as a separate
  commit. Not done in this wave pending the author's decision.
* F5 (cosmetic, left as-is): moments.ser_inv contains a dead first
  loop superseded by its own "redo properly" pass; harmless, kept to
  preserve the author's file byte-pristine.
