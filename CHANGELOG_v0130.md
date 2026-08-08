# v0.13.0 — the ancestry wave (2026-08-08)

New module: combinatorial (the 2024 lineage as certified tools).
Tests: tests/test_combinatorial.py (73, all fast, ~1 s).

* Figurate/Faulhaber engine, exact Fractions: S_p = C(p) T^q - R_p with
  C(p) = 2^q/(p+1); jump-rule recurrence re-certified against the closed
  form; R(1) = C(p)-1 anchored. Two routes share no steps (direct
  summation vs Bernoulli B+ closed form).
* The sigma involution n -> -1-n: Q[n]^sigma = Q[T]; odd sector
  (2n+1)Q[T] with (2n+1)^2 = 8T+1. Odd-p power sums sigma-even, even-p
  sigma-odd — the figurate principle is a parity selection rule
  (structural companion to the du03 eta-parity, recorded as structural
  only). Parity certified by two disjoint routes (downward recurrence at
  negative integers vs coefficient substitution).
* s-gonal absorption defect closed form: [A_s P_s^q - S_p]_{n^p}
  = -(s-3)/(s-2), p-independent, zero iff s = 3; triangular is the
  unique sigma-invariant s-gonal family. The 2024 honest negative now
  carries an exact certificate.
* Graph uncertainty: [D,L]_{uv} = A_{uv}(d_v - d_u) two-route exact;
  zero iff regular; Robertson certified on complex states (real states
  make <[X,Y]> vanish identically for real-symmetric pairs — noted).
* Number-phase: Maassen-Uffink H_K + H_Theta >= log(n+1) certified;
  saturating basis states hit the bound to 1e-9. ARCHIVED-2024
  REGRESSION RESOLVED: the index-units pass missed the archived
  (1.3712, 1.5331, 0.8869) by exactly 2*pi/6, pinning the 2024
  convention as angular (spectrum 2*pi*j/d); under it all four archived
  numbers reproduce (Cert 5e-4). Discovery route preserved in
  number_phase_regression().
* q-thermodynamics: Gaussian binomials by q-Pascal vs product-formula
  exact division; Galois anchors q=2 (1,2,5,16,67,374); Gibbs chain
  logZ -> U, Var, S, C, F, Fisher each certified by two routes
  (ensemble moments vs finite differences); multiplicity-as-energy
  Z(0) = n+1, Z(1) = 2^n / Galois numbers exact.
* THE BRIDGE: <w>_beta = -zeta'(beta+1) on the primon ensemble (Pr,
  one line from sum w_n n^{-s} = -zeta(s) zeta'(s+1); identity itself
  certified sieve-vs-mpmath with explicit tail majorant, gap 4.7e-6
  <= bound 8.9e-4 at s=2, 1.3e-10 <= 1.9e-8 at s=3). Endpoint
  beta -> 1+: <w> -> -zeta'(2) = cold-gas alpha = 2*T_INF (the stored
  constant is -zeta'(2)/2; agreement 8e-11, the stored precision).
  Primon Fisher metric = log-zeta curvature, gap 2.1e-8. Sieve's
  independent route: mtft.weight_array, agreement 1e-12 over n <= 300.

## Addendum — crypto wave landed (2026-08-08, K3 staging)

The last INTEGRATION_v0120 open item is closed: mtft.crypto is now a
package (tier-5f module precedent — import mtft.crypto, not flattened).

* mtft.crypto.jacobian_order: |A_i(F_q)| = N_{K_i/Q}(q+1-a_q) for the
  three newform orbits (dims 1, 4, 6) on J0(143); 93-prime table at
  src/mtft/crypto/_data/jacobian_orders_N143.csv (279 records, all
  good primes p <= 499); PARI/GP companion _compute_orders.gp.
* K3 independent certification (E2, two routes): A_1 all 93/93 primes
  exact vs direct point counting on 143a1 (y^2+y = x^3-x^2-x-2);
  orbit field polys y^4-4y^2-y+1 (disc 1957) and
  y^6-10y^4-2y^3+24y^2+7y-12 (disc 194616205) are the same fields as
  the repo's T2 charpoly factors (irreducibility checked); A_2/A_3
  20/20 spot checks |A_i(F_q)| = F_{i,q}(q+1) from the engine's
  integer Hecke charpoly factors; embedding degrees with minimality.
* K3 gate tests/test_jacobian_order.py (6, all fast): table integrity
  (Hasse-Weil, security ratio), the three exact routes above, and
  Atkin-Lehner eigenvalues {1:(+1,+1), 2:(-1,+1), 3:(+1,-1)}.

## Auditor disclosures (K3)

* BU-2 (staged, disclosed): tests/test_jacobian_order.py is an
  auditor-written gate, not source-session material; it certifies the
  CSV rather than reproducing its generation.
* Provenance flag: the GP companion header lists the 28-prime pilot;
  the CSV is the 93-prime production run. The CSV is authoritative.
* Gate arithmetic: 484 (v0.12.0) + 73 (combinatorial) + 6 (crypto)
  = 563 tests, 1 skipped (opt-in slow tier unchanged).
* BU-3 (staged, disclosed): the new mtft.crypto package shadowed the
  pre-existing flat crypto.py module (package wins over module),
  breaking 6 legacy tests (ArithmeticHash / BurningShipPRNG /
  ModularKeyExchange / ArithmeticLattice imports). Resolution, by the
  auditor: crypto.py merged INTO the package as its __init__.py — every
  legacy public name re-exported unchanged, jacobian_order added as a
  submodule. Full suite re-gated green after the merge (563/1).
