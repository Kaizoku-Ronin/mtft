# studies/crypto_2026aug — the "Modular curves and encryption" session (2026-08-07)

What the v0.11.x X0(143) structures say about cryptography, in both
directions. All three scripts run on the repo's own Manin-symbol engine
(level 143, weight 2, 29-dim) with studies/ on PYTHONPATH.

Offense — Hecke-CDH (eprint 2025/1681) does not stand:
- hecke_ke_check: the eigenbasis attack the paper dismisses in Sec 3.3,
  run honestly. Eve recovers the shared secret to 2.6e-15 with one
  eigendecomposition and d divisions; NP-hardness of factoring the N_i is
  irrelevant to CDH.
- hecke_ke_exact: the stronger break. No eigendecomposition, no floats, no
  semisimplicity assumption — commutativity alone. Eve solves X in T with
  X(f) = H_A(f) over Q exactly (one gauss_jordan_solve); then X(H_B f) is
  the shared secret, exactly.

Defense — what the curve offers:
- crypto_probe: Q1 every isotypic block has multiplicity 2 (rank-2 Hecke
  module, the Kyber shape); Q2 the integral Hecke order is GLUED —
  disc ratio 7234523136 = 2^12 * 3^2 * 443^2, congruence primes {2,3,443}
  (resultants Res(F1,F3) = -12, Res(F2,F3) = 7088 = 2^4*443); Q3 the Hecke
  correspondence graph on P^1(Z/143) is Ramanujan at every good prime
  p <= 43 measured.

Not landed here (referenced by INTEGRATION_v0120, still open):
mtft.crypto jacobian_order module + 93-prime CSV; periods.gp variants.
