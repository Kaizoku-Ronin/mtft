# CC-33 — the self-dual tensor's p2 coefficient in R2C-03 (2026-09-23; EXACT; applied in v0.33.0)

## The error (found while writing the proof of compendium theorem IV.3)
R2C-03 quoted the p2 coefficient of a six-dimensional self-dual tensor as +1/360.  That number is the p1^2 coefficient of the tensor's
index density -L8/8, with L8 = (7 p2 - p1^2)/45; the p2 coefficient is -7/360.  In units of the complex Weyl fermion (A-roof_8 =
(7 p1^2 - 4 p2)/5760, p2 coefficient -1/1440) the three fields count 1 : 28 : 245 (Weyl : tensor : gravitino, the gravitino from
[A-roof (tr_V e^R - 1)]_8 = (275 p1^2 - 980 p2)/5760).

## Three independent witnesses (`gravitational_anomaly.tensor_coefficient_witnesses`)
1. The Chern-root genera above, expanded directly.
2. The (1,0) supergravity count H - V + 29 T = 273: the gravity multiplet contributes 245 (gravitino) + 28 (self-dual tensor) = 273 and
   a tensor multiplet 28 + 1 = 29.  With 4 units a tensor multiplet would count 5, not 29.
3. The (2,0) tensor multiplet (one self-dual tensor, two Weyl): p2 magnitude 1/48 = (28 + 2)/1440.

## What changes
- Cancelling the p2 term of a fermion spectrum with n_grav signed Weyl dimensions needs n_grav/28 net anti-self-dual tensors: an integer
  only if n_grav = 0 mod 28 (was: n_grav/4, n_grav = 0 mod 4).
- None of M1's four ledger-preserving assignments (n_grav = 24, 22, 18, 16) is completed by chiral tensors, with or without a gravitino
  (269/28).  R2C-03's survivors (+,+) -> 6 tensors and (-,-) -> 4 tensors are withdrawn.
- R2C-04's pure-gaugino rule becomes dim(adj) = 0 mod 28: E6 (78), E7 (133) and E8 (248) all fail; the E8 survivor is withdrawn.
- Gate `six_d_anomaly_lift` now carries the full table as witness (no survivors); its verdict (FAIL) is unchanged, since the colour-cubic
  obstruction already closed the gate.
- Unchanged: the Weyl and gravitino coefficients, "no gravitino with tensors alone", H - V + 29 T = 273, kappa = 1, and every statement of
  R2C-03 that does not use the tensor coefficient.

## Package (v0.33.0)
`gravitational_anomaly.P2_SELF_DUAL_TENSOR = -7/360` (with `P1SQ_SELF_DUAL_TENSOR = 1/360` recorded), `tensor_integrality` (rule field),
`m1_tensor_survivors` (empty), `unified_parent.gaugino_p2_integrality` (/28), `pipeline.six_d_anomaly_lift` (full-table witness).
Tests re-derive the constant from the genera rather than re-asserting it: `tests/test_route2_r2c03.py::test_tensor_coefficient_cc33`,
`tests/test_v0330.py::test_cc33_tensor_units_and_rule`; `tests/test_route2_r2c04.py` updated for E8.  Legend: `cc33_tensor_coefficient`.
Provenance: theorem compendium Part D (IV.3 and App. P, P-D1/P-D2), written 2026-09-23.
