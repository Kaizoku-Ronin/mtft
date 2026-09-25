# mtft v0.19.0 — the integral-model arc becomes package machinery

Everything certified in the 2026-08-24 arc (certificates v1-v9: small-prime
counts, saturation, codifferent theorem, CI-A factorization, product-lattice
chain, Atkin-Lehner splitting) now ships as modules, frozen data, and
call-time gates.  Claude built; Sol's gates adopted; Kimi audit pending.

## New modules

* **mtft.integral_lattice** — exact integer-lattice toolkit: mod-p
  rank/kernel, p-saturation with unimodular pivots, multi-prime
  `saturate`, Hermite form, Smith invariants (HNF ping-pong — see
  Fixes), `quotient_invariants`, HNF membership and `class_order`,
  `rational_kernel`, and `operator_matrix` via the exact Gram route.
* **mtft.canonical.integral** — integral models of X0(143):
  `saturated_qexpansions` (step ledger {2:25, 3:8, 5,7,13,19,103,
  5560463: 1}), `count_points_modp(p, model)` exposing the instructive
  model triple at p = 2 (packaged s2 -> 7, adapted+saturated ideal -> 3,
  fully saturated -> 4 = the four cusps), `points_modp`,
  `cusp_reductions` (exact normal-equations route), `ci_a_codifferent`
  (a = -637 recomputed live), `quadratic_saturation_obstruction` (the
  Q(L) invariant), `al_splitting` ((Z/2)^6 (+) Z/52 and the sequential
  2^6 * 2 * 26 filtration), `al_denominator` (W11 -> 1, W13 -> 13).
* **mtft.canonical.integral_gates** — certificates v1-v9 as live gates;
  nothing asserted from memory.
* **mtft.codifferent** — Canonical Codifferent Theorem at instance:
  frozen exact eigenvalue tables for both orbits, gamma tables,
  pure-Python Newton-trace verification of every packaged coefficient
  (n <= 140), orbit indices 576 * 1957 and 2304 * 194616205.
* **mtft.exception_spectrum** — the Exception-Spacing Curvature Law:
  exact cumulant/Brioschi curvature, the two-exception constant and
  phase trichotomy in rho = 2 e1/e2, line-atom and flat-base variation
  lemmas (lambda-dressed), the signed defect spectrum with exact
  e_j = e1 * m cancellation, and the RELATIVE truncation rule.
* **mtft.quadratic_forms** — the Gauss-Legendre three-squares seed
  (L1-L7) promoted to a module with its E2 gate.

## New data (src/mtft/canonical/_data/)

* `X0_143_I2_adapted_saturated.txt` — 55 quadrics in the ADAPTED frame,
  integrally saturated, with a mixed-model warning in the header.
* `X0_143_f2_eigen_an.txt`, `X0_143_f3_eigen_an.txt` — exact power-basis
  a_n coordinates, n <= 140 (PARI mfeigenbasis, monogenic via a_2).

## Frame clarification (the v0.18.0 pitfall, certificate v6)

`ideal_basis()` is and was s2-frame — its data-file header always said
so; the pitfall was API discoverability next to COORDINATE_LABELS.
v0.19.0 adds `IDEAL_BASIS_FRAME = "s2"`, `ideal_basis_adapted()`, a
frame note in the module docstring, and `gate_frame()` asserting the
residual behavior of both files.  Lesson filed: read data headers
before interpreting arrays.

## Fixes and lessons during the build (append-only honesty)

* Two module-name collisions with existing exports (`mtft.lattice` is
  the gauge module; `spectrum` is a riemann-toolkit function) — new
  modules named `integral_lattice` and `exception_spectrum`.
* Newton-identities k-range bug in the trace table (caught by the
  codifferent gate returning False; fixed and re-verified).
* `cusp_reductions` solved coordinates in an hnf-permuted basis and
  then with a transposed change of basis — both caught against the
  certified route; final implementation mirrors it exactly (Gram
  normal equations).
* Classical two-sided Smith elimination exploded on 13 x 13 quotient
  matrices; replaced by HNF ping-pong with a divisibility repair.
* `operator_matrix` on the solve_in_lattice path returned a wrong
  denominator for W11; rewritten on the Gram route and re-verified
  against the PARI-certified values (1 and 13).  The HNF membership
  path (`class_order`) is unchanged and self-asserting in gates.
* `Q(L)` is computed as HNF -> saturate -> small-matrix Smith, never
  Smith of the raw 141-row product matrix (165 s+ -> ~5 s).

## Tests

13 fast + 5 slow new tests (`test_integral_v0190`, `test_codifferent_
v0190`, `test_spectrum_v0190`, `test_quadratic_forms_v0190`); slow tier
~33 s total.  Existing suite re-run green (piecewise, zero failures).
