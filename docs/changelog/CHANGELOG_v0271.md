# mtft v0.27.1 — bimodule gate set completed (2026-09-09)

From Astra's weighted-Dirac audit: `surface.bimodule` checked order-zero, first-order and
one-forms but NOT the reality condition J D = D J, self-adjointness, or the grading relations,
so it could pass a non-real model.  `RealTriple` now carries the full finite real-spectral-triple
gate set (self-adjoint, reality, D odd, algebra even, JΓ = −ΓJ, order-zero, first-order, one-form
size) plus a Poincaré-pairing report (Q_ij = Tr Γ π(e_i) e_j°, unit-in-radical flag) and the
differential rank.  `symmetric_dirac_block(W, h) = W h + h W⁻¹` is the Dirac rule that keeps all
gates for real self-adjoint h (my earlier M = W h fails reality; Astra's correction verified).
Exact two-state cyclic control (`two_state_cyclic_control`) reproduces the addendum's table.
Helpers: `cyclic_permutation`, `kron_power`, `tensor_alphabet` (dense, size-guarded).  Recorded
correction to the three-factor audit note: rank preservation needs invertible h (T2 has a zero
character); K_h spectra mix weights (λ³ − (x²+y²+z²)λ − 2xyz on a 3-cycle).
Tests: `tests/test_surface_bimodule_axioms.py` (3).  Release pin bumped four-way.
