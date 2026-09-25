# mtft v0.27.3 — bimodule representation gates (2026-09-10)

From Astra's rational-descent release audit: `Doubling.left` used W a Wᵀ, correct for orthogonal
Atkin–Lehner twists but non-unital for a general twist (probe: W = 2I, alphabet {I}, D = 0 passed
every gate).  Now W a W⁻¹, plus `representation_gates` (π unital, multiplicative on alphabet
products, opposite action anti-multiplicative, twist orthogonal in the working frame) folded into
`axiom_gates`; the probe now fails on `representation`.  `differential_rank` docstring states that
it ranks generator commutators, not the generated algebra.  All AF-09 and cyclic-control results
unchanged (AL twists are orthogonal).
Base: live 0.27.1 tarball + the 0.27.2 periods-frame changes; rebase onto the published 0.27.2
tree if it diverged.  Tests: `tests/test_surface_bimodule_probes.py` (2).  Release pin four-way.
