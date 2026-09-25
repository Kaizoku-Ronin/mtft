# mtft v0.26.0 — first QFT on X0(N), frozen certified data, dynamics with controls (2026-09-03)

## `mtft.surface.ising` — the Ising model on the dual Manin graph (EXACT)
Spins on the Manin triangles, couplings on the edges, cusps as faces.  Fisher
triangle decoration + the Cimasoni-Reshetikhin spin-structure Pfaffian sum
(4^g Kasteleyn orientations from a mod-2 face system and H^1(Z/2) twists;
q_eps(C) = n_K(C~)+1 on lifted dual cycles; Arf on a mod-2 symplectic basis
built from the exact intersection form; reference-matching normalisation).
Two-route gate against brute-force enumeration of all 2^F configurations:
N = 6, 11, 15, 35, 55 (genus 0-5) agree to <= 2e-15.  N = 77 (4^7 terms) runs
in 28 s.  N = 143: all structural gates pass; the 4^13 sum is a checkpointed
CPU-day (`full_sum_job`); `sample` gives per-spin-structure Pfaffians
(r = 0 baseline: even/odd |Pf| ratio 1.001, i.e. no parity separation yet —
the continuum prediction that odd-parity Pfaffians vanish is the refinement
target, DIAGNOSTIC).  Critical coupling of the honeycomb refinement family
tanh(beta_c) = 1/sqrt(3).

## `mtft.surface.frozen` — GP-free verification of the certified layers
N = 143 cycle basis, intersection form, T2, T3, W11, W13, W143, period matrix
and J_true frozen under `surface/_data/` with the generating GP script's
SHA-1.  Twelve gates re-verified at call time (integrality, involutions,
commutation, self-adjointness, Riemann relations, J commutators, and that the
live tree/cotree basis still equals the frozen one).

## `mtft.surface.dynamics` — Waves 11-13 core with mandatory genericity controls
Hamiltonian flows and Lie closure (after Sol).  Audit finding recorded: full
closure 351 = dim sp(26,R) is the generic outcome (random Hamiltonians: 351);
the gate discriminates structure (Hecke-commuting: 4; block-diagonal: 127 =
dim sp(2)+sp(4)+sp(8)+sp(12)).  Closure residual gates are frame-sensitive at
double precision; `Stage.orthonormal()` conjugates by G^{1/2} so residuals
measure Hodge-metric size (raw cycle frame: AMBIGUOUS; orthonormal: separations
>= 1e4).  `off_block_fraction` is the informative quantity.  Falsifier F8:
a closure dimension that random generators also reach is not evidence.

## Packaging
sympy already declared (verified); `MANIFEST.in` includes `surface/_data`.
Tests: `tests/test_surface_qft.py` (7 tests, no GP needed).

## Added before push (2026-09-04)

### `mtft.surface.bimodule` — doubled-space census
Real spectral-triple instrument on H = V ⊕ V: left action diag(a, W a W⁻¹),
real structure = swap, D_M off-diagonal.  Measures order-zero, the first-order
solution space and the one-forms it generates — all with absolute scales
(a relative SVD threshold counts numerical dust; that mistake is recorded).
Exact result at N = 143 with A = R[T2, T3, U13]: untwisted and W11-twisted
doublings fail order-zero (U13 is not normal, ‖[U13,U13ᵀ]‖/‖U13‖² = 0.89);
W13 and W143 twists restore it to 7e-15 through U13* = W13 U13 W13; in every
configuration the one-forms vanish (≤ 2e-14): involutive twists give
[D, a] ∝ λ_{σ²(i)} − λ_i = 0.  Sectors (2,12,10,2) on V and (252,88,88,248) on
V⊗V are Hom-space dimensions.  No non-involutive twist exists over Q:
Aut(K4) = Aut(K6) = {1} (Galois closures S4, S6).  `random_alphabet` is the
F8 control.

### Frozen data extended (`surface/_data/x0143_certified.json`)
U11, U13 transported (integral, route-B checked; U11: (x−1)¹²(x+1)¹⁴,
U13: (x−1)¹²(x+1)¹⁰(x²−4x+13)²) and the saturated Hecke-block lattices
H_1(Z) ∩ V_block (PARI matkerint).  New gates (19 total, GP-free): U commute
with good Hecke and with J; AL adjoint identities; blocks annihilated and full
rank; covolume² = polarization determinant.

### EXACT block invariants (`frozen.block_invariants`)
Elementary divisors of the intersection form on the saturated Hecke blocks:
ell (4,4); ghost (2,2,18,18); q4 (2⁶,18,18); q6 (1,1,2⁸,4,4).  The Hodge
covolumes are the square roots of the polarization determinants (4, 36, 144,
64) because det J|_B = 1 on a J-stable sublattice — so they carry no
transcendental Petersson information (correction to the AF-gate audit note);
the transcendental block content is the block period shape (e.g. j(143a1)).
Tests: `tests/test_surface_bimodule.py` (3 tests).
