# Unreleased — after v1.1.0

Everything after v1.0.0 shipped in v1.1.0 ([release notes](CHANGELOG_v110.md)). Since then:

## Wave 7 (2026-11-19): twists and knots

- `mtft.matter.twist_costs` (TC-01): the exact Hessian at the Pati–Salam pattern for any winding weight, the twist table of the spin
  structures of 143a1 (tree Yukawas, V = 0, M_A = 0 and the twelve's flatness go together on all nine pairs), the continuous twist's
  KK spectrum, and `hold_pull`, the pull of a one-loop energy on τ against the γ-flux's well. From records FF-12 and FF-13.
- `mtft.lines.genus` (KN-01): genus characters and the trace-4 law Ψ(A_p) = χ₋₃(p)Ψ(A), which is why the systole of X₀(143) winds
  round the divisor of (η(τ)/η(11τ))¹² and not round the other cusps' units. From record KN-1.
- Ledger entries `docs/tools/TC01_TWIST_COSTS.md` and `docs/tools/KN01_GENUS_SIGNS.md` (kept off the puzzle map, which waits for results); tests `test_twist_costs.py`, `test_lines_genus.py`.

## Wave 8 (2026-11-20): the twist space

- `mtft.lines.twists` (ZT-01): the twists ζ ∈ H¹(X₀(143), ℤ/2) in the records' 26-bit encoding; the involutions on twists; the mod-2
  cup form by route A (FF-17's C row for row, and the signed crossings mod 2 on 716² pairs); 143a1's pullbacks φ*(a, b); the transfer
  N(ζ) = φ*φ_!ζ and the flux charge q. From records FF-14, FF-15, FF-17 and FF-18.
- `mtft.lines.intersection` (IF-01): the integral intersection form on the 26 Hecke-block coordinates from signed crossings
  (block-diagonal; 143a1 block [[0, 1/4], [−1/4, 0]], so deg φ = 4); the homology lattice (index 2¹⁴3⁴, unimodular) and the flux-free
  lattice ker φ_* (rank 24, det 16, rank 22 mod 2: no handlebody carries the γ-flux). From record FF-19. Planned as
  `mtft.lines.forms`, a name already taken twice.
- `mtft.lines.walls` (WL-01): closed geodesics by class mod 2, the exact wall ball and the wall distance between twists, a second
  route by sums of at most three geodesics, the systole band, and the bottleneck closure of one-way costs. From records FF-14, FF-17
  and FF-20; it agrees with FF-20 class by class on all 120,827 classes within 16.
- Ledger entries `docs/tools/ZT01_TWISTS.md`, `docs/tools/IF01_INTERSECTION.md` and `docs/tools/WL01_WALLS.md` (off the puzzle map);
  tests `test_zt01_twists.py`, `test_if01_intersection.py`, `test_wl01_walls.py`.

## Wave 9 (2026-11-20): the twist landscape

- `mtft.lines.twist_energy` (TE-01): the relative heat trace D_ζ of a twist from the Selberg side (a lower bound from the
  geodesics up to length 14, and a central value with the tail); the per-class weights of the one-loop energy at the pattern
  (FF-15, the w = 0 winding included) and of the twelve's pins (FF-14); the energy and pins of one twist by direct sums, and of
  all 2²⁶ twists by the Walsh–Hadamard transform; the depth bands by the number of systoles twisted (FF-20); the twelve's
  curvature eigenvalues and the hold of τ from a D curve. The classes are checked against `mtft.maass.trace`'s frozen weights
  per trace. Not the glueball fit's `mtft.maass.glueballs.twist_energy` (GB-02).
- Ledger entry `docs/tools/TE01_TWIST_ENERGY.md` (off the puzzle map); tests `test_te01_twist_energy.py`.
