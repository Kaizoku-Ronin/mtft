# TE-01 — the twist landscape: the relative heat trace, and the energy and pin of every twist (2026-11-20)

Module `mtft.lines.twist_energy`; tests `tests/test_te01_twist_energy.py`. Ported from the project records FF-14 (the relative heat
trace and the pin) and FF-15 (the energy and the transform over all twists), with FF-20's depth bands, 19–20 November 2026, filed
under `claude/experiments/` in the MTFT Python project. This is not the glueball fit's "constant twist energy"
(`mtft.maass.glueballs.twist_energy`, GB-02), which is a different quantity.

A twist ζ ∈ H¹(X₀(143), ℤ/2) (`mtft.lines.twists`) lifts one fermion zero mode of every adjoint component on the curve
(MATTER-04). On the Selberg side, the bosons' heat trace minus the twisted fermions' is a sum over the primitive oriented closed
geodesics that ζ makes antiperiodic:

D_ζ(t̃) = Σ_γ₀ Σ_{k odd} 2[ζ(γ₀) = −1] ℓ₀ e^{−t̃/4} e^{−(kℓ₀)²/4t̃} / (√(4πt̃) · 2 sinh(kℓ₀/2)).

Every term is ≥ 0, so a truncation is a lower bound. The central value adds a tail for the geodesics beyond L_max, with ζ = −1 on
half of them: ½ erfc((L_max − t̃)/2√t̃). Through VE-01's torus sums over the windings of the ledger's torus, D_ζ gives two things:
- the one-loop energy at the pattern, V(ζ) = −32 (4π/1536 r R⁶) ∫ t̃⁻⁴ (S₀ + 1) D_ζ dt̃ (FF-15; the "+ 1" is the w = 0 winding);
- the twelve's curvature on each cycle, pin_a(ζ) = 64π² (4π/1536 r R⁶) ∫ t̃⁻⁴ S_a D_ζ dt̃ (FF-14, in FF-12's convention).

Both are linear in which classes are twisted, so each class carries a weight. The value of every twist at once is a
Walsh–Hadamard transform over F₂²⁶.

What it holds:
- `classes`: the 93,720 geodesics up to length 14 with their class bits; `per_trace_check` is the second route (below).
- `heat_terms`, `tail`, `relative_heat_trace`, `mean_relative_heat_trace` and `J`: the heat trace, as a lower bound and a central
  value.
- `lattice_moments`, `switch_check`, `energy_weights` and `pin_weights`: the torus sums, and the per-class weights with their tails.
- `energy` and `pin`: one twist, by direct sums. `fwht`, `landscape` and `summary`: every twist, by the transform (2²⁶ doubles,
  512 MB).
- `bands`: the depth V/V_min by the number of systoles twisted.
- `pin_matrix` and `hold`, from a D curve: the twelve's curvature eigenvalues, and the hold of τ (V, its gradient in τ and sinh d).

Status: EXACT within VE-01's product approximation; the central values are ESTIMATES (the tail model). Gates:
- **The classes.** 93,720 classes up to length 14, of which 248 are never twisted (the shortest has length 9.897), as FF-14
  found.
- **Second route to the classes.** Per trace, Σ ℓ₀ over the classes and their powers equals `mtft.maass.trace`'s frozen weights
  at level 143 on all 321 traces with classes, to 3.1 × 10⁻¹⁴. 773 traces are empty in both routes.
- **The records' numbers**, to 10⁻⁹:
  - FF-15's energies of its eleven named twists at r = 230 and 280;
  - FF-14's J for its weakest, hardest and mean twists;
  - FF-14's pins, the twelve's eigenvalues and the hold for its hardest twist;
  - FF-14's printed V for its weakest twist.
- **Two routes inside the tool.** The per-class weights summed over a twist, plus the tail, equal the integral over its D curve,
  for the pins and the energy.
- **Outside the suite.** The landscape over all 2²⁶ twists reproduces FF-15's energy minimum (45796554), mean and maximum at both
  radii, FF-14's pin lower bounds, and FF-20's five bands (wave 9's record, `WAVE9_Plan.md`).
- **Ties.** The energy and the pins depend only on the lengths of the twisted geodesics, so a twist and its images under W₁₁, W₁₃,
  W₁₄₃ and the mirror have equal values. The extremes therefore come in pairs: the minimum is {45796554, 48201642} and the maximum
  {394080, 6533120 = 394080 ⊕ φ*η}. `summary` lists the tied twists instead of letting the order of summation pick one.

Conventions:
- The energy includes the w = 0 winding, as FF-15 did.
- FF-14's hold printed V without it, because `twist_costs.windings` has no (0, 0). `hold` returns V both ways; the gradient and
  sinh d do not depend on it.
- The t̃ grid stops at 600. FF-14 found that the omitted part only raises the values.

Not here: the twelve's masses. They need FF-13's Wilson-line masses, which are not in the package.
