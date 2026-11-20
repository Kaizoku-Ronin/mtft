# ZT-01 — the twists of X₀(143): the records' encoding, the cup form, 143a1's pullbacks and the flux charge (2026-11-20)

Module `mtft.lines.twists`; tests `tests/test_zt01_twists.py`. Ported from the project records FF-14 (the encoding and the
pullbacks), FF-15 (the involutions on twists), FF-17 (the cup form by route A) and FF-18 (the transfer and the flux charge),
18–20 November 2026, filed under `claude/experiments/` in the MTFT Python project.

A twist is MATTER-04's ζ: a 2-torsion line bundle on X₀(143), one sign for each class of closed curve,
ζ ∈ H¹(X₀(143), ℤ/2) = J₀(143)[2] = (ℤ/2)²⁶. Fermions see it; bosons do not. The module writes a twist as the records' 26-bit
integer. Bit i is the coefficient of `basis()[i]`, FF-14's basis of the 29-bit characters of the punctured curve that are even on
every cusp loop, found from the face cycles of `mtft.lines.darts` by elimination.

What it holds:
- the encoding: `basis`, `to_twist`, `to_character`, `class_bits` (the 26 bits Y of a class, with ζ(γ) = parity(ζ & Y)), `value`;
- W₁₁, W₁₃, W₁₄₃ and the mirror on twists (`act`), from the class matrices of `mtft.lines.involutions`;
- the cup form by route A (`cup`, `cup_form`, `form`). Each cusp face is coned from a centre and each edge halved; ordering
  centres < midpoints < vertices and applying Alexander–Whitney gives (α ∪ β)[X̄] = Σ_faces Σ_{j<i} α(e_j)β(e_i) + Σ_e α(e)β(e);
- Poincaré duality mod 2: a twist δ is dual to the class Cδ (`dual_class`), classes meet as (C⁻¹Y₁)·Y₂ (`class_form`,
  `intersection`), and a wall along a curve of class Y changes ζ by C⁻¹Y (`wall_twist`);
- 143a1's twists on the curve, φ*(a, b)(γ) = a n₁ + b n₂ mod 2, with (n₁, n₂) the class's 143a1 coordinates from
  `mtft.lines.periods` (`ell_matrix`, `pullback`, `W`);
- the transfer N(ζ) = ζ + W₁₁ζ + W₁₃ζ + W₁₄₃ζ = φ*(q₂, q₁) and the flux charge q(ζ) = (ζ ∪ φ*(1,0), ζ ∪ φ*(0,1))
  (`transfer`, `flux_charge`, `check_transfer`, `mirror_on_flux_classes`).

Status: EXACT (finite linear algebra over F₂). Gates:
- the records' integers: φ*(1,0) = 6665056 (φ*η), φ*(0,1) = 14025081, φ*(1,1) = 11776537. W₁₃ sends the minimum 45796554 to
  48201642, and W₁₁ fixes it.
- The cup form equals FF-17's C row for row. It is alternating of rank 26, invariant under W₁₁, W₁₃, W₁₄₃ and the mirror, and
  independent of where each face word starts.
- Second route: route A agrees mod 2 with the signed crossings of all 716² ordered pairs of geodesics up to level-1 trace 24
  (`mtft.lines.intersection.cup_against_crossings`). 211,894 pairs cross an odd number of times, and none disagree.
- The transfer identity holds on all 26 basis vectors; both sides are linear, so it holds on all 2²⁶ twists. W is isotropic.
  FF-18's flux charges of the named twists and of its three orbit starts are reproduced. The mirror fixes the flux class (0,1) and
  swaps (1,0) with (1,1), as FF-18 found.

Not here: the twist energy (FF-14/15's heat trace over all 2²⁶ twists), and which twist a vacuum takes. The encoding works at any
level `darts` handles. The involutions, the pullbacks, the transfer and the flux charge are for N = 143 only.
