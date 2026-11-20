# IF-01 — the intersection form of X₀(143) on its Hecke blocks, deg φ = 4 and the flux-free lattice (2026-11-20)

Module `mtft.lines.intersection`; tests `tests/test_if01_intersection.py`. Ported from the project record FF-19 (no bubble of
nothing), 19 November 2026, filed under `claude/experiments/` in the MTFT Python project. The plan named the module
`mtft.lines.forms`, but that name is already used twice (`mtft.forms`, `mtft.periods.forms`).

`mtft.lines.periods` gives every class 26 integer coordinates on the Hecke blocks: 143a1, the level-11 oldspace, f2 and f3
(2 + 4 + 8 + 12). `mtft.lines.crossings` counts, with signs, how two closed geodesics cross. The intersection form in block
coordinates, ι(γ, δ) = b(γ)ᵀ J b(δ), is solved from 26 geodesics with independent coordinates and checked on every pair:
- J is block-diagonal: the Hecke blocks are orthogonal.
- On the 143a1 block J = [[0, 1/4], [−1/4, 0]], so deg φ = |(J⁻¹)₁₂| = 4 (`modular_degree`).
- On Λ = H₁(X₀(143), ℤ) the form is integral and unimodular (`homology_lattice`). Λ has index 1,327,104 = 2¹⁴3⁴ in ℤ²⁶.
- The flux-free classes ker φ_* form a rank-24 lattice (`flux_free_lattice`). Its form has determinant 16 = (deg φ)² and rank 22
  mod 2, and its isotropic sublattices have rank at most 12.
- Consequence: a handlebody bounded by the curve kills 13 independent meridians, which would all have to be flux-free. Since 12 < 13,
  no handlebody carries the γ-flux (FF-19's seal). Mod 2 the obstruction is invisible, because deg φ is even.

Status: EXACT. The form reproduces all 512,656 signed crossings of the 716 geodesics up to level-1 trace 24. FF-19 used 1,426
geodesics and 2,033,476 pairs; trace 16 already spans the 26 coordinates. 143a1's twists vanish on a basis of the flux-free classes
of ℤ²⁹ (`pullbacks_on_flux_free`; rank 27 with the cusp loops). The cup form of `mtft.lines.twists` agrees with the crossings mod 2
on every pair (`cup_against_crossings`).

What this is not: FF-19's Stokes argument that there is no bubble of nothing stays in the record; the module supplies its lattice
facts. Timing: about 17 s, nearly all of it `mtft.lines.crossings` on the 716 geodesics.
