r"""mtft.lines — the closed lines of Γ₀(143)\PSL(2, ℝ): their thirty integers, their crossings and linking numbers (v1.1.0).

The unit tangent bundle of X₀(143) is a 3-fold, M = Γ₀(143)\PSL(2, ℝ), and the geodesic flow on it
is the frozen frame of the lab's instrument 12.7.  Its closed orbits are the *closed lines*.  At
level 1 a closed line has one integer, Rademacher's Ψ, its linking number with the trefoil (Ghys).
At level 143 the first homology of M has rank 30, and a closed line has thirty.

Layers, with their claim classes:

  darts        EXACT      P¹(ℤ/N), the graph on its darts, the class of an element in ℤ^b₁
  rademacher   EXACT      Dedekind sums, Ψ, the windings of η-quotients round the cusps
  closed       EXACT      the closed lines over a motion: ℓ, γ = g A^ℓ g⁻¹, (Ψ; cls); the transfer
  periods      EXACT      the 26 periods by Hecke block (2 + 4 + 8 + 12) by two routes;
               CERTIFIED  the period lattice of 143a1 (1e-8), also by two routes
  involutions  EXACT      the Atkin–Lehner involutions and the mirror on ℤ²⁹; a recorded miss
  cutting      EXACT      a closed geodesic cut along the fundamental domain, in rational arithmetic
  crossings    EXACT      every crossing of two closed geodesics on the curve; Lk(K_γ, K_δ) = −I(γ, δ)
  lorenz       EXACT      Ghys's Lorenz template: words and linking numbers at level 1; the transfer
  table        frozen     the lab's record, with SHA-1 provenance, for the gates to compare with

Provenance: the Modular Plane lab, records MP12 to MP14 (12 November 2026); the lab's scripts and
number files are preserved under ``studies/lines_2026nov/``.  What this is and is not: mathematics
of the curve's unit tangent bundle.  Nothing here is a statement about gravity; whether any of the
thirty integers is a charge of a frozen field is not known.
"""
from __future__ import annotations

from . import closed, crossings, cutting, darts, involutions, lorenz, periods, rademacher, table

__all__ = ["closed", "crossings", "cutting", "darts", "involutions", "lorenz", "periods", "rademacher", "table"]
