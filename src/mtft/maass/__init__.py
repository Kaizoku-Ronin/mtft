"""mtft.maass — the Laplace spectrum of X₀(143): Maass cusp forms by Hejhal's method (v1.1.0).

Layers, with their claim classes:

  kbessel     CERTIFIED  K_{iR}(x): a fast trapezoid rule (R ≤ 3) checked against mpmath
  domain      EXACT      the reduction into the Ford domain of Γ₀*(N), Atkin–Lehner characters
                         and the η-twist's sign
  hejhal      CERTIFIED  the solver: normalised system, two-horocycle refinement, pole-free detector
  spectrum    CERTIFIED  the frozen table of Γ₀(143): 41 untwisted forms below R = 2.0 and 42
                         η-twisted forms below R = 1.8, with coefficient vectors
  scattering  EXACT      the scattering matrix of the four cusps, Φ(s) = φ(s) M₁₁(s) ⊗ M₁₃(s),
                         Φ(½) = −I, and the Eisenstein series
  ford        CERTIFIED  the Ford domain of Γ₀*(143) with quadrature: norms, overlaps,
                         participation ratios

Provenance: the Modular Plane lab, records MP14 to MP20 (12 to 14 November 2026).  The lab's own
scripts and number files are preserved under ``studies/maass_2026nov/``.  Nothing here is a
statement about physics; the lattice comparison lives apart, with its forecasts and misses.
"""
from __future__ import annotations

from . import domain, ford, hejhal, kbessel, scattering, spectrum

__all__ = ["domain", "ford", "hejhal", "kbessel", "scattering", "spectrum"]
