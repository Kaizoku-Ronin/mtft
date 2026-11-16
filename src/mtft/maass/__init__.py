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
  trace       CERTIFIED  the Selberg trace formula by sector on frozen class-number weights
                         (lengths to 17); the positivity bound: no form missing below R = 1.95
  ford        CERTIFIED  the Ford domain of Γ₀*(143) with quadrature: norms, overlaps,
                         participation ratios
  glueballs   PHENO      H24 against the lattice: the dictionary, every registered forecast with
                         its outcome, the fits, and two sealed forecasts that are still open

Provenance: the Modular Plane lab, records MP14 to MP20 (12 to 14 November 2026).  The lab's own
scripts and number files are preserved under ``studies/maass_2026nov/``.  Nothing here is a
statement about physics except ``glueballs``, which is a comparison under a stated dictionary
and records that the hypothesis it tests fails.
"""
from __future__ import annotations

from . import domain, ford, glueballs, hejhal, kbessel, scattering, spectrum, trace

__all__ = ["domain", "ford", "glueballs", "hejhal", "kbessel", "scattering", "spectrum", "trace"]
