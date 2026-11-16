r"""mtft.plane — the Modular Plane lab's kernel: one matrix read eight ways, and the lattice of quantities (v1.1.0).

  kernel      EXACT   the eight readings of a real 2 × 2 matrix as a motion or as a generator
                      (``read``); its numerical range, spectral norm and Jacobi monodromy; exp(tL)
  quantities  EXACT   the lattice ℤ⁴ of dimensions, the natural and Planck sublattices with their
                      Smith forms, the Planck decomposition and the prospector (``site``)

The lab's page and its matrix glossary are teaching material and stay out of the package; this is
the part of its core that other tools can call.
"""
from __future__ import annotations

from . import kernel, quantities
from .kernel import LABELS, READINGS, read

__all__ = ["LABELS", "READINGS", "kernel", "quantities", "read"]
