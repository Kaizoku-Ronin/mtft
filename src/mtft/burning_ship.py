"""
Burning Ship Fractal: The Fermion Vacuum
=========================================

The Burning Ship iteration:
    z_{n+1} = (|Re(z_n)| + i|Im(z_n)|)² + c

Unlike the Mandelbrot z² + c (analytic), the absolute values create
fold lines where the map is non-analytic. This governs the fermion
sector of MTFT:

    Mandelbrot (analytic)      → gauge sector  → Yang-Mills mass gap
    Burning Ship (non-analytic) → fermion sector → Koide formula

Phase space has four sectors from sgn(x), sgn(y):
    S₊₊, S₋₊, S₊₋, S₋₋

with Jacobian:
    J = [[2x, −2y], [2 sgn(x)sgn(y)|y|, 2 sgn(x)sgn(y)|x|]]
    det(J) = 4|z|²  (continuous)
    tr(J)  = 2x + 2 sgn(x)sgn(y)|x|  (discontinuous at x=0)

Anisotropic Feigenbaum constants:
    δₓ ≈ 7.30 (attractor fold, quark direction)
    δᵧ ≈ 3.38 (barrier fold, lepton direction)
    β_x = +0.29, β_y = −0.21, |β_x| + |β_y| = 1/2

Reference: Paper 4, Paper 6, Chapter 10–11.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from mtft.constants import FEIGENBAUM_DELTA, DELTA_X, DELTA_Y, BETA_X, BETA_Y


# ═══════════════════════════════════════════════════════════════
#  Burning Ship Iteration
# ═══════════════════════════════════════════════════════════════

def burning_ship_iterate(c: complex, max_iter: int = 200, escape: float = 100.0) -> int:
    """
    Iterate z_{n+1} = (|Re z| + i|Im z|)² + c.

    Returns escape time (max_iter if bounded).
    """
    z = 0.0 + 0.0j
    for n in range(max_iter):
        x, y = abs(z.real), abs(z.imag)
        z = complex(x, y) ** 2 + c
        if abs(z) > escape:
            return n
    return max_iter


def burning_ship_array(
    re_range: tuple = (-2.5, 1.5),
    im_range: tuple = (-2.0, 1.0),
    resolution: int = 500,
    max_iter: int = 200,
) -> np.ndarray:
    """
    Compute escape-time array for the Burning Ship fractal.

    Returns (resolution × resolution) integer array.
    """
    re = np.linspace(re_range[0], re_range[1], resolution)
    im = np.linspace(im_range[0], im_range[1], resolution)
    result = np.zeros((resolution, resolution), dtype=int)
    for i, y in enumerate(im):
        for j, x in enumerate(re):
            result[i, j] = burning_ship_iterate(complex(x, y), max_iter)
    return result


# ═══════════════════════════════════════════════════════════════
#  Phase Space Sectors and Jacobian
# ═══════════════════════════════════════════════════════════════

def sector(z: complex) -> str:
    """Identify the phase-space sector of z."""
    sx = "+" if z.real >= 0 else "-"
    sy = "+" if z.imag >= 0 else "-"
    return f"S_{sx}{sy}"


def jacobian(z: complex) -> np.ndarray:
    """
    Jacobian of the Burning Ship map at z = x + iy:

        J = [[2x, −2y],
             [2 sgn(x)sgn(y)|y|, 2 sgn(x)sgn(y)|x|]]

    det(J) = 4|z|² (continuous across folds)
    tr(J) = 2x + 2 sgn(x)sgn(y)|x| (discontinuous at x=0)
    """
    x, y = z.real, z.imag
    sx = 1.0 if x >= 0 else -1.0
    sy = 1.0 if y >= 0 else -1.0
    return np.array([
        [2 * x, -2 * y],
        [2 * sx * sy * abs(y), 2 * sx * sy * abs(x)],
    ])


def jacobian_determinant(z: complex) -> float:
    """det(J) = 4|z|² — continuous across all folds."""
    return 4.0 * abs(z) ** 2


def jacobian_trace(z: complex) -> float:
    """tr(J) = 2x + 2 sgn(x)sgn(y)|x| — discontinuous at x=0."""
    x, y = z.real, z.imag
    sx = 1.0 if x >= 0 else -1.0
    sy = 1.0 if y >= 0 else -1.0
    return 2.0 * x + 2.0 * sx * sy * abs(x)


# ═══════════════════════════════════════════════════════════════
#  Anisotropic Feigenbaum Constants
# ═══════════════════════════════════════════════════════════════

@dataclass(frozen=True)
class AnisotropicScaling:
    """
    Anisotropic universality of the Burning Ship (Paper 4).

    The standard Feigenbaum δ splits into direction-dependent constants:
        δₓ ≈ δ^1.29 ≈ 7.30  (attractor fold / quark direction)
        δᵧ ≈ δ^0.79 ≈ 3.38  (barrier fold / lepton direction)

    Fold approach exponents:
        β_x = +0.29  (attractor approach rate)
        β_y = −0.21  (barrier approach rate)
        |β_x| + |β_y| = 1/2  (fold budget constraint)

    Power laws:
        δ ∝ ρ_x^{+0.29}  (x-fold)
        δ ∝ ρ_y^{−0.21}  (y-fold)
    """
    delta: float = FEIGENBAUM_DELTA
    delta_x: float = DELTA_X
    delta_y: float = DELTA_Y
    beta_x: float = BETA_X
    beta_y: float = BETA_Y

    @property
    def fold_budget(self) -> float:
        """Should be ≈ 0.5."""
        return abs(self.beta_x) + abs(self.beta_y)

    def generation_mass_ratio(self, direction: str = "lepton") -> float:
        """
        Mass ratio between adjacent generations:
            m_{n+1}/m_n ~ δ_eff²

        For leptons: δ_y² ≈ 11.4
        For quarks:  δ_x² ≈ 53.3
        """
        d = self.delta_y if direction == "lepton" else self.delta_x
        return d ** 2


ANISOTROPIC = AnisotropicScaling()


# ═══════════════════════════════════════════════════════════════
#  Fold Proximity Data (Paper 4, Table 2)
# ═══════════════════════════════════════════════════════════════

FOLD_PROXIMITY = {
    # period: (min_|x|, min_|y|, x_crossings, y_crossings)
    3:   (2.51e-1, 4.31e-1, 0, 0),
    6:   (4.37e-2, 1.32e-1, 3, 0),
    12:  (2.82e-2, 4.72e-2, 3, 0),
    24:  (2.22e-2, 8.35e-4, 7, 4),
    48:  (7.02e-3, 1.92e-3, 15, 4),
    384: (9.02e-4, 2.60e-6, 255, 128),
    768: (1.11e-5, 3.70e-6, 511, 255),
}


# ═══════════════════════════════════════════════════════════════
#  Three Generations from Self-Similar Satellites (Paper 6 §5)
# ═══════════════════════════════════════════════════════════════

def satellite_mass_scale(n: int, delta_eff: float = None) -> float:
    """
    Mass scale of the n-th generation satellite relative to 1st:

        m_n / m_1 ~ δ_eff^{−2(n−1)}

    Higher-order satellites (n ≥ 4) fall below the EW scale
    and are not realised as stable particles.
    """
    if delta_eff is None:
        delta_eff = DELTA_Y
    return delta_eff ** (-2 * (n - 1))
