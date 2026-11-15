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

    CC-72 (2026-11-15), caveat: β_x and β_y were fitted to points that are not one
    cascade; see ``cc72_cascade_check``.  The values stay exposed and unconfirmed.
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


# ═══════════════════════════════════════════════════════════════
#  CC-72 (2026-11-15): the fold exponents wait on a cascade search
# ═══════════════════════════════════════════════════════════════
#
# Caveat, no numeric change.  β_x = +0.29 and β_y = −0.21 were fitted to the five points of
# Paper 4's off-axis table (periods 48, 96, 192, 384, 768).  The four gaps between those points
# do not behave like one period-doubling cascade: a cascade's gaps shrink by Feigenbaum's δ at
# each doubling, δ³ = 101.8 over three, while these shrink 1.81 times over three doublings and
# the third gap is larger than the second.  So the five points are not consecutive doublings of
# one cascade, and exponents fitted to them are unconfirmed.  β_x, β_y, DELTA_X and DELTA_Y stay
# exposed (append-only).  What would settle them is an off-axis search that confirms each
# doubling before any fit.  On the real axis the Burning Ship's map is x² + c, so its real
# cascade has Feigenbaum's δ exactly: ``real_axis_cascade`` is that control.
#
# Found on 2026-11-06 in the Modular Plane lab (record MP05, Room 11).  The gaps are as
# transcribed there from Paper 4's table; the transcription is to be checked against the paper
# at audit.

PAPER4_OFFAXIS_GAPS = (2.42e-4, 1.81e-4, 2.08e-4, 1.34e-4)   # |c_n − c_{n−1}|, periods 48 → 768

CC72_DEPENDENTS = (
    "mtft.constants.BETA_X", "mtft.constants.BETA_Y", "mtft.constants.DELTA_X", "mtft.constants.DELTA_Y",
    "mtft.constants.MassRatios.m_tau_over_m_mu",
    "mtft.burning_ship.AnisotropicScaling", "mtft.burning_ship.satellite_mass_scale",
)


def real_axis_cascade(n_max: int = 10) -> dict:
    """Superstable parameters of x ↦ x² + c with period 2ⁿ, n = 0 … n_max, by Newton's method.

    c_n solves f_c^{2ⁿ}(0) = 0; the derivative in c follows the orbit, x′ ↦ 2x·x′ + 1.  The gaps
    c_{n−1} − c_n shrink by a ratio that tends to Feigenbaum's δ, and c_n tends to the Feigenbaum
    point −1.401155189092….  In double precision this is reliable to n_max = 11 (period 2048).

    Returns ``{"c": [...], "gaps": [...], "ratios": [...], "c_inf": float}``.
    """
    if not 2 <= n_max <= 11:
        raise ValueError("n_max must lie between 2 and 11 (double precision)")
    cs = [0.0, -1.0]
    for n in range(2, n_max + 1):
        c = cs[-1] + (cs[-1] - cs[-2]) / FEIGENBAUM_DELTA
        for _ in range(80):
            x, dx = 0.0, 0.0
            for _ in range(2 ** n):
                x, dx = x * x + c, 2.0 * x * dx + 1.0
            step = x / dx
            c -= step
            if abs(step) < 1e-15:
                break
        cs.append(c)
    gaps = [cs[i - 1] - cs[i] for i in range(1, len(cs))]
    ratios = [gaps[i - 1] / gaps[i] for i in range(1, len(gaps))]
    return {"c": cs, "gaps": gaps, "ratios": ratios, "c_inf": cs[-1] + (cs[-1] - cs[-2]) / (ratios[-1] - 1.0)}


def cc72_cascade_check() -> dict:
    """CC-72: are the five points of Paper 4's off-axis table one period-doubling cascade?  No.

    Compares the table's four gaps with what a cascade does.  ``ratios`` are the successive gap
    ratios (a cascade's tend to δ = 4.669…); ``shrink`` is the first gap over the last (a cascade
    over three doublings: δ³ = 101.8).  ``one_cascade`` asks that every ratio lie within a quarter
    of δ, a window far wider than a cascade needs at these periods.
    """
    g = PAPER4_OFFAXIS_GAPS
    ratios = [g[i] / g[i + 1] for i in range(len(g) - 1)]
    monotone = all(r > 1.0 for r in ratios)
    return {
        "gaps": g,
        "ratios": ratios,
        "shrink": g[0] / g[-1],
        "cascade_shrink": FEIGENBAUM_DELTA ** (len(g) - 1),
        "monotone": monotone,
        "one_cascade": monotone and all(abs(r / FEIGENBAUM_DELTA - 1.0) < 0.25 for r in ratios),
        "status": "caveat: β_x, β_y unconfirmed until an off-axis cascade search confirms each doubling",
        "dependents": CC72_DEPENDENTS,
    }
