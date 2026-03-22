"""
Information geometry layer: the Counting → Geometry → Dynamics bridge.

The logistic map x_{n+1} = r x_n (1−x_n) generates period-m orbits
whose statistical ensemble defines a Fisher-Rao metric on parameter
space.  The Ricci scalar R_core at the Feigenbaum accumulation point
is a universal curvature invariant that enters MTFT as:

    m_τ² = c_m |R_core|

connecting microscopic chaos to macroscopic τ-geometry.

The spectral determinant identity
    det(1 − q^{1/m} P_m) = η(τ)^{−1/m} θ₃(0,τ)^{1/m}
is the Rosetta Stone linking:
    - Fredholm determinants (chaos / RG)
    - Modular forms (string partition functions)
    - CFT (chiral bosons at R² = m)
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np


# ── Logistic map tools ───────────────────────────────────────────

def logistic_iterate(r: float, x0: float, n_transient: int = 1000, n_sample: int = 500) -> np.ndarray:
    """
    Iterate the logistic map x_{n+1} = r x_n (1−x_n).

    Returns the last n_sample iterates after discarding transients.
    """
    x = x0
    for _ in range(n_transient):
        x = r * x * (1.0 - x)
    samples = np.empty(n_sample)
    for i in range(n_sample):
        x = r * x * (1.0 - x)
        samples[i] = x
    return samples


def lyapunov_exponent(r: float, x0: float = 0.3, n_iter: int = 10000) -> float:
    """
    Compute the Lyapunov exponent λ = ⟨ln|r(1−2x)|⟩ for the logistic map.
    """
    x = x0
    total = 0.0
    for _ in range(200):  # transient
        x = r * x * (1.0 - x)
    for _ in range(n_iter):
        deriv = abs(r * (1.0 - 2.0 * x))
        if deriv > 0:
            total += math.log(deriv)
        x = r * x * (1.0 - x)
    return total / n_iter


# Feigenbaum constants
FEIGENBAUM_DELTA = 4.669201609102990
FEIGENBAUM_ALPHA = 2.502907875095892
R_ACCUMULATION = 3.5699456  # r_∞ (onset of chaos)


# ── Fisher-Rao metric ───────────────────────────────────────────

def _histogram_distribution(samples: np.ndarray, n_bins: int = 100) -> np.ndarray:
    """Normalised histogram as a probability distribution."""
    counts, _ = np.histogram(samples, bins=n_bins, range=(0, 1))
    total = counts.sum()
    if total == 0:
        return np.ones(n_bins) / n_bins
    return counts / total


def fisher_rao_metric(
    r: float,
    dr: float = 1e-4,
    x0: float = 0.3,
    n_bins: int = 100,
) -> float:
    """
    Fisher information (scalar metric component) g_{rr} at parameter r.

        g_{rr} = ∑_i (∂_r √p_i)² / p_i ≈ ∑_i (∂_r p_i)² / p_i

    Computed via central-difference of the invariant density.
    """
    p_plus = _histogram_distribution(
        logistic_iterate(r + dr, x0), n_bins
    )
    p_minus = _histogram_distribution(
        logistic_iterate(r - dr, x0), n_bins
    )
    p_center = _histogram_distribution(
        logistic_iterate(r, x0), n_bins
    )

    dp = (p_plus - p_minus) / (2.0 * dr)

    # Fisher metric: g = ∑ (dp)² / p, avoiding division by zero
    mask = p_center > 1e-15
    g = np.sum(dp[mask] ** 2 / p_center[mask])
    return float(g)


def ricci_scalar_logistic(
    r: float,
    dr: float = 5e-4,
    x0: float = 0.3,
    n_bins: int = 100,
) -> float:
    """
    Ricci scalar R of the 1D Fisher-Rao manifold at parameter r.

    For a 1D Riemannian manifold with metric g(r):

        R = −(1/g) d²(ln g)/dr²  + (1/(2g²)) (dg/dr)²

    This is a curvature diagnostic — it diverges at Feigenbaum bifurcation
    points (R → −∞), signalling the chaos–geometry bridge.
    """
    g0 = fisher_rao_metric(r, dr=dr / 2, x0=x0, n_bins=n_bins)
    g_p = fisher_rao_metric(r + dr, dr=dr / 2, x0=x0, n_bins=n_bins)
    g_m = fisher_rao_metric(r - dr, dr=dr / 2, x0=x0, n_bins=n_bins)

    if g0 < 1e-20:
        return 0.0

    dg = (g_p - g_m) / (2.0 * dr)
    ddg = (g_p + g_m - 2.0 * g0) / (dr ** 2)

    ln_g = math.log(max(g0, 1e-30))
    ln_gp = math.log(max(g_p, 1e-30))
    ln_gm = math.log(max(g_m, 1e-30))
    dd_ln_g = (ln_gp + ln_gm - 2.0 * ln_g) / (dr ** 2)

    R = -dd_ln_g / g0 + dg ** 2 / (2.0 * g0 ** 2)
    return float(R)


# ── R_core: universal curvature at Feigenbaum point ──────────────

def R_core(dr: float = 1e-4, n_bins: int = 150) -> float:
    """
    Compute the Ricci curvature R_core at the Feigenbaum accumulation
    point r_∞ ≈ 3.5699456.

    This is the universal invariant bridging chaos → geometry → dynamics:
        m_τ² = c_m |R_core|
    """
    return ricci_scalar_logistic(R_ACCUMULATION, dr=dr, n_bins=n_bins)


# ── Zamolodchikov gradient flow ──────────────────────────────────

@dataclass
class ZamolodchikovFlow:
    """
    RG flow on the information manifold, identified with modular time flow:

        μ d/dμ = (1/κ)(Im τ) ∂_{Im τ}

    The β-function is the gradient of a c-function (Zamolodchikov's theorem)
    with respect to the Fisher-Rao metric.
    """

    kappa: float = 1.0

    def beta(self, r: float, dr: float = 1e-4) -> float:
        """
        Approximate β-function from the c-function gradient:

            β(r) ≈ −dg/dr / (2g)

        where g is the Fisher metric.
        """
        g_p = fisher_rao_metric(r + dr)
        g_m = fisher_rao_metric(r - dr)
        g_c = fisher_rao_metric(r)
        dg = (g_p - g_m) / (2.0 * dr)
        if g_c < 1e-20:
            return 0.0
        return -dg / (2.0 * g_c)
