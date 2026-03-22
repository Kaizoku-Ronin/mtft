"""
MTFT Arithmetic Weights and Holonomy Stiffness
================================================

The arithmetic weights wₙ = Σ_{d|n} (log d)/d are the atoms of MTFT.
Their Dirichlet series satisfies W(s) = −ζ(s)ζ′(s+1).

Key series at modular depth y:
    aₙ(y)  = wₙ e^{−2πyn}      (damped weights)
    S₁(y)  = Σ aₙ(y)           (partition function / stiffness action)
    S(y)   = Σ n² aₙ(y)        (holonomy stiffness)
    C_N(y) = Σ n² aₙ(y) cos(2πn/N)  (center-projected stiffness)
    μ_N(y) = min_m Σ n² aₙ(y) (1 − cos(2πnm/N))  (mass gap)

The fine-structure constant is α⁻¹ = 2/S₁(y_spec).
The confinement lock is C_3(y_conf) = 0  (Hessian isotropy).

Reference: Papers 5, 7, 9, 13, 24; Dictionary IV–V.
"""

from __future__ import annotations

import math
from functools import lru_cache

import numpy as np


# ── Arithmetic weight function ───────────────────────────────

def weight(n: int) -> float:
    """
    Compute wₙ = Σ_{d|n} (log d)/d.

    These are rigid number-theoretic objects — not fitted.
    w₁ = 0, w₂ = (log 2)/2, w₃ = (log 3)/3, w₆ = w₂ + w₃ + (log 6)/6, ...
    """
    if n < 1:
        raise ValueError("n must be ≥ 1")
    total = 0.0
    for d in range(1, n + 1):
        if n % d == 0:
            if d > 1:
                total += math.log(d) / d
    return total


@lru_cache(maxsize=8192)
def _weight_cached(n: int) -> float:
    return weight(n)


def weight_array(n_max: int) -> np.ndarray:
    """Compute wₙ for n = 1, ..., n_max as a NumPy array."""
    return np.array([_weight_cached(n) for n in range(1, n_max + 1)])


# ── Euler-shifted weights (Paper 20) ─────────────────────────

def weight_euler(n: int) -> float:
    """
    Euler-shifted weight: w_n^(e) = w_n + σ_{-1}(n).

    The substitution log k → log(ek) activates w₁ = 0 → w₁^(e) = 1,
    creating the electron as the vacuum's first excitation.
    """
    sigma_minus1 = sum(1.0 / d for d in range(1, n + 1) if n % d == 0)
    return weight(n) + sigma_minus1


# ── Damped coefficients ──────────────────────────────────────

def damped_weight(n: int, y: float) -> float:
    """aₙ(y) = wₙ exp(−2πyn)."""
    return _weight_cached(n) * math.exp(-2.0 * math.pi * y * n)


def damped_weights(y: float, n_max: int = 500) -> np.ndarray:
    """Array of aₙ(y) for n = 1..n_max."""
    ns = np.arange(1, n_max + 1)
    ws = weight_array(n_max)
    return ws * np.exp(-2.0 * math.pi * y * ns)


# ── Stiffness series ─────────────────────────────────────────

def S1(y: float, n_max: int = 500) -> float:
    """
    Partition function / stiffness action:
        S₁(y) = Σ wₙ e^{−2πyn}

    The fine-structure constant is α⁻¹ = 2/S₁(y_spec).
    """
    return float(np.sum(damped_weights(y, n_max)))


def stiffness_S(y: float, n_max: int = 500) -> float:
    """
    Holonomy stiffness:
        S(y) = Σ n² wₙ e^{−2πyn}
    """
    ns = np.arange(1, n_max + 1)
    ws = weight_array(n_max)
    return float(np.sum(ns**2 * ws * np.exp(-2 * math.pi * y * ns)))


def center_stiffness(y: float, N: int = 3, n_max: int = 500) -> float:
    """
    Z_N center-projected stiffness:
        C_N(y) = Σ n² wₙ e^{−2πyn} cos(2πn/N)

    C_3(y_conf) = 0 is the confinement lock.
    """
    ns = np.arange(1, n_max + 1)
    ws = weight_array(n_max)
    cos_terms = np.cos(2.0 * math.pi * ns / N)
    return float(np.sum(ns**2 * ws * np.exp(-2 * math.pi * y * ns) * cos_terms))


def mass_gap_stiffness(y: float, N: int = 3, n_max: int = 500) -> float:
    """
    SU(N) mass gap stiffness (minimum over modes m):
        μ_N(y) = min_{m=1..N-1} Σ n² wₙ e^{−2πyn} (1 − cos(2πnm/N))

    Theorem: μ_N(y) > 0 unconditionally for all N≥2, y>0.
    """
    ns = np.arange(1, n_max + 1)
    ws = weight_array(n_max)
    damped = ns**2 * ws * np.exp(-2 * math.pi * y * ns)

    mu_min = float("inf")
    for m in range(1, N):
        gap_terms = 1.0 - np.cos(2 * math.pi * ns * m / N)
        mu_m = float(np.sum(damped * gap_terms))
        mu_min = min(mu_min, mu_m)
    return mu_min


# ── SU(3) Hessian eigenvalues (Paper 5 §24.4, Bridge §4.3) ──

def su3_hessian_eigenvalues(y: float, n_max: int = 500) -> tuple[float, float]:
    """
    SU(3) holonomy Hessian eigenvalues:
        λ₁ = (S_N + 2C_N) / 3
        λ₂ = (S_N − C_N) / 3

    At y_conf ≈ 0.18174:  λ₁/λ₂ = 1.000001 (the spectral lock).
    """
    S = stiffness_S(y, n_max)
    C = center_stiffness(y, N=3, n_max=n_max)
    lam1 = (S + 2 * C) / 3.0
    lam2 = (S - C) / 3.0
    return lam1, lam2


def find_confinement_depth(
    y_min: float = 0.10, y_max: float = 0.25, tol: float = 1e-8, n_max: int = 500
) -> float:
    """
    Find y_conf where C_3(y) = 0 (SU(3) Hessian isotropy).

    Uses bisection on C_3(y).
    """
    # C_3 starts positive and crosses zero
    c_lo = center_stiffness(y_min, N=3, n_max=n_max)
    c_hi = center_stiffness(y_max, N=3, n_max=n_max)

    if c_lo * c_hi > 0:
        # Try a wider scan
        for y_test in np.linspace(y_min, y_max, 200):
            c_test = center_stiffness(y_test, N=3, n_max=n_max)
            if c_test * c_lo < 0:
                y_max = y_test
                c_hi = c_test
                break

    lo, hi = y_min, y_max
    for _ in range(100):
        mid = (lo + hi) / 2.0
        c_mid = center_stiffness(mid, N=3, n_max=n_max)
        if abs(c_mid) < tol:
            return mid
        if c_mid * c_lo > 0:
            lo = mid
            c_lo = c_mid
        else:
            hi = mid
    return (lo + hi) / 2.0


# ── Lambert series ───────────────────────────────────────────

def lambert_series(q: complex, n_max: int = 200) -> complex:
    """
    MTFT Lambert series: W(q) = Σ_{d≥1} (log d)/d · q^d/(1−q^d).

    This is the q-series feeding the modular holonomy potential.
    """
    total = 0.0 + 0j
    for d in range(2, n_max + 1):
        qd = q ** d
        if abs(1.0 - qd) < 1e-30:
            continue
        total += (math.log(d) / d) * qd / (1.0 - qd)
    return total


# ── Torque flux convergence (Paper 7, §9) ────────────────────

def torque_partial_sum(N: int) -> float:
    """
    Partial sum of the normalised torque: (1/N²) Σ_{n=1}^{N} n·wₙ.

    Converges to T∞ = −ζ′(2)/2 ≈ 0.4688.
    """
    total = 0.0
    for n in range(1, N + 1):
        total += n * _weight_cached(n)
    return total / (N ** 2)


# ── Dirichlet series identity check ──────────────────────────

def dirichlet_F(s: complex, n_max: int = 500) -> complex:
    """
    Compute F(s) = Σ f(n)/n^s where f = σ * Λ.

    Should equal −ζ(s−1)ζ′(s) (the master Dirichlet identity).
    """
    total = 0.0 + 0j
    for n in range(1, n_max + 1):
        f_n = sum(math.log(d) for d in range(1, n + 1) if n % d == 0)
        total += f_n / n ** s
    return total
