"""
Modular forms relevant to MTFT.

Implements Dedekind η, Jacobi θ₃, Eisenstein E_{2k}, the modular
discriminant Δ(τ) = η(τ)²⁴, and the spectral determinant identity:

    det(1 − q^{1/m} P_m) = η(τ)^{−1/m} θ₃(0,τ)^{1/m}

which bridges Fredholm determinants (chaos/RG) to CFT partition functions.
"""

from __future__ import annotations

import cmath
import math

import numpy as np


# ── nome ─────────────────────────────────────────────────────────

def nome(tau: complex) -> complex:
    """Compute the nome q = exp(2πiτ)."""
    return cmath.exp(2j * math.pi * tau)


def nome_array(tau: np.ndarray) -> np.ndarray:
    """Vectorised nome."""
    return np.exp(2j * np.pi * tau)


# ── Dedekind eta ─────────────────────────────────────────────────

def dedekind_eta(tau: complex, N: int = 80) -> complex:
    """
    Dedekind eta function:

        η(τ) = q^{1/24} ∏_{n=1}^{N} (1 − qⁿ)

    where q = e^{2πiτ}.
    """
    q = nome(tau)
    q24 = cmath.exp(2j * math.pi * tau / 24.0)
    prod = 1.0 + 0j
    qn = 1.0 + 0j
    for n in range(1, N + 1):
        qn *= q
        prod *= (1.0 - qn)
    return q24 * prod


def dedekind_eta_array(tau: np.ndarray, N: int = 80) -> np.ndarray:
    """Vectorised Dedekind η."""
    return np.vectorize(lambda t: dedekind_eta(complex(t), N))(tau)


# ── Jacobi theta ─────────────────────────────────────────────────

def jacobi_theta3(z: complex, tau: complex, N: int = 60) -> complex:
    """
    Jacobi theta function θ₃(z, τ):

        θ₃(z,τ) = ∑_{n=−N}^{N} q^{n²} e^{2πinz}

    where q = e^{πiτ}.
    """
    q = cmath.exp(1j * math.pi * tau)
    total = 0.0 + 0j
    for n in range(-N, N + 1):
        total += q ** (n * n) * cmath.exp(2j * math.pi * n * z)
    return total


# ── Eisenstein series ────────────────────────────────────────────

def eisenstein_E2k(tau: complex, k: int = 2, N: int = 200) -> complex:
    """
    Normalised Eisenstein series E_{2k}(τ) via q-expansion.

    E_{2k}(τ) = 1 − (4k / B_{2k}) ∑_{n=1}^{N} σ_{2k−1}(n) qⁿ

    where σ_k(n) = ∑_{d|n} d^k and B_{2k} are Bernoulli numbers.

    Parameters
    ----------
    k : int
        Half-weight (k=1 → E₂, k=2 → E₄, k=3 → E₆).
    """
    bernoulli = {1: 1 / 6, 2: -1 / 30, 3: 1 / 42, 4: -1 / 30, 5: 5 / 66}
    if k not in bernoulli:
        raise ValueError(f"k must be in 1..5, got {k}")
    B2k = bernoulli[k]
    weight = 2 * k
    coeff = -2.0 * weight / B2k

    q = nome(tau)
    total = 0.0 + 0j
    qn = 1.0 + 0j
    for n in range(1, N + 1):
        qn *= q
        sigma = sum(d ** (weight - 1) for d in range(1, n + 1) if n % d == 0)
        total += sigma * qn
    return 1.0 + coeff * total


# ── Modular discriminant ─────────────────────────────────────────

def modular_discriminant(tau: complex) -> complex:
    """Δ(τ) = η(τ)²⁴  (weight-12 cusp form)."""
    eta = dedekind_eta(tau)
    return eta ** 24


# ── Klein j-invariant ────────────────────────────────────────────

def j_invariant(tau: complex) -> complex:
    """
    j(τ) = E₄(τ)³ / Δ(τ) = 1728 E₄³ / (E₄³ − E₆²)
    """
    E4 = eisenstein_E2k(tau, k=2)
    E6 = eisenstein_E2k(tau, k=3)
    E4_cubed = E4 ** 3
    denom = E4_cubed - E6 ** 2
    if abs(denom) < 1e-30:
        return complex("inf")
    return 1728.0 * E4_cubed / denom


def _j_invariant_array(tau: np.ndarray) -> np.ndarray:
    """Vectorised j-invariant (internal helper)."""
    return np.vectorize(lambda t: j_invariant(complex(t)))(tau)


# ── Spectral determinant identity ────────────────────────────────

def spectral_determinant_lhs(tau: complex, m: int, N: int = 40) -> complex:
    """
    Left-hand side of the spectral determinant identity:

        det(1 − q^{1/m} P_m) ≈ ∏_{n=1}^{N} (1 − q^{(2n−1)/m})

    This is a truncated Fredholm determinant for the period-m
    transfer operator P_m with eigenvalues ≈ q^{(2n−1)/m}.
    """
    q = nome(tau)
    prod = 1.0 + 0j
    for n in range(1, N + 1):
        prod *= (1.0 - q ** ((2 * n - 1) / m))
    return prod


def spectral_determinant_rhs(tau: complex, m: int) -> complex:
    """
    Right-hand side:

        η(τ)^{−1/m} · θ₃(0, τ)^{1/m}

    This is the CFT partition function of a chiral boson at radius R² = m.
    """
    eta = dedekind_eta(tau)
    theta = jacobi_theta3(0, tau)
    return eta ** (-1.0 / m) * theta ** (1.0 / m)


def verify_spectral_identity(tau: complex, m: int = 2, N: int = 40) -> dict:
    """
    Verify the spectral determinant identity for given τ and m.

    Returns dict with 'lhs', 'rhs', 'ratio', 'relative_error'.
    """
    lhs = spectral_determinant_lhs(tau, m, N)
    rhs = spectral_determinant_rhs(tau, m)
    ratio = lhs / rhs if abs(rhs) > 1e-30 else complex("nan")
    rel_err = abs(ratio - 1.0)
    return {"lhs": lhs, "rhs": rhs, "ratio": ratio, "relative_error": rel_err}
