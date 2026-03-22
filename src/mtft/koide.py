"""
Geometric Koide Theorem
========================

The Koide formula K = (Σmᵢ)/(Σ√mᵢ)² = 2/3 is not a numerical coincidence.
It is a geometric constraint from mass eigenstate stability in the MTFT
electroweak fold.

Key structures:
    - Mass-root space: w = (√m₁, √m₂, √m₃) ∈ ℝ³₊
    - Koide ratio: Q = |w|² / (w·1)²
    - Koide manifold: K ≅ S¹ × ℝ₊  (circle × scale)
    - Koide phase: w(φ) = r[1/√3 · 1 + √(2/3) · n(φ)]
    - 45° bisector stability under RG → K = 2/3 exactly

The value 2/3 arises from:
    K = 1/(3 cos²(45°)) = 1/(3 · 1/2) = 2/3

Reference: Paper 6, Chapter 16, Paper 26.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from mtft.constants import LEPTONS, QUARKS


# ═══════════════════════════════════════════════════════════════
#  Koide Ratio
# ═══════════════════════════════════════════════════════════════

def koide_ratio(m1: float, m2: float, m3: float) -> float:
    """
    Compute the Koide ratio:

        K = (m₁ + m₂ + m₃) / (√m₁ + √m₂ + √m₃)²

    For charged leptons: K = 0.666661 ≈ 2/3 to 10⁻⁵.
    """
    mass_sum = m1 + m2 + m3
    sqrt_sum = math.sqrt(m1) + math.sqrt(m2) + math.sqrt(m3)
    if sqrt_sum == 0:
        return float("nan")
    return mass_sum / (sqrt_sum ** 2)


def koide_leptons() -> float:
    """K for charged leptons (e, μ, τ)."""
    return koide_ratio(LEPTONS.e, LEPTONS.mu, LEPTONS.tau)


def koide_up_quarks() -> float:
    """K for up-type quarks (u, c, t)."""
    return koide_ratio(QUARKS.u, QUARKS.c, QUARKS.t)


def koide_down_quarks() -> float:
    """K for down-type quarks (d, s, b)."""
    return koide_ratio(QUARKS.d, QUARKS.s, QUARKS.b)


# ═══════════════════════════════════════════════════════════════
#  Koide Angle (Phase on the Koide Manifold)
# ═══════════════════════════════════════════════════════════════

def koide_angle(m1: float, m2: float, m3: float) -> float:
    """
    Extract the Koide phase φ from three masses.

    The parametrisation is:
        √mᵢ = M(1 + √2 cos(φ + 2πi/3)),  i = 0, 1, 2

    where M = (√m₁ + √m₂ + √m₃)/3.

    Returns φ in radians.
    """
    sqrts = [math.sqrt(m) for m in [m1, m2, m3]]
    M = sum(sqrts) / 3.0
    if M == 0:
        return float("nan")
    cos_phi = (sqrts[0] / M - 1.0) / math.sqrt(2.0)
    cos_phi = max(-1.0, min(1.0, cos_phi))
    return math.acos(cos_phi)


def koide_angle_leptons() -> float:
    """Koide phase for charged leptons. Should be ≈ 0.2222 rad ≈ 12.73°."""
    # NB: the "angle form" convention differs from Tano formula convention
    # This gives the angle in the Koide parametrisation
    return koide_angle(LEPTONS.e, LEPTONS.mu, LEPTONS.tau)


# ═══════════════════════════════════════════════════════════════
#  Mass-Root Geometry
# ═══════════════════════════════════════════════════════════════

def mass_root_vector(m1: float, m2: float, m3: float) -> np.ndarray:
    """Mass-root vector w = (√m₁, √m₂, √m₃)."""
    return np.array([math.sqrt(m1), math.sqrt(m2), math.sqrt(m3)])


def koide_as_angle(m1: float, m2: float, m3: float) -> float:
    """
    The Koide ratio as the angle between w and 1 = (1,1,1):

        Q = |w|²/(w·1)² = 1/(3 cos²θ)

    Returns θ in radians.  Q = 2/3 ↔ θ ≈ 54.74° (the "magic angle").
    """
    w = mass_root_vector(m1, m2, m3)
    ones = np.ones(3)
    cos_theta = np.dot(w, ones) / (np.linalg.norm(w) * np.linalg.norm(ones))
    return math.acos(cos_theta)


MAGIC_ANGLE = math.acos(1.0 / math.sqrt(3.0))  # ≈ 54.74° = arccos(1/√3)


# ═══════════════════════════════════════════════════════════════
#  Koide Manifold: K ≅ S¹ × ℝ₊
# ═══════════════════════════════════════════════════════════════

def koide_manifold_point(
    phi: float, r: float = 1.0
) -> np.ndarray:
    """
    Point on the Koide manifold K parametrised by (φ, r):

        w(φ, r) = r · [1/√3 · 1 + √(2/3) · n(φ)]

    where n(φ) = (cos φ, cos(φ+2π/3), cos(φ+4π/3)).
    All points satisfy Q = 2/3 exactly.

    Parameters
    ----------
    phi : float — Koide phase in [0, 2π)
    r : float   — overall scale |w|

    Returns mass-root vector w ∈ ℝ³.
    """
    n = np.array([
        math.cos(phi),
        math.cos(phi + 2 * math.pi / 3),
        math.cos(phi + 4 * math.pi / 3),
    ])
    w = r * (np.ones(3) / math.sqrt(3) + math.sqrt(2.0 / 3.0) * n)
    return w


def masses_from_koide(phi: float, r: float) -> np.ndarray:
    """Compute (m₁, m₂, m₃) from Koide manifold coordinates (φ, r)."""
    w = koide_manifold_point(phi, r)
    return w ** 2


def verify_koide_on_manifold(phi: float, r: float = 1.0) -> float:
    """Verify that Q = 2/3 for any point on the Koide manifold."""
    m = masses_from_koide(phi, r)
    return koide_ratio(m[0], m[1], m[2])


# ═══════════════════════════════════════════════════════════════
#  Tau Mass Prediction (Chapter 16, §7.1)
# ═══════════════════════════════════════════════════════════════

def predict_tau_mass(
    m_e_MeV: float = 0.51099895,
    m_mu_MeV: float = 105.6583755,
) -> dict:
    """
    Predict m_τ from the Koide constraint K = 2/3 + (m_e, m_μ).

    This gives m_τ to 0.006% precision.

    Returns dict with 'predicted_MeV', 'experimental_MeV', 'error_percent'.
    """
    sqrt_e = math.sqrt(m_e_MeV)
    sqrt_mu = math.sqrt(m_mu_MeV)

    # From K = 2/3:  m_τ + m_e + m_μ = (2/3)(√m_e + √m_μ + √m_τ)²
    # This is a quadratic in √m_τ.  Let x = √m_τ:
    # x² + (m_e + m_μ) = (2/3)(sqrt_e + sqrt_mu + x)²
    # x² + S_m = (2/3)(S_s + x)²  where S_m = m_e+m_μ, S_s = √m_e+√m_μ
    S_m = m_e_MeV + m_mu_MeV
    S_s = sqrt_e + sqrt_mu

    # Expand: x² + S_m = (2/3)(S_s² + 2S_s x + x²)
    # x² + S_m = (2/3)S_s² + (4/3)S_s x + (2/3)x²
    # (1/3)x² − (4/3)S_s x + S_m − (2/3)S_s² = 0
    # x² − 4 S_s x + 3 S_m − 2 S_s² = 0
    a_coeff = 1.0
    b_coeff = -4.0 * S_s
    c_coeff = 3.0 * S_m - 2.0 * S_s ** 2

    disc = b_coeff ** 2 - 4 * a_coeff * c_coeff
    if disc < 0:
        return {"error": "No real solution"}

    x_plus = (-b_coeff + math.sqrt(disc)) / (2 * a_coeff)
    x_minus = (-b_coeff - math.sqrt(disc)) / (2 * a_coeff)

    # Take the larger root (τ is heavy)
    sqrt_tau = max(x_plus, x_minus)
    m_tau_pred = sqrt_tau ** 2

    m_tau_exp = 1776.86
    error = abs(m_tau_pred - m_tau_exp) / m_tau_exp * 100

    return {
        "predicted_MeV": m_tau_pred,
        "experimental_MeV": m_tau_exp,
        "error_percent": error,
    }


# ═══════════════════════════════════════════════════════════════
#  Bisector Stability Argument (Paper 6 §6.3)
# ═══════════════════════════════════════════════════════════════

def bisector_stability_K() -> float:
    """
    Why K = 2/3 from the Burning Ship fold geometry:

        K = 1 / (3 cos²(45°)) = 1 / (3 · 1/2) = 2/3

    The 45° bisector is the unique stable fixed point of RG flow
    in the fold quadrant with δₓ ≠ δᵧ.
    """
    bisector_angle = math.pi / 4  # 45°
    return 1.0 / (3.0 * math.cos(bisector_angle) ** 2)


# ═══════════════════════════════════════════════════════════════
#  Comprehensive Koide Report
# ═══════════════════════════════════════════════════════════════

def koide_report() -> dict:
    """Full Koide analysis across leptons and quarks."""
    return {
        "leptons": {
            "K": koide_leptons(),
            "K_exact": 2 / 3,
            "error": abs(koide_leptons() - 2 / 3),
            "angle_rad": koide_angle_leptons(),
            "angle_deg": math.degrees(koide_angle_leptons()),
            "magic_angle_deg": math.degrees(MAGIC_ANGLE),
        },
        "up_quarks": {"K": koide_up_quarks()},
        "down_quarks": {"K": koide_down_quarks()},
        "bisector_prediction": bisector_stability_K(),
        "tau_prediction": predict_tau_mass(),
    }
