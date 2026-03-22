"""
Modular time field and SL(2,ℤ) geometry on the upper half-plane ℍ.

The central object is the modular time field:

    τ(x) ≡ t_R(x) / t_U(x) = x(x) + i·y(x),   y > 0

mapping spacetime M^{1,3} → ℍ, with the Poincaré metric
ds² = (dx² + dy²) / y² giving constant curvature K = −1.

Modular symmetry:  τ ↦ (aτ + b)/(cτ + d),  [[a,b],[c,d]] ∈ SL(2,ℤ)
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Tuple

import numpy as np


# ── SL(2,ℤ) generators ──────────────────────────────────────────

S_MATRIX = np.array([[0, -1], [1, 0]], dtype=int)   # S: τ → −1/τ
T_MATRIX = np.array([[1, 1], [0, 1]], dtype=int)     # T: τ → τ + 1


def sl2z_transform(tau: complex, matrix: np.ndarray) -> complex:
    """
    Apply an SL(2,ℤ) transformation to τ ∈ ℍ.

    Parameters
    ----------
    tau : complex
        Point in the upper half-plane (Im τ > 0).
    matrix : array-like, shape (2, 2)
        Integer matrix with det = 1.

    Returns
    -------
    complex
        Transformed τ' = (aτ + b) / (cτ + d).

    Raises
    ------
    ValueError
        If τ is not in ℍ or matrix is not in SL(2,ℤ).
    """
    if tau.imag <= 0:
        raise ValueError(f"τ must be in ℍ (Im τ > 0), got Im τ = {tau.imag}")
    m = np.asarray(matrix, dtype=int)
    if m.shape != (2, 2) or int(np.linalg.det(m).round()) != 1:
        raise ValueError("Matrix must be 2×2 with determinant 1")
    a, b, c, d = m[0, 0], m[0, 1], m[1, 0], m[1, 1]
    return (a * tau + b) / (c * tau + d)


def reduce_to_fundamental(tau: complex, max_iter: int = 200) -> Tuple[complex, np.ndarray]:
    """
    Reduce τ to the SL(2,ℤ) fundamental domain F = {|τ|≥1, |Re τ|≤1/2}.

    Returns (τ_reduced, accumulated_matrix) so that
    τ_reduced = (aτ + b)/(cτ + d) with [[a,b],[c,d]] = accumulated_matrix.
    """
    acc = np.eye(2, dtype=int)
    z = tau
    for _ in range(max_iter):
        # T-shift: bring Re τ into (−1/2, 1/2]
        n = round(z.real)
        if n != 0:
            z -= n
            shift = np.array([[1, -n], [0, 1]], dtype=int)
            acc = shift @ acc
        # S-inversion if |τ| < 1
        if abs(z) < 1.0 - 1e-12:
            z = -1.0 / z
            acc = S_MATRIX @ acc
        else:
            break
    return z, acc


# ── Hyperbolic geometry ──────────────────────────────────────────

def hyperbolic_distance(tau1: complex, tau2: complex) -> float:
    """
    Geodesic distance on ℍ between τ₁ and τ₂.

    d(τ₁, τ₂) = arccosh(1 + |τ₁ − τ₂|² / (2 · Im τ₁ · Im τ₂))
    """
    if tau1.imag <= 0 or tau2.imag <= 0:
        raise ValueError("Both points must be in ℍ")
    num = abs(tau1 - tau2) ** 2
    denom = 2.0 * tau1.imag * tau2.imag
    return math.acosh(1.0 + num / denom)


def poincare_metric(tau: complex) -> float:
    """
    Conformal factor of the Poincaré metric at τ: ds² = (dx²+dy²)/y².

    Returns 1/y² where y = Im τ.
    """
    return 1.0 / (tau.imag ** 2)


# ── τ-field class ────────────────────────────────────────────────

@dataclass
class TauField:
    """
    A τ-field configuration on a spatial grid.

    Represents τ(x) : M → ℍ discretised on an N-point grid.

    Attributes
    ----------
    values : np.ndarray, shape (..., ) of complex
        τ values at each grid point, all with Im > 0.
    """

    values: np.ndarray

    def __post_init__(self):
        self.values = np.asarray(self.values, dtype=complex)
        if np.any(self.values.imag <= 0):
            raise ValueError("All τ values must lie in ℍ (Im τ > 0)")

    # ── accessors ────────────────────────────────────────────
    @property
    def x(self) -> np.ndarray:
        """Re(τ) — related to phase of rotational time."""
        return self.values.real

    @property
    def y(self) -> np.ndarray:
        """Im(τ) — related to ratio |t_R|/t_U."""
        return self.values.imag

    @property
    def j_invariant(self) -> np.ndarray:
        """Klein j-invariant j(τ) (numerical, via q-expansion to 20 terms)."""
        from mtft.forms import _j_invariant_array
        return _j_invariant_array(self.values)

    # ── energy density ───────────────────────────────────────
    def energy_density(self, grad_tau: np.ndarray, V: np.ndarray | None = None) -> np.ndarray:
        """
        τ-field energy density:  ρ_τ = ½|∇τ|² + V(τ)

        Parameters
        ----------
        grad_tau : array
            |∇τ|² at each grid point.
        V : array, optional
            Potential V(τ) at each grid point.  Defaults to 0.
        """
        rho = 0.5 * np.asarray(grad_tau)
        if V is not None:
            rho = rho + np.asarray(V)
        return rho

    # ── SL(2,ℤ) action ──────────────────────────────────────
    def transform(self, matrix: np.ndarray) -> "TauField":
        """Apply an SL(2,ℤ) matrix element-wise."""
        m = np.asarray(matrix, dtype=int)
        a, b, c, d = m[0, 0], m[0, 1], m[1, 0], m[1, 1]
        new_vals = (a * self.values + b) / (c * self.values + d)
        return TauField(new_vals)

    # ── static constructors ──────────────────────────────────
    @staticmethod
    def logarithmic_vortex(
        r: np.ndarray,
        A: float,
        r0: float = 1.0,
        y0: float = 1.0,
    ) -> "TauField":
        """
        Static spherical τ-vortex ansatz:

            τ(r) = A ln(r/r₀) + i·y₀

        Gives ρ_τ = A²/(2r²) → flat rotation curves v²∞ = 2πGA².
        """
        real_part = A * np.log(np.maximum(r, 1e-30) / r0)
        imag_part = np.full_like(r, y0, dtype=float)
        return TauField(real_part + 1j * imag_part)

    @staticmethod
    def constant(tau0: complex, shape: tuple) -> "TauField":
        """Uniform τ-field (cosmological background)."""
        return TauField(np.full(shape, tau0, dtype=complex))
