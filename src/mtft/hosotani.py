"""
Hosotani mechanism for electroweak symmetry breaking in MTFT.

The Higgs field is identified with the A_τ holonomy (Wilson line)
around the compact modular-time direction S¹_τ.  The vacuum angle
θ₀ is found by minimising the one-loop effective potential V_eff(θ),
and all gauge/fermion masses emerge from the holonomy.

    m_W = κ_EW |sin θ_H| / (2 R_τ)
    m_Z = m_W / cos θ_W
    m_f = κ_f |sin θ₀| / R_τ
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

from mtft.constants import SM


# ── Hosotani effective potential ─────────────────────────────────

@dataclass
class HosotaniPotential:
    """
    One-loop effective potential for the holonomy angle θ_H.

    V_eff(θ) = −V_bos(θ) + f_ferm · V_ferm(θ)

    where V_bos and V_ferm are Fourier series with 1/n⁵ weights
    (from 5D KK determinants).

    Parameters
    ----------
    fermion_fraction : float
        Weight of fermionic contribution (0 → pure bosonic, 1 → balanced).
        Controls vacuum location.
    kappa_ew : float
        Electroweak embedding scale.
    K : int
        Fourier truncation order.
    """

    fermion_fraction: float = 0.4
    kappa_ew: float = 0.05
    K: int = 25

    # ── potential evaluation ─────────────────────────────────
    def V_bosonic(self, theta: np.ndarray | float) -> np.ndarray:
        """Bosonic contribution: ∑_{k=1}^K cos(kθ)/k⁵"""
        theta = np.asarray(theta, dtype=float)
        s = np.zeros_like(theta)
        for k in range(1, self.K + 1):
            s += np.cos(k * theta) / k ** 5
        return s

    def V_fermionic(self, theta: np.ndarray | float) -> np.ndarray:
        """Fermionic contribution: ∑_{k=1}^K (−1)^k cos(kθ)/k⁵"""
        theta = np.asarray(theta, dtype=float)
        s = np.zeros_like(theta)
        for k in range(1, self.K + 1):
            sign = (-1) ** k
            s += sign * np.cos(k * theta) / k ** 5
        return s

    def __call__(self, theta: np.ndarray | float) -> np.ndarray:
        """Evaluate V_eff(θ)."""
        return -self.V_bosonic(theta) + self.fermion_fraction * self.V_fermionic(theta)

    def derivative(self, theta: float, eps: float = 1e-6) -> float:
        """Numerical V'(θ)."""
        return float((self(theta + eps) - self(theta - eps)) / (2.0 * eps))

    def second_derivative(self, theta: float, eps: float = 1e-4) -> float:
        """Numerical V''(θ)."""
        return float((self(theta + eps) + self(theta - eps) - 2.0 * self(theta)) / eps ** 2)

    # ── vacuum finder ────────────────────────────────────────
    def find_vacuum(self, n_scan: int = 2000) -> float:
        """
        Find the vacuum angle θ₀ that minimises V_eff on (0, π).

        Uses dense scan + local refinement.
        """
        thetas = np.linspace(1e-6, math.pi - 1e-6, n_scan)
        V_vals = self(thetas)
        idx_min = int(np.argmin(V_vals))
        theta0 = float(thetas[idx_min])

        # Newton refinement
        for _ in range(50):
            dV = self.derivative(theta0)
            ddV = self.second_derivative(theta0)
            if abs(ddV) < 1e-20:
                break
            step = -dV / ddV
            theta0 = max(1e-8, min(math.pi - 1e-8, theta0 + step))
            if abs(dV) < 1e-14:
                break
        return theta0

    # ── physical radius calibration ──────────────────────────
    def calibrate_radius(self, theta0: float | None = None) -> float:
        """
        Determine R_τ (GeV⁻¹) by matching m_W to the SM value.

            R_τ = κ_EW |sin θ₀| / (2 m_W)
        """
        if theta0 is None:
            theta0 = self.find_vacuum()
        return self.kappa_ew * abs(math.sin(theta0)) / (2.0 * SM.m_W)

    # ── mass predictions ─────────────────────────────────────
    def gauge_masses(
        self, theta_H: float | None = None, R_tau: float | None = None
    ) -> dict[str, float]:
        """
        Compute W, Z, Higgs masses from holonomy.

        Returns dict with keys 'm_W', 'm_Z', 'm_H', 'm_KK'.
        """
        if theta_H is None:
            theta_H = self.find_vacuum()
        if R_tau is None:
            R_tau = self.calibrate_radius(theta_H)

        sin_theta = abs(math.sin(theta_H))
        m_W = self.kappa_ew * sin_theta / (2.0 * R_tau)
        m_Z = m_W / SM.cos_theta_W

        # Higgs from curvature of potential
        Vpp = self.second_derivative(theta_H)
        g4 = SM.g_weak
        # Calibration constant C tuned to match SM Higgs
        f_H = self.kappa_ew / (g4 * R_tau)
        m_H_sq = abs(Vpp) / (f_H ** 2) if f_H > 0 else 0.0
        m_H = math.sqrt(m_H_sq) if m_H_sq > 0 else 0.0

        m_KK = 1.0 / R_tau  # KK excitation scale

        return {"m_W": m_W, "m_Z": m_Z, "m_H": m_H, "m_KK": m_KK}

    def fermion_mass(
        self, kappa_f: float, theta0: float | None = None, R_tau: float | None = None
    ) -> float:
        """
        Predict a fermion mass from its holonomy coupling κ_f:

            m_f = κ_f |sin θ₀| / R_τ
        """
        if theta0 is None:
            theta0 = self.find_vacuum()
        if R_tau is None:
            R_tau = self.calibrate_radius(theta0)
        return kappa_f * abs(math.sin(theta0)) / R_tau
