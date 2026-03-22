"""
MTFT cosmology: modified Friedmann equations with τ-field dark sector.

The τ-field contributes an energy density ρ_τ and pressure p_τ to the
Friedmann equations.  Key features:

    - Oscillatory Newton's constant G(t_U) from temporal axion θ_T
    - τ-field perturbation theory in Newtonian gauge
    - Green-Schwarz anomaly cancellation: S_GS = ∫ B ∧ F_UR
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

from mtft.constants import PhysicalConstants as PC


# ── Background cosmology ─────────────────────────────────────────

@dataclass
class FriedmannMTFT:
    """
    Modified Friedmann equations with τ-field dark sector.

    The Hubble rate is:

        H² = (8πG/3)(ρ_matter + ρ_radiation + ρ_τ + ρ_Λ)

    where ρ_τ = ½ τ̇² + V(τ) is the τ-field energy density, acting as
    a dynamical dark component.

    Parameters
    ----------
    H0 : float
        Hubble constant in km/s/Mpc.
    Omega_m : float
        Present matter density parameter.
    Omega_r : float
        Present radiation density parameter.
    Omega_tau : float
        Present τ-field density parameter.
    Omega_Lambda : float
        Cosmological constant (residual after τ-field).
    tau_dot_0 : float
        Present τ̇ (kinetic energy contribution).
    V_tau_0 : float
        Present V(τ) (potential energy contribution).
    theta_T_amplitude : float
        Amplitude of temporal axion oscillation in G(t).
    theta_T_omega : float
        Angular frequency of G-oscillation (Gyr⁻¹).
    """

    H0: float = 67.4             # km/s/Mpc
    Omega_m: float = 0.315
    Omega_r: float = 9.1e-5
    Omega_tau: float = 0.25       # τ-field replaces part of CDM
    Omega_Lambda: float = 0.685 - 0.25  # residual Λ
    tau_dot_0: float = 0.0
    V_tau_0: float = 0.0
    theta_T_amplitude: float = 1e-12   # dimensionless
    theta_T_omega: float = 0.1         # Gyr⁻¹

    @property
    def H0_natural(self) -> float:
        """H₀ in GeV (for natural-unit calculations)."""
        # 1 km/s/Mpc ≈ 3.24e-20 s⁻¹ ≈ 2.13e-44 GeV
        return self.H0 * 2.133e-44

    @property
    def rho_crit(self) -> float:
        """Critical density ρ_c = 3H₀²/(8πG) in GeV⁴."""
        return 3.0 * self.H0_natural ** 2 / (8.0 * math.pi * PC.G_N)

    # ── Hubble rate ──────────────────────────────────────────
    def H_squared(self, a: float) -> float:
        """
        H²(a) in units of H₀².

        Parameters
        ----------
        a : float
            Scale factor (a=1 today).
        """
        return (
            self.Omega_r / a ** 4
            + self.Omega_m / a ** 3
            + self.Omega_tau * self._tau_scaling(a)
            + self.Omega_Lambda
        )

    def _tau_scaling(self, a: float) -> float:
        """
        τ-field density scaling with scale factor.

        For kinetic-dominated τ: ρ_τ ∝ a⁻⁶ (stiff fluid)
        For potential-dominated: ρ_τ ∝ a⁰  (cosmological constant-like)
        We interpolate with w_τ = (τ̇² − 2V) / (τ̇² + 2V)
        """
        # Simplified: assume tracking behaviour ρ_τ ∝ a^{-3(1+w)}
        w_tau = self._equation_of_state()
        return a ** (-3.0 * (1.0 + w_tau))

    def _equation_of_state(self) -> float:
        """Effective equation of state w_τ = (½τ̇² − V) / (½τ̇² + V)."""
        kinetic = 0.5 * self.tau_dot_0 ** 2
        denom = kinetic + self.V_tau_0
        if denom < 1e-30:
            return -1.0  # vacuum-like
        return (kinetic - self.V_tau_0) / denom

    # ── Oscillatory G(t) ─────────────────────────────────────
    def G_effective(self, t_Gyr: float) -> float:
        """
        Effective Newton's constant with temporal-axion modulation:

            G_eff(t) = G_N [1 + δ_T sin(ω_T t)]

        where δ_T = θ_T amplitude.

        Parameters
        ----------
        t_Gyr : float
            Cosmic time in Gyr.
        """
        modulation = self.theta_T_amplitude * math.sin(self.theta_T_omega * t_Gyr)
        return PC.G_N * (1.0 + modulation)

    # ── Perturbation theory (Newtonian gauge) ────────────────
    def perturbation_eom(
        self,
        k: float,
        a: float,
        delta_tau: float,
        delta_tau_dot: float,
        Phi: float,
        Phi_dot: float,
        V_prime: float = 0.0,
        V_double_prime: float = 0.0,
        alpha_coupling: float = 0.0,
        delta_rho_b: float = 0.0,
    ) -> float:
        """
        Linearised τ-perturbation equation in Fourier space:

            δτ̈_k + 3H δτ̇_k + (k²/a² + V'') δτ_k
                = 4 τ̇₀ Φ̇_k − 2V'Φ_k + α δρ_{b,k}

        Returns δτ̈_k.
        """
        H = math.sqrt(max(self.H_squared(a), 0)) * self.H0_natural

        ddelta_tau = (
            -3.0 * H * delta_tau_dot
            - (k ** 2 / a ** 2 + V_double_prime) * delta_tau
            + 4.0 * self.tau_dot_0 * Phi_dot
            - 2.0 * V_prime * Phi
            + alpha_coupling * delta_rho_b
        )
        return ddelta_tau

    def poisson_equation(
        self,
        k: float,
        a: float,
        delta_rho_b: float,
        delta_tau: float,
        delta_tau_dot: float,
        Phi: float,
        V_prime: float = 0.0,
    ) -> float:
        """
        Linearised Poisson equation:

            k² Φ_k = 4π G a² (δρ_b + δρ_τ)

        where δρ_τ = τ̇₀ δτ̇ + V' δτ − τ̇₀² Φ.

        Returns Φ_k.
        """
        delta_rho_tau = (
            self.tau_dot_0 * delta_tau_dot
            + V_prime * delta_tau
            - self.tau_dot_0 ** 2 * Phi
        )
        rhs = 4.0 * math.pi * PC.G_N * a ** 2 * (delta_rho_b + delta_rho_tau)
        if abs(k) < 1e-30:
            return 0.0
        return rhs / (k ** 2)

    # ── Convenience: expansion history ───────────────────────
    def expansion_history(
        self, a_min: float = 1e-4, a_max: float = 2.0, n_points: int = 500
    ) -> dict[str, np.ndarray]:
        """
        Compute H(a)/H₀ over a range of scale factors.

        Returns dict with 'a', 'H_over_H0', 'q' (deceleration parameter).
        """
        a_arr = np.linspace(a_min, a_max, n_points)
        H2 = np.array([self.H_squared(ai) for ai in a_arr])
        H_ratio = np.sqrt(np.maximum(H2, 0))

        # Deceleration parameter q ≈ −1 − d ln H / d ln a
        ln_H = np.log(np.maximum(H_ratio, 1e-30))
        ln_a = np.log(a_arr)
        dln_H = np.gradient(ln_H, ln_a)
        q = -(1.0 + dln_H)

        return {"a": a_arr, "H_over_H0": H_ratio, "q": q}
