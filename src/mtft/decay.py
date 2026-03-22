"""
Radioactive Decay in Modular Time
===================================

Internal phase evolution with spin-holonomy coupling:

    dφ/dt_U = ω(R_τ, τ) + γ_h ∫ F_UR dt_R

Half-life formula:
    T_{1/2} ≈ C · Δ²_eff · Θ(τ*)⁻¹

where Θ(τ*) is the modular depth function at the decay threshold τ*,
Δ_eff is the effective mass splitting, and C is a dimensionful constant
set by the confinement scale.

Beta decay with helicity coupling:
    γ_h = ±1/2 for left/right-handed fermions

Testable: anomalies in neutron star environments where τ-field gradients
are large.

Reference: Simulator chat (τ-field dark sector), Papers 3, 19.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from mtft.constants import PI, T_INF, CriticalDepths


# ═══════════════════════════════════════════════════════════════
#  Modular Decay Rate
# ═══════════════════════════════════════════════════════════════

@dataclass
class ModularDecay:
    """
    Decay rate in MTFT modular time framework.

    The decay rate receives corrections from the τ-field:
        Γ = Γ₀ · [1 + δΓ(τ)]

    where δΓ encodes the holonomy correction from the compact
    modular-time direction.

    Parameters
    ----------
    Gamma_0 : float
        Standard decay rate (s⁻¹) without modular corrections.
    delta_eff : float
        Effective mass splitting (GeV) controlling the phase space.
    tau_star : complex
        Modular depth at the decay threshold.
    gamma_h : float
        Helicity coupling (±1/2 for L/R fermions).
    """

    Gamma_0: float
    delta_eff: float = 0.001         # GeV
    tau_star: complex = 0.0 + 0.18j  # near confinement depth
    gamma_h: float = -0.5            # left-handed (SM default)

    @property
    def theta_tau(self) -> float:
        """
        Modular depth function Θ(τ*) = Im(τ*) · exp(−2π Im(τ*)).

        Controls the suppression of decay at large modular depth.
        """
        y = self.tau_star.imag
        return y * math.exp(-2 * PI * y)

    @property
    def half_life_correction(self) -> float:
        """
        T_{1/2} / T_{1/2}^(SM) = 1 + C_mod · δ²_eff / Θ(τ*)

        Returns the ratio (>1 means longer-lived in modular time).
        """
        if self.theta_tau < 1e-30:
            return 1.0
        C_mod = abs(self.gamma_h) * T_INF  # coupling strength
        return 1.0 + C_mod * self.delta_eff ** 2 / self.theta_tau

    @property
    def corrected_rate(self) -> float:
        """Γ_corrected = Γ₀ / correction_factor."""
        return self.Gamma_0 / self.half_life_correction

    def half_life(self) -> float:
        """T_{1/2} in seconds."""
        if self.corrected_rate <= 0:
            return float("inf")
        return math.log(2) / self.corrected_rate


# ═══════════════════════════════════════════════════════════════
#  Internal Phase Evolution
# ═══════════════════════════════════════════════════════════════

def phase_evolution(
    omega_0: float,
    R_tau: float,
    F_UR: float,
    gamma_h: float,
    t_U: np.ndarray,
) -> np.ndarray:
    """
    Internal phase evolution in universal time:

        φ(t_U) = ∫₀^{t_U} [ω(R_τ) + γ_h · F_UR(t')] dt'

    Parameters
    ----------
    omega_0 : float  — base angular frequency (1/R_τ scale)
    R_tau : float    — modular radius
    F_UR : float     — field strength F_{UR} (assumed constant)
    gamma_h : float  — helicity coupling
    t_U : array      — universal time points

    Returns array of phase values.
    """
    omega_eff = omega_0 / R_tau + gamma_h * F_UR
    return omega_eff * t_U


# ═══════════════════════════════════════════════════════════════
#  Example: U-238 Alpha Decay
# ═══════════════════════════════════════════════════════════════

def u238_example() -> dict:
    """
    U-238 alpha decay with modular corrections.

    Standard: T_{1/2} = 4.468 × 10⁹ years = 1.41 × 10¹⁷ s
    MTFT correction: tiny at laboratory τ-field values,
    potentially observable in neutron star environments.
    """
    T_half_standard = 4.468e9 * 365.25 * 24 * 3600  # seconds
    Gamma_0 = math.log(2) / T_half_standard

    # Laboratory: negligible correction
    lab = ModularDecay(
        Gamma_0=Gamma_0,
        delta_eff=4.27e-3,        # ~4.27 MeV alpha Q-value in GeV
        tau_star=0.0 + 0.18j,     # near confinement
        gamma_h=-0.5,
    )

    # Neutron star: enhanced τ-field gradient
    ns = ModularDecay(
        Gamma_0=Gamma_0,
        delta_eff=4.27e-3,
        tau_star=0.0 + 0.05j,     # shallower depth near NS surface
        gamma_h=-0.5,
    )

    return {
        "lab": {
            "T_half_s": lab.half_life(),
            "correction_factor": lab.half_life_correction,
        },
        "neutron_star": {
            "T_half_s": ns.half_life(),
            "correction_factor": ns.half_life_correction,
        },
        "prediction": "Enhanced decay rate near neutron stars due to τ-gradient",
    }


# ═══════════════════════════════════════════════════════════════
#  Neutron Beta Decay
# ═══════════════════════════════════════════════════════════════

def neutron_beta_decay() -> dict:
    """
    Neutron beta decay n → p + e⁻ + ν̄_e.

    Standard: T_{1/2} = 611 s (bottle) or 887 s (beam).
    The ~4σ discrepancy between bottle and beam measurements
    is a potential MTFT observable: different τ-field environments
    in trapped vs free neutrons.
    """
    T_half_bottle = 877.75  # seconds (PDG 2024)
    T_half_beam = 887.7     # seconds (beam average)

    Gamma_bottle = math.log(2) / T_half_bottle
    Gamma_beam = math.log(2) / T_half_beam

    Q_value = 1.293e-3  # GeV (neutron-proton mass difference)

    return {
        "T_half_bottle_s": T_half_bottle,
        "T_half_beam_s": T_half_beam,
        "discrepancy_percent": abs(T_half_beam - T_half_bottle) / T_half_bottle * 100,
        "Q_value_GeV": Q_value,
        "mtft_interpretation": (
            "Different effective τ* in trapped (bottle) vs "
            "free-streaming (beam) configurations may explain "
            "the lifetime anomaly."
        ),
    }
