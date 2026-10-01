"""Independent massive Born e- e+ -> gamma gamma reference formulas.

This module contains no gamma matrices, spinors, or imports from the engine.
It is an analytic oracle, not an MTFT derivation of the electron mass or alpha.
Natural units hbar=c=1; energies GeV; cross sections GeV**-2.
All cross sections count EVENTS and include identical final photons' 1/2!.
"""
from __future__ import annotations

import math

MASS_GEV = 0.00051099895069
ALPHA = 1.0 / 137.035999177
GEV2_TO_NB = 389379.3721
BETA_GRID = (0.05, 0.2, 0.5, 0.8, 0.95, 0.99)
COSINE_GRID = (-0.999, -0.95, -0.5, 0.0, 0.37, 0.8, 0.99)


def _inputs(beta: float, mass: float, alpha: float) -> None:
    if not math.isfinite(beta) or not 0.0 <= beta < 1.0:
        raise ValueError("CM beta must satisfy 0 <= beta < 1")
    if not math.isfinite(mass) or mass <= 0.0:
        raise ValueError("mass must be positive and finite")
    if not math.isfinite(alpha) or alpha <= 0.0:
        raise ValueError("alpha must be positive and finite")


def mandelstam(beta: float, cosine: float, mass: float = MASS_GEV) -> tuple[float, float, float]:
    """s,t,u from the electron-to-first-photon CM polar angle."""
    _inputs(beta, mass, ALPHA)
    if not math.isfinite(cosine) or not -1.0 <= cosine <= 1.0:
        raise ValueError("cos(theta) must satisfy -1 <= cosine <= 1")
    s = 4.0 * mass**2 / (1.0 - beta**2)
    t = mass**2 - 0.5 * s * (1.0 - beta * cosine)
    u = mass**2 - 0.5 * s * (1.0 + beta * cosine)
    return s, t, u


def msq_mandelstam(s: float, t: float, u: float, mass: float = MASS_GEV, alpha: float = ALPHA) -> float:
    """Spin averaged and physical-polarization summed Born matrix element."""
    if not all(math.isfinite(x) for x in (s, t, u, mass, alpha)) or mass <= 0.0 or alpha <= 0.0 or s < 4.0 * mass**2:
        raise ValueError("require mass>0, alpha>0, s>=4m^2")
    scale = max(abs(s), abs(t), abs(u), mass**2)
    if abs(s + t + u - 2.0 * mass**2) > 1e-12 * scale:
        raise ValueError("on-shell identity s+t+u=2m^2 is not satisfied")
    shifted_t, shifted_u = t - mass**2, u - mass**2
    if shifted_t >= 0.0 or shifted_u >= 0.0:
        raise ValueError("physical annihilation requires t-m^2<0 and u-m^2<0")
    e4 = (4.0 * math.pi * alpha)**2
    return 2.0 * e4 * (
        shifted_t / shifted_u + shifted_u / shifted_t
        + 4.0 * mass**2 * s / (shifted_t * shifted_u)
        - 4.0 * mass**4 * (1.0 / shifted_t + 1.0 / shifted_u)**2
    )


def msq_cm(beta: float, cosine: float, alpha: float = ALPHA) -> float:
    """Separate CM expression: does not call the Mandelstam implementation."""
    _inputs(beta, MASS_GEV, alpha)
    if not math.isfinite(cosine) or not -1.0 <= cosine <= 1.0:
        raise ValueError("cos(theta) must satisfy -1 <= cosine <= 1")
    b2 = beta**2
    denominator = 1.0 - b2 * cosine**2
    e4 = (4.0 * math.pi * alpha)**2
    return 4.0 * e4 * (
        (3.0 - 2.0 * b2 + b2 * cosine**2) / denominator
        - 2.0 * (1.0 - b2)**2 / denominator**2
    )


def dsigma_domega(beta: float, cosine: float, mass: float = MASS_GEV, alpha: float = ALPHA) -> float:
    """Full-sphere event dσ/dΩ; 1/2! is already included."""
    _inputs(beta, mass, alpha)
    if beta == 0.0:
        raise ValueError("cross section is singular at zero incident flux; use sigma_v_moller limit")
    s = 4.0 * mass**2 / (1.0 - beta**2)
    return msq_cm(beta, cosine, alpha) / (128.0 * math.pi**2 * s * beta)


def dsigma_dcosine(beta: float, cosine: float, mass: float = MASS_GEV, alpha: float = ALPHA) -> float:
    return 2.0 * math.pi * dsigma_domega(beta, cosine, mass, alpha)


def _total_bracket(beta: float) -> float:
    """Stable B=(3-beta^4)atanh(beta)/beta-(2-beta^2)."""
    if beta < 0.01:
        x = beta**2
        # The exact Taylor coefficients through beta^12, evaluated by Horner.
        coefficients = (1.0, 2.0, -2.0/5.0, 2.0/21.0, 2.0/15.0, 10.0/77.0, 14.0/117.0)
        value = coefficients[-1]
        for coefficient in reversed(coefficients[:-1]):
            value = coefficient + x * value
        return value
    return (3.0 - beta**4) * math.atanh(beta) / beta - (2.0 - beta**2)


def sigma_total(beta: float, mass: float = MASS_GEV, alpha: float = ALPHA) -> float:
    """Analytic event total over the full photon direction sphere."""
    _inputs(beta, mass, alpha)
    if beta == 0.0:
        return math.inf
    return math.pi * alpha**2 * (1.0-beta**2) / (2.0 * mass**2 * beta) * _total_bracket(beta)


def sigma_v_moller(beta: float, mass: float = MASS_GEV, alpha: float = ALPHA) -> float:
    """σ v_Møller in CM, v_Møller=2β; finite at threshold."""
    _inputs(beta, mass, alpha)
    return math.pi * alpha**2 * (1.0-beta**2) / mass**2 * _total_bracket(beta)


def sigma_total_heitler(beta: float, mass: float = MASS_GEV, alpha: float = ALPHA) -> float:
    """Independent target-rest Heitler formula evaluated at the same invariant s.

    Geant4 Physics Reference Manual, Positron-Electron Annihilation, per electron.
    Gamma_lab=(1+beta_CM^2)/(1-beta_CM^2), r_e=alpha/m in natural units.
    This expression is deliberately not rewritten through the CM total bracket.
    """
    _inputs(beta, mass, alpha)
    if beta == 0.0:
        return math.inf
    gamma = (1.0 + beta**2) / (1.0 - beta**2)
    radical = math.sqrt(gamma**2 - 1.0)
    return math.pi * (alpha / mass)**2 / (gamma + 1.0) * (
        (gamma**2 + 4.0*gamma + 1.0) / (gamma**2 - 1.0) * math.acosh(gamma)
        - (gamma + 3.0) / radical
    )


def sigma_high_energy(beta: float, mass: float = MASS_GEV, alpha: float = ALPHA) -> float:
    """Leading s>>m² asymptotic event total; not a finite-energy replacement."""
    _inputs(beta, mass, alpha)
    s = 4.0 * mass**2 / (1.0 - beta**2)
    return 2.0 * math.pi * alpha**2 / s * (math.log(s/mass**2) - 1.0)
