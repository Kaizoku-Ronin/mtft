"""
The Dimensional Bridge: Electron Mass from η(τ)
=================================================

The formula (Paper 22, Theorem):

    m_e / m_P = |η(τ_c^(e))|^{2α⁻¹/π}

connects the Dedekind eta function at the Euler-shifted confinement
depth to the electron-to-Planck mass ratio.

Equivalently:
    ln(m_P/m_e) = (2/πα) |ln η(τ_c)|

The chain: -1/12 → η → Δ → m_e:
    1. Ramanujan: 1+2+3+... = ζ(-1) = -1/12
    2. Casimir: E₀ = 1/24
    3. η(τ) = q^{1/24} Π(1-qⁿ)
    4. Δ(τ) = η(τ)²⁴
    5. m_e/m_P = |η(τ_c)|^{2α⁻¹/π}

Four Feigenbaum-electromagnetism identities:
    (i)   e = √(2/δ)         — charge from Feigenbaum (0.02%)
    (ii)  α⁻¹ = 2πδ²         — fine structure (leading, 0.04%)
    (iii) α_F = exp(40πα)     — spatial Feigenbaum (0.04%)
    (iv)  (δ/α_F)² ≈ 7/2     — modular impedance (0.57%)

Reference: Paper 22, Paper 20 (Euler shift).
"""

from __future__ import annotations

import cmath
import math

from mtft.constants import (
    FEIGENBAUM_DELTA, FEIGENBAUM_ALPHA, GAUGE, PhysicalConstants as PC,
    PI, CriticalDepths,
)
from mtft.forms import dedekind_eta


# ═══════════════════════════════════════════════════════════════
#  Feigenbaum–Electromagnetism Identities (Paper 22 §2)
# ═══════════════════════════════════════════════════════════════

def charge_from_feigenbaum() -> dict:
    """
    Identity (i): e = √(2/δ) in Gaussian natural units.

    Verification: √(4πα) vs √(2)/δ — should match to 0.02%.
    """
    e_exact = math.sqrt(4 * PI * GAUGE.alpha)
    e_feig = math.sqrt(2.0) / FEIGENBAUM_DELTA
    return {
        "e_exact": e_exact,
        "e_feigenbaum": e_feig,
        "error_percent": abs(e_exact - e_feig) / e_exact * 100,
    }


def alpha_F_from_impedance() -> dict:
    """
    Identity (iii): α_F = exp(40πα).

    The spatial Feigenbaum constant is the exponential of the
    vacuum impedance summed over 10 metric components.
    """
    alpha_F_pred = math.exp(40 * PI * GAUGE.alpha)
    return {
        "alpha_F_predicted": alpha_F_pred,
        "alpha_F_actual": FEIGENBAUM_ALPHA,
        "error_percent": abs(alpha_F_pred - FEIGENBAUM_ALPHA) / FEIGENBAUM_ALPHA * 100,
    }


def feigenbaum_product_lock() -> dict:
    """
    Conjecture: δ² · ln(α_F) = 20 (0.009% deviation).

    If exact, knowing either Feigenbaum constant determines the other.
    """
    product = FEIGENBAUM_DELTA ** 2 * math.log(FEIGENBAUM_ALPHA)
    return {
        "delta2_ln_alphaF": product,
        "target": 20.0,
        "error_percent": abs(product - 20.0) / 20.0 * 100,
    }


def modular_impedance() -> dict:
    """
    Identity (iv): (δ/α_F)² ≈ 7/2 (0.57%).

    The ratio of temporal to spatial Feigenbaum scaling.
    """
    ratio_sq = (FEIGENBAUM_DELTA / FEIGENBAUM_ALPHA) ** 2
    return {
        "ratio_squared": ratio_sq,
        "target": 3.5,
        "error_percent": abs(ratio_sq - 3.5) / 3.5 * 100,
    }


# ═══════════════════════════════════════════════════════════════
#  Electron Mass Formula (Paper 22, main result)
# ═══════════════════════════════════════════════════════════════

def electron_mass_from_eta(
    y_c: float | None = None,
    alpha_inv: float | None = None,
) -> dict:
    """
    The dimensional bridge formula:

        m_e / m_P = |η(i·y_c^(e))|^{2α⁻¹/π}

    Equivalently:
        ln(m_P/m_e) = (2α⁻¹/π) · |ln|η(τ_c)||

    Parameters
    ----------
    y_c : float
        Euler-shifted confinement depth (default: 0.18174)
    alpha_inv : float
        Fine-structure constant inverse (default: MTFT value)

    Returns dict with prediction vs measurement.
    """
    if y_c is None:
        y_c = CriticalDepths.y_conf
    if alpha_inv is None:
        alpha_inv = GAUGE.alpha_inv

    tau_c = 1j * y_c
    eta_val = dedekind_eta(tau_c)
    abs_eta = abs(eta_val)
    ln_abs_eta = math.log(abs_eta) if abs_eta > 0 else float("-inf")

    exponent = 2.0 * alpha_inv / PI
    ratio_predicted = abs_eta ** exponent   # m_e / m_P

    ln_ratio_predicted = exponent * ln_abs_eta
    ln_ratio_measured = -math.log(PC.M_Pl / (0.000511 * 1e-3 * 5.609588e26))
    # More directly: m_e = 0.511 MeV, m_P = 1.22089e19 GeV
    ln_mP_over_me = math.log(PC.M_Pl / 0.000511)  # ≈ 51.53

    ln_predicted = exponent * abs(ln_abs_eta)

    return {
        "y_c": y_c,
        "eta_at_tau_c": eta_val,
        "abs_eta": abs_eta,
        "exponent_2alpha_inv_over_pi": exponent,
        "ln_mP_over_me_predicted": ln_predicted,
        "ln_mP_over_me_measured": ln_mP_over_me,
        "error_percent": abs(ln_predicted - ln_mP_over_me) / ln_mP_over_me * 100,
    }


# ═══════════════════════════════════════════════════════════════
#  Verification at Multiple Critical Depths (Paper 22 §end)
# ═══════════════════════════════════════════════════════════════

def verify_at_all_depths() -> list[dict]:
    """
    Check the eta formula at all three MTFT critical depths:
        y_s1 = 0.1236:  |Δ(y)|² ≈ m_e/m_P  (exponent 48)
        y_c  = 0.1817:  |η|^{2α⁻¹/π}       (main formula)
        y_s2 = 0.2106:  |η|^{111}           (integer exponent)
    """
    depths = [
        ("y_s1 (skeleton 1)", CriticalDepths.y_s1, 48.0),
        ("y_c (confinement)", CriticalDepths.y_conf, 2 * GAUGE.alpha_inv / PI),
        ("y_s2 (skeleton 2)", CriticalDepths.y_s2, 111.0),
    ]
    results = []
    for label, y, exp in depths:
        tau = 1j * y
        eta_val = dedekind_eta(tau)
        abs_eta = abs(eta_val)
        ratio = abs_eta ** exp
        me_over_mP = 0.000511 / PC.M_Pl
        results.append({
            "depth": label,
            "y": y,
            "exponent": exp,
            "|eta|": abs_eta,
            "predicted_ratio": ratio,
            "measured_ratio": me_over_mP,
            "log_ratio_pred": exp * math.log(abs_eta),
            "log_ratio_meas": math.log(me_over_mP),
        })
    return results
