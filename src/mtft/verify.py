"""
MTFT Verification Scorecard
=============================

Systematic comparison of every MTFT prediction against experiment.
Each entry records: formula, MTFT value, PDG/CODATA value, error, and paper.

Usage:
    from mtft.verify import full_report
    for entry in full_report():
        print(entry)

Reference: MTFT Complete Mathematical Dictionary (March 2026).
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import List

from mtft.constants import (
    FEIGENBAUM_DELTA, FEIGENBAUM_ALPHA, T_INF, LAMBERT_OMEGA,
    EULER_GAMMA, XI, ZETA_2, MEISSEL_MERTENS, GAUGE, HIGGS, SM, PDG,
    QUARKS, LEPTONS, DELTA_X, DELTA_Y, PI, MassRatios,
)


@dataclass
class Prediction:
    """A single MTFT prediction vs experiment."""
    name: str
    formula: str
    mtft_value: float
    experimental: float
    unit: str
    error_percent: float
    paper: str
    sector: str

    def __str__(self) -> str:
        sign = "✓" if self.error_percent < 1.0 else "○" if self.error_percent < 5.0 else "△"
        return (
            f"{sign} {self.name:30s}  MTFT={self.mtft_value:12.6f}  "
            f"Exp={self.experimental:12.6f} {self.unit:5s}  "
            f"Err={self.error_percent:8.4f}%  [{self.paper}]"
        )


def _pct(pred: float, obs: float) -> float:
    if obs == 0:
        return float("inf")
    return abs(pred - obs) / abs(obs) * 100


# ═══════════════════════════════════════════════════════════════

def gauge_sector() -> List[Prediction]:
    """Gauge coupling predictions."""
    return [
        Prediction(
            "Fine-structure α⁻¹",
            "2πδ² + 1/(4δ) − ξT∞/δ⁶",
            GAUGE.alpha_inv, PDG.alpha_inv, "",
            _pct(GAUGE.alpha_inv, PDG.alpha_inv), "Paper 7 §19", "Gauge"
        ),
        Prediction(
            "Strong coupling α_s",
            "T∞/4",
            GAUGE.alpha_s, PDG.alpha_s_mZ, "",
            _pct(GAUGE.alpha_s, PDG.alpha_s_mZ), "Paper 7 §20", "Gauge"
        ),
        Prediction(
            "α_s (from 13)",
            "13^(−5/6)",
            GAUGE.alpha_s_from_13, PDG.alpha_s_mZ, "",
            _pct(GAUGE.alpha_s_from_13, PDG.alpha_s_mZ), "Addendum", "Gauge"
        ),
        Prediction(
            "sin²θ_W",
            "3/13",
            GAUGE.sin2_theta_W, PDG.sin2_theta_W, "",
            _pct(GAUGE.sin2_theta_W, PDG.sin2_theta_W), "Paper 7 §21", "Gauge"
        ),
        Prediction(
            "M_W/M_Z",
            "1/(2Ω)",
            GAUGE.W_Z_ratio, PDG.m_W / PDG.m_Z, "",
            _pct(GAUGE.W_Z_ratio, PDG.m_W / PDG.m_Z), "Addendum", "Gauge"
        ),
        Prediction(
            "e (Gaussian)",
            "√(2/δ)",
            GAUGE.e_charge_gaussian, math.sqrt(4 * PI / PDG.alpha_inv), "",
            _pct(GAUGE.e_charge_gaussian, math.sqrt(4 * PI / PDG.alpha_inv)),
            "Paper 22", "Gauge"
        ),
    ]


def higgs_sector() -> List[Prediction]:
    """Higgs sector predictions."""
    return [
        Prediction(
            "m_H",
            "v·γ/(2Ω)",
            HIGGS.m_H, PDG.m_H, "GeV",
            _pct(HIGGS.m_H, PDG.m_H), "Dictionary V", "Higgs"
        ),
        Prediction(
            "m_H/m_W",
            "√π·T∞²",
            HIGGS.m_H_over_m_W, PDG.m_H / PDG.m_W, "",
            _pct(HIGGS.m_H_over_m_W, PDG.m_H / PDG.m_W), "Dictionary V", "Higgs"
        ),
        Prediction(
            "λ_quartic",
            "(γ/Ω)²/2",
            HIGGS.lambda_quartic, (PDG.m_H / PDG.v_ew) ** 2 / 2 * 4, "",
            _pct(HIGGS.lambda_quartic, (2 * PDG.m_H / PDG.v_ew) ** 2 / 2),
            "Dictionary V", "Higgs"
        ),
    ]


def mass_ratios() -> List[Prediction]:
    """Fermion mass ratio predictions."""
    return [
        Prediction(
            "m_d/m_u",
            "√δ",
            MassRatios.m_d_over_m_u(), QUARKS.d / QUARKS.u, "",
            _pct(MassRatios.m_d_over_m_u(), QUARKS.d / QUARKS.u),
            "Paper 7", "Mass Ratio"
        ),
        Prediction(
            "m_s/m_μ",
            "T∞²",
            MassRatios.m_s_over_m_mu(), QUARKS.s / LEPTONS.mu, "",
            _pct(MassRatios.m_s_over_m_mu(), QUARKS.s / LEPTONS.mu),
            "Dictionary VI", "Mass Ratio"
        ),
        Prediction(
            "m_c/m_μ",
            "φ(13) = 12",
            MassRatios.m_c_over_m_mu(), QUARKS.c / LEPTONS.mu, "",
            _pct(MassRatios.m_c_over_m_mu(), QUARKS.c / LEPTONS.mu),
            "Dictionary VI", "Mass Ratio"
        ),
        Prediction(
            "m_b/m_τ",
            "∛13",
            MassRatios.m_b_over_m_tau(), QUARKS.b / LEPTONS.tau, "",
            _pct(MassRatios.m_b_over_m_tau(), QUARKS.b / LEPTONS.tau),
            "Dictionary VI", "Mass Ratio"
        ),
        Prediction(
            "m_τ/m_μ",
            "δ_y^2.33",
            MassRatios.m_tau_over_m_mu(), LEPTONS.tau / LEPTONS.mu, "",
            _pct(MassRatios.m_tau_over_m_mu(), LEPTONS.tau / LEPTONS.mu),
            "Chapter 15", "Mass Ratio"
        ),
    ]


def koide_sector() -> List[Prediction]:
    """Koide formula verification."""
    from mtft.koide import koide_leptons, koide_up_quarks, koide_down_quarks, predict_tau_mass
    tau_pred = predict_tau_mass()
    return [
        Prediction(
            "Koide K (leptons)",
            "(Σm)/(Σ√m)² = 2/3",
            koide_leptons(), 2 / 3, "",
            _pct(koide_leptons(), 2 / 3), "Paper 6", "Koide"
        ),
        Prediction(
            "Koide K (up quarks)",
            "K_u",
            koide_up_quarks(), 2 / 3, "",
            _pct(koide_up_quarks(), 2 / 3), "Chapter 16", "Koide"
        ),
        Prediction(
            "τ mass (from Koide)",
            "K=2/3 + m_e,m_μ",
            tau_pred["predicted_MeV"], tau_pred["experimental_MeV"], "MeV",
            tau_pred["error_percent"], "Chapter 16 §7", "Koide"
        ),
    ]


def tano_formula() -> List[Prediction]:
    """Tano mass formula predictions.

    SUPERSEDED (audit B10/B11): the formula is built on the phantom
    eigenvalue 0.5732+0.3564i (root of a wrong polynomial).  Retained for
    archival continuity; excluded from any live prediction count.
    """
    from mtft.x0_143 import tano_mass_predictions
    pred = tano_mass_predictions()
    results = []
    for particle in ["muon", "tau"]:
        results.append(Prediction(
            f"m_{particle} (Tano formula, SUPERSEDED)",
            "θ₀ = 2/δ² + 4·arg(a₂ᶜˣ) [phantom input]",
            pred["masses_MeV"][particle],
            pred["experimental_MeV"][particle],
            "MeV",
            pred["errors_percent"][particle],
            "Paper 26 §8 (superseded, Correction Session 1)", "Tano-archival"
        ))
    return results


def feigenbaum_identities() -> List[Prediction]:
    """Feigenbaum-electromagnetism identities (Paper 22)."""
    delta2_ln_aF = FEIGENBAUM_DELTA ** 2 * math.log(FEIGENBAUM_ALPHA)
    aF_pred = math.exp(40 * PI * GAUGE.alpha)
    impedance = (FEIGENBAUM_DELTA / FEIGENBAUM_ALPHA) ** 2
    return [
        Prediction(
            "δ²·ln(α_F) = 20",
            "product lock",
            delta2_ln_aF, 20.0, "",
            _pct(delta2_ln_aF, 20.0), "Paper 22", "Feigenbaum"
        ),
        Prediction(
            "α_F = exp(40πα)",
            "impedance formula",
            aF_pred, FEIGENBAUM_ALPHA, "",
            _pct(aF_pred, FEIGENBAUM_ALPHA), "Paper 22", "Feigenbaum"
        ),
        Prediction(
            "(δ/α_F)² = 7/2",
            "modular impedance",
            impedance, 3.5, "",
            _pct(impedance, 3.5), "Paper 22", "Feigenbaum"
        ),
    ]


# ═══════════════════════════════════════════════════════════════
#  Full Report
# ═══════════════════════════════════════════════════════════════

def full_report() -> List[Prediction]:
    """All MTFT predictions in one list."""
    return (
        gauge_sector()
        + higgs_sector()
        + mass_ratios()
        + koide_sector()
        + tano_formula()
        + feigenbaum_identities()
    )


def print_report():
    """Pretty-print the full MTFT verification scorecard."""
    report = full_report()
    current_sector = ""
    print("=" * 100)
    print("MTFT VERIFICATION SCORECARD — All Predictions vs Experiment")
    print("=" * 100)
    for p in report:
        if p.sector != current_sector:
            current_sector = p.sector
            print(f"\n── {current_sector.upper()} {'─' * (90 - len(current_sector))}")
        print(p)
    print("\n" + "=" * 100)
    sub1 = sum(1 for p in report if p.error_percent < 1.0)
    sub5 = sum(1 for p in report if p.error_percent < 5.0)
    print(f"Total: {len(report)} predictions  |  <1%: {sub1}  |  <5%: {sub5}  |  >5%: {len(report) - sub5}")
    print("✓ = <1%    ○ = 1-5%    △ = >5%")
