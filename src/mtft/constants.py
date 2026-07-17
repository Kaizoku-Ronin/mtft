"""
MTFT Constants: The Arithmetic Alphabet
========================================

Every prediction is built from computable mathematical constants.
None require physical measurement — all derivable from the integers.

Tier 1 (Pillars):  π, e, γ
Tier 2 (Feigenbaum): δ, α_F
Tier 3 (Derived):  T∞, M, Ω, ζ-values

Physical observables are *derived* via MTFT arithmetic relations.
PDG values kept for comparison but NOT used as inputs.

Reference: MTFT Complete Mathematical Dictionary (March 2026)
"""
from __future__ import annotations
import math
from dataclasses import dataclass
from typing import ClassVar

PI = math.pi
E  = math.e
EULER_GAMMA = 0.5772156649015328606065120900824024310421

# ── Tier 2: Feigenbaum ───────────────────────────────────────
FEIGENBAUM_DELTA = 4.669201609102990671853203821578
FEIGENBAUM_ALPHA = 2.502907875095892822283902873218
DELTA_X = FEIGENBAUM_DELTA ** 1.29          # ≈7.30 attractor fold (DEFINED)
DELTA_Y = FEIGENBAUM_DELTA ** 0.79          # ≈3.38 barrier fold (DEFINED)
# Measured Burning Ship values (Paper 16) — distinct from the defined
# exponentials above (audit B9: keep the two provenances separate):
DELTA_X_MEASURED = 7.2422                   # attractor fold, measured
DELTA_Y_MEASURED = 3.3699                   # barrier fold, measured
BETA_X  = +0.29
BETA_Y  = -0.21
XI      = math.sqrt(7.0 / 5.0)             # sector metric ≈1.183

# ── Tier 3: Number-theoretic ─────────────────────────────────
ZETA_PRIME_2    = -0.9375482543988
# Vacuum-torque conventions (audit B3 — the repo's canonical resolution):
#   TORQUE_FULL = −ζ′(2) = 0.9375 is the PROVED Cesàro limit
#       lim_{N→∞} (1/N) Σ_{n≤N} w_n,  w_n = Σ_{d|n} (log d)/d
#     (regression-locked in tests/test_x0_143_verified.py).
#   T_INF = TORQUE_FULL/2 is the DEFINITION used by the physics chain
#     (α_s = T_INF/4 = −ζ′(2)/8).  AG Df 0.5.5 conflates the two;
#     the factor of 2 lives here, deliberately, in one place.
T_INF           = -ZETA_PRIME_2 / 2.0      # vacuum torque ≈0.4688 (defined)
TORQUE_FULL    = -ZETA_PRIME_2                  # full torque = −ζ'(2) ≈ 0.9375 (Cesàro limit)
# ln|M| of the Monster group (audit B2: AG Df 10.2.1 misprints 124.01348).
LN_MONSTER      = 124.126423366            # ln(8.08017...×10⁵³)
MEISSEL_MERTENS = 0.2614972128476427837554268386086958590516
LAMBERT_OMEGA   = 0.5671432904097838729999686622103555497538
ZETA_2 = PI**2 / 6.0
ZETA_3 = 1.2020569031595942853997381
ZETA_4 = PI**4 / 90.0
GLAISHER_A = 1.2824271291006226368753425688697917
R_ACCUMULATION = 3.5699456718709449018420051513864989

# ── MTFT Gauge Sector (Paper 7, Dictionary IV) ───────────────
class _MTFTGauge:
    """Gauge couplings — zero free parameters."""
    @property
    def alpha_inv(self) -> float:
        """α⁻¹ = 2πδ² + 1/(4δ) − ξT∞/δ⁶"""
        d = FEIGENBAUM_DELTA
        return 2*PI*d**2 + 1/(4*d) - XI*T_INF/d**6

    @property
    def alpha_inv_4term(self) -> float:
        """α⁻¹ = 2πδ² + 1/(4δ) − ξT∞/δ⁶ − y_c/δ⁹  (137.03599917)."""
        d = FEIGENBAUM_DELTA
        return self.alpha_inv - CriticalDepths.y_conf / d**9

    @property
    def alpha_inv_monster(self) -> float:
        """α⁻¹ = ln|M| + genus − 1/dim S₂^new = 137.0355 (3.5 ppm, parameter-free)."""
        return LN_MONSTER + 13.0 - 1.0/11.0

    @property
    def alpha(self) -> float:
        return 1.0 / self.alpha_inv

    @property
    def e_charge_gaussian(self) -> float:
        """e = √(2/δ) in Gaussian units (0.02%)"""
        return math.sqrt(2.0) / FEIGENBAUM_DELTA

    @property
    def alpha_s(self) -> float:
        """α_s = T∞/4 (0.7σ from PDG)"""
        return T_INF / 4.0

    @property
    def alpha_s_from_13(self) -> float:
        """α_s = 13^(−5/6) (0.046%)"""
        return 13.0 ** (-5.0/6.0)

    @property
    def sin2_theta_W(self) -> float:
        """sin²θ_W = 3/13 = (1/4)·φ(13)/13"""
        return 3.0 / 13.0

    @property
    def cos2_theta_W(self) -> float:
        return 10.0 / 13.0

    @property
    def cos_theta_W(self) -> float:
        return math.sqrt(10.0/13.0)

    @property
    def sin_theta_W(self) -> float:
        return math.sqrt(3.0/13.0)

    @property
    def W_Z_ratio(self) -> float:
        """M_W/M_Z = 1/(2Ω) (0.024%)"""
        return 1.0 / (2.0 * LAMBERT_OMEGA)

GAUGE = _MTFTGauge()

# ── MTFT Higgs Sector (Dictionary V) ─────────────────────────
class _MTFTHiggs:
    """Higgs from γ, Ω, T∞."""
    V_EW = 246.22  # GeV — single dimensional anchor
    @property
    def m_H(self) -> float:
        """m_H = v·γ/(2Ω) (0.037%)"""
        return self.V_EW * EULER_GAMMA / (2*LAMBERT_OMEGA)

    @property
    def m_H_over_m_W(self) -> float:
        """m_H/m_W = √π·T∞² (0.019%)"""
        return math.sqrt(PI) * TORQUE_FULL**2

    @property
    def lambda_quartic(self) -> float:
        """λ = (γ/Ω)²/2 (0.074%)"""
        return (EULER_GAMMA / LAMBERT_OMEGA)**2 / 2.0

HIGGS = _MTFTHiggs()

# ── SM from MTFT arithmetic + VEV anchor ─────────────────────
class _MTFTSM:
    v_ew: ClassVar[float] = 246.22

    @property
    def m_W(self) -> float:
        e_c = math.sqrt(4*PI*GAUGE.alpha)
        g2 = e_c / GAUGE.sin_theta_W
        return self.v_ew * g2 / 2.0

    @property
    def m_Z(self) -> float:
        return self.m_W / GAUGE.cos_theta_W

    @property
    def m_H(self) -> float:
        return HIGGS.m_H

    @property
    def G_F(self) -> float:
        return 1.0 / (math.sqrt(2) * self.v_ew**2)

    @property
    def g_weak(self) -> float:
        return 2.0 * self.m_W / self.v_ew

    @property
    def g_prime(self) -> float:
        return self.g_weak * GAUGE.sin_theta_W / GAUGE.cos_theta_W

    @property
    def cos_theta_W(self) -> float:
        """cos θ_W = √(10/13) from MTFT gauge sector."""
        return GAUGE.cos_theta_W

    @property
    def sin_theta_W(self) -> float:
        """sin θ_W = √(3/13) from MTFT gauge sector."""
        return GAUGE.sin_theta_W

    @property
    def sin2_theta_W(self) -> float:
        """sin²θ_W = 3/13 from MTFT gauge sector."""
        return GAUGE.sin2_theta_W

SM = _MTFTSM()

# ── PDG Reference (verification only, NOT inputs) ────────────
@dataclass(frozen=True)
class _PDGReference:
    m_W: float = 80.3692
    m_Z: float = 91.1876
    m_H: float = 125.25
    v_ew: float = 246.22
    sin2_theta_W: float = 0.23122
    alpha_inv: float = 137.035999084
    alpha_s_mZ: float = 0.1180
    G_F: float = 1.1663787e-5

PDG = _PDGReference()

# ── Physical constants (unit conversions) ─────────────────────
class PhysicalConstants:
    G_N: ClassVar[float] = 6.70883e-39
    M_Pl: ClassVar[float] = 1.22089e19
    M_Pl_reduced: ClassVar[float] = 2.435e18
    k_B: ClassVar[float] = 8.617333e-14
    c_SI: ClassVar[float] = 2.99792458e8
    hbar_SI: ClassVar[float] = 6.582119569e-25
    GeV_inv_to_m: ClassVar[float] = 1.9732698e-16
    M_sun_GeV: ClassVar[float] = 1.116e57
    kpc_GeV_inv: ClassVar[float] = 1.5637e35

# ── Critical Depths (Dictionary XIII) ─────────────────────────
class CriticalDepths:
    y_s1: ClassVar[float] = 0.1236
    y_conf: ClassVar[float] = 0.18174
    y_s2: ClassVar[float] = 0.2106
    y_spec: ClassVar[float] = 0.2710
    @staticmethod
    def y_conf_from_arithmetic() -> float:
        """y_c ≈ (γ + Ω)/(2π) (0.072%)"""
        return (EULER_GAMMA + LAMBERT_OMEGA) / (2*PI)

# ── Fermion Masses (PDG + MTFT scaling notes) ─────────────────
@dataclass(frozen=True)
class QuarkMasses:
    u: float = 0.00216
    d: float = 0.00467
    s: float = 0.0934
    c: float = 1.27
    b: float = 4.18
    t: float = 172.69

@dataclass(frozen=True)
class LeptonMasses:
    e: float = 0.000511
    mu: float = 0.10566
    tau: float = 1.7768
    delta_m21_sq: float = 7.53e-5
    delta_m32_sq: float = 2.453e-3

QUARKS = QuarkMasses()
LEPTONS = LeptonMasses()

# ── Mass Ratio Predictions (Dictionary VI) ────────────────────
class MassRatios:
    @staticmethod
    def m_d_over_m_u() -> float:
        """√δ ≈ 2.161"""
        return math.sqrt(FEIGENBAUM_DELTA)
    @staticmethod
    def m_s_over_m_mu() -> float:
        """T∞² ≈ 0.884"""
        return TORQUE_FULL**2
    @staticmethod
    def m_c_over_m_mu() -> float:
        """φ(13) = 12"""
        return 12.0
    @staticmethod
    def m_b_over_m_tau() -> float:
        """∛13 ≈ 2.351 (0.048%)"""
        return 13.0**(1/3)
    @staticmethod
    def m_tau_over_m_mu() -> float:
        """δ_y^2.33 ≈ 16.8"""
        return DELTA_Y**2.33
    @staticmethod
    def m_H_over_m_W() -> float:
        """√π·T∞² (0.019%)"""
        return math.sqrt(PI) * TORQUE_FULL**2
