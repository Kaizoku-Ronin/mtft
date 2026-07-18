"""
The Modular Curve X₀(143): Arithmetic Stage of MTFT
=====================================================

X₀(143) is the compactification of ℍ/Γ₀(143), where 143 = 11 × 13.
It is the self-referential fixed point of MTFT:

    genus = 13             (cusp form dimension)
    index = 168 = |PSL(2,7)| = dim SU(13)
    newform dim = 11       (S₂^new)
    Galois orbits = 3      with dimensions [1, 4, 6]

The Atkin-Lehner involutions W₁₁, W₁₃ form ℤ₂ × ℤ₂, identical to the
parity-time orbifold mapping Mandelbrot → Burning Ship.

SUPERSESSION NOTICE (independent audit + Correction Sessions 1–6, July 2026)
----------------------------------------------------------------------------
The Tano Mass Formula (Paper 26, Theorem 8.1),

    θ₀ = 2/δ² + dim(f₂) · arg(a₂ᶜˣ(f₃)),

was built on the polynomial x⁶−x⁵−9x⁴+11x³+13x²−20x+8, which is NOT the
characteristic polynomial of T₂ on f₃.  The correct charpoly (PARI
mfeigenbasis; independently verified by orbit traces at a₂, a₄, a₉) is

    x⁶ − 10x⁴ + 2x³ + 24x² − 7x − 12,

whose 6 roots are ALL REAL (as they must be: trivial nebentypus ⇒ every
Hecke eigenvalue on S₂^new is totally real).  The "complex eigenvalue"
a₂ᶜˣ = 0.5732 + 0.3564i is a root of the wrong polynomial only — a phantom.
The formula is retained below, clearly marked SUPERSEDED, for historical
continuity; do not use it in new work.

All other data LMFDB-verified; per-orbit traces independently verified by
curve point-counts and companion-matrix evaluation of the PARI polmods.
Reference: Papers 25, 26, 29, 30; Correction Sessions 1–6 (2026).
"""

from __future__ import annotations

import cmath
import math
from dataclasses import dataclass
from typing import List, Tuple

import numpy as np

from mtft.constants import FEIGENBAUM_DELTA, LEPTONS


# ═══════════════════════════════════════════════════════════════
#  Structural Constants (LMFDB-verified, Table 1 of Paper 26)
# ═══════════════════════════════════════════════════════════════

LEVEL = 143                          # = 11 × 13
GENUS = 13                           # cusp form dimension
INDEX = 168                          # = |PSL(2,7)| = dim SU(13)
NEW_DIMENSION = 11                   # dim S₂^new(Γ₀(143))
EISENSTEIN_DIMENSION = 3
NUM_NEWFORMS = 3                     # Galois orbits
ORBIT_DIMENSIONS = (1, 4, 6)         # f₁, f₂, f₃

# Atkin-Lehner eigenspace decomposition (Table 2)
# (w₁₁, w₁₃) → dimension
AL_DECOMPOSITION = {
    (+1, +1): 1,    # f₁ — electron sector
    (-1, +1): 4,    # f₂ — muon sector
    (+1, -1): 6,    # f₃ — tau sector
    (-1, -1): 0,    # FORBIDDEN — why exactly 3 generations
}

# Root numbers
ROOT_NUMBERS = {
    "f1": -1,    # ε₁ = −1 → L(f₁, 1) = 0 (forced vanishing)
    "f2": +1,    # ε₂ = +1 → L(f₂, 1) ≠ 0
    "f3": +1,    # ε₃ = +1 → L(f₃, 1) ≠ 0
}


# ═══════════════════════════════════════════════════════════════
#  Elliptic Curve 143a1 (rational newform f₁)
# ═══════════════════════════════════════════════════════════════

@dataclass(frozen=True)
class EllipticCurve143a1:
    """The elliptic curve corresponding to f₁ = 143.2.a.a."""
    label: str = "143a1"
    a_invariants: tuple = (0, -1, 1, -1, -2)     # y² + y = x³ − x² − x − 2
    conductor: int = 143
    j_invariant: str = "-1/15"
    analytic_rank: int = 1     # ε = −1 forces odd rank; L′(1) ≠ 0 → exactly 1 (PARI MW rank 1)
    root_number: int = -1

    # First Hecke eigenvalues aₚ(f₁) — verified by direct point counting
    # on y² + y = x³ − x² − x − 2 AND against the PARI q-expansion
    # (multiplicativity checked to a₅₀).  [Audit fix: the previous table
    # was wrong at p = 37, 41, 43, 47, 53, 59, 61.]
    hecke_eigenvalues: dict = None

    def __post_init__(self):
        if self.hecke_eigenvalues is None:
            object.__setattr__(self, "hecke_eigenvalues", {
                2: 0, 3: -1, 5: -1, 7: -2, 11: -1, 13: -1,
                17: -4, 19: 2, 23: 7, 29: -2, 31: -3, 37: -11,
                41: 10, 43: -4, 47: -4, 53: 2, 59: -1, 61: -2,
            })


CURVE_143A1 = EllipticCurve143a1()


# ═══════════════════════════════════════════════════════════════
#  Hecke Polynomials (Paper 26 §6, Paper 30 §2.3)
# ═══════════════════════════════════════════════════════════════

def hecke_polynomial_f2_T2() -> np.ndarray:
    """
    Characteristic polynomial of T₂ on f₂ (dimension 4):
        P₄(x) = x⁴ − 3x³ − x² + 5x + 1

    All four roots are real and satisfy |a₂| < 2√2.
    """
    return np.array([1, -3, -1, 5, 1], dtype=float)


def hecke_polynomial_f3_T2() -> np.ndarray:
    """
    Characteristic polynomial of T₂ on f₃ (dimension 6):
        P₆(x) = x⁶ − 10x⁴ + 2x³ + 24x² − 7x − 12

    ALL SIX ROOTS REAL, |a₂| ≤ 2.71 < 2√2 (Ramanujan).  Verified via
    PARI mfeigenbasis and by orbit traces Tr_f₃(a₂) = 0, Tr_f₃(a₄) = 8,
    Tr_f₃(a₉) = 13 (Newton-identity moments).

    [Audit fix B11: the polynomial previously hardcoded here,
    x⁶−x⁵−9x⁴+11x³+13x²−20x+8, is NOT a Hecke charpoly of this space.
    Its complex pair 0.5732±0.3564i was the "phantom eigenvalue" on
    which the superseded Tano formula was built.]
    """
    return np.array([1, 0, -10, 2, 24, -7, -12], dtype=float)


def hecke_polynomial_f2_T3() -> np.ndarray:
    """
    Characteristic polynomial of T₃ on f₂ (dimension 4):
        x⁴ − 7x² + 4x + 1
    All roots real, |a₃| ≤ 1.84 < 2√3.  (PARI-verified.)
    """
    return np.array([1, 0, -7, 4, 1], dtype=float)


def hecke_polynomial_f3_T3() -> np.ndarray:
    """
    Characteristic polynomial of T₃ on f₃ (dimension 6):
        x⁶ − 3x⁵ − 11x⁴ + 33x³ + 25x² − 91x + 28
    All roots real, |a₃| ≤ 3.31 < 2√3.  (PARI-verified.)
    """
    return np.array([1, -3, -11, 33, 25, -91, 28], dtype=float)


def hecke_roots_f2() -> np.ndarray:
    """All 4 roots of P₄(x) — real Hecke eigenvalues of T₂ on f₂."""
    return np.roots(hecke_polynomial_f2_T2())


def hecke_roots_f3() -> np.ndarray:
    """All 6 roots of P₆(x) — ALL REAL (trivial nebentypus)."""
    return np.roots(hecke_polynomial_f3_T2())


# ── SUPERSEDED: the "complex Hecke eigenvalue" ────────────────────────────
# The value below is a root of the WRONG polynomial (see hecke_polynomial_f3_T2
# docstring).  It does not exist in the arithmetic of X₀(143): every Hecke
# eigenvalue on S₂^new(Γ₀(143)) is totally real.  Retained ONLY so that
# historical code (Paper 26 §8, examples/03) still runs; the Tano Mass
# Formula built on it is superseded.  Do not use in new work.
A2_COMPLEX = 0.5732 + 0.3564j     # PHANTOM — root of a wrong polynomial
A2_COMPLEX_CONJ = A2_COMPLEX.conjugate()
_TANO_SUPERSEDED = True


def verify_complex_eigenvalue() -> dict:
    """
    Audit replacement for the old phantom check.

    Verifies the TRUTH: the correct T₂|f₃ charpoly has six real roots,
    all inside the Ramanujan bound, and the historical value
    0.5732 + 0.3564i is NOT among them.
    """
    roots = hecke_roots_f3()
    max_imag = float(max(abs(r.imag) for r in roots))
    max_abs = float(max(abs(r) for r in roots))
    phantom_dist = float(min(abs(r - A2_COMPLEX) for r in roots))
    return {
        "all_roots_real": max_imag < 1e-8,
        "max_imaginary_part": max_imag,
        "max_abs_root": max_abs,
        "ramanujan_bound": 2.0 * math.sqrt(2.0),
        "ramanujan_satisfied": max_abs < 2.0 * math.sqrt(2.0) + 1e-8,
        "phantom_value": A2_COMPLEX,
        "phantom_min_distance_to_true_root": phantom_dist,
        "phantom_is_root": phantom_dist < 1e-6,
        "status": "SUPERSEDED — the complex eigenvalue was an artifact of a wrong polynomial",
    }


# ═══════════════════════════════════════════════════════════════
#  Atkin-Lehner Spectral Separation (Paper 26 §6)
# ═══════════════════════════════════════════════════════════════

def atkin_lehner_separate(
    tr_Tm: float, tr_T11m: float, am_f1: float
) -> Tuple[float, float]:
    """
    Separate Hecke traces using Atkin-Lehner twists:

        Tr(aₘ(f₂)) = (Tr(Tₘ) + Tr(T₁₁ₘ)) / 2
        Tr(aₘ(f₃)) = (Tr(Tₘ) − Tr(T₁₁ₘ) − 2aₘ(f₁)) / 2

    Parameters
    ----------
    tr_Tm : float   — Tr(Tₘ) on full 11-dim new space
    tr_T11m : float — Tr(T₁₁ₘ) on full space
    am_f1 : float   — aₘ(f₁) from elliptic curve

    Returns (Tr_f2, Tr_f3)
    """
    tr_f2 = (tr_Tm + tr_T11m) / 2.0
    tr_f3 = (tr_Tm - tr_T11m - 2.0 * am_f1) / 2.0
    return tr_f2, tr_f3


# ═══════════════════════════════════════════════════════════════
#  Generation Counting Theorem (Paper 26, Theorem 4.1)
# ═══════════════════════════════════════════════════════════════

def generation_count() -> int:
    """
    Theorem: The number of fermion generations equals the number of
    nonempty irreducible representations of ℤ₂ × ℤ₂ on S₂^new(Γ₀(143)).

    Since ℤ₂ × ℤ₂ has 4 irreps and dim(−,−) = 0, exactly 3 exist.
    """
    return sum(1 for dim in AL_DECOMPOSITION.values() if dim > 0)


# ═══════════════════════════════════════════════════════════════
#  The Tano Mass Formula (Paper 26, Theorem 8.1)
# ═══════════════════════════════════════════════════════════════

def koide_angle_tano() -> float:
    """
    SUPERSEDED (audit B10/B11; Correction Session 1).

    The Tano Mass Formula for the Koide angle:

        θ₀ = 2/δ² + dim(f₂) · arg(a₂ᶜˣ(f₃))

    is built on the phantom value a₂ᶜˣ = 0.5732 + 0.3564i, which is NOT a
    Hecke eigenvalue of X₀(143) (all such eigenvalues are real; see
    hecke_polynomial_f3_T2).  The formula is therefore not evaluable as
    arithmetic of the modular curve.  This function still returns the
    historical number (θ₀ ≈ 2.31687 rad) for continuity checks only.
    """
    feigenbaum_term = 2.0 / FEIGENBAUM_DELTA ** 2     # ≈ 0.09174 rad (5.26°)
    arithmetic_angle = ORBIT_DIMENSIONS[1] * cmath.phase(A2_COMPLEX)  # 4 × 0.55628
    return feigenbaum_term + arithmetic_angle


def koide_angle_experimental() -> float:
    """Extract the Koide angle from experimental lepton masses."""
    me, mmu, mtau = LEPTONS.e * 1e3, LEPTONS.mu * 1e3, LEPTONS.tau * 1e3  # MeV
    sqrt_me, sqrt_mmu, sqrt_mtau = math.sqrt(me), math.sqrt(mmu), math.sqrt(mtau)
    M = (sqrt_me + sqrt_mmu + sqrt_mtau) / 3.0

    # Solve √mᵢ = M(1 + √2 cos(θ₀ + 2πi/3)) for θ₀
    # From the electron (i=0): cos(θ₀) = (√mₑ/M − 1)/√2
    cos_theta = (sqrt_me / M - 1.0) / math.sqrt(2.0)
    cos_theta = max(-1.0, min(1.0, cos_theta))
    return math.acos(cos_theta)


def tano_mass_predictions(m_e_MeV: float = 0.51099895) -> dict:
    """
    SUPERSEDED (audit B10/B11; Correction Session 1).

    Historical evaluation of the Tano Mass Formula.  The formula depends
    on the phantom complex eigenvalue and is not part of the verified
    arithmetic of X₀(143).  Retained for continuity; the dict gains a
    'superseded' flag.  Do not cite as a prediction.
    """
    theta0 = koide_angle_tano()

    # Koide parametrization: √mᵢ = M(1 + √2 cos(θ₀ + 2πi/3))
    sqrt_me = math.sqrt(m_e_MeV)

    # From i=0 (electron): √mₑ = M(1 + √2 cos θ₀) → M = √mₑ / (1 + √2 cos θ₀)
    M = sqrt_me / (1.0 + math.sqrt(2) * math.cos(theta0))

    masses_MeV = []
    for i in range(3):
        sqrt_mi = M * (1.0 + math.sqrt(2) * math.cos(theta0 + 2 * math.pi * i / 3.0))
        masses_MeV.append(sqrt_mi ** 2)

    exp_MeV = [m_e_MeV, LEPTONS.mu * 1e3, LEPTONS.tau * 1e3]
    errors = [
        abs(masses_MeV[i] - exp_MeV[i]) / exp_MeV[i] * 100
        for i in range(3)
    ]

    return {
        "theta0_rad": theta0,
        "theta0_deg": math.degrees(theta0),
        "M_MeV": M,
        "masses_MeV": {"electron": masses_MeV[0], "muon": masses_MeV[1], "tau": masses_MeV[2]},
        "experimental_MeV": {"electron": exp_MeV[0], "muon": exp_MeV[1], "tau": exp_MeV[2]},
        "errors_percent": {"electron": errors[0], "muon": errors[1], "tau": errors[2]},
        "superseded": True,
        "supersession_note": "Built on the phantom eigenvalue 0.5732+0.3564i "
                             "(root of a wrong polynomial). Not a verified prediction.",
    }


# ═══════════════════════════════════════════════════════════════
#  L-values (Paper 26 §7)
# ═══════════════════════════════════════════════════════════════

# Numerically computed central L-values
L_VALUES = {
    "L'(f1, 1)": 0.945696,     # derivative (forced vanishing of L, ε = −1); AFE + PARI ellL1 agree
    "sum_L(f2, 1)": 6.023,     # sum over 4 conjugates
    "sum_L(f3, 1)": 8.210,     # sum over 6 conjugates
}


# ═══════════════════════════════════════════════════════════════
#  Jacobian Stiffness Matrix (Paper 30)
# ═══════════════════════════════════════════════════════════════

@dataclass
class JacobianStiffness:
    """
    3×3 Hermitian matrix on the Jacobian J₀(143) = A₁ × A₂ × A₃.

    Eigenvalue ratios encode the Feigenbaum constants:
        λ₃/λ₂ ≈ 7.07 ≈ δₓ = 7.24  (2.3% match)
        λ₂/λ₁ ≈ 4.21 ≈ δ_F = 4.67 (9.8% match)

    Structural results (Paper 30):
        - Electron decouples from muon (coupling 0.018)
        - Muon-tau strongly mixed (coupling 0.72)
        - 32° rotation from Atkin-Lehner basis to mass basis
    """
    lambda1: float = 1.0       # electron eigenvalue (normalised)
    lambda2: float = 4.21      # muon eigenvalue
    lambda3: float = 29.79     # tau eigenvalue
    coupling_12: float = 0.018 # electron-muon coupling
    coupling_23: float = 0.72  # muon-tau coupling
    rotation_deg: float = 32.0 # AL basis → mass basis

    @property
    def ratio_32(self) -> float:
        """λ₃/λ₂ ≈ δₓ = 7.24"""
        return self.lambda3 / self.lambda2

    @property
    def ratio_21(self) -> float:
        """λ₂/λ₁ ≈ δ_F = 4.67"""
        return self.lambda2 / self.lambda1

    # Petersson norm ratios (Paper 30 §7)
    petersson_21: float = 1.72   # ⟨f₂,f₂⟩/⟨f₁,f₁⟩
    petersson_31: float = 2.61   # ⟨f₃,f₃⟩/⟨f₁,f₁⟩


JACOBIAN = JacobianStiffness()


# ═══════════════════════════════════════════════════════════════
#  Verified Hecke Data (audit + Correction Sessions, July 2026)
# ═══════════════════════════════════════════════════════════════
#
# Provenance: PARI mfeigenbasis polmods (Session 1), converted to exact
# integer traces via companion matrices of the coefficient-field polys.
# Cross-checks: f₁ from direct point-counts on 143.a1; per-orbit sums
# reproduce the trace form Tr(T_n) at every n ≤ 50 (51/51); q-expansion
# multiplicativity a_mn = a_m·a_n (336 relations) passes; T₃ charpoly
# moments Tr(a₄), Tr(a₉) match Newton identities.

# Per-orbit traces at primes: p → (a_p(f₁), Tr_f₂(a_p), Tr_f₃(a_p))
ORBIT_TRACES_VERIFIED = {
    2:  ( 0,  3,   0),   3:  (-1,  0,   3),   5:  (-1,  0,   1),
    7:  (-2,  6,   4),  11:  (-1,  4,  -6),  13: (-1, -4,   6),
    17: (-4,  6,   0),  19: ( 2,  8, -10),  23: ( 7, -4,  11),
}

# Full per-orbit traces Tr_i(a_n), n = 1..50 (index 0 ↔ n = 1)
ORBIT_TRACE_F1 = [1, 0, -1, -2, -1, 0, -2, 0, -2, 0, -1, 2, -1, 0, 1, 4, -4, 0,
                  2, 2, 2, 0, 7, 0, -4, 0, 5, 4, -2, 0, -3, 0, 1, 0, 2, 4, -11,
                  0, 1, 0, 10, 0, -4, 2, 2, 0, -4, -4, -3, 0]
ORBIT_TRACE_F2 = [4, 3, 0, 3, 0, -1, 6, 9, 2, -8, 4, 4, -4, -4, -10, 5, 6, -15,
                  8, -24, -2, 3, -4, 2, 12, -3, -12, -1, -10, 8, 2, -4, 0, 6,
                  -6, -28, 12, -5, 0, -30, 8, -13, 26, 3, 26, 6, -18, -21, 6, 29]
ORBIT_TRACE_F3 = [6, 0, 3, 8, 1, -3, 4, -6, 13, 6, -6, -6, 6, -12, 3, 8, 0, 6,
                  -10, 4, -12, 0, 11, -38, 23, 0, 9, 9, 2, -56, -9, -17, -3,
                  -40, -24, 11, 15, -9, 3, 16, -4, 19, -2, -8, -26, -6, 6, 19,
                  20, -4]
# Trace form Tr(T_n) on S₂^new, n = 1..50 — equals F1+F2+F3 at every n
TRACE_TOTALS_50 = [11, 3, 2, 9, 0, -4, 8, 3, 13, -2, -3, 0, 1, -16, -6, 17, 2,
                   -9, 0, -18, -12, 3, 14, -36, 31, -3, 2, 12, -10, -48, -10,
                   -21, -2, -34, -28, -13, 16, -14, 4, -14, 14, 6, 20, -3, 2,
                   0, -16, -6, 23, 25]

# Coefficient fields (PARI mfparams; Galois structure by Chebotarev census)
FIELD_POLY_F2 = [1, 0, -4, -1, 1]          # y⁴ − 4y² − y + 1; disc = 1957 = 19·103
FIELD_POLY_F3 = [1, 0, -10, -2, 24, 7, -12]  # y⁶−10y⁴−2y³+24y²+7y−12; disc = 194616205 = 5·7·5560463
FIELD_DISCRIMINANTS = {"f2": 1957, "f3": 194616205}
GALOIS_GROUPS = {"f2": "S4", "f3": "S6"}   # both totally real, unramified at 11, 13

# Root numbers re-derived without PARI: for p ‖ N the U_p eigenvalue is
# rational ±1 (visible: Tr(a_p) = ±dim), ε = −w₁₁·w₁₃:
#   f₁: (a₁₁,a₁₃) = (−1,−1) → ε = −1; f₂: (+1,−1) → +1; f₃: (−1,+1) → +1
ROOT_NUMBERS_LIST = (-1, +1, +1)

# Frob₁₁ on the f₃ coefficient field: P₃ ≡ (y²−2y−5)(y⁴+2y³−y²−5y−2) (mod 11),
# cycle type [2,4]; the quadratic factor has Hensel roots 1+3296i, 1+11345i
# (mod 11⁴, i²=−1).  NOTE (canonicity): Gal = S₆ acts transitively on the 15
# pair-partitions, so no arithmetic structure canonically selects a Q₂ pair;
# Paper 32's Q₂/Q₄ split is an analytic datum (period phases), not algebraic.
FROB_11_CYCLE_TYPE_F3 = (2, 4)

# Rankin–Selberg coupling matrix at s = 3 (Correction Session 4 definitions):
#   D_ij = Σ_{n≤N} Tr_i(a_n)Tr_j(a_n)/n³,  ε_ij = D_ij/(d_i d_j) − 1,
#   Q = ⟨ε, R⟩_F with R = ε⃗ε⃗ᵀ,  Q_corr = Q/(‖ε‖_F‖R‖_F).
# Independently reproduced from the raw traces above:
RS_COUPLING = {
    "s": 3,
    "eps_matrix_N50": [
        [+0.149686, -0.037330, -0.071683],
        [-0.037330, +0.118562, +0.009885],
        [-0.071683, +0.009885, +0.064248],
    ],
    "Q_N50": +0.570292,           # → +0.580690 (N=500 exact subset)
    "Q_corr_N50": +0.819186,      # → +0.816328 (N=500 exact subset)
    "Q_N1500": +0.587882,         # Session 4 value (NMAX = 1500)
    "Q_corr_N1500": +0.813978,
    "strict_sign_rule": True,     # sign(ε_ij) = ε_i·ε_j for all i ≤ j (6/6)
}


def rankin_selberg_epsilon(n_max: int = 50) -> list:
    """
    Recompute the ε-matrix at s = 3 from the verified trace tables
    (self-certifying: no external data needed).  Valid for n_max ≤ 50.
    """
    if n_max > 50:
        raise ValueError("verified tables extend to n = 50 only")
    dims = ORBIT_DIMENSIONS
    traces = (ORBIT_TRACE_F1, ORBIT_TRACE_F2, ORBIT_TRACE_F3)
    eps = [[0.0] * 3 for _ in range(3)]
    for i in range(3):
        for j in range(3):
            D = sum(traces[i][n - 1] * traces[j][n - 1] / n ** 3
                    for n in range(1, n_max + 1))
            eps[i][j] = D / (dims[i] * dims[j]) - 1.0
    return eps


def rankin_selberg_Q(n_max: int = 50) -> dict:
    """Q and Q_corr at s = 3 with R = outer(ROOT_NUMBERS_LIST)."""
    eps = rankin_selberg_epsilon(n_max)
    R = [[ROOT_NUMBERS_LIST[i] * ROOT_NUMBERS_LIST[j] for j in range(3)]
         for i in range(3)]
    Q = sum(eps[i][j] * R[i][j] for i in range(3) for j in range(3))
    eps_norm = math.sqrt(sum(eps[i][j] ** 2 for i in range(3) for j in range(3)))
    R_norm = 3.0
    strict = all(
        (eps[i][j] > 0) == (R[i][j] > 0)
        for i in range(3) for j in range(i, 3)
    )
    return {"Q": Q, "Q_corr": Q / (eps_norm * R_norm), "strict": strict,
            "eps_matrix": eps, "n_max": n_max, "s": 3}
