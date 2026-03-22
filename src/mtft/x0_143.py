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

The Tano Mass Formula (Paper 26, Theorem 8.1):

    θ₀ = 2/δ² + dim(f₂) · arg(a₂ᶜˣ(f₃))

predicts muon and tau masses from the Hecke eigenvalue structure
with zero free parameters beyond the electron mass anchor.

All data LMFDB-verified.  Reference: Papers 25, 26, 29, 30.
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
    analytic_rank: int = 0
    root_number: int = -1

    # First Hecke eigenvalues aₚ(f₁) — verified by point counting
    hecke_eigenvalues: dict = None

    def __post_init__(self):
        if self.hecke_eigenvalues is None:
            object.__setattr__(self, "hecke_eigenvalues", {
                2: 0, 3: -1, 5: -1, 7: -2, 11: -1, 13: -1,
                17: -4, 19: 2, 23: 7, 29: -2, 31: -3, 37: 3,
                41: -10, 43: 4, 47: 0, 53: -2, 59: 12, 61: 2,
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
        P₆(x) = x⁶ − x⁵ − 9x⁴ + 11x³ + 13x² − 20x + 8

    Has 4 real roots and ONE complex conjugate pair:
        a₂ᶜˣ = 0.5732 ± 0.3564i

    This is the ONLY non-real datum in the entire newform
    decomposition — its phase determines the lepton masses.
    """
    return np.array([1, -1, -9, 11, 13, -20, 8], dtype=float)


def hecke_roots_f2() -> np.ndarray:
    """All 4 roots of P₄(x) — real Hecke eigenvalues of T₂ on f₂."""
    return np.roots(hecke_polynomial_f2_T2())


def hecke_roots_f3() -> np.ndarray:
    """All 6 roots of P₆(x) — includes the complex conjugate pair."""
    return np.roots(hecke_polynomial_f3_T2())


# The critical complex eigenvalue
A2_COMPLEX = 0.5732 + 0.3564j     # a₂ᶜˣ(f₃)
A2_COMPLEX_CONJ = A2_COMPLEX.conjugate()


def verify_complex_eigenvalue() -> dict:
    """Verify a₂ᶜˣ against the roots of P₆."""
    roots = hecke_roots_f3()
    complex_roots = [r for r in roots if abs(r.imag) > 0.01]
    if len(complex_roots) >= 2:
        # Find the one closest to our stored value
        best = min(complex_roots, key=lambda r: abs(r - A2_COMPLEX))
        return {
            "stored": A2_COMPLEX,
            "computed": best,
            "error": abs(best - A2_COMPLEX),
            "phase_rad": cmath.phase(best),
            "phase_deg": math.degrees(cmath.phase(best)),
        }
    return {"error": "No complex roots found"}


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
    The Tano Mass Formula for the Koide angle:

        θ₀ = 2/δ² + dim(f₂) · arg(a₂ᶜˣ(f₃))

    where:
        δ = 4.669201609...   (Feigenbaum)
        dim(f₂) = 4          (muon-sector orbit dimension)
        a₂ᶜˣ(f₃) = 0.5732 + 0.3564i  (complex Hecke eigenvalue)

    Returns θ₀ in radians.
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
    Predict muon and tau masses from the Tano Mass Formula.

    Uses only m_e as dimensional anchor + arithmetic of X₀(143).

    Returns dict with 'theta0_rad', 'theta0_deg', 'masses_MeV',
    'experimental_MeV', 'errors_percent'.
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
    }


# ═══════════════════════════════════════════════════════════════
#  L-values (Paper 26 §7)
# ═══════════════════════════════════════════════════════════════

# Numerically computed central L-values
L_VALUES = {
    "L'(f1, 1)": 0.791,        # derivative (forced vanishing, ε = −1)
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
