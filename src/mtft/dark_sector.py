"""
τ-field dark sector: vortex halos, flat rotation curves, Tully-Fisher.

A static spherical τ-vortex τ(r) = A ln(r/r₀) produces:

    ρ_τ(r) = A² / (2r²)
    M_τ(<r) = 2π A² r
    v²(r)   = G M_eff / r = 2π G A² = v²_∞ = const

giving flat rotation curves with no dark matter particles.

Tully-Fisher relation emerges from BH entropy scaling:
    A₀ ∝ M_BH^{1/4}  →  v⁴_∞ ∝ M_BH ∝ M_baryonic
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from mtft.constants import PhysicalConstants as PC


# ── τ-vortex halo ───────────────────────────────────────────────

@dataclass
class TauVortexHalo:
    """
    Static spherical τ-vortex dark halo.

    Parameters
    ----------
    A : float
        Vortex strength (dimensionless in geometric units, or GeV in natural).
    r0 : float
        Core radius (GeV⁻¹).  ρ_τ diverges as r→0; physically capped near r₀.
    epsilon : float
        Anisotropy parameter from BH spin: ρ_τ ∝ [1 + 2ε P₂(cos θ)].
        0 = spherical, >0 = prolate, <0 = oblate.
    """

    A: float
    r0: float = 1.0
    epsilon: float = 0.0

    # ── density ──────────────────────────────────────────────
    def rho(self, r: np.ndarray | float) -> np.ndarray:
        """
        τ-field energy density:

            ρ_τ(r) = A² / (2 r²)

        with soft core at r₀.
        """
        r = np.asarray(r, dtype=float)
        r_eff = np.maximum(r, self.r0)
        return self.A ** 2 / (2.0 * r_eff ** 2)

    def rho_anisotropic(self, r: np.ndarray, cos_theta: np.ndarray) -> np.ndarray:
        """
        Anisotropic density from BH spin:

            ρ_τ(r, θ) = ρ₀(r) [1 + 2ε P₂(cos θ)]

        where P₂(x) = (3x² − 1)/2.
        """
        P2 = (3.0 * cos_theta ** 2 - 1.0) / 2.0
        return self.rho(r) * (1.0 + 2.0 * self.epsilon * P2)

    # ── enclosed mass ────────────────────────────────────────
    def enclosed_mass(self, r: np.ndarray | float) -> np.ndarray:
        """
        Enclosed τ-mass:

            M_τ(<r) = 4π ∫₀ʳ r'² ρ_τ(r') dr' = 2π A² r   (for r ≫ r₀)
        """
        r = np.asarray(r, dtype=float)
        # Exact integral with core softening:
        r_eff = np.maximum(r, self.r0)
        return 2.0 * math.pi * self.A ** 2 * (r_eff - self.r0) + (
            2.0 * math.pi * self.A ** 2 * self.r0 * np.minimum(r / self.r0, 1.0)
        )

    # ── rotation curve ───────────────────────────────────────
    def v_circular(self, r: np.ndarray | float, M_baryon: float = 0.0) -> np.ndarray:
        """
        Circular velocity:

            v²(r) = G [M_baryon + M_τ(<r)] / r

        Returns v in natural units (dimensionless, v/c).
        """
        r = np.asarray(r, dtype=float)
        M_total = M_baryon + self.enclosed_mass(r)
        r_safe = np.maximum(r, 1e-30)
        return np.sqrt(PC.G_N * M_total / r_safe)

    @property
    def v_infinity(self) -> float:
        """
        Asymptotic flat velocity:

            v²_∞ = 2π G A²
        """
        return math.sqrt(2.0 * math.pi * PC.G_N * self.A ** 2)


# ── Convenience functions ────────────────────────────────────────

def rotation_curve(
    r: np.ndarray,
    A: float,
    M_baryon: float = 0.0,
    r0: float = 1.0,
) -> np.ndarray:
    """
    Compute rotation curve v(r) for a τ-vortex halo.

    Parameters
    ----------
    r : array
        Radial distances.
    A : float
        Vortex amplitude.
    M_baryon : float
        Central baryonic mass.
    r0 : float
        Core radius.

    Returns
    -------
    v : array
        Circular velocity at each r.
    """
    halo = TauVortexHalo(A=A, r0=r0)
    return halo.v_circular(r, M_baryon=M_baryon)


def tully_fisher(v_flat: float) -> float:
    """
    Tully-Fisher relation from BH entropy scaling.

    In MTFT: A ∝ M_BH^{1/4}, so v⁴_∞ = (2πG)² A⁴ ∝ M_BH ∝ M_baryonic.

    Parameters
    ----------
    v_flat : float
        Observed flat rotation velocity (natural units, v/c).

    Returns
    -------
    M_baryonic_estimate : float
        Estimated baryonic mass (GeV).
    """
    # Calibration: v⁴ / (2πG)² = A⁴ → M_BH = C·A⁴
    # Using empirical TF calibration constant
    C_TF = 50.0  # Calibration to be refined against observations
    return C_TF * v_flat ** 4 / (2.0 * math.pi * PC.G_N) ** 2


def vortex_from_bh_mass(M_bh: float) -> float:
    """
    Compute vortex amplitude A from black hole mass via entropy scaling:

        A = (S_BH / (4π))^{1/4} = (M_BH² / (4 M_Pl²))^{1/4}

    (using Bekenstein-Hawking entropy S = A_horizon / 4 = 4πG²M² / 4 = πG²M²)
    """
    S_bh = math.pi * PC.G_N ** 2 * M_bh ** 2  # in natural units
    return (S_bh / (4.0 * math.pi)) ** 0.25
