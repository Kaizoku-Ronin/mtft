"""
MTFT Materials Science Layer: The Tano Metric
===============================================

The Tano metric T_a = S_m / M_a (entropy per unit mass) quantifies
elemental 'etherealness'. Combined with the Geometry Index G, it
predicts superconducting T_c via:

    λ_eff = λ₀ + a·ΔT + b·G
    T_c ∝ Θ · exp(−1/λ_eff)

where ΔT = T_light − T_heavy captures the heavy-cage + light-injector
principle for superconductivity.

The Tano-Seebeck Bridge (Paper 19) connects entropy-per-mass to
entropy-per-charge: Ψ = S_Seebeck · M / (T_a · e).

Reference: Paper 1 §37, Paper 19.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Dict, List, Optional

from mtft.constants import CriticalDepths, FEIGENBAUM_DELTA


# ═══════════════════════════════════════════════════════════════
#  Element Database
# ═══════════════════════════════════════════════════════════════

@dataclass(frozen=True)
class Element:
    """Elemental data for Tano metric computation."""
    symbol: str
    name: str
    Z: int
    molar_mass: float      # g/mol
    molar_entropy: float   # J/(mol·K) at 298 K
    debye_temp: float      # K

    @property
    def tano_metric(self) -> float:
        """T_a = S_m / M_a — entropy per unit mass."""
        return self.molar_entropy / self.molar_mass


# Standard thermodynamic data (CRC Handbook, NIST)
_ELEMENT_DATA = {
    "H":  Element("H",  "Hydrogen",   1,   1.008, 130.68, 110),
    "He": Element("He", "Helium",     2,   4.003,  126.15, 25),
    "Li": Element("Li", "Lithium",    3,   6.941,   29.12, 344),
    "Be": Element("Be", "Beryllium",  4,   9.012,    9.50, 1440),
    "B":  Element("B",  "Boron",      5,  10.81,     5.90, 1250),
    "C":  Element("C",  "Carbon",     6,  12.011,    5.74, 2230),
    "N":  Element("N",  "Nitrogen",   7,  14.007, 191.61, 70),
    "O":  Element("O",  "Oxygen",     8,  15.999, 205.14, 91),
    "Na": Element("Na", "Sodium",    11,  22.990,   51.30, 158),
    "Mg": Element("Mg", "Magnesium", 12,  24.305,   32.67, 400),
    "Al": Element("Al", "Aluminium", 13,  26.982,   28.33, 428),
    "Si": Element("Si", "Silicon",   14,  28.086,   18.83, 645),
    "S":  Element("S",  "Sulfur",    16,  32.06,    32.05, 250),
    "K":  Element("K",  "Potassium", 19,  39.098,   64.68, 91),
    "Ca": Element("Ca", "Calcium",   20,  40.078,   41.59, 230),
    "Ti": Element("Ti", "Titanium",  22,  47.867,   30.72, 420),
    "V":  Element("V",  "Vanadium",  23,  50.942,   28.91, 380),
    "Cr": Element("Cr", "Chromium",  24,  51.996,   23.77, 630),
    "Fe": Element("Fe", "Iron",      26,  55.845,   27.28, 470),
    "Ni": Element("Ni", "Nickel",    28,  58.693,   29.87, 450),
    "Cu": Element("Cu", "Copper",    29,  63.546,   33.15, 343),
    "Zn": Element("Zn", "Zinc",      30,  65.38,    41.63, 327),
    "Ga": Element("Ga", "Gallium",   31,  69.723,   40.88, 320),
    "Nb": Element("Nb", "Niobium",   41,  92.906,   36.40, 275),
    "Mo": Element("Mo", "Molybdenum",42,  95.95,    28.66, 450),
    "Ag": Element("Ag", "Silver",    47, 107.868,   42.55, 225),
    "Sn": Element("Sn", "Tin",       50, 118.71,    51.18, 200),
    "Ba": Element("Ba", "Barium",    56, 137.327,   62.42, 110),
    "La": Element("La", "Lanthanum", 57, 138.905,   56.90, 142),
    "Y":  Element("Y",  "Yttrium",   39,  88.906,   44.43, 280),
    "Pb": Element("Pb", "Lead",      82, 207.2,     64.81, 105),
    "Bi": Element("Bi", "Bismuth",   83, 208.98,    56.74, 119),
    "Au": Element("Au", "Gold",      79, 196.967,   47.49, 165),
    "Pt": Element("Pt", "Platinum",  78, 195.084,   41.63, 240),
    "W":  Element("W",  "Tungsten",  74, 183.84,    32.64, 400),
    "Hg": Element("Hg", "Mercury",   80, 200.59,    75.90, 72),
}

ELEMENTS: Dict[str, Element] = _ELEMENT_DATA


def get_element(symbol: str) -> Element:
    """Look up an element by symbol."""
    if symbol not in ELEMENTS:
        raise KeyError(f"Element '{symbol}' not in database. "
                       f"Available: {sorted(ELEMENTS.keys())}")
    return ELEMENTS[symbol]


# ═══════════════════════════════════════════════════════════════
#  Tano Metric Core
# ═══════════════════════════════════════════════════════════════

def tano_contrast(light: str, heavy: str) -> float:
    """
    ΔT = T_light − T_heavy.

    Large ΔT → strong effective coupling → high T_c.
    The 'heavy cage + light injector' principle.
    """
    return get_element(light).tano_metric - get_element(heavy).tano_metric


GEOMETRY_INDEX_MAP = {
    "simple_metal": 0.2,
    "covalent": 0.4,
    "layered": 0.6,
    "moire": 0.8,
    "interface": 0.8,
    "clathrate": 1.0,
    "cage": 1.0,
}


def geometry_index(structure: str) -> float:
    """
    G ∈ [0.2, 1.0] from crystal structure type.

    Maps modular curvature |R(τ)| to a discrete classification:
      simple_metal → 0.2,  covalent → 0.4,  layered → 0.6,
      moire/interface → 0.8,  clathrate/cage → 1.0
    """
    key = structure.lower().replace(" ", "_").replace("-", "_")
    if key not in GEOMETRY_INDEX_MAP:
        raise ValueError(
            f"Unknown structure '{structure}'. "
            f"Valid: {list(GEOMETRY_INDEX_MAP.keys())}"
        )
    return GEOMETRY_INDEX_MAP[key]


def effective_coupling(
    delta_T: float,
    G: float,
    lambda_0: float = 0.3,
    a: float = 0.01,
    b: float = 0.5,
) -> float:
    """λ_eff = λ₀ + a·ΔT + b·G  [Eq. 115]"""
    return lambda_0 + a * delta_T + b * G


def predict_Tc(
    light: str,
    heavy: str,
    structure: str,
    Theta: float = 300.0,
    lambda_0: float = 0.3,
    a: float = 0.01,
    b: float = 0.5,
) -> dict:
    """
    Predict superconducting T_c from Tano metric + geometry.

    T_c = Θ · exp(−1/λ_eff)  [Eq. 116]

    Parameters
    ----------
    light : str
        Light element symbol (the 'injector').
    heavy : str
        Heavy element symbol (the 'cage').
    structure : str
        Crystal structure type.
    Theta : float
        Characteristic phonon temperature (K).
    """
    dT = tano_contrast(light, heavy)
    G = geometry_index(structure)
    lam = effective_coupling(dT, G, lambda_0, a, b)

    if lam <= 0:
        Tc = 0.0
    else:
        Tc = Theta * math.exp(-1.0 / lam)

    return {
        "T_c_K": Tc,
        "lambda_eff": lam,
        "delta_T": dT,
        "G": G,
        "Theta_K": Theta,
        "light": light,
        "heavy": heavy,
        "structure": structure,
    }


# ═══════════════════════════════════════════════════════════════
#  Materials Screening (Paper 1, Table 9)
# ═══════════════════════════════════════════════════════════════

_DEFAULT_CANDIDATES = [
    # (light, heavy, structure, Theta, observed_Tc)
    ("H",  "La", "cage",         1500, 260.0),   # LaH₁₀
    ("H",  "S",  "cage",         1200, 203.0),   # H₃S
    ("O",  "Ba", "layered",       500,  92.0),   # YBCO (Y-Ba-Cu-O)
    ("B",  "Mg", "covalent",      900,  39.0),   # MgB₂
    ("Al", "Al", "simple_metal",  428,   1.2),   # Al
    ("Pb", "Pb", "simple_metal",  105,   7.2),   # Pb
]


def materials_screening(
    candidates: Optional[List[dict]] = None,
) -> List[dict]:
    """
    Screen materials for superconductivity using the Tano metric.

    Default reproduces Paper 1, Table 9.

    Each candidate dict: {light, heavy, structure, Theta, observed_Tc}
    """
    if candidates is None:
        rows = []
        for light, heavy, struct, theta, obs_Tc in _DEFAULT_CANDIDATES:
            pred = predict_Tc(light, heavy, struct, Theta=theta)
            pred["observed_Tc_K"] = obs_Tc
            pred["error_percent"] = (
                abs(pred["T_c_K"] - obs_Tc) / obs_Tc * 100
                if obs_Tc > 0 else 0.0
            )
            rows.append(pred)
        return rows
    else:
        rows = []
        for c in candidates:
            pred = predict_Tc(
                c["light"], c["heavy"], c["structure"],
                Theta=c.get("Theta", 300),
            )
            if "observed_Tc" in c:
                pred["observed_Tc_K"] = c["observed_Tc"]
                pred["error_percent"] = (
                    abs(pred["T_c_K"] - c["observed_Tc"]) / c["observed_Tc"] * 100
                    if c["observed_Tc"] > 0 else 0.0
                )
            rows.append(pred)
        return rows


# ═══════════════════════════════════════════════════════════════
#  Tano-Seebeck Bridge (Paper 19)
# ═══════════════════════════════════════════════════════════════

E_CHARGE_SI = 1.602176634e-19  # C


def seebeck_diagnostic(
    symbol: str,
    S_seebeck_uVK: float,
) -> dict:
    """
    Ψ = S_Seebeck · M / (T_a · e) — dimensionless diagnostic.

    Measures how efficiently a material's entropic 'etherealness'
    converts to thermoelectric voltage.

    Parameters
    ----------
    symbol : str
        Element symbol.
    S_seebeck_uVK : float
        Seebeck coefficient in μV/K.
    """
    el = get_element(symbol)
    Ta = el.tano_metric
    # Convert μV/K → V/K
    S_SI = S_seebeck_uVK * 1e-6
    # Ψ = S · M / (T_a · e)  [dimensionless in SI]
    # T_a in J/(mol·K·(g/mol)) = J/(g·K), M in g/mol
    # Need consistent units: use S in V/K, M in kg/mol, T_a in J/(kg·K)
    M_kg = el.molar_mass * 1e-3  # g/mol → kg/mol
    Ta_SI = el.molar_entropy / M_kg  # J/(kg·K)
    psi = S_SI * M_kg / (Ta_SI * E_CHARGE_SI) if Ta_SI > 0 else 0.0

    return {
        "element": symbol,
        "psi": psi,
        "S_seebeck_uVK": S_seebeck_uVK,
        "tano_metric": Ta,
        "molar_mass": el.molar_mass,
        "interpretation": (
            "Large |Ψ| → efficient entropy-to-voltage conversion. "
            "Ψ saturates when scattering asymmetry is maximised, "
            "corresponding to y ≈ y_c in MTFT."
        ),
    }


# ═══════════════════════════════════════════════════════════════
#  Josephson Holonomy (Paper 1 §40.2)
# ═══════════════════════════════════════════════════════════════

def josephson_holonomy() -> dict:
    """
    MTFT prediction: universal holonomy phase shift in Josephson junctions.

    Φ_H/Φ₀ ≈ −2.3% (measured in Sn-InSb, Paper 1 §36.2)

    Should be:
      1. Independent of junction material
      2. Independent of junction geometry
      3. Present in all superconductors

    BCS predicts Φ_H = 0.
    """
    y_c = CriticalDepths.y_conf
    phi_ratio = -y_c / (2 * math.pi) * (1 + 1 / FEIGENBAUM_DELTA)
    return {
        "phi_H_over_phi_0_predicted": phi_ratio,
        "phi_H_over_phi_0_measured": -0.023,
        "percent_predicted": phi_ratio * 100,
        "percent_measured": -2.3,
        "BCS_prediction": 0.0,
        "material_independent": True,
        "falsification": "If different materials show different Φ_H, MTFT is falsified",
    }
