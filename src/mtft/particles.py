"""
Standard Model particle database with MTFT modular-time embeddings.

Each particle has a holonomy coupling κ_f such that its mass is:

    m_f = κ_f |sin θ₀| / R_τ

The κ hierarchy across generations encodes how each fermion's
wavefunction overlaps with the compact modular-time direction.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from enum import Enum
from typing import List, Optional

import numpy as np


class ParticleType(Enum):
    GAUGE = "gauge"
    FERMION = "fermion"
    SCALAR = "scalar"


@dataclass(frozen=True)
class Particle:
    """
    A Standard Model particle with its MTFT modular-time embedding.

    Attributes
    ----------
    name : str
        Human-readable name.
    symbol : str
        LaTeX-style symbol (e.g. 'e⁻', 'W±').
    ptype : ParticleType
        gauge, fermion, or scalar.
    mass_GeV : float
        PDG mass in GeV.
    charge_em : float
        Electromagnetic charge in units of e.
    color : int
        Number of colour degrees of freedom (1 for leptons, 3 for quarks).
    representation : str
        SU(3)×SU(2)×U(1)_Y quantum numbers, e.g. '(3,2,+1/6)'.
    generation : int or None
        1, 2, or 3 for fermions; None for bosons.
    kappa : float
        Holonomy coupling κ_f.
    spin : float
        Spin quantum number.
    """

    name: str
    symbol: str
    ptype: ParticleType
    mass_GeV: float
    charge_em: float
    color: int
    representation: str
    generation: Optional[int]
    kappa: float
    spin: float

    def predicted_mass(self, theta0: float, R_tau: float) -> float:
        """Mass from MTFT: m = κ |sin θ₀| / R_τ."""
        return self.kappa * abs(math.sin(theta0)) / R_tau

    def compute_kappa(self, theta0: float, R_tau: float) -> float:
        """Infer κ from observed mass: κ = m R_τ / |sin θ₀|."""
        sin_theta = abs(math.sin(theta0))
        if sin_theta < 1e-30:
            return float("inf")
        return self.mass_GeV * R_tau / sin_theta


# ── Standard Model catalog ───────────────────────────────────────

def _build_catalog() -> List[Particle]:
    """Construct the full SM particle list with PDG masses and κ values."""
    P = Particle
    G = ParticleType.GAUGE
    F = ParticleType.FERMION
    S = ParticleType.SCALAR

    return [
        # ── Gauge bosons ─────────────────────────────────────
        P("Photon",  "γ",   G,   0.0,       0.0,  1, "(1,1,0)",     None, 0.0,       1.0),
        P("W boson", "W±",  G,  80.3692,     1.0,  1, "(1,3,0)",     None, 0.05,      1.0),
        P("Z boson", "Z⁰",  G,  91.1876,    0.0,  1, "(1,3,0)",     None, 0.05,      1.0),
        P("Gluon",   "g",   G,   0.0,       0.0,  8, "(8,1,0)",     None, 0.0,       1.0),

        # ── Higgs ────────────────────────────────────────────
        P("Higgs",   "H⁰",  S, 125.20,     0.0,  1, "(1,2,+1/2)",  None, 1.0,       0.0),

        # ── Generation 1 quarks ──────────────────────────────
        P("Up",      "u",   F,   0.00216,   2/3,  3, "(3,2,+1/6)",  1,    9.5e-6,    0.5),
        P("Down",    "d",   F,   0.00467,  -1/3,  3, "(3,2,+1/6)",  1,    1.97e-5,   0.5),

        # ── Generation 2 quarks ──────────────────────────────
        P("Charm",   "c",   F,   1.27,      2/3,  3, "(3,2,+1/6)",  2,    5.16e-3,   0.5),
        P("Strange", "s",   F,   0.0934,   -1/3,  3, "(3,2,+1/6)",  2,    3.79e-4,   0.5),

        # ── Generation 3 quarks ──────────────────────────────
        P("Top",     "t",   F, 172.69,      2/3,  3, "(3,2,+1/6)",  3,    0.702,     0.5),
        P("Bottom",  "b",   F,   4.18,     -1/3,  3, "(3,2,+1/6)",  3,    0.0183,    0.5),

        # ── Generation 1 leptons ─────────────────────────────
        P("Electron",    "e⁻",  F, 0.000511,  -1.0, 1, "(1,1,-1)",    1, 2.068e-6,  0.5),
        P("Electron ν",  "νe",  F, 1e-9,       0.0, 1, "(1,2,-1/2)",  1, 1e-11,     0.5),

        # ── Generation 2 leptons ─────────────────────────────
        P("Muon",        "μ⁻",  F, 0.10566,   -1.0, 1, "(1,1,-1)",    2, 4.29e-4,   0.5),
        P("Muon ν",      "νμ",  F, 1e-9,       0.0, 1, "(1,2,-1/2)",  2, 1e-11,     0.5),

        # ── Generation 3 leptons ─────────────────────────────
        P("Tau",         "τ⁻",  F, 1.7768,    -1.0, 1, "(1,1,-1)",    3, 7.21e-3,   0.5),
        P("Tau ν",       "ντ",  F, 1e-9,       0.0, 1, "(1,2,-1/2)",  3, 1e-11,     0.5),
    ]


class StandardModel:
    """
    Container for the full SM particle catalog with query helpers.

    Usage
    -----
    >>> sm = StandardModel()
    >>> sm.by_name("Top").mass_GeV
    172.69
    >>> sm.quarks(generation=3)
    [Particle(Top), Particle(Bottom)]
    """

    def __init__(self):
        self._particles = _build_catalog()
        self._by_name = {p.name.lower(): p for p in self._particles}

    @property
    def all(self) -> List[Particle]:
        return list(self._particles)

    def by_name(self, name: str) -> Particle:
        key = name.lower()
        if key not in self._by_name:
            raise KeyError(f"Unknown particle: {name}")
        return self._by_name[key]

    def fermions(self, generation: Optional[int] = None) -> List[Particle]:
        out = [p for p in self._particles if p.ptype == ParticleType.FERMION]
        if generation is not None:
            out = [p for p in out if p.generation == generation]
        return out

    def quarks(self, generation: Optional[int] = None) -> List[Particle]:
        out = [p for p in self.fermions(generation) if p.color == 3]
        return out

    def leptons(self, generation: Optional[int] = None) -> List[Particle]:
        out = [p for p in self.fermions(generation) if p.color == 1]
        return out

    def bosons(self) -> List[Particle]:
        return [p for p in self._particles if p.ptype in (ParticleType.GAUGE, ParticleType.SCALAR)]

    def kappa_hierarchy(self) -> List[tuple[str, float]]:
        """Return fermions sorted by κ (ascending) — the generation hierarchy."""
        ferm = self.fermions()
        ferm_sorted = sorted(ferm, key=lambda p: p.kappa)
        return [(p.name, p.kappa) for p in ferm_sorted]

    def mass_table(self, theta0: float, R_tau: float) -> list[dict]:
        """
        Compare PDG masses with MTFT predictions for all fermions.

        Returns list of dicts with 'name', 'mass_pdg', 'mass_mtft', 'ratio'.
        """
        rows = []
        for p in self.fermions():
            pred = p.predicted_mass(theta0, R_tau)
            ratio = pred / p.mass_GeV if p.mass_GeV > 0 else float("nan")
            rows.append({
                "name": p.name,
                "generation": p.generation,
                "mass_pdg_GeV": p.mass_GeV,
                "mass_mtft_GeV": pred,
                "ratio": ratio,
                "kappa": p.kappa,
            })
        return rows
