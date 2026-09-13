"""mtft.surface.condensation — tachyon condensation as bundle extension (COND-01, v0.29.0).

U(2) flux model on X0(143): the off-diagonal scalar is the extension class e in Ext^1(K^3, O) = H^1(K^-3);
condensation forms 0 -> O -> V -> K^3 -> 0.  Exact statements (constant curvature -1/R^2, area 48 pi R^2):
  tachyon mass^2 = -2 pi |d| / A = -|d| / (24 R^2)      (Nielsen–Olesen; d = deg K^-3 = -72 -> -3/R^2)
  Landau degeneracy |d| + g - 1 = 84 = h^1(K^-3)
  Yang–Mills energies (E = (1/2g²)∫tr B²): split (0,72) = 216 pi/(g²R²), balanced 108 pi/(g²R²), drop 108 pi/(g²R²)
  long exact sequence: 0 -> H^0(O) -> H^0(V) -> H^0(K^3) --delta--> H^1(O) -> H^1(V) -> 0, delta = e cup . = the Yukawa map;
  generic e: rank 13, h^0(V) = 48, h^1(V) = 0, index 48; AL-invariant e: sector ranks (1,6,5,1), remainder (11,12,12,12)+1.
Stability of the non-split extension is generic (no line subbundle of degree >= 36) and recorded as a condition.
"""
from __future__ import annotations

import numpy as np

from .yukawa import yukawa_matrix, _SECTORS


def tachyon_mass2(deg: int, g: int = 13, R: float = 1.0) -> float:
    return -abs(deg) / (2 * (g - 1) * R * R)


def condensation_energy(deg: int, g: int = 13, R: float = 1.0, gauge_coupling: float = 1.0, degrees=None) -> dict:
    """Yang–Mills energies in the convention E = (1/2 g_YM^2) ∫ tr(B^2) dA on the compact curvature −1/R² model
    (A = 4 pi (g−1) R²), for the split bundle L1 ⊕ L2 of degrees (d1, d2) (default (0, deg)) and the balanced
    (polystable) bound at the same total degree.  KK05 correction: the earlier scalar return value was the SPLIT
    energy 216 pi/(g² R²) for (0,72); the drop to the balanced bound is 108 pi/(g² R²).  Attaining the bound on a
    given holomorphic bundle requires polystability (Atiyah–Bott / Harder–Narasimhan)."""
    A = 4 * np.pi * (g - 1) * R * R
    d1, d2 = degrees if degrees is not None else (0, -deg if deg < 0 else deg)
    split = (2 * np.pi) ** 2 * (d1 * d1 + d2 * d2) / (2 * A * gauge_coupling ** 2)
    balanced = (2 * np.pi) ** 2 * (d1 + d2) ** 2 / (4 * A * gauge_coupling ** 2)
    return {"split": split, "balanced_bound": balanced, "drop": split - balanced, "convention": "E = (1/2 g^2) int tr B^2, A = 48 pi R^2"}


def extension_cohomology(Y: dict, e, g: int = 13) -> dict:
    T = Y["tensor"]; M = yukawa_matrix(T, e)
    d, nb = M.shape; r = int(np.linalg.matrix_rank(M, tol=1e-9 * max(np.abs(M).max(), 1e-300)))
    sec = {}
    for s in _SECTORS:
        rows = [a for a in range(d) if tuple(Y["omega_sectors"][a]) == s]; cols = [b for b in range(nb) if tuple(Y["cubic_sectors"][b]) == s]
        sub = M[np.ix_(rows, cols)]
        rk = int(np.linalg.matrix_rank(sub, tol=1e-9 * max(np.abs(M).max(), 1e-300))) if sub.size else 0
        sec[s] = {"H0_K3": len(cols), "H1_O": len(rows), "rank": rk, "massless_remainder": len(cols) - rk}
    return {"rank_delta": r, "h0_V": 1 + nb - r, "h1_V": d - r, "index": 1 + nb - d, "sectors": sec}
