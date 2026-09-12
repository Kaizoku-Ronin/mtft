"""mtft.surface.condensation — tachyon condensation as bundle extension (COND-01, v0.29.0).

U(2) flux model on X0(143): the off-diagonal scalar is the extension class e in Ext^1(K^3, O) = H^1(K^-3);
condensation forms 0 -> O -> V -> K^3 -> 0.  Exact statements (constant curvature -1/R^2, area 48 pi R^2):
  tachyon mass^2 = -2 pi |d| / A = -|d| / (24 R^2)      (Nielsen–Olesen; d = deg K^-3 = -72 -> -3/R^2)
  Landau degeneracy |d| + g - 1 = 84 = h^1(K^-3)
  condensation energy (Atiyah–Bott) = (2 pi)^2 d^2 / (2 A) = 216 pi / (g^2 R^2)
  long exact sequence: 0 -> H^0(O) -> H^0(V) -> H^0(K^3) --delta--> H^1(O) -> H^1(V) -> 0, delta = e cup . = the Yukawa map;
  generic e: rank 13, h^0(V) = 48, h^1(V) = 0, index 48; AL-invariant e: sector ranks (1,6,5,1), remainder (11,12,12,12)+1.
Stability of the non-split extension is generic (no line subbundle of degree >= 36) and recorded as a condition.
"""
from __future__ import annotations

import numpy as np

from .yukawa import yukawa_matrix, _SECTORS


def tachyon_mass2(deg: int, g: int = 13, R: float = 1.0) -> float:
    return -abs(deg) / (2 * (g - 1) * R * R)


def condensation_energy(deg: int, g: int = 13, R: float = 1.0, gauge_coupling: float = 1.0) -> float:
    A = 4 * np.pi * (g - 1) * R * R
    return (2 * np.pi) ** 2 * deg ** 2 / (2 * A) / gauge_coupling ** 2


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
