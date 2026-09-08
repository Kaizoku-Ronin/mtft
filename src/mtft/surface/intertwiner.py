"""mtft.surface.intertwiner — the exact map between the surface cycle frame and the canonical
homology frame of :mod:`mtft.homology`, and frame-independent Hodge cross-checks.

Pi = L^T C is the pairing matrix between the canonical mslattice symbol basis L (29 x 26,
regenerated deterministically by PARI and equal to the frozen ``mtft.homology`` data) and the
Wave 10.1 cycle basis C (in Manin-generator coordinates).  Facts, all exact integers:

    Pi in GL(26, Z)  (det -1):  the canonical symbol lattice is the Z-dual of H_1(Z) in the cycle frame
    Pi W_cyc = W_L^T Pi          for W11, W13 (AL self-adjoint for the pairing)
    P_canonical = - Pi Jint^-1 Pi^T   (Poincare duality between the two intersection forms)

The periods frame of ``mtft.periods`` (frozen tau0) was rounded from float and carries no integer
map to the canonical frame, so it is compared through frame-independent invariants instead:
the j-invariants of the two 2-dimensional Atkin-Lehner sectors.  Result at N = 143:
(+,+) -> j(143a1) = -141.0134480904 in both frames (3e-14); (-,-) -> j = -4096/11 = j(11a3)
in both frames (2e-14): the (-,-) sector of H_1(X0(143), Z) is the curve 11a3.
"""
from __future__ import annotations

from typing import Dict, Optional, Tuple

import mpmath as mp
import numpy as np


def sector_j(J: np.ndarray, E: np.ndarray, W11: np.ndarray, W13: np.ndarray, a: int, b: int,
             dps: int = 30) -> Tuple[Optional[complex], Optional[complex], int]:
    """(tau, j, rank) of the saturated (a, b) Atkin-Lehner sector lattice with J restricted; rank must be 2."""
    from .hodge_structure import _integer_kernel_basis
    n = J.shape[0]
    I = np.eye(n, dtype=np.int64)
    K = _integer_kernel_basis(np.vstack([I - a * W11, I - b * W13]).astype(np.int64))
    if K.shape[1] != 2:
        return None, None, int(K.shape[1])
    v1, v2 = K[:, 0], K[:, 1]
    if float(v1 @ E @ v2) < 0:
        v1, v2 = v2, v1
    coef, *_ = np.linalg.lstsq(np.column_stack([v1, v2]), J @ v1, rcond=None)
    tau = (1j - coef[0]) / coef[1]
    if tau.imag < 0:
        tau = tau.conjugate()
    with mp.workdps(dps):
        j = complex(1728 * mp.kleinj(mp.mpc(tau.real, tau.imag)))
    return tau, j, 2


J_143A1 = -141.01344809037116
J_11A3 = -4096 / 11


def canonical_frame_gates() -> Dict[str, bool]:
    """Exact gates linking the frozen cycle frame to the frozen mtft.homology canonical frame."""
    import mtft.homology as H
    from .frozen import x0143
    d = x0143()
    Pi = d["Pi_canonical"].astype(np.int64)
    Jint = d["intersection_cycles"].astype(np.int64)
    m = H.matrices()
    can = {k: np.array(m[k], dtype=object).astype(np.int64) for k in ("W11", "W13", "P")}
    Ji = np.rint(np.linalg.inv(Jint.astype(float))).astype(np.int64)
    return {
        "Pi_unimodular": abs(int(round(np.linalg.det(Pi.astype(float))))) == 1,
        "Pi_intertwines_W11": bool(np.array_equal(Pi @ d["W11"].astype(np.int64), can["W11"].T @ Pi)),
        "Pi_intertwines_W13": bool(np.array_equal(Pi @ d["W13"].astype(np.int64), can["W13"].T @ Pi)),
        "poincare_duality_P_equals_minus_Pi_Jint_inv_PiT": bool(np.array_equal(-Pi @ Ji @ Pi.T, can["P"])),
    }


def cross_frame_hodge_check(dps: int = 40) -> Dict:
    """j-invariants of the (+,+) and (-,-) AL sectors in the cycle frame (J_true from periods) and in
    the mtft.periods frame (J from the frozen Riemann matrix tau0)."""
    import mtft.homology as H
    import mtft.periods as P
    from .frozen import x0143
    d = x0143()
    po = {k: np.array(v, dtype=object).astype(np.int64) for k, v in H.periods_frame_ops().items()}
    Jm = P.hodge_complex_structure(dps)
    Jp = np.array([[float(Jm[i, j]) for j in range(26)] for i in range(26)])
    Ep = np.array(H.standard_J(), dtype=object).astype(np.int64)
    out = {}
    for (a, b), label, ref in (((1, 1), "(+,+)", J_143A1), ((-1, -1), "(-,-)", J_11A3)):
        tc, jc, rc = sector_j(d["J_true"], d["intersection_cycles"].astype(np.int64), d["W11"].astype(np.int64), d["W13"].astype(np.int64), a, b)
        tp, jp, rp = sector_j(Jp, Ep, po["W11"], po["W13"], a, b)
        out[label] = {"rank": (rc, rp), "tau_cycle": tc, "tau_periods": tp, "j_cycle": jc, "j_periods": jp,
                      "cross_frame_rel": abs(jc - jp) / abs(ref) if jc is not None and jp is not None else None,
                      "reference_j": ref, "curve": "143a1" if label == "(+,+)" else "11a3"}
    return out
