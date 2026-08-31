"""mtft.al_morphology — what Atkin-Lehner does to shape (v0.24.0).

The Atkin-Lehner involutions are usually used as *operators* (eigenvalue
bookkeeping on newforms).  This module reads them instead as *maps of
surfaces*: each W_Q is a holomorphic involution of X_0(143), so its trace on
H_1 determines a fixed-point count by Lefschetz, and that count determines
the genus of the quotient by Riemann-Hurwitz.

For X_0(143) the resulting cascade is

    g = 13  ->  g(X/W_11) = 7, g(X/W_13) = 6, g(X/W_143) = 2  ->  g(X*) = 1

computed here from the certified integral matrices, with the full-quotient
genus obtained two independent ways (Riemann-Hurwitz on the Klein group, and
the character formula dim Omega^G) as a gate.

Two structural facts fall out and are exposed as gates:

  * W_11 acts *freely* on X_0(143) (trace +2, zero fixed points), so
    X_0(143) -> X_0(143)/W_11 is an unramified double cover.
  * all three involutions act freely on the four cusps, which therefore form
    a single free (Z/2)^2-torsor.  Consequently the boundary lattice is the
    augmentation ideal of Z[G], and the 26 -> 29 puncture extension is the
    equivariant extension whose class lives in H^2(G, L_26).  This is the
    structural reason the puncture can carry an obstruction at all.

Fixed points are counted, not identified: the CM decomposition of the 20
fixed points of W_143 into class numbers is left open and flagged for audit.
"""
from __future__ import annotations

from functools import lru_cache

import numpy as np

from . import homology as HO
from . import periods as PP
from . import hecke as HK

__all__ = ["al_traces", "morphology", "cusp_torsor"]

_LABELS = (143, 13, 11, 1)   # model order; CC-15 corrected divisor labels


def _frame_ops():
    """AL involutions on H_1 in the hecke/K frame (the cuspidal_hecke frame)."""
    C = np.array(PP.hecke_to_symplectic_change(), dtype=object)
    Ci = np.array(HO.int_inverse(C), dtype=object)
    ops = HO.periods_frame_ops()
    W11 = Ci @ np.array(ops["W11"], dtype=object) @ C
    W13 = Ci @ np.array(ops["W13"], dtype=object) @ C
    return W11, W13, W11 @ W13


@lru_cache(maxsize=1)
def al_traces() -> dict:
    """Traces of W_11, W_13, W_143 on H_1(X_0(143), Z).  EXACT."""
    W11, W13, W143 = _frame_ops()
    return {11: int(np.trace(W11)), 13: int(np.trace(W13)),
            143: int(np.trace(W143))}


@lru_cache(maxsize=1)
def morphology() -> dict:
    """Fixed-point counts and quotient genera for every AL involution.

    Lefschetz for a holomorphic involution of a curve of genus g:
        #Fix(w) = 2 - tr(w | H_1),
    and Riemann-Hurwitz for the degree-2 quotient:
        g(X/w) = (2g + 2 - #Fix(w)) / 4.
    """
    g = 13
    tr = al_traces()
    fixed = {q: 2 - t for q, t in tr.items()}
    genera = {}
    for q, f in fixed.items():
        assert f >= 0 and f % 2 == 0, f"implausible fixed-point count {f}"
        num = 2 * g + 2 - f
        assert num % 4 == 0, f"Riemann-Hurwitz non-integral for W_{q}"
        genera[q] = num // 4
    total_fixed = sum(fixed.values())
    # Riemann-Hurwitz for the Klein group: 2g - 2 = |G|(2g* - 2) + sum_P (e_P - 1).
    # Each point fixed by exactly one involution has e_P = 2 and contributes 1.
    lhs = 2 * g - 2 - total_fixed          # = 4(2g* - 2) = 8g* - 8
    assert (lhs + 8) % 8 == 0, "Klein-quotient Riemann-Hurwitz non-integral"
    g_star_rh = (lhs + 8) // 8
    # Character formula: dim H^0(Omega)^G = (1/|G|) sum_h tr(h | H^0(Omega)),
    # and tr(h | H_1) = 2 tr(h | H^0(Omega)) for these real involutions.
    tr_sum = sum(tr.values())
    assert tr_sum % 2 == 0, "odd trace sum"
    inv4 = g + tr_sum // 2
    assert inv4 % 4 == 0, "character formula non-integral"
    g_star_chi = inv4 // 4
    return dict(
        genus=g, traces=tr, fixed_points=fixed, quotient_genera=genera,
        free_involutions=[q for q, f in fixed.items() if f == 0],
        total_fixed=total_fixed,
        genus_full_quotient=g_star_chi,
        genus_full_quotient_riemann_hurwitz=g_star_rh,
        methods_agree=(g_star_chi == g_star_rh),
        cascade=f"13 -> {{{genera[11]}, {genera[13]}, {genera[143]}}} "
                f"-> {g_star_chi}",
        open_items=["CM decomposition of the 20 fixed points of W_143 "
                    "into class numbers is not certified here"],
    )


@lru_cache(maxsize=1)
def cusp_torsor() -> dict:
    """Action of the AL group on the four cusps.

    Verifies that all three involutions act without fixed cusps, so the
    cusp set is a free (Z/2)^2-torsor and the boundary lattice is the
    augmentation ideal of Z[G].
    """
    def toggle(d, Q):
        dq = 1
        for p in (11, 13):
            if Q % p == 0 and d % p == 0:
                dq *= p
        return d // dq * (Q // dq)

    perms = {}
    for Q in (11, 13, 143):
        perms[Q] = {d: toggle(d, Q) for d in _LABELS}
    free = {Q: all(perms[Q][d] != d for d in _LABELS) for Q in perms}
    orbit = {perms[143][_LABELS[0]]}
    reach = {_LABELS[0]}
    for Q in (11, 13, 143):
        reach.add(perms[Q][_LABELS[0]])
    return dict(labels=list(_LABELS), permutations=perms, free=free,
                all_free=all(free.values()),
                single_orbit=(len(reach) == 4),
                boundary_module="augmentation ideal I_G of Z[G]")
