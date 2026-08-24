"""Gates for `mtft.canonical.integral` — the 2026-08-24 arc, recomputed live.

Every headline number of certificates v1-v9 is re-derived from the frozen
q-expansions at call time.  Nothing is asserted from memory.
"""
from __future__ import annotations

from math import gcd

import numpy as np

from .. import integral_lattice as lat
from . import MONOMIALS, adapted_qexpansions, ideal_basis, ideal_basis_adapted
from .integral import (LEDGER, SATURATION_PRIMES, adapted_matrix,
                       ci_a_codifferent, count_points_modp, cusp_reductions,
                       points_modp, quadratic_saturation_obstruction,
                       al_splitting, al_denominator, saturated_qexpansions,
                       sector_columns)


def _residual(A, qvec):
    res = np.zeros(282, dtype=object)
    for m in range(91):
        c = int(qvec[m])
        if c:
            i, j = MONOMIALS[m]
            conv = np.convolve(A[:, i], A[:, j])
            res[:len(conv)] += c * conv
    return max(abs(int(x)) for x in res[:141])


def gate_frame():
    """ideal_basis is s2-frame; ideal_basis_adapted is adapted-frame."""
    Aad = adapted_matrix()
    import numpy as _np
    from . import s2_qexpansions
    As2 = _np.array([[int(v) for v in row]
                     for row in _np.array(s2_qexpansions(), dtype=object)],
                    dtype=object)
    if As2.shape[0] < As2.shape[1]:
        As2 = As2.T
    Qs2 = np.array(ideal_basis(), dtype=object)
    Qad = np.array(ideal_basis_adapted(), dtype=object)
    assert _residual(As2, Qs2[:, 0]) == 0 and _residual(As2, Qs2[:, 54]) == 0
    assert _residual(Aad, Qs2[:, 0]) != 0, "frame mismatch marker"
    assert _residual(Aad, Qad[:, 0]) == 0 and _residual(Aad, Qad[:, 54]) == 0
    assert lat.rank_modp(Qad, 2) == 55 and lat.rank_modp(Qad, 3) == 55
    return {"ideal_basis_frame": "s2", "ideal_basis_adapted_frame": "adapted"}


def gate_integral_model():
    """The Integral Model Gate: rank_{F_p} = g after saturation, with the
    certified step ledger."""
    B, steps = saturated_qexpansions()
    assert steps == LEDGER["full_saturation_steps"], steps
    for p in SATURATION_PRIMES:
        assert lat.rank_modp(B, p) == 13
    return {"saturation_steps": steps}


def gate_counts_mod2():
    got = {m: count_points_modp(2, m)
           for m in ("saturated", "packaged_s2", "adapted_mixed")}
    assert got == LEDGER["counts_mod2"], got
    return got


def gate_counts_mod3():
    n = count_points_modp(3, "saturated")
    assert n == LEDGER["counts_mod3"]["saturated"], n
    return {"saturated": n}


def gate_cusps(p):
    cusps = cusp_reductions(p)
    pts = set(points_modp(p))
    norm = set()
    for t in cusps.values():
        norm.add(t)
    covered = all(any(tuple((s * np.array(t)) % p) in pts
                      for s in range(1, p)) or t in pts for t in norm)
    assert covered and len(norm) == len(pts) == 4, (cusps, pts)
    return {"cusps": cusps, "bijection": True}


def gate_ci_a():
    a, quad = ci_a_codifferent()
    assert a == LEDGER["ci_a_codifferent_a"] == -(7 ** 2) * 13
    A = adapted_matrix()
    pairs = [(0, 0)] + [(0, j) for j in range(8, 12)] \
        + [(i, j) for i in range(8, 12) for j in range(i, 12)]
    cols = [[int(x) for x in np.convolve(A[:, i], A[:, j])[:141]]
            for (i, j) in pairs]
    ker = lat.rational_kernel(np.array(cols, dtype=object).T)
    assert len(ker) == 1
    v = ker[0]
    den = 1
    for e in v:
        den = den * e.denominator // gcd(den, e.denominator)
    vi = [int(e * den) for e in v]
    g = 0
    for e in vi:
        g = gcd(g, abs(e))
    vi = [e // g for e in vi]
    if vi[0] > 0:
        vi = [-e for e in vi]
    assert vi[0] == LEDGER["ci_a_packaged_a"] == -(7 ** 2) * 13 * 1957 ** 2
    return {"a_packaged": vi[0], "a_codifferent": a}


def gate_product_chain():
    Lmp = sector_columns("(-,+)")
    inv_full = quadratic_saturation_obstruction(Lmp, (7, 13))
    assert inv_full == LEDGER["Q_full_sector"], inv_full
    A = adapted_matrix()
    F2, _ = lat.saturate(A[:, [8, 9, 10, 11]], (2, 3, 19, 103))
    inv_f2 = quadratic_saturation_obstruction(F2, (7, 13))
    assert inv_f2 == LEDGER["Q_f2_codifferent"], inv_f2
    f1 = A[:, 0]
    g = 0
    for v in f1:
        g = gcd(g, abs(int(v)))
    f1p = np.array([int(v) // g for v in f1], dtype=object)
    f1sq = [int(x) for x in np.convolve(f1p, f1p)[:141]]
    pairs5 = [(i, j) for i in range(5) for j in range(i, 5)]
    M = np.array([[int(x) for x in np.convolve(Lmp[:, i], Lmp[:, j])[:141]]
                  for (i, j) in pairs5], dtype=object).T
    H = lat.hnf(M)
    assert lat.class_order(H, f1sq, [1, 7, 13]) == 13
    pairs4 = [(i, j) for i in range(4) for j in range(i, 4)]
    M4 = np.array([[int(x) for x in np.convolve(F2[:, i], F2[:, j])[:141]]
                   for (i, j) in pairs4], dtype=object).T
    H4 = lat.hnf(M4)
    assert lat.class_order(H4, f1sq, [1, 7, 13, 49, 91, 637]) == 637
    assert 637 // 13 == LEDGER["product_index_full_over_f2"] == 49
    from .. import eisenstein
    assert eisenstein.eisenstein_modulus("f2_quartic") == 49
    return {"Q_full": inv_full, "Q_f2": inv_f2, "index": 49}


def gate_al_splitting():
    q = al_splitting()
    assert q["four_sector"] == LEDGER["al_four_sector_snf"], q
    assert q["W11_split"] == LEDGER["al_W11_split_snf"], q
    assert q["W13_in_W11plus"] == LEDGER["al_W13_in_W11plus"], q
    assert q["W13_in_W11minus"] == LEDGER["al_W13_in_W11minus"], q
    assert al_denominator(11) == 1 and al_denominator(13) == 13
    return q
