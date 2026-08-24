"""mtft.codifferent — the Canonical Codifferent Theorem at X0(143) (v0.19.0).

Certificates v3-v7: for each nontrivial newform orbit f of level 143 the
packaged rational orbit basis is, via the trace pairing
coefficient-of-q^n <-> Tr_{K/Q}(gamma a_n), an integral submodule
M subset O_K with pure {2,3} index, the saturated orbit lattice is the
codifferent D^-1, and the orbit saturation defect is

    [L_sat : L_pkg] = [O_K : M] * |Delta_K|.

Both eigenvalue fields are monogenic via a_2 (poldisc = nfdisc), shipped
here as frozen data with the exact power-basis coordinates of a_n.
Everything below recomputes from that data in pure Python.
"""
from __future__ import annotations

from math import gcd

import numpy as np

from .canonical import adapted_qexpansions, data_path

__all__ = ["ORBITS", "field_trace_table", "eigen_an", "gamma_table",
           "verify_orbit", "orbit_indices"]

ORBITS = {
    "f2": {
        "poly_low": [1, -1, -4, 0, 1],       # a_0..a_4, monic: y^4 - 4y^2 - y + 1
        "degree": 4,
        "disc": 1957,
        "columns": [8, 9, 10, 11],
        "an_file": "X0_143_f2_eigen_an.txt",
        "gamma": [[2, 0, 0, 0], [-4, -26, -2, 8],
                  [12, -18, -6, 6], [2, 28, 0, -6]],
        "index_OK": 576,                     # 2^6 3^2
    },
    "f3": {
        "poly_low": [-12, 7, 24, -2, -10, 0, 1],
        "degree": 6,
        "disc": 194616205,
        "columns": [1, 2, 3, 4, 5, 6],
        "an_file": "X0_143_f3_eigen_an.txt",
        "gamma": [[-28, 49, 21, -27, -3, 3], [56, -76, -54, 50, 8, -6],
                  [0, -15, -5, 9, 1, -1], [-16, 15, 7, -9, -1, 1],
                  [-40, 43, 35, -25, -5, 3], [-12, -20, 4, 16, 0, -2]],
        "index_OK": 2304,                    # 2^8 3^2
    },
}


def _newton_traces(poly_low, upto):
    """t_m = Tr(alpha^m), m = 0..upto, for monic poly (low-to-high coeffs)."""
    d = len(poly_low) - 1
    a = poly_low                              # a[k] coeff of x^k, a[d] = 1
    t = [d]
    for m in range(1, upto + 1):
        kmax = min(m - 1, d)
        s = sum(a[d - k] * t[m - k] for k in range(1, kmax + 1))
        if m <= d:
            s += m * a[d - m]
        t.append(-s)
    return t


def field_trace_table(orbit):
    d = ORBITS[orbit]["degree"]
    return _newton_traces(ORBITS[orbit]["poly_low"], d - 1)


def _polymulmod(u, v, poly_low):
    d = len(poly_low) - 1
    w = [0] * (len(u) + len(v) - 1)
    for i, ui in enumerate(u):
        if ui:
            for j, vj in enumerate(v):
                w[i + j] += ui * vj
    for k in range(len(w) - 1, d - 1, -1):
        c = w[k]
        if c:
            for j in range(d + 1):
                w[k - d + j] -= c * poly_low[j]
    return w[:d]


def eigen_an(orbit):
    """[a_0, ..., a_140] as power-basis integer coordinate lists."""
    out = []
    with open(data_path(ORBITS[orbit]["an_file"])) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            vals = [int(t) for t in line.strip().split(",")]
            out.append(vals[1:])
    return out


def gamma_table(orbit):
    return [list(g) for g in ORBITS[orbit]["gamma"]]


def verify_orbit(orbit):
    """Every packaged orbit column equals n -> Tr(gamma_i a_n), n <= 140."""
    info = ORBITS[orbit]
    A = np.array([[int(v) for v in row]
                  for row in np.array(adapted_qexpansions(), dtype=object)],
                 dtype=object)
    if A.shape[0] < A.shape[1]:
        A = A.T
    an = eigen_an(orbit)
    tr = field_trace_table(orbit)
    ok = True
    for gi, col in zip(info["gamma"], info["columns"]):
        for n in range(1, 141):
            prod = _polymulmod(gi, an[n], info["poly_low"])
            val = sum(c * t for c, t in zip(prod, tr))
            if val != int(A[n, col]):
                ok = False
    return ok


def orbit_indices(orbit):
    """{[O_K : M], [D^-1 : M]} from the gamma determinant and |Delta|."""
    info = ORBITS[orbit]
    d = info["degree"]
    G = [[info["gamma"][j][i] for j in range(d)] for i in range(d)]
    det = _int_det(G)
    return {"index_OK": abs(det),
            "index_codiff": abs(det) * info["disc"]}


def _int_det(M):
    from fractions import Fraction
    n = len(M)
    A = [[Fraction(x) for x in row] for row in M]
    det = Fraction(1)
    for c in range(n):
        pr = next((i for i in range(c, n) if A[i][c] != 0), None)
        if pr is None:
            return 0
        if pr != c:
            A[c], A[pr] = A[pr], A[c]
            det = -det
        det *= A[c][c]
        inv = 1 / A[c][c]
        A[c] = [x * inv for x in A[c]]
        for i in range(c + 1, n):
            if A[i][c] != 0:
                f = A[i][c]
                A[i] = [x - f * y for x, y in zip(A[i], A[c])]
    assert det.denominator == 1
    return int(det)
