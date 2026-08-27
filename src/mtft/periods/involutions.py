"""Involutions on the promoted X_0(143) homology stage.

Atkin-Lehner W_11, W_13, W_143 as frozen exact integer matrices on the 26D
``mtft.hecke`` cuspidal basis, the star involution transported to the
symplectic charge basis, the exact oldspace spectral projector, and the
sign/census decode connecting the homology stage to the canonical
differential sector pinning ((+,+),(+,-),(-,+),(-,-)) = (1,6,5,1).

Epistemic classes
-----------------
* frozen W matrices: integer-recognized from 80-dps periods (probe A1,
  distances ~4e-50), then EXACT-certified at call time by ``gates``:
  W^2=I, commuting, det +1, symplectic for E_H, [W,T_p]=0 (p=2,3,5,7),
  eigenspace dims (14,12,4), U_11=-W_11, im[U_13,W_13]=oldspace.
* sign decode / census / Route-2 intersections / star identities: EXACT
  Fraction linear algebra at call time.
"""
from __future__ import annotations

from fractions import Fraction as Fr
from functools import lru_cache
import json

from mtft import hecke as H
from .core import data_path, intersection_form
from .bridge import cuspidal_basis_change, hecke_to_symplectic_change

N = 26
BLOCK_DIMS = {"ell": 2, "ghost": 4, "q4": 8, "q6": 12}

# ------------------------------------------------------------ Fraction linalg
def _mul(A, B):
    k, m2 = len(B), len(B[0])
    return [[sum(A[i][t] * B[t][j] for t in range(k)) for j in range(m2)]
            for i in range(len(A))]


def _add(A, B, s=1):
    return [[A[i][j] + s * B[i][j] for j in range(len(A[0]))]
            for i in range(len(A))]


def _eye(n=N):
    return [[Fr(1 if i == j else 0) for j in range(n)] for i in range(n)]


def _rref(rows):
    M = [list(r) for r in rows]
    nr, nc = len(M), len(M[0]); piv = []; r = 0
    for c in range(nc):
        p = next((i for i in range(r, nr) if M[i][c]), None)
        if p is None:
            continue
        M[r], M[p] = M[p], M[r]
        q = M[r][c]; M[r] = [x / q for x in M[r]]
        for i in range(nr):
            if i != r and M[i][c]:
                q = M[i][c]
                M[i] = [a - q * b for a, b in zip(M[i], M[r])]
        piv.append(c); r += 1
        if r == nr:
            break
    return M, piv


def _rank(A):
    return len(_rref(A)[1])


def _nullspace(A):
    R, piv = _rref(A); nc = len(A[0])
    out = []
    for f in (c for c in range(nc) if c not in piv):
        v = [Fr(0)] * nc; v[f] = Fr(1)
        for r_, c_ in enumerate(piv):
            v[c_] = -R[r_][f]
        out.append(v)
    return out


def _inverse(A):
    n = len(A)
    M = [[Fr(x) for x in row] + [Fr(1 if i == j else 0) for j in range(n)]
         for i, row in enumerate(A)]
    R, piv = _rref(M)
    if piv[:n] != list(range(n)):
        raise ValueError("singular")
    return [r[n:] for r in R]


def _cols(vecs):
    return [[vecs[j][i] for j in range(len(vecs))] for i in range(N)]


def _fr(rows):
    return [[Fr(x) for x in r] for r in rows]


# ------------------------------------------------------------ frozen AL data
@lru_cache(maxsize=1)
def _al_record():
    return json.loads(data_path("X0_143_atkin_lehner_v022.json").read_text())


@lru_cache(maxsize=4)
def al_matrix(q: int):
    """W_q on the promoted 26D Hecke basis, EXACT integers; q in {11,13,143}."""
    rec = _al_record()
    if q == 11:
        return tuple(tuple(int(x) for x in r) for r in rec["W11"])
    if q == 13:
        return tuple(tuple(int(x) for x in r) for r in rec["W13"])
    if q == 143:
        P = _mul(_fr(rec["W11"]), _fr(rec["W13"]))
        return tuple(tuple(int(x) for x in r) for r in P)
    raise ValueError("q must be 11, 13, or 143")


@lru_cache(maxsize=1)
def transported_intersection():
    """E_H = C^T E C on the promoted Hecke basis, EXACT unimodular."""
    C = cuspidal_basis_change()
    E = intersection_form()
    Ct = [[Fr(C[i][j]) for i in range(N)] for j in range(N)]
    EH = _mul(_mul(Ct, _fr(E)), _fr(C))
    return tuple(tuple(int(x) for x in r) for r in EH)


# ------------------------------------------------------------ Hecke blocks
@lru_cache(maxsize=1)
def hecke_blocks():
    """Exact Fraction bases of the four T_2-charpoly blocks (dims 2,4,8,12)."""
    T2 = _fr([[Fr(x.numerator, x.denominator) for x in r]
              for r in H.cuspidal_hecke(2)])

    def poly(coeffs):
        out = [[Fr(0)] * N for _ in range(N)]
        for c in reversed(coeffs):
            out = _mul(out, T2)
            for i in range(N):
                out[i][i] += c
        return out

    blocks = {
        "ell": _nullspace(T2),
        "ghost": _nullspace(_add(T2, _eye(), 2)),
        "q4": _nullspace(poly(list(H.G4))),
        "q6": _nullspace(poly(list(H.H6))),
    }
    assert {k: len(v) for k, v in blocks.items()} == BLOCK_DIMS
    return blocks


@lru_cache(maxsize=1)
def al_signs():
    """Per-block eigenspace split of each W_q, plus the newform sign vectors.

    Returns ``{"split": {Wq: {block: (dim+, dim-)}},
               "eps":   {Wq: {"ell": s, "q4": s, "q6": s}},
               "ghost_W11": sign}``.
    EXACT.  eps(W11)=(+,-,+) and eps(W13)=(+,+,-) on (143a1, f2, f3);
    the ghost sits at W11=-1 and splits 2+2 under W13.
    """
    blocks = hecke_blocks()
    Ws = {q: _fr(al_matrix(q)) for q in (11, 13, 143)}
    split = {}
    for q, W in Ws.items():
        row = {}
        for b, vecs in blocks.items():
            B = _cols(vecs)
            plus = _rank(_mul(_add(W, _eye(), 1), B))
            minus = _rank(_mul(_add(W, _eye(), -1), B))
            assert plus + minus == BLOCK_DIMS[b]
            row[b] = (plus, minus)
        split[q] = row
    sgn = lambda pm, d: "+" if pm[0] == d else ("-" if pm[1] == d else "?")
    eps = {q: {b: sgn(split[q][b], BLOCK_DIMS[b]) for b in ("ell", "q4", "q6")}
           for q in (11, 13)}
    gp, gm = split[11]["ghost"]
    ghost11 = "+" if gp == 4 else ("-" if gm == 4 else "?")
    return {"split": split, "eps": eps, "ghost_W11": ghost11}


def sector_census():
    """Reconstruct the differential sector census from homology AL signs.

    Returns the tuple for ((+,+),(+,-),(-,+),(-,-)); EXACT (1, 6, 5, 1),
    matching the canonical quadric-side pinning and resolving the 5-sector
    as f2(4) + old_plus(1).
    """
    d = al_signs()
    cnt = {(a, b): 0 for a in "+-" for b in "+-"}
    cnt[(d["eps"][11]["ell"], d["eps"][13]["ell"])] += 1
    cnt[(d["eps"][11]["q4"], d["eps"][13]["q4"])] += 4
    cnt[(d["eps"][11]["q6"], d["eps"][13]["q6"])] += 6
    gp, gm = d["split"][13]["ghost"]
    assert (gp, gm) == (2, 2)
    cnt[(d["ghost_W11"], "+")] += 1
    cnt[(d["ghost_W11"], "-")] += 1
    return (cnt[("+", "+")], cnt[("+", "-")], cnt[("-", "+")], cnt[("-", "-")])


def route2_fixed_intersections():
    """dim(fix(W_143) ∩ block) per block: (ell,ghost,q4,q6)=(2,2,0,0), EXACT.

    Homology-side confirmation of Jac(X/W_143) ~ 143a1 x 11a1 (Route 2).
    """
    blocks = hecke_blocks()
    fix = _nullspace(_add(_fr(al_matrix(143)), _eye(), -1))
    assert len(fix) == 4
    F = _cols(fix)
    out = {}
    for b, vecs in blocks.items():
        B = _cols(vecs)
        joint = [[*B[i], *F[i]] for i in range(N)]
        out[b] = len(vecs) + 4 - _rank(joint)
    return out


# ------------------------------------------------------------ oldspace
@lru_cache(maxsize=1)
def oldspace_projector():
    """P_old = (U13^2 - I)(I - 2 U13)/90, EXACT Fractions, idempotent rank 4."""
    U = _fr([[Fr(x.numerator, x.denominator) for x in r]
             for r in H.cuspidal_hecke(13)])
    U2 = _mul(U, U)
    A = _add(U2, _eye(), -1)
    B = _add(_eye(), U, -2)
    P = [[x / 90 for x in row] for row in _mul(A, B)]
    assert _mul(P, P) == P and _rank(P) == 4
    return tuple(tuple(x for x in r) for r in P)


# ------------------------------------------------------------ star / charges
@lru_cache(maxsize=1)
def star_symplectic():
    """Star involution in symplectic charge coordinates, EXACT integers.

    Satisfies star^2 = I and star^T J star = -J (anti-symplectic), hence
    together with {star, J_H}=0 (CERT) the charge energy is star-invariant.
    """
    P = _fr(hecke_to_symplectic_change())
    S = _fr([[int(v) for v in r] for r in H.star_involution()])
    Ss = _mul(_mul(P, S), _inverse(P))
    assert all(x.denominator == 1 for r in Ss for x in r)
    return tuple(tuple(int(x) for x in r) for r in Ss)


def star_charge_orbit():
    """The exact 4-element star orbit of the minimal-shell charges.

    star(n3) = -(n5 - n6) and star(n5 - n6) = -n3 (1-based n-indices), so
    E(n3) = E(n5-n6) is symmetry-protected.  Returns the orbit as integer
    26-vectors in symplectic coordinates.
    """
    S = star_symplectic()
    n3 = [0] * N; n3[2] = 1
    n56 = [0] * N; n56[4] = 1; n56[5] = -1
    img = [sum(S[i][j] * n3[j] for j in range(N)) for i in range(N)]
    assert img == [-x for x in n56]
    img2 = [sum(S[i][j] * n56[j] for j in range(N)) for i in range(N)]
    assert img2 == [-x for x in n3]
    return (tuple(n3), tuple(-x for x in n56), tuple(-x for x in n3),
            tuple(n56))


__all__ = ["al_matrix", "transported_intersection", "hecke_blocks",
           "al_signs", "sector_census", "route2_fixed_intersections",
           "oldspace_projector", "star_symplectic", "star_charge_orbit",
           "BLOCK_DIMS"]
