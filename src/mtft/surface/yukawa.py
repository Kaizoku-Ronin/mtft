"""mtft.surface.yukawa — exact Yukawa tensors on X0(143) as canonical-ring multiplication (YUK-01, v0.29.0).

For a scalar in H^1(K^-3) coupling H^0(K^3) (charge +1 modes) to H^1(O) (charge −3 modes), Serre
duality makes the Yukawa tensor the multiplication map  H^0(K) (x) H^0(K^3) -> H^0(K^4)  of the
canonical ring.  Everything is computed GP-free from the frozen weight-2 basis
(``_data/x0143_weight2_basis.json``): AL-eigen basis by exact projectors, monomials by convolution,
ranks modulo 2^31−1, coordinates exactly over Q.  Gates: dims 60/84 (Riemann–Roch), AL sectors
(12,18,17,13)/(24,18,19,23) (Lefschetz), zero selection-rule violations, generic rank 13.
"""
from __future__ import annotations

import itertools, json
from math import lcm
from pathlib import Path
from typing import Dict, List

import numpy as np
from sympy import Matrix, Rational, QQ
from sympy.polys.matrices import DomainMatrix

_P = 2**31 - 1
_SECTORS = ((1, 1), (1, -1), (-1, 1), (-1, -1))


def load_basis() -> Dict:
    return json.loads((Path(__file__).parent / "_data" / "x0143_weight2_basis.json").read_text())


def al_eigenbasis(data: Dict, prec: int = 130):
    """Exact simultaneous W11/W13 eigenbasis of S_2(143): integer coefficient arrays, sector labels, scalings."""
    W11 = Matrix([[Rational(x) for x in r] for r in data["W11"]]); W13 = Matrix([[Rational(x) for x in r] for r in data["W13"]])
    F = [np.array(c[: prec + 1], dtype=object) for c in data["coefficients"]]; d = len(F); I = Matrix.eye(d)
    forms, labels, dens, vecs = [], [], [], []
    for sa, sb in _SECTORS:
        for v in ((I + sa * W11) * (I + sb * W13) / 4).columnspace():
            den = lcm(*[Rational(x).q for x in v]); vv = [int(Rational(x) * den) for x in v]
            forms.append(sum(vv[k] * F[k] for k in range(d))); labels.append((sa, sb)); dens.append(den); vecs.append(vv)
    return forms, labels, dens, vecs


def conv(a, b):
    n = len(a); out = [0] * n
    for i in range(1, n):
        if a[i] == 0: continue
        for j in range(1, n - i): out[i + j] += a[i] * b[j]
    return np.array(out, dtype=object)


def rank_mod(cols, nrows):
    M = np.array([[int(c[m]) % _P for c in cols] for m in range(1, nrows + 1)], dtype=np.int64)
    r, piv = 0, []
    for c in range(M.shape[1]):
        pr = next((i for i in range(r, M.shape[0]) if M[i, c] % _P), None)
        if pr is None: continue
        M[[r, pr]] = M[[pr, r]]; inv = pow(int(M[r, c]), _P - 2, _P); M[r] = (M[r] * inv) % _P
        for i in range(M.shape[0]):
            if i != r and M[i, c]: M[i] = (M[i] - M[i, c] * M[r]) % _P
        piv.append(c); r += 1
        if r == M.shape[0]: break
    return r, piv


def yukawa_tensor(verbose: bool = False) -> Dict:
    data = load_basis(); forms, labels, dens, vecs = al_eigenbasis(data); d = len(forms); PREC = len(forms[0]) - 1
    cub, cubl = [], []
    for i, j, k in itertools.combinations_with_replacement(range(d), 3):
        cub.append(conv(conv(forms[i], forms[j]), forms[k])); cubl.append((labels[i][0] * labels[j][0] * labels[k][0], labels[i][1] * labels[j][1] * labels[k][1]))
    rc, pivc = rank_mod(cub, PREC); cubbasis = [cub[c] for c in pivc]; cubsec = [cubl[c] for c in pivc]
    cub_triples = [list(itertools.combinations_with_replacement(range(d), 3))[c] for c in pivc]
    quart, quartl = [], []
    for a in range(d):
        for b in range(len(cubbasis)):
            quart.append(conv(forms[a], cubbasis[b])); quartl.append((labels[a][0] * cubsec[b][0], labels[a][1] * cubsec[b][1]))
    rq, pivq = rank_mod(quart, PREC); basis = [quart[c] for c in pivq]; bsec = [quartl[c] for c in pivq]
    sec_c = [rank_mod([cub[c] for c in range(len(cub)) if cubl[c] == s], PREC)[0] for s in _SECTORS]
    sec_q = [rank_mod([quart[c] for c in range(len(quart)) if quartl[c] == s], PREC)[0] for s in _SECTORS]
    # independent coefficient rows of the quartic basis, then exact solve
    Mt = np.array([[int(b[m]) % _P for m in range(1, PREC + 1)] for b in basis], dtype=np.int64)
    rows = rank_mod([np.array([0] + [int(b[m]) for m in range(1, PREC + 1)], dtype=object).T for b in []] or [], 1)[1] if False else None
    # (row selection: pivot columns of the transposed basis matrix)
    _, rows = rank_mod([np.array([0] + list(Mt[:, m - 1]), dtype=object) for m in range(1, PREC + 1)], len(basis))
    A = DomainMatrix([[QQ(int(basis[j][m + 1])) for j in range(len(basis))] for m in rows], (len(rows), len(basis)), QQ)
    B = DomainMatrix([[QQ(int(quart[c][m + 1])) for c in range(len(quart))] for m in rows], (len(rows), len(quart)), QQ)
    X = A.lu_solve(B).to_Matrix()
    nb, nq = len(cubbasis), len(basis)
    T = np.zeros((d, nb, nq), dtype=object); viol = 0; nz = 0
    for a in range(d):
        for b in range(nb):
            for j in range(nq):
                v = X[j, a * nb + b]; T[a, b, j] = v
                if v != 0:
                    nz += 1
                    if (labels[a][0] * cubsec[b][0], labels[a][1] * cubsec[b][1]) != bsec[j]: viol += 1
    return {"tensor": T, "omega_sectors": labels, "cubic_sectors": cubsec, "quartic_sectors": bsec, "cubic_triples": cub_triples,
            "quartic_pairs": [(c // nb, c % nb) for c in pivq], "forms": forms, "scalings": dens, "eigenvectors": vecs,
            "gates": {"cubic_rank": rc, "quartic_rank": rq, "cubic_sectors": sec_c, "quartic_sectors": sec_q, "selection_rule_violations": viol, "nonzero": nz}}


def yukawa_matrix(T, e) -> np.ndarray:
    """13 x 60 Yukawa matrix (= connecting map of the extension) for a scalar direction e in H^0(K^4)^* (length 84)."""
    Tf = np.vectorize(float)(T); return np.einsum("abj,j->ab", Tf, np.asarray(e, float))
