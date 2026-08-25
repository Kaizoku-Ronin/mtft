"""mtft.origami.instances — the two certified instances.

galashin_24()   Example 2.18 / Figure 5: type (2,4), the smallest nontrivial
                origami.  Carries an explicit Kenyon-Smirnov t-embedding.
hexagonal_prism_36()  Figure 14 / Example 9.21: type (3,6), C3-symmetric,
                the instance with four distinct perfect t-embeddings.
"""
from __future__ import annotations

import numpy as np
import sympy as sp

from .dimer import DimerGraph

__all__ = [
    "galashin_24", "prism_36", "PRISM_C", "PRISM_LAMBDA0",
    "t_embedding_24", "mandelstams_24", "closed_curvature_B",
]

# --------------------------------------------------------------- (2,4) case
_POS24 = {"u1": (0.0, -2.0), "u2": (-2.0, 0.0), "u3": (0.0, 2.0),
          "u4": (2.0, 0.0), "bb": (0.0, -1.0), "wx": (-1.0, 0.0),
          "bt": (0.0, 1.0), "wy": (1.0, 0.0)}
_COLOR24 = {"u1": "w", "u2": "b", "u3": "w", "u4": "b",
            "bb": "b", "wx": "w", "bt": "b", "wy": "w"}
_FPOS24 = {"f1": (-1.4, -1.4), "f2": (-1.4, 1.4), "f3": (1.4, 1.4),
           "f4": (1.4, -1.4), "g": (0.0, 0.0)}
_FACES24 = {"wx": ("f1", "g", "f2"), "bt": ("f2", "g", "f3"),
            "wy": ("f3", "g", "f4"), "bb": ("f4", "g", "f1")}


def _edges24(p, q, r, s):
    return [("u1", "bb", 1, ("f1", "f4")),
            ("wx", "u2", 1, ("f1", "f2")),
            ("u3", "bt", 1, ("f2", "f3")),
            ("wy", "u4", 1, ("f3", "f4")),
            ("wx", "bb", p, ("f1", "g")),
            ("wx", "bt", q, ("f2", "g")),
            ("wy", "bt", r, ("f3", "g")),
            ("wy", "bb", -s, ("f4", "g"))]


def galashin_24(p=None, q=None, r=None, s=None):
    """The (2,4) DimerGraph with ORDINARY (unsigned) weights."""
    if p is None:
        p, q, r, s = sp.symbols("p q r s", positive=True)
    edges = [(u, v, abs(w) if not hasattr(w, "free_symbols") else w)
             for (u, v, w, _f) in _edges24(p, q, r, s)]
    edges = [(u, v, sp.Abs(w) if w == -s else w) for (u, v, w) in edges]
    edges = [(u, v, s if (u, v) == ("wy", "bb") else w) for (u, v, w) in edges]
    return DimerGraph(_COLOR24, edges, ["u1", "u2", "u3", "u4"])


def _ang(u, v):
    return float(np.arctan2(u[0] * v[1] - u[1] * v[0],
                            u[0] * v[0] + u[1] * v[1]))


def t_embedding_24(p, q, r, s):
    """Kenyon-Smirnov primitive: dict face -> (x_hat, x_check), plus (lam, lamt).

    Returns (emb, lam, lamt, residual).  ``residual`` is the cycle-closure
    residual; it must vanish for the primitive to be well defined (local
    closure => global integrability).
    """
    Fo = {"u1": -p + s * 1j, "wx": 1.0 + 0j, "u3": -q - r * 1j, "wy": 1j}
    Ft = {"u2": p + q * 1j, "bb": -1.0 + 0j, "u4": -s + r * 1j, "bt": -1j}
    lam = np.array([[p, 1, q, 0], [-s, 0, r, 1]], float)
    lamt = np.array([[1, -p, 0, s], [0, q, -1, r]], float)
    eqs = []
    for w, b, K, (fa, fb) in _edges24(p, q, r, s):
        mw, mb = np.array(_POS24[w]), np.array(_POS24[b])
        d = mb - mw
        mid = (mw + mb) / 2
        v = np.array(_FPOS24[fa]) - mid
        L, R = (fa, fb) if (d[0] * v[1] - d[1] * v[0]) > 0 else (fb, fa)
        eqs.append((L, R, Fo[w] * K * Ft[b], np.conj(Fo[w]) * K * Ft[b]))
    faces = ["f1", "f2", "f3", "f4", "g"]
    idx = {f: i for i, f in enumerate(faces)}
    A = np.zeros((len(eqs) + 1, 5))
    bh = np.zeros(len(eqs) + 1, complex)
    bc = np.zeros(len(eqs) + 1, complex)
    for k, (L, R, dh, dc) in enumerate(eqs):
        A[k, idx[L]], A[k, idx[R]] = 1, -1
        bh[k], bc[k] = dh, dc
    A[-1, idx["f1"]] = 1
    xh, *_ = np.linalg.lstsq(A, bh, rcond=None)
    xc, *_ = np.linalg.lstsq(A, bc, rcond=None)
    res = float(np.linalg.norm(A @ xh - bh) + np.linalg.norm(A @ xc - bc))
    return {f: (xh[idx[f]], xc[idx[f]]) for f in faces}, lam, lamt, res


def mandelstams_24(p, q, r, s):
    """Planar Mandelstams via Lemma 1.10, with the metric-defect check.

    Returns dict (i,j) -> S(i,j).  At this instance S(1,3) = pr and
    S(2,4) = qs EXACTLY: the two terms of the Gr(2,4) Plucker exchange
    relation Delta_13 Delta_24 = Delta_12 Delta_34 + Delta_14 Delta_23 are
    individually the two scattering channels.
    """
    emb, lam, lamt, res = t_embedding_24(p, q, r, s)
    assert res < 1e-9, f"Kenyon-Smirnov closure failed: {res}"

    def br(M, i, j):
        return float(np.linalg.det(np.stack([M[:, i], M[:, j]], 1)))

    out = {}
    for (i, j), (pp, qq) in {(1, 3): (1, 2), (2, 4): (2, 3)}.items():
        S = br(lam, pp, qq) * br(lamt, pp, qq)
        lhs = (abs(emb[f"f{i}"][0] - emb[f"f{j}"][0]) ** 2
               - abs(emb[f"f{i}"][1] - emb[f"f{j}"][1]) ** 2)
        assert abs(lhs - 4 * S) < 1e-9, (i, j, lhs, 4 * S)
        out[(i, j)] = S
    return out


def closed_curvature_B():
    """Certified closed forms for the (2,4) section-B Fisher geometry.

    Returns (X, Y, c, det_g, K).  With c = (1+q)(1+s) the odds ratio and
    X = p, Y = r:

        det g = XY(XY + cX + cY + c) / (c + X + Y + XY)^3
        K     = -(c-1)(c-X^2)(c-Y^2) / (4(XY + cX + cY + c)^2)

    Zero-curvature walls X = sqrt(c), Y = sqrt(c) are the zero-effective-field
    loci of the equivalent two-spin Ising model (J = (1/4) log c); the
    complement involution (X,Y) -> (c/X, c/Y) is the global spin flip.
    """
    X, Y, c = sp.symbols("X Y c", positive=True)
    D = X * Y + c * X + c * Y + c
    det_g = X * Y * D / (c + X + Y + X * Y) ** 3
    K = -(c - 1) * (c - X ** 2) * (c - Y ** 2) / (4 * D ** 2)
    return X, Y, c, det_g, K


# --------------------------------------------------------------- (3,6) case
def prism_36(heavy=2, offset=0):
    """Hexagonal prism of type (3,6); ``heavy`` on an alternating inner triple.

    18 vertices (6 boundary spokes, 6 outer, 6 inner), 24 edges.  With
    heavy = 2 this is Galashin Figure 14 / Example 9.21, which admits four
    distinct perfect t-embeddings.
    """
    heavy = sp.sympify(heavy)
    color, edges = {}, []
    for i in range(6):
        color[("u", i)] = "w" if i % 2 == 0 else "b"
        color[("o", i)] = "b" if i % 2 == 0 else "w"
        color[("n", i)] = "w" if i % 2 == 0 else "b"
    for i in range(6):
        edges.append((("u", i), ("o", i), sp.Integer(1)))
        edges.append((("o", i), ("o", (i + 1) % 6), sp.Integer(1)))
        edges.append((("o", i), ("n", i), sp.Integer(1)))
        w = heavy if i % 2 == offset else sp.Integer(1)
        edges.append((("n", i), ("n", (i + 1) % 6), w))
    return DimerGraph(color, edges, [("u", i) for i in range(6)])


#: certified rational representative of Meas(prism) in Gr_{>0}(3,6),
#: scaled by Delta_{012} = 14.  All 20 maximal minors positive.
PRISM_C = np.array([[14, 0, 0, 2, 21, 13],
                    [0, 14, 0, -14, -49, -21],
                    [0, 0, 14, 10, 14, 2]], dtype=float) / 14.0

#: the exact C3-fixed perfect branch (integer representative).  Certified:
#: Theta(lam0) C^T = 0 exactly; brackets 20,25,20,25,20,25; t_i alternate
#: -3/100, -7/100; wind(lam0) = 2 pi, wind(Theta lam0) = 4 pi.
PRISM_LAMBDA0 = np.array([[7, 5, 0, -4, -7, -1],
                          [-4, 0, 5, 3, -1, -3]], dtype=float)
