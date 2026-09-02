"""mtft.surface.hodge — discrete Hodge layer on the Manin complex.

CERTIFIED_NUMERICAL(tol): harmonic nullity 2g, intersection form recovered
from period-dual Whitney representatives (integral, skew, det +-1, later
cross-certified by exact Hecke/AL transport in :mod:`transport`).

DIAGNOSTIC: positive Laplace spectra.  The Whitney metric is piecewise
equilateral on the Manin triangles.  By Voevodsky-Shabat that metric is in
the conformal class of X0(N) (the Manin triangulation is the dessin of the
Belyi map X0(N) -> X(1)), so Hodge-theoretic objects (harmonic classes,
1-form star, periods) of the refinement family converge to the true ones,
while the Laplace spectrum is NOT that of the hyperbolic metric.  The
exact/coexact bracket is the non-vacuous convergence criterion for the
spectrum (both branches share one continuum limit on a closed surface).

The polar complex structure ``polar_star`` is a compatible complex structure
for (metric, wedge).  ``star^2 = -I`` is a tautology of the construction and
is not a gate; distance to ``hodge_structure.J_true`` is.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Dict, List, Tuple

import numpy as np
from scipy.linalg import eigh

from .cycles import CycleBasis, bareiss_det
from .manin import ManinComplex


# ------------------------------------------------------------ unweighted
def unweighted_hodge(cx: ManinComplex, cb: CycleBasis, tol: float = 1e-10) -> Dict:
    b1 = cx.boundary_1.astype(float)
    b2 = cx.boundary_2.astype(float)
    L1 = b1.T @ b1 + b2 @ b2.T
    w, v = np.linalg.eigh(L1)
    H = v[:, np.abs(w) <= tol]
    n = cb.rank
    if H.shape[1] != n:
        raise ArithmeticError(f"N={cx.N}: nullity {H.shape[1]} != 2g={n}")
    C = cb.basis_matrix.astype(float)
    Pi = H.T @ C
    D = H @ np.linalg.inv(Pi.T)                      # periods on cycles = identity
    cup = D.T @ cx.wedge_1 @ D
    cup = (cup - cup.T) / 2
    Icup = np.rint(cup).astype(np.int64)
    rounding = float(np.max(np.abs(cup - Icup)))
    det = bareiss_det(Icup.tolist())
    pos = w[np.abs(w) > tol]
    return {
        "class": "CERTIFIED_NUMERICAL",
        "tolerance": tol,
        "nullity": int(H.shape[1]),
        "harmonic": H,
        "period_dual": D,
        "intersection_dual": Icup,          # on the cohomology basis dual to the cycles
        "intersection_cycles": np.rint(np.linalg.inv(Icup)).astype(np.int64),
        "intersection_rounding": rounding,
        "intersection_det": int(det),
        "first_positive_eigenvalue_DIAGNOSTIC": float(pos[0]) if len(pos) else None,
        "gates": [
            {"name": "harmonic_nullity_is_2g", "status": "PASS"},
            {"name": "intersection_skew_unimodular", "status": "PASS" if (np.array_equal(Icup.T, -Icup) and abs(det) == 1) else "FAIL"},
            {"name": "intersection_rounding_below_tol", "status": "PASS" if rounding < tol else "FAIL"},
        ],
    }


# --------------------------------------------------------------- Whitney
@dataclass
class Mesh:
    vertex_count: int
    edges: List[Tuple[int, int]]
    faces: List[List[Tuple[int, int]]]
    level: int = 0


def base_mesh(cx: ManinComplex) -> Mesh:
    edges = [(cx.vertex_of[d], cx.vertex_of[cx.S(d)]) for d in cx.edges]
    faces = [[cx.edge_of[d] for d in o] for o in cx.faces]
    return Mesh(len(cx.vertices), edges, faces, 0)


def refine(mesh: Mesh) -> Mesh:
    """Split every edge, quarter every triangle; quotient-aware (loops, multi-edges)."""
    mid = [mesh.vertex_count + e for e in range(len(mesh.edges))]
    edges, halves = [], []
    for e, (a, b) in enumerate(mesh.edges):
        k = len(edges)
        edges.extend(((a, mid[e]), (mid[e], b)))
        halves.append((k, k + 1))
    faces = []
    for face in mesh.faces:
        paths = []
        for e, s in face:
            f, sc = halves[e]
            paths.append(([(f, 1), (sc, 1)] if s == 1 else [(sc, -1), (f, -1)], mid[e]))
        (p01, m01), (p12, m12), (p20, m20) = paths
        i01 = len(edges); edges.append((m01, m12))
        i12 = len(edges); edges.append((m12, m20))
        i20 = len(edges); edges.append((m20, m01))
        faces.extend(([p01[0], (i20, -1), p20[1]], [p12[0], (i01, -1), p01[1]],
                      [p20[0], (i12, -1), p12[1]], [(i01, 1), (i12, 1), (i20, 1)]))
    return Mesh(mesh.vertex_count + len(mesh.edges), edges, faces, mesh.level + 1)


def lift_cycles(C: np.ndarray, from_edges: int, to_edges: int) -> np.ndarray:
    """Cycle coefficients on a mesh -> its refinement (both halves inherit, inner edges 0)."""
    out = np.zeros((to_edges, C.shape[1]))
    for e in range(from_edges):
        out[2 * e] = C[e]
        out[2 * e + 1] = C[e]
    return out


def wedge_matrix(mesh: Mesh) -> np.ndarray:
    n = len(mesh.edges)
    W = np.zeros((n, n))
    for face in mesh.faces:
        for i, j in ((0, 1), (1, 2), (2, 0)):
            (ei, si), (ej, sj) = face[i], face[j]
            W[ei, ej] += si * sj / 6
            W[ej, ei] -= si * sj / 6
    return W


def _reference_masses(scale: float):
    area = math.sqrt(3) * scale * scale / 4
    M0 = area / 12 * np.array([[2, 1, 1], [1, 2, 1], [1, 1, 2]], float)
    grads = np.array([[-1, -1 / math.sqrt(3)], [1, -1 / math.sqrt(3)], [0, 2 / math.sqrt(3)]]) / scale
    coef = []
    for a, b in ((0, 1), (1, 2), (2, 0)):
        c = np.zeros((3, 2))
        c[a] += grads[b]
        c[b] -= grads[a]
        coef.append(c)
    bp = np.full((3, 3), area / 12)
    np.fill_diagonal(bp, area / 6)
    M1 = np.array([[sum(np.dot(coef[i][a], coef[j][b]) * bp[a, b] for a in range(3) for b in range(3))
                    for j in range(3)] for i in range(3)])
    return M0, M1, 1 / area


def assemble(mesh: Mesh) -> Dict[str, np.ndarray]:
    V, E, F = mesh.vertex_count, len(mesh.edges), len(mesh.faces)
    b1 = np.zeros((V, E)); b2 = np.zeros((E, F))
    M0 = np.zeros((V, V)); M1 = np.zeros((E, E)); M2d = np.zeros(F)
    l0, l1, l2 = _reference_masses(0.5 ** mesh.level)
    for e, (a, b) in enumerate(mesh.edges):
        b1[a, e] -= 1; b1[b, e] += 1
    for f, face in enumerate(mesh.faces):
        verts = []
        for e, s in face:
            a, b = mesh.edges[e]
            verts.append(a if s == 1 else b)
            b2[e, f] += s
        for i, vi in enumerate(verts):
            for j, vj in enumerate(verts):
                M0[vi, vj] += l0[i, j]
        for i, (ei, si) in enumerate(face):
            for j, (ej, sj) in enumerate(face):
                M1[ei, ej] += si * sj * l1[i, j]
        M2d[f] = l2
    d0, d1 = b1.T, b2.T
    K = M1 @ d0 @ np.linalg.solve(M0, d0.T @ M1) + (d1.T * M2d) @ d1
    return {"b1": b1, "b2": b2, "M0": M0, "M1": M1, "M2d": M2d, "K1": (K + K.T) / 2}


def weighted_harmonics(A: Dict[str, np.ndarray], n: int, tol: float = 1e-8) -> Tuple[np.ndarray, np.ndarray]:
    """M1-orthonormal harmonic basis (E x n) and the positive spectrum head."""
    vals, vecs = eigh(A["K1"], A["M1"], subset_by_index=(0, min(len(A["M1"]) - 1, n + 20)), driver="gvx")
    H = vecs[:, np.abs(vals) <= tol]
    if H.shape[1] != n:
        raise ArithmeticError(f"weighted nullity {H.shape[1]} != {n}")
    R = np.linalg.cholesky(H.T @ A["M1"] @ H)
    H = H @ np.linalg.inv(R.T)
    return H, vals[np.abs(vals) > tol]


def branch_split(A: Dict[str, np.ndarray], k: int = 4) -> Dict[str, List[float]]:
    """Exact (P1 scalar, conforming upper bound) vs coexact (dual) first eigenvalues."""
    d0, d1, M0, M1 = A["b1"].T, A["b2"].T, A["M0"], A["M1"]
    l0 = eigh(d0.T @ M1 @ d0, M0, eigvals_only=True)
    l0 = l0[l0 > 1e-8][:k]
    l2 = eigh(d1 @ np.linalg.solve(M1, d1.T), np.diag(1 / A["M2d"]), eigvals_only=True)
    l2 = l2[l2 > 1e-8][:k]
    return {"exact": [float(x) for x in l0], "coexact": [float(x) for x in l2],
            "bracket_ratio": float(max(l0[0], l2[0]) / min(l0[0], l2[0]))}


# ----------------------------------------------------------- polar star
def polar_star_on_cycles(H: np.ndarray, wedge: np.ndarray, C: np.ndarray, Jint: np.ndarray):
    """Compatible complex structure of (metric for which H is orthonormal, wedge),
    transported to cycle coordinates via the period pairing, sign fixed by Riemann
    positivity against the cycle intersection form ``Jint``.
    Returns (J_cycles, G = Jint J symmetric positive definite, wedge singular values)."""
    W = H.T @ wedge @ H
    W = (W - W.T) / 2
    pos, vec = np.linalg.eigh(-W @ W)
    if np.min(pos) <= 0:
        raise ArithmeticError("degenerate harmonic wedge form")
    J = W @ (vec @ np.diag(1 / np.sqrt(pos)) @ vec.T)
    Pi = H.T @ C
    Jc = np.linalg.solve(Pi, J.T @ Pi)
    G = Jint @ Jc
    G = (G + G.T) / 2
    if np.linalg.eigvalsh(G)[0] < 0:
        Jc, G = -Jc, -G
    return Jc, G, np.sqrt(np.sort(pos))


def siegel_distance(G1: np.ndarray, G2: np.ndarray) -> Tuple[float, float]:
    """Invariant distance between two complex structures compatible with one
    symplectic form, via their metrics G_i = Jint J_i: eigenvalues of G1^-1 G2
    pair as (l, 1/l); returns (sqrt(sum log^2 l), max |log l|)."""
    lam = eigh(G2, G1, eigvals_only=True)
    logs = np.log(lam)
    return float(np.sqrt(np.sum(logs ** 2))), float(np.max(np.abs(logs)))


def refinement_family(cx: ManinComplex, cb: CycleBasis, Jint: np.ndarray, max_level: int = 2) -> List[Dict]:
    """Polar stars of the Whitney family r = 0..max_level in cycle coordinates."""
    n = cb.rank
    mesh = base_mesh(cx)
    C = cb.basis_matrix.astype(float)
    out = []
    for r in range(max_level + 1):
        if r:
            new = refine(mesh)
            C = lift_cycles(C, len(mesh.edges), len(new.edges))
            mesh = new
        A = assemble(mesh)
        H, pos = weighted_harmonics(A, n)
        Jc, G, sv = polar_star_on_cycles(H, wedge_matrix(mesh), C, Jint)
        chi = mesh.vertex_count - len(mesh.edges) + len(mesh.faces)
        out.append({"level": r, "cells": (mesh.vertex_count, len(mesh.edges), len(mesh.faces)),
                    "euler": chi, "J": Jc, "G": G, "wedge_singular_values": sv,
                    "first_positive_eigenvalue_DIAGNOSTIC": float(pos[0]),
                    "branch": branch_split(A)})
    return out
