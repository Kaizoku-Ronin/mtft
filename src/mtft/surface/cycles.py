"""mtft.surface.cycles — EXACT layer: deterministic integral basis of H_1(X0(N), Z).

Lexicographic BFS primal tree on the cusp graph, BFS dual tree on the face
graph; the 2g leftover edges close to fundamental cycles.  Gate: the chord
matrix [omitted-face boundaries | leftover cycles] has Bareiss determinant
+-1 in the cycle lattice, so the basis is integral (not a rational
complement).  The basis is reproducible from N but not intrinsic.
"""
from __future__ import annotations

from collections import deque
from dataclasses import dataclass
from typing import List

import numpy as np

from .manin import ManinComplex


def bareiss_det(matrix: List[List[int]]) -> int:
    if not matrix:
        return 1
    a = [row[:] for row in matrix]
    n = len(a)
    if any(len(r) != n for r in a):
        raise ValueError("square matrix required")
    sign, prev = 1, 1
    for k in range(n - 1):
        if a[k][k] == 0:
            sw = next((r for r in range(k + 1, n) if a[r][k] != 0), None)
            if sw is None:
                return 0
            a[k], a[sw] = a[sw], a[k]
            sign = -sign
        piv = a[k][k]
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                num = a[i][j] * piv - a[i][k] * a[k][j]
                if num % prev:
                    raise ArithmeticError("non-exact Bareiss division")
                a[i][j] = num // prev
            a[i][k] = 0
        prev = piv
    return sign * a[-1][-1]


@dataclass
class CycleBasis:
    primal_tree: List[int]
    dual_tree: List[int]
    leftover_edges: List[int]
    basis_matrix: np.ndarray      # E x 2g integer
    unimodular_determinant: int
    gates: List[dict]

    @property
    def rank(self) -> int:
        return self.basis_matrix.shape[1]


def tree_cotree(cx: ManinComplex) -> CycleBasis:
    V, E, F = cx.counts
    adj: List[List[tuple]] = [[] for _ in range(V)]
    for eid, d in enumerate(cx.edges):
        a, b = cx.vertex_of[d], cx.vertex_of[cx.S(d)]
        adj[a].append((b, eid, 1))
        adj[b].append((a, eid, -1))
    for row in adj:
        row.sort(key=lambda t: (t[1], t[0], t[2]))

    primal, seen, q = set(), {0}, deque([0])
    while q:
        v = q.popleft()
        for nb, eid, _ in adj[v]:
            if nb not in seen:
                seen.add(nb)
                primal.add(eid)
                q.append(nb)

    dadj: List[List[tuple]] = [[] for _ in range(F)]
    for eid, d in enumerate(cx.edges):
        if eid in primal:
            continue
        l, r = cx.face_of[d], cx.face_of[cx.S(d)]
        dadj[l].append((r, eid))
        dadj[r].append((l, eid))
    for row in dadj:
        row.sort(key=lambda t: (t[1], t[0]))
    dual, seenf, q = set(), {0}, deque([0])
    while q:
        f = q.popleft()
        for nb, eid in dadj[f]:
            if nb not in seenf:
                seenf.add(nb)
                dual.add(eid)
                q.append(nb)

    leftovers = sorted(set(range(E)) - primal - dual)
    tadj = [[t for t in row if t[1] in primal] for row in adj]

    def path(start: int, target: int) -> np.ndarray:
        par = {start: None}
        q = deque([start])
        while q:
            v = q.popleft()
            if v == target:
                break
            for nb, eid, s in tadj[v]:
                if nb not in par:
                    par[nb] = (v, eid, s)
                    q.append(nb)
        vec = np.zeros(E, dtype=np.int64)
        v = target
        while v != start:
            pv, eid, s = par[v]
            vec[eid] += s
            v = pv
        return vec

    B = np.zeros((E, len(leftovers)), dtype=np.int64)
    for col, eid in enumerate(leftovers):
        d = cx.edges[eid]
        a, b = cx.vertex_of[d], cx.vertex_of[cx.S(d)]
        cyc = np.zeros(E, dtype=np.int64)
        cyc[eid] = 1
        cyc += path(b, a)
        B[:, col] = cyc

    chords = sorted(set(range(E)) - primal)
    square = []
    for eid in chords:
        row = [int(cx.boundary_2[eid, f]) for f in range(1, F)]
        row.extend(1 if eid == lo else 0 for lo in leftovers)
        square.append(row)
    det = bareiss_det(square)
    g2 = cx.inv.b1_compact
    gates = [
        {"name": "primal_tree_spans", "class": "EXACT", "status": "PASS" if len(seen) == V else "FAIL"},
        {"name": "dual_tree_spans", "class": "EXACT", "status": "PASS" if len(seenf) == F else "FAIL"},
        {"name": "tree_cotree_rank_is_2g", "class": "EXACT", "status": "PASS" if len(leftovers) == g2 else "FAIL"},
        {"name": "cycles_are_closed", "class": "EXACT", "status": "PASS" if not np.any(cx.boundary_1 @ B) else "FAIL"},
        {"name": "cycle_lattice_change_is_unimodular", "class": "EXACT", "status": "PASS" if abs(det) == 1 else "FAIL"},
    ]
    return CycleBasis(sorted(primal), sorted(dual), leftovers, B, det, gates)
