#!/usr/bin/env python3
"""
x0143_graph_uncertainty.py — where uncertainty lives on the house graph
============================================================================

MIT License — Copyright (c) 2026 Roger Tano

The v0.13.0 ancestry theorem  [D, L]_{uv} = A_{uv} (d_v - d_u)  (mtft.
combinatorial, Pr/EXACT: noncommutativity is generated exactly by degree
gradients across edges) evaluated on the house's own dual graph — the
Farey skeleton of X0(143).

CONSTRUCTION (all exact integer / permutation work)
---------------------------------------------------
Flags = P1(Z/143), |P1| = (11+1)(13+1) = 168 = [PSL2(Z) : Gamma0(143)],
with the right actions  S: (c:d) -> (d:-c),  T: (c:d) -> (c:c+d),
R = ST: (c:d) -> (d:d-c).  S and R act freely (eps_2 = eps_3 = 0).
Triangles of the Farey tessellation on X0(143) = R-orbits (56 of size
3); edge gluings = S-orbits (84 of size 2); cusps = T-orbits, with
orbit size = cusp width: {1, 11, 13, 143} (cusps oo, 1/13, 1/11, 0).
Dual graph: vertex per triangle, edge per S-orbit joining two
triangles; an S-orbit inside one triangle is a self-loop.

RESULTS (pre-registered prediction -> confirmed)
------------------------------------------------
R1 (EXACT).  The skeleton has V = 56, 83 proper edges + exactly ONE
  self-loop, no multi-edges, connected, b1 = 84 - 56 + 1 = 29, and
  loop-dropped tr L = 166 — reproducing the du-corpus anchors.

R2 (EXACT).  Kirchhoff via fraction-free Bareiss on the 55x55 reduced
  Laplacian: kappa = 3 518 081 582 959 364 640, matching the archived
  spanning-tree count exactly.  This certifies the object built here IS
  the corpus dual graph (independent route: P1-flag construction vs the
  original du-session construction).

R3 (Pr/EXACT).  As a ribbon MULTIGRAPH (self-loop counted, degree 3
  everywhere) the skeleton is 3-regular, so by the ancestry theorem
  [D, L] = 0 identically: the skeleton is commutative.

R4 (Pr/EXACT).  Under the standard simple-Laplacian convention the
  loop is dropped; the loop triangle becomes the UNIQUE degree defect
  (degree 1 vs 3), and the ENTIRE noncommutativity of the house graph
  concentrates on the single edge from the loop triangle to its
  neighbor, with entries [D,L] = +/-2.  Two disjoint routes
  (matrix products vs the degree-gradient closed form) agree.

R5 (EXACT).  Geometric identification: the self-loop is the fan of the
  WIDTH-1 CUSP AT INFINITY closing on itself in one step.  The loop
  triangle is the base triangle {0, 1, oo}, flags
  {(0:1), (1:1), (1:0)}, corner widths [1, 143, 143]; its neighbor has
  corner widths [143, 143, 143].  Graph uncertainty on X0(143) is the
  shadow of the cusp at infinity.

R6 (Pr + Cert).  ROBERTSON CAPACITY = 1.  For any unit state,
  <[D,L]> = 4i Im(conj(psi_a) psi_b) on the defect edge (a, b), so the
  Robertson lower bound is 2|Im(conj(psi_a) psi_b)| <= 1, with
  equality at |psi_a| = |psi_b| = 1/sqrt(2), quarter phase.  The
  maximal uncertainty the modular skeleton can enforce is exactly 1.
  Certified numerically (maximizer hits 1.000000000000; 500 random
  states never exceed it; all Robertson margins nonnegative).

STRUCTURAL RHYME, FILED AS RHYME ONLY (AG-D5): du03 found cusp wells
annihilating the harmonic stage; here the cusp at infinity is the sole
source of skeleton noncommutativity.  Same locus, different theorems;
no derivation link is claimed.

GATES
-----
G1 flag space: |P1| = 168; S, R free; T-widths {1, 11, 13, 143}
G2 skeleton anchors: V=56, E=83+1 loop, simple, connected, b1=29,
   tr L = 166
G3 spanning trees (Bareiss, exact int) == 3518081582959364640
G4 commutator two routes agree; support = exactly one edge at the loop
   vertex; entries +/-2
G5 multigraph 3-regular and commutator identically zero
G6 loop triangle = {(0:1),(1:1),(1:0)}, corner widths [1,143,143]
G7 Robertson: maximizer bound = 1 (1e-12), random margins >= -1e-9,
   random bounds <= 1 + 1e-9

Run:  python studies/x0143_graph_uncertainty.py     (~15 s)
Writes x0143_graph_uncertainty_ledger.json next to itself.
"""

from __future__ import annotations

import json
import os
import sys
from math import gcd

import numpy as np

from mtft import combinatorial as C

N = 143
ARCHIVED_TREES = 3518081582959364640


# ── exact flag space and actions ────────────────────────────────────

def build_flags():
    units = [u for u in range(1, N) if gcd(u, N) == 1]

    def valid(c, d):
        return not (c % 11 == 0 and d % 11 == 0) and \
               not (c % 13 == 0 and d % 13 == 0)

    def canon(c, d):
        return min(((c * u) % N, (d * u) % N) for u in units)

    P1 = sorted({canon(c, d) for c in range(N) for d in range(N)
                 if valid(c, d)})
    idx = {p: i for i, p in enumerate(P1)}

    def perm(f):
        return [idx[canon(*f(*p))] for p in P1]

    sS = perm(lambda c, d: (d % N, (-c) % N))
    sT = perm(lambda c, d: (c, (c + d) % N))
    sR = perm(lambda c, d: (d % N, (d - c) % N))
    return P1, sS, sT, sR


def orbits(sigma):
    seen, out = set(), []
    for i in range(len(sigma)):
        if i in seen:
            continue
        o, j = [], i
        while j not in seen:
            seen.add(j)
            o.append(j)
            j = sigma[j]
        out.append(o)
    return out


def bareiss_det(M):
    """Exact integer determinant, fraction-free Bareiss."""
    M = [row[:] for row in M]
    n, prev, sign = len(M), 1, 1
    for k in range(n - 1):
        if M[k][k] == 0:
            for r in range(k + 1, n):
                if M[r][k]:
                    M[k], M[r] = M[r], M[k]
                    sign = -sign
                    break
            else:
                return 0
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                M[i][j] = (M[i][j] * M[k][k] - M[i][k] * M[k][j]) // prev
        prev = M[k][k]
    return sign * M[n - 1][n - 1]


def main() -> int:
    ledger = {"study": "x0143_graph_uncertainty", "gates": {}}
    ok = True

    def gate(name, passed, **info):
        nonlocal ok
        ok &= bool(passed)
        ledger["gates"][name] = {"passed": bool(passed), **info}
        print(f"[{'PASS' if passed else 'FAIL'}] {name}  " +
              "  ".join(f"{k}={v}" for k, v in info.items()))

    # G1 flag space
    P1, sS, sT, sR = build_flags()
    widths_all = orbits(sT)
    widths = sorted(len(o) for o in widths_all)
    gate("G1_flags", len(P1) == 168
         and all(sS[i] != i for i in range(168))
         and all(sR[i] != i for i in range(168))
         and widths == [1, 11, 13, 143],
         n_flags=len(P1), widths=widths)
    cusp_w = {}
    for o in widths_all:
        for f in o:
            cusp_w[f] = len(o)

    # skeleton
    tris = orbits(sR)
    tri_of = {}
    for ti, t in enumerate(tris):
        for f in t:
            tri_of[f] = ti
    V = len(tris)
    A = np.zeros((V, V), dtype=np.int64)
    loops, proper = [], 0
    for o in orbits(sS):
        a, b = tri_of[o[0]], tri_of[sS[o[0]]]
        if a == b:
            loops.append(a)
        else:
            A[a, b] += 1
            A[b, a] += 1
            proper += 1
    deg = A.sum(1)
    seen, stack = {0}, [0]
    while stack:
        v = stack.pop()
        for u in np.nonzero(A[v])[0]:
            if int(u) not in seen:
                seen.add(int(u))
                stack.append(int(u))

    # G2 anchors
    gate("G2_skeleton", V == 56 and proper == 83 and len(loops) == 1
         and int(A.max()) == 1 and len(seen) == V
         and proper + len(loops) - V + 1 == 29
         and int(deg.sum()) == 166,
         V=V, proper_edges=proper, loops=len(loops),
         b1=proper + len(loops) - V + 1, trL=int(deg.sum()))

    # G3 spanning trees, exact
    L = np.diag(deg) - A
    kappa = bareiss_det([[int(L[i, j]) for j in range(1, V)]
                         for i in range(1, V)])
    gate("G3_spanning_trees", kappa == ARCHIVED_TREES,
         kappa=kappa, archived=ARCHIVED_TREES)
    ledger["spanning_trees"] = kappa

    # G4 the theorem, loop-dropped convention
    comm = C.commutator_DL(A)
    two_routes = np.array_equal(comm, C.degree_gradient_matrix(A))
    nz = [(int(u), int(v), int(comm[u, v]))
          for u, v in np.argwhere(comm != 0)]
    lv = loops[0]
    support_ok = (len(nz) == 2 and sorted(abs(x[2]) for x in nz) == [2, 2]
                  and all(lv in (u, v) for u, v, _ in nz))
    gate("G4_defect_edge", two_routes and support_ok,
         nonzero=nz, loop_vertex=lv, loop_degree=int(deg[lv]))
    nb = nz[0][1] if nz[0][0] == lv else nz[0][0]

    # G5 multigraph regularity
    A2 = A.copy()
    A2[lv, lv] = 2
    gate("G5_multigraph_commutes",
         bool(np.all(A2.sum(1) == 3)) and not C.commutator_DL(A2).any())

    # G6 geometric identification
    flags = sorted(P1[f] for f in tris[lv])
    corners = sorted(cusp_w[f] for f in tris[lv])
    gate("G6_cusp_at_infinity",
         flags == [(0, 1), (1, 0), (1, 1)] and corners == [1, 143, 143],
         loop_triangle_flags=flags, corner_widths=corners,
         neighbor_corner_widths=sorted(cusp_w[f] for f in tris[nb]))

    # G7 Robertson capacity
    D = np.diag(deg).astype(float)
    Lf = D - A
    psi = np.zeros(V, complex)
    psi[lv], psi[nb] = 1 / np.sqrt(2), 1j / np.sqrt(2)
    r = C.robertson_margin(D, Lf, psi)
    rng = np.random.default_rng(N)
    margins, bounds = [], []
    for _ in range(500):
        p = rng.standard_normal(V) + 1j * rng.standard_normal(V)
        rr = C.robertson_margin(D, Lf, p)
        margins.append(rr["margin"])
        bounds.append(rr["bound"])
    gate("G7_capacity_one", abs(r["bound"] - 1.0) < 1e-12
         and min(margins) >= -1e-9 and max(bounds) <= 1 + 1e-9,
         maximizer_bound=f"{r['bound']:.12f}",
         worst_margin=f"{min(margins):.2e}",
         max_random_bound=f"{max(bounds):.4f}")

    ledger["all_passed"] = ok
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "x0143_graph_uncertainty_ledger.json")
    with open(out, "w") as f:
        json.dump(ledger, f, indent=2, default=str)
    print(f"\nledger -> {out}")
    print("ALL GATES PASS" if ok else "GATE FAILURE")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
