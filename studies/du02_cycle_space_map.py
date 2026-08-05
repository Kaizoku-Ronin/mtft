#!/usr/bin/env python3
"""
du02_cycle_space_map.py — dynamical units, session 2:
                          the cycle-space map (drawn-loop mathematics)
======================================================================
Roger Tano / MTFT Research Program — built with Claude, August 2026

WHY THIS STUDY. du01 froze the two clocks' unit-free content and left
the anchor count at 2: no certified internal process couples graph time
to Hecke time. The spine identity b1(dual graph) = 29 = 2g + (c-1)
located the only geometric channel: the graph clock's cycle space
contains the Hecke clock's arena. This study COMPUTES that channel.

THE THREE SPACES (all derived from P^1(Z/143), zero parameters):
  cycle space  Z = ker(d1) in dual-edge space, dim 29
                 = H1 of the open surface X0(143) minus its 4 cusps
  symbols      M = Manin space / (sigma, tau relations), dim 29
                 = H1(X, cusps)  (relative homology; Hecke acts here)
  homology     H1(X) = Z / (cusp links)  =  cuspidal subspace of M,
                 dim 26 (both presentations)

STAGES
  A  The crossing pairing. Each dual edge (i -> j across symbol x)
     meets the oriented geodesic of x once; pairing = +1 on x, -1 on
     x.sigma. Certificates (all EXACT, integer arithmetic):
       A1  pairing of every basis cycle kills every sigma- and
           tau-relation (validates the global orientation convention);
       A2  LEFSCHETZ DUALITY: the 29 x 29 pairing matrix D between an
           integer cycle basis and the Manin quotient basis has
           det D = +-1 (perfect, unimodular). Poincare-Lefschetz
           duality for X0(143) computed from raw combinatorics.
  B  Cusp links. For each cusp c the unique cycle b_c with
     <b_c, s> = (coefficient of c in the boundary of s) for all s.
     Certificates: b_c integral (EXACT, forced by det D = +-1);
     sum_c b_c = 0; rank of the link lattice = 3; the width-1 cusp's
     link IS the self-loop edge (EXACT identification — explains du01's
     self-loop); empirical width census over all 168 cosets.
  C  Transported Hecke. T2 exact on M (engine, Fractions); transport
     T2* = D T2^T D^{-1} onto cycle space. Certificates: T2* is
     INTEGRAL (Hecke preserves the integer lattice of cycles); T2*
     preserves the link lattice; induced action on Z/links has the
     SAME characteristic polynomial as cuspidal T2 (E2: the engine's
     exact restrict_to_cuspidal vs the transported route); the 3-dim
     Eisenstein block's char poly is reported.
  D  Hodge decomposition of edge space. d1 = incidence (56 x 84),
     d2 = link matrix (84 x 4). Delta_1 = d1^T d1 + d2 d2^T.
     Certificates: dim ker Delta_1 = 26; the nonzero spectrum is
     EXACTLY {55 nonzero graph-clock levels} u {3 link-Gram levels}
     (E2 vs du01's spectrum); trace re-derivation 166 = 2*83.
  E  The structural theorem for dynamical units:
       On the shared 26-dim stage (harmonic homology), the FREE graph
       clock is frozen (kernel of Delta_1, EXACT) while the Hecke
       clock runs all 13 lines x2. Therefore NO free-evolution
       exchange rate chi_H/chi_g exists — the du01/H4 honest negative
       is upgraded from "not found" to "structurally excluded at free
       level". An internal exchange rate requires an INTERACTION that
       lifts the harmonic degeneracy (the du03 dispersion program).
  F  Interface export: everything the drawn-loop React stage needs
     (layout, edges, pairing rows, boundary, links, Hecke, spectral
     signature basis with Hecke-time phases), written as JSON.

Deferred honestly: the loop-x-loop Goldman intersection number on the
ribbon graph (needs its own sign audit — session 3, with Kimi on the
conventions). The loop-x-symbol pairing shipped here IS the
intersection pairing of a drawn loop with the 29 geodesic arcs.

Epistemic classes on every number: EXACT / Cert / DIAGNOSTIC.
"""

from __future__ import annotations
from fractions import Fraction
import json
import os
import time

import numpy as np
import sympy as sp

from x0143_particle_box import (P1, SIGMA, TAU, tessellation, dual_graph,
                                ModularSymbols, cusp_divisor)

_HERE = os.path.dirname(os.path.abspath(__file__))
LEDGER: dict = {"study": "du02_cycle_space_map", "version_context": "0.11.3"}
N = 143


# ----------------------------------------------------------------------
# exact integer linear algebra (Bareiss; shared with du01)
# ----------------------------------------------------------------------

def bareiss(M):
    M = [row[:] for row in M]
    n, m = len(M), (len(M[0]) if M else 0)
    prev, rank, sign = 1, 0, 1
    for col in range(m):
        piv = None
        for r in range(rank, n):
            if M[r][col] != 0:
                piv = r
                break
        if piv is None:
            continue
        if piv != rank:
            M[rank], M[piv] = M[piv], M[rank]
            sign = -sign
        p = M[rank][col]
        for r in range(rank + 1, n):
            for c in range(col + 1, m):
                M[r][c] = (M[r][c] * p - M[r][col] * M[rank][c]) // prev
            M[r][col] = 0
        prev = p
        rank += 1
        if rank == n:
            break
    det = sign * prev if (n == m and rank == n) else 0
    return rank, det


def int_solve_unimodular(D, rhs):
    """Solve D^T a = rhs exactly (D integer unimodular, rhs integer);
    returns integer vector. Uses Fractions then certifies integrality."""
    n = len(D)
    A = [[Fraction(D[j][i]) for j in range(n)] for i in range(n)]  # D^T
    b = [Fraction(x) for x in rhs]
    # Gaussian elimination over Q
    for col in range(n):
        piv = next(r for r in range(col, n) if A[r][col] != 0)
        A[col], A[piv] = A[piv], A[col]
        b[col], b[piv] = b[piv], b[col]
        pv = A[col][col]
        for r in range(n):
            if r != col and A[r][col] != 0:
                f = A[r][col] / pv
                for c in range(col, n):
                    A[r][c] -= f * A[col][c]
                b[r] -= f * b[col]
    out = [b[i] / A[i][i] for i in range(n)]
    assert all(v.denominator == 1 for v in out), "non-integral solve"
    return [int(v) for v in out]


# ----------------------------------------------------------------------
# Stage A — crossing pairing and Lefschetz duality
# ----------------------------------------------------------------------

def build_all():
    p1, tris, edges, cert = tessellation(N)
    A, edge_list, tri_of = dual_graph(p1, tris, edges)
    ms = ModularSymbols(p1)
    assert ms.dim == 29 and ms.cuspidal_dim == 26
    return p1, tris, edges, A, edge_list, tri_of, ms


def crossing_matrix(p1, edge_list, sign=+1):
    """84 x 168 signed crossings: oriented dual edge (i->j across x)
    pairs +sign with x, -sign with x.sigma."""
    X = [[0] * len(p1.reps) for _ in range(len(edge_list))]
    for e, (i, j, x) in enumerate(edge_list):
        xs = p1.canon[tuple(np.array(p1.act(x, SIGMA)))] \
            if False else p1.act(x, SIGMA)
        X[e][p1.index[x]] += sign
        X[e][p1.index[xs]] -= sign
    return X


def incidence(edge_list, n_tri=56):
    """Signed vertex-edge incidence d1: column e has -1 at tail i,
    +1 at head j; self-loop column is zero."""
    B = [[0] * len(edge_list) for _ in range(n_tri)]
    for e, (i, j, _x) in enumerate(edge_list):
        if i == j:
            continue
        B[i][e] -= 1
        B[j][e] += 1
    return B


def cycle_basis(edge_list, n_tri=56):
    """Integer fundamental-cycle basis from a spanning tree of the
    multigraph. Returns K (29 x 84 int) and the tree edge set."""
    adj = {v: [] for v in range(n_tri)}
    for e, (i, j, _x) in enumerate(edge_list):
        if i != j:
            adj[i].append((j, e))
            adj[j].append((i, e))
    parent = {0: (None, None)}
    order = [0]
    tree = set()
    for v in order:
        for (w, e) in adj[v]:
            if w not in parent:
                parent[w] = (v, e)
                tree.add(e)
                order.append(w)
    assert len(parent) == n_tri, "dual graph not connected"

    def tree_path(u):  # signed edge chain from root 0 to u
        chain = {}
        while parent[u][0] is not None:
            v, e = parent[u]
            i, j, _x = edge_list[e]
            # edge oriented i->j; traversing v->u
            s = +1 if (i, j) == (v, u) else -1
            chain[e] = chain.get(e, 0) + s
            u = v
        return chain

    K = []
    for e, (i, j, _x) in enumerate(edge_list):
        if e in tree:
            continue
        row = [0] * len(edge_list)
        row[e] = 1                      # traverse e as i -> j
        if i != j:
            for ee, s in tree_path(i).items():
                row[ee] += s            # root -> i
            for ee, s in tree_path(j).items():
                row[ee] -= s            # j -> root
        K.append(row)
    assert len(K) == 29
    return K, tree


def stage_A(p1, edge_list, ms):
    print("=" * 70)
    print("STAGE A — crossing pairing, relations, LEFSCHETZ DUALITY")
    print("=" * 70)

    for sign in (+1, -1):
        X = crossing_matrix(p1, edge_list, sign)
        K, tree = cycle_basis(edge_list)
        # A1: every basis cycle must kill every relation, exactly.
        ok = True
        rel_rows = []
        seen = set()
        n = len(p1.reps)
        for x in p1.reps:
            i = p1.index[x]
            j = p1.index[p1.act(x, SIGMA)]
            key = ("s", tuple(sorted((i, j))))
            if key not in seen:
                seen.add(key)
                r = [0] * n
                r[i] += 1
                r[j] += 1
                rel_rows.append(r)
            a = p1.index[p1.act(x, TAU)]
            b = p1.index[p1.act(p1.act(x, TAU), TAU)]
            key = ("t", tuple(sorted((i, a, b))))
            if key not in seen:
                seen.add(key)
                r = [0] * n
                r[i] += 1
                r[a] += 1
                r[b] += 1
                rel_rows.append(r)
        F = [[sum(K[c][e] * X[e][m] for e in range(84)) for m in range(n)]
             for c in range(29)]
        worst = max(abs(sum(F[c][m] * r[m] for m in range(n)))
                    for c in range(29) for r in rel_rows)
        if worst == 0:
            print(f"  A1 orientation convention sign={sign:+d}: all "
                  f"{len(rel_rows)} relations killed by all 29 basis "
                  f"cycles, worst residue 0  PASS (EXACT)")
            ok = True
            break
        else:
            print(f"  A1 sign={sign:+d}: worst relation residue {worst} "
                  f"— flipping global orientation")
            ok = False
    assert ok, "no global orientation works — convention bug"
    LEDGER["A1 relations_killed (EXACT)"] = True
    LEDGER["A1 crossing_sign"] = sign

    # A2: the 29 x 29 pairing in quotient coordinates; det = +-1.
    free = ms.free
    D = [[F[c][free[k]] for k in range(29)] for c in range(29)]
    _r, det = bareiss(D)
    print(f"  A2 Lefschetz pairing D (cycles x symbols, 29 x 29): "
          f"det D = {det}  "
          f"{'PASS' if abs(det) == 1 else 'FAIL'} (EXACT — perfect, "
          f"unimodular duality)")
    assert abs(det) == 1
    LEDGER["A2 det_D (EXACT)"] = det
    return X, K, tree, F, D


# ----------------------------------------------------------------------
# Stage B — cusp links
# ----------------------------------------------------------------------

def stage_B(p1, edge_list, ms, K, D):
    print("=" * 70)
    print("STAGE B — cusp links: the 3-dim kernel, realized")
    print("=" * 70)
    # boundary matrix rows in quotient coords (4 x 29), exact integers
    Bd = [[int(ms.boundary[i, j]) for j in range(29)]
          for i in range(len(ms.divisors))]
    assert all(all(Fraction(ms.boundary[i, j]).denominator == 1
                   for j in range(29)) for i in range(len(ms.divisors)))

    # width census over all 168 cosets: terminal cusp of each symbol
    tally = {d: 0 for d in ms.divisors}
    for x in p1.reps:
        (a, b), (c0, d0) = p1.lift(x)
        tally[cusp_divisor(N, a, c0)] += 1
    print(f"  B1 terminal-end census over 168 cosets by divisor: "
          f"{tally}  (EXACT; the four cusp widths)")
    LEDGER["B1 width_census (EXACT)"] = {str(k): v for k, v in tally.items()}

    links = {}
    for ci, d in enumerate(ms.divisors):
        a = int_solve_unimodular(D, Bd[ci])          # coeffs in K-basis
        vec = [sum(a[c] * K[c][e] for c in range(29)) for e in range(84)]
        links[d] = vec
    ssum = [sum(links[d][e] for d in ms.divisors) for e in range(84)]
    print(f"  B2 links integral (forced by det D = +-1): PASS (EXACT); "
          f"sum of the 4 links = 0: "
          f"{'PASS' if all(v == 0 for v in ssum) else 'FAIL'} (EXACT)")
    assert all(v == 0 for v in ssum)
    Lmat = [links[d] for d in ms.divisors]
    rank_links, _ = bareiss([row[:] for row in Lmat])
    print(f"  B3 link lattice rank = {rank_links}  "
          f"{'PASS' if rank_links == 3 else 'FAIL'} (EXACT) — "
          f"H1 = cycles/links has dim 26")
    assert rank_links == 3

    loop_e = next(e for e, (i, j, _x) in enumerate(edge_list) if i == j)
    ident = None
    for d in ms.divisors:
        v = links[d]
        if all((abs(v[e]) == (1 if e == loop_e else 0)) for e in range(84)):
            ident = d
    print(f"  B4 the SELF-LOOP edge is exactly the link of the divisor-"
          f"{ident} cusp (width {tally[ident]})  "
          f"{'PASS' if ident is not None and tally[ident] == 1 else 'FAIL'}"
          f" (EXACT) — du01's self-loop explained: one triangle wraps "
          f"the width-1 cusp")
    assert ident is not None and tally[ident] == 1
    LEDGER["B4 selfloop_is_width1_link (EXACT)"] = str(ident)
    supp = {str(d): sum(1 for v in links[d] if v != 0) for d in ms.divisors}
    print(f"  B5 link supports (#edges used): {supp}  (EXACT, reported)")
    LEDGER["B5 link_supports (EXACT)"] = supp
    return Bd, links, Lmat


# ----------------------------------------------------------------------
# Stage C — transported Hecke
# ----------------------------------------------------------------------

def stage_C(ms, D, Lmat):
    print("=" * 70)
    print("STAGE C — the Hecke clock carried onto the graph side")
    print("=" * 70)
    t0 = time.time()
    T2 = ms.hecke_on_quotient(2)                    # exact 29 x 29
    print(f"  C1 exact T2 on the 29-dim symbol space: {time.time()-t0:.1f}s")
    T2i = [[Fraction(sp.Rational(T2[i, j]).p, sp.Rational(T2[i, j]).q)
            for j in range(29)] for i in range(29)]

    # geometric transport through the perfect pairing <v, s> = v^T D s:
    # the self-transpose Hecke correspondence acts on cycles by the
    # pairing adjoint  T2* = D^{-T} T2^T D^T
    Dm = sp.Matrix(D)
    DT = Dm.T
    T2star = DT.inv() * T2.T * DT
    ints = all(sp.Rational(T2star[i, j]).q == 1
               for i in range(29) for j in range(29))
    print(f"  C2 T2* integral on the cycle lattice: "
          f"{'PASS' if ints else 'FAIL'} (EXACT — Hecke preserves the "
          f"integer cycles)")
    assert ints
    LEDGER["C2 T2star_integral (EXACT)"] = True

    # link lattice invariance: T2* rows of links stay in link span
    Lsp = sp.Matrix([[Fraction(v) for v in row] for row in Lmat])
    # coordinates of links in K-basis: solve  a D-route? links were built
    # in edge coords; convert: link = sum a_c K_c with a from stage B —
    # recompute coefficients through D^T a = Bd row is equivalent; here
    # verify invariance in the 29-dim K-coordinate system directly:
    # a_link rows:
    return T2, T2star


def stage_C2(ms, D, Bd, T2, T2star):
    # link coefficient vectors in K-basis
    A_link = [int_solve_unimodular(D, Bd[ci]) for ci in range(4)]
    Asp = sp.Matrix([[Fraction(v) for v in row] for row in A_link])
    M = sp.Matrix(T2star)
    img = Asp * M.T          # images of link vectors under T2*
    aug = Asp.col_join(img)
    inv_ok = (Asp.rank() == aug.rank() == 3)
    print(f"  C3 T2* preserves the link lattice (Eisenstein sector): "
          f"{'PASS' if inv_ok else 'FAIL'} (EXACT)")
    assert inv_ok

    # Eisenstein block: T2* on span(links) in the basis of the first
    # three links (the fourth is minus their sum). Exact via normal
    # equations on the unimodular lattice coordinates.
    basis = Asp[0:3, :]                       # 3 x 29
    G = basis * basis.T                       # 3 x 3, invertible
    E3 = (G.LUsolve(basis * img[0:3, :].T)).T # rows: images in basis
    cp_e = sp.factor(E3.charpoly(sp.Symbol('x')).as_expr())
    print(f"  C4 Eisenstein 3-block char poly: {cp_e}  (EXACT) "
          f"[prediction: (x - 3)^3 — both degeneracy images of a cusp "
          f"share its divisor class, multiplicity 1 + p]")
    LEDGER["C4 eisenstein_charpoly"] = str(cp_e)

    # E2: char poly of T2* on cycles/links  vs engine cuspidal T2
    t0 = time.time()
    x = sp.Symbol('x')
    cp_full = sp.factor(sp.Matrix(T2star).charpoly(x).as_expr())
    A26, _B = ms.restrict_to_cuspidal(T2)
    cp_cusp = sp.factor(A26.charpoly(x).as_expr())
    quot = sp.cancel(cp_full / (cp_cusp * cp_e))
    ok = sp.simplify(quot - 1) == 0 or sp.degree(quot, x) == 0
    print(f"  C5 char poly factorization ({time.time()-t0:.1f}s):")
    print(f"     charpoly(T2*) = charpoly(cuspidal T2) * "
          f"charpoly(Eis 3-block)")
    print(f"     residual factor: {sp.simplify(quot)}  "
          f"{'PASS' if ok else 'FAIL'} (EXACT — the transported clock "
          f"carries the SAME 13 lines; E2 route pair: duality transport "
          f"vs engine restriction)")
    assert ok
    LEDGER["C5 charpoly_match (EXACT)"] = True
    LEDGER["C5 cuspidal_charpoly"] = str(cp_cusp)
    return cp_cusp, cp_e


# ----------------------------------------------------------------------
# Stage D — Hodge decomposition of edge space
# ----------------------------------------------------------------------

def stage_D(edge_list, links, ms):
    print("=" * 70)
    print("STAGE D — Hodge decomposition: the shared 26-dim stage")
    print("=" * 70)
    B1 = np.array(incidence(edge_list), dtype=float)          # 56 x 84
    B2 = np.array([links[d] for d in ms.divisors], float).T   # 84 x 4
    D1 = B1.T @ B1 + B2 @ B2.T
    ev = np.linalg.eigvalsh(D1)
    nker = int(np.sum(np.abs(ev) < 1e-9))
    print(f"  D1 dim ker Delta_1 = {nker}  "
          f"{'PASS' if nker == 26 else 'FAIL'} (Cert) — the harmonic "
          f"homology; 26 = 2g")
    assert nker == 26

    L0 = B1 @ B1.T
    ev0 = np.linalg.eigvalsh(L0)
    gram = B2.T @ B2
    evg = np.linalg.eigvalsh(gram)
    target = np.sort(np.concatenate([ev0[np.abs(ev0) > 1e-9],
                                     evg[np.abs(evg) > 1e-9]]))
    got = np.sort(ev[np.abs(ev) > 1e-9])
    err = float(np.max(np.abs(got - target)))
    print(f"  D2 nonzero spec(Delta_1) = {{55 graph-clock levels}} u "
          f"{{3 link-Gram levels}}: max err {err:.1e}  "
          f"{'PASS' if err < 1e-8 else 'FAIL'} (Cert; E2 vs du01 "
          f"spectrum route)")
    assert err < 1e-8
    tr = int(round(np.trace(L0)))
    print(f"  D3 trace(B1 B1^T) = {tr} = 2 * 83 non-loop edges — the "
          f"du01 trace-166 fact re-derived by a second route  "
          f"{'PASS' if tr == 166 else 'FAIL'} (EXACT)")
    assert tr == 166
    print(f"  D4 link-Gram eigenvalues (the 3 Eisenstein stiffnesses of "
          f"the edge clock): {np.round(evg[np.abs(evg) > 1e-9], 6)}  "
          f"(Cert, reported)")
    LEDGER["D4 link_gram_levels (Cert)"] = [round(float(v), 8)
                                            for v in evg[np.abs(evg) > 1e-9]]
    return D1, nker


# ----------------------------------------------------------------------
# Stage E — the structural theorem
# ----------------------------------------------------------------------

def stage_E():
    print("=" * 70)
    print("STAGE E — structural theorem for dynamical units")
    print("=" * 70)
    print("  E1 On the shared 26-dim stage: the free graph clock acts "
          "as ZERO (harmonic kernel, EXACT from D1) while the Hecke "
          "clock acts with the full 13-line spectrum x2 (C5, EXACT).")
    print("  E2 THEOREM (free-level obstruction): no ratio of free-"
          "evolution rates chi_H/chi_g can be formed on homology — one "
          "of the two frequencies is identically zero there. The H4 "
          "lifetime non-correspondence (v0.11.0) is thereby explained, "
          "not merely recorded: it was structurally forced.")
    print("  E3 COROLLARY (the du03 program): an internal exchange rate "
          "requires an interaction lifting the harmonic degeneracy — "
          "e.g. the cusp wells transported to edge space, or curvature/"
          "holonomy terms. The SPLITTING PATTERN of the 26 harmonic "
          "modes under such a term, compared line-by-line against the "
          "13 Hecke lines, is the box's dispersion relation — the "
          "clock-to-clock map the Smith-chart reading predicted.")
    LEDGER["E theorem"] = ("free-level exchange rate structurally "
                           "excluded; interaction-induced splitting of "
                           "the 26 harmonic modes = the dispersion "
                           "relation (du03)")


# ----------------------------------------------------------------------
# Stage F — interface export
# ----------------------------------------------------------------------

def stage_F(p1, edge_list, ms, X, K, D, Bd, links, T2, cp_cusp):
    print("=" * 70)
    print("STAGE F — drawn-loop interface data export")
    print("=" * 70)
    free = ms.free
    freepos = {c: k for k, c in enumerate(free)}

    # per-oriented-edge functional on the 29-dim quotient (only valid
    # summed over CLOSED loops; the interface enforces closure)
    edge_rows = []
    for e, (i, j, x) in enumerate(edge_list):
        xs = p1.act(x, SIGMA)
        row = []
        for idx, s in ((p1.index[x], +1), (p1.index[xs], -1)):
            if idx in freepos:
                row.append([freepos[idx], s])
        edge_rows.append(row)

    # layout: Hall (Fiedler-plane) + light force relaxation
    A = np.zeros((56, 56))
    for (i, j, _x) in edge_list:
        if i != j:
            A[i, j] += 1
            A[j, i] += 1
    L0 = np.diag(A.sum(1)) - A
    w, V = np.linalg.eigh(L0)
    pos = V[:, 1:3].copy()
    pos /= np.abs(pos).max()
    rng = np.random.default_rng(143)
    pos += 0.01 * rng.standard_normal(pos.shape)
    for _ in range(300):
        disp = np.zeros_like(pos)
        for a in range(56):
            d = pos[a] - pos
            r2 = (d ** 2).sum(1) + 1e-4
            disp[a] += (d / r2[:, None]).sum(0) * 0.0015
        for (i, j, _x) in edge_list:
            if i == j:
                continue
            d = pos[j] - pos[i]
            disp[i] += 0.06 * d
            disp[j] -= 0.06 * d
        pos += np.clip(disp, -0.03, 0.03)
        pos /= np.abs(pos).max()

    # spectral signature basis: left eigen-structure of exact T2 (float),
    # orthonormalized WITHIN each eigenvalue group (QR) so degenerate
    # blocks (Eis x3, old x4, each f-line x2) stay well-conditioned.
    T2f = np.array([[float(T2[i, j]) for j in range(29)]
                    for i in range(29)])
    wl, Wl = np.linalg.eig(T2f.T)
    wl = wl.real
    groups: dict = {}
    for k in range(29):
        key = round(float(wl[k]), 6)
        groups.setdefault(key, []).append(k)
    a2_f2 = [-1.126757, -0.197126, 1.747468, 2.576415]
    cols, labels, avals = [], [], []
    for lam in sorted(groups):
        idx = groups[lam]
        block = Wl[:, idx].real
        Q, _ = np.linalg.qr(block)
        if abs(lam - 3.0) < 1e-6:
            lab = "Eis"
        elif abs(lam + 2.0) < 1e-6:
            lab = "old"
        elif abs(lam) < 1e-6:
            lab = "f1"
        elif min(abs(lam - a) for a in a2_f2) < 1e-4:
            lab = "f2"
        else:
            lab = "f3"
        for c in range(Q.shape[1]):
            cols.append(Q[:, c])
            labels.append(lab)
            avals.append(round(float(lam), 8))
    Wsig = np.array(cols).T                     # 29 x 29
    Winv = np.linalg.inv(Wsig)
    cond = float(np.linalg.cond(Wsig))
    resid = float(np.max(np.abs(Wsig @ Winv - np.eye(29))))
    counts = {l: labels.count(l) for l in ("Eis", "old", "f1", "f2", "f3")}
    print(f"  F1 signature basis (left T2 eigenvectors, groupwise QR; "
          f"DIAGNOSTIC Euclidean normalization): counts {counts} "
          f"[expect Eis:3 old:4 f1:2 f2:8 f3:12]; cond = {cond:.2f}, "
          f"inv residual {resid:.1e}  "
          f"{'PASS' if cond < 1e3 and resid < 1e-10 else 'FAIL'} (Cert)")
    ok = counts == {"Eis": 3, "old": 4, "f1": 2, "f2": 8, "f3": 12}
    assert ok and cond < 1e3 and resid < 1e-10
    LEDGER["F1 signature_counts (Cert)"] = counts
    LEDGER["F1 basis_cond (Cert)"] = round(cond, 3)

    data = {
        "meta": {"level": 143, "genus": 13, "triangles": 56, "edges": 84,
                 "b1": 29, "H1": 26, "study": "du02 v0.11.3",
                 "classes": {"pairing": "EXACT", "hecke": "EXACT",
                             "signature_weights": "DIAGNOSTIC",
                             "loop_x_loop": "session 3 (Goldman)"}},
        "positions": [[round(float(x), 5) for x in p] for p in pos],
        "edges": [{"i": i, "j": j, "sym": p1.index[x]}
                  for (i, j, x) in edge_list],
        "edge_rows": edge_rows,
        "boundary": Bd,
        "links": {str(d): links[d] for d in ms.divisors},
        "T2": [[round(float(T2[i, j]), 10) for j in range(29)]
               for i in range(29)],
        "Wsig": [[round(float(Wsig[i, k]), 8) for k in range(29)]
                 for i in range(29)],
        "Winv": [[round(float(Winv[i, k]), 10) for k in range(29)]
                 for i in range(29)],
        "labels": labels,
        "avals": avals,
        "cuspidal_charpoly": str(cp_cusp),
    }
    out = os.path.join(_HERE, "du02_interface_data.json")
    with open(out, "w") as f:
        json.dump(data, f)
    print(f"  F2 interface data written: {out} "
          f"({os.path.getsize(out)//1024} KB)")
    return out


# ----------------------------------------------------------------------

def main():
    p1, tris, edges, A, edge_list, tri_of, ms = build_all()
    X, K, tree, F, D = stage_A(p1, edge_list, ms)
    Bd, links, Lmat = stage_B(p1, edge_list, ms, K, D)
    T2, T2star = stage_C(ms, D, Lmat)
    cp_cusp, cp_e = stage_C2(ms, D, Bd, T2, T2star)
    stage_D(edge_list, links, ms)
    stage_E()
    stage_F(p1, edge_list, ms, X, K, D, Bd, links, T2, cp_cusp)
    out = os.path.join(_HERE, "du02_cycle_space_map.json")
    with open(out, "w") as f:
        json.dump(LEDGER, f, indent=1)
    print("=" * 70)
    print(f"ledger written: {out}")
    print("du02 complete — the cycle-space map exists, is perfect "
          "(det +-1), carries Hecke integrally with the same 13 lines, "
          "and the free-level exchange rate is now a THEOREM-level "
          "exclusion. Interface data staged for the React loop.")


if __name__ == "__main__":
    main()
