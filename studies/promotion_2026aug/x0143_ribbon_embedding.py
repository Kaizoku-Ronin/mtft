#!/usr/bin/env python3
"""
x0143_ribbon_embedding.py — the surface embedding of the house skeleton
============================================================================

MIT License — Copyright (c) 2026 Roger Tano

Sequel to x0143_graph_uncertainty.py.  The R-cyclic orderings make the
dual skeleton a RIBBON GRAPH (darts = P1(Z/143), vertex rotation R,
edge involution S), and the ribbon structure embeds it in a closed
oriented surface — which is X0(143) itself.  This study certifies the
embedding, works both sides of the face/vertex duality, derives the
homological origin of the number 26, and extracts a certified
cusp-type quotient spectrum.

RESULTS (classes as tagged; pre-registrations noted)
----------------------------------------------------
R1 (Pr/EXACT, pre-registered -> confirmed).  The face permutation of
  the ribbon graph (apply S, then R) equals T pointwise on all 168
  darts — the identity S·R = S·(S·T) = T in PSL2(Z), realized.  Faces
  = cusp fans, sizes {1, 11, 13, 143}.

R2 (EXACT).  Euler characteristic V - E + F = 56 - 84 + 4 = -24 =
  2 - 2·13: the ribbon structure knows the genus of X0(143).

R3 (EXACT, pre-registered -> confirmed).  Face/vertex duality: the
  primal graph is a 4-vertex multigraph on the cusps with degrees =
  widths (1, 11, 13, 143).  Exact multiplicity structure:
     m(oo, 0) = 1        (the unique edge at infinity — the SAME
                          S-orbit {(0:1),(1:0)} that is the dual
                          self-loop: one object, both defects)
     m(w11, w13) = 1     (the two middle cusps meet exactly once)
     m(w11, 0) = 10,  m(w13, 0) = 12,  loops at cusp 0: 60,
     loops elsewhere: 0.        Total 1+1+10+12+60 = 84 edges.

R4 (Pr/EXACT + Cert).  The width-gradient commutator on the primal
  graph: [D,L]_{ij} = m_ij (h_j - h_i), entries {142, 2, 1320, 1560}
  up to sign, two disjoint routes agreeing.  Primal Robertson capacity
  (exact algebraic number via the 4x4 Pfaffian):
     capacity* = (1/2) sqrt( (P + sqrt(P^2 - 4 Pf^2)) / 2 ),
     P = 4 196 168,  Pf = 284,   capacity* = 1024.22751144...
  against dual capacity exactly 1.  The SAME theorem on the SAME
  surface: the dual side concentrates uncertainty into one unit at one
  edge; the primal side exposes three orders of magnitude of it across
  the width gradients.

R5 (EXACT — the headline).  HOMOLOGICAL ORIGIN OF 26.  The boundary
  map d2: Z^4(faces) -> Z^84(edges) of the CW structure has
  SNF = (1, 1, 1): rank 3, torsion-free, sum of all face boundaries
  = 0 the single relation.  Hence
     H1(X0(143)) = Z^(29-3) = Z^26 = Z^(2·13),
     b1(skeleton) = 29 = 2g + (#cusps - 1) = 26 + 3.
  The three directions of the skeleton's cycle space killed on the
  surface are exactly the three independent cusp-fan boundaries.  The
  du03 stage dimension 26 thus has a ribbon-theoretic derivation on
  the skeleton side; identification with du03's dynamical mechanism
  (cusp wells annihilating the harmonic stage) is deferred — the
  bookkeeping identity is certified, the mechanism match is not
  claimed (AG-D5).

R6 (EXACT).  FAN GEOMETRY.  The width-11 and width-13 fans are SIMPLE
  cycles in the dual graph (11 and 13 distinct triangles in 11 and 13
  steps — embedded circles on the surface).  The width-143 fan is a
  closed walk visiting ALL 56 triangles; indeed every triangle has at
  least one corner at cusp 0.  Triangle type census (corner widths):
  (1,143,143) x1 [base], (143,143,143) x33, (13,143,143) x11,
  (11,143,143) x9, (11,13,143) x2 — corner budgets 1/11/13/143 all
  exactly consumed.

R7 (EXACT; pre-registration FALSIFIED, preserved).  Pre-registered
  guess: corner-type seeding refines to the discrete partition (the
  skeleton is WL-asymmetric).  FALSE.  Color refinement seeded by the
  corner types stabilizes at 30 equitable classes on 56 vertices: the
  skeleton carries strictly more regularity than its cusp types force.
  The 30x30 Laplacian quotient has an exact integer characteristic
  polynomial that DIVIDES the exact characteristic polynomial of the
  full 56x56 Laplacian over Z (division certificate, remainder 0): a
  certified 30-harmonic cusp-type spectrum inside the skeleton's 56.

R8 (EXACT).  THIRD ROUTE TO THE TREE COUNT.  The exact charpoly of L
  (Bareiss interpolation at 57 integer points, monic, integral) gives
  kappa = -[x^1]/56 = 3 518 081 582 959 364 640 — agreeing with the
  Bareiss-minor route and the archived du-session value.  Three
  routes, no shared steps.

R9 (EXACT census; AG-D5 filing).  Primal spanning trees kappa* = 142
  = m(w11,0)·m(w13,0) + m(w11,0)·m(w11,w13) + m(w13,0)·m(w11,w13)
  = 120 + 10 + 12 (the unique oo-edge is forced).  FILED, DISMISSED:
  142 also appears as the width gradient h(0) - h(oo) = 143 - 1 on
  the oo-edge; proximity of two different 142s, no mechanism, not
  evidence.  gcd(kappa, kappa*) = 2; 71 does not divide kappa.

R10 (Pr/EXACT — mechanism identified, proximity dismissed).  The WL
  profile 4 + 26x2 demanded a mechanism before any filing against
  2g = 26.  Candidate: complex conjugation z -> -conj(z), which
  preserves the Farey tessellation and normalizes Gamma0(143), acting
  on flags by iota(c:d) = (-c:d), with iota^2 = 1, iota S = S iota,
  iota T iota = T^{-1}.  Right-action conjugation sends the R-orbit
  partition to its T^{-1}-translate, so the triangle involution is
  psi = iota followed by (.)·T.  (First attempt psi' = iota·T^{-1}
  FAILED to permute triangles; the failure diagnosed the handedness of
  right-action conjugation and is preserved as the discovery route.)
  Certified: psi^2 = 1, psi permutes R-orbits, the induced tau is an
  automorphism of the dual graph fixing the loop vertex, with EXACTLY
  4 fixed triangles and 26 two-cycles, and the WL partition equals the
  tau-orbit partition class-for-class (singletons = fixed set, pairs =
  swaps, both directions).  Fixed triangles: the base (1,143,143);
  {(1:2),(1:72),(1:142)} of type (143,143,143); one (13,143,143); one
  (11,143,143) — the (11,13,143) type contains no fixed triangle (its
  two members are a tau-pair).  No unoriented edge is psi-fixed.
  AG-D5 RESOLUTION: the pair count 26 = (56 - 4)/2 is real-involution
  arithmetic, NOT the homological 26 = 2g of R5.  Two different 26s;
  the proximity is dismissed, and the mechanism hunt is what earned
  the dismissal.

GATES: G1-G10 mirror R1-R10; all exact integer/permutation work except
the capacity SVD cross-check and quotient-eigenvalue reporting (Cert
1e-9).  Run: python studies/x0143_ribbon_embedding.py  (~60 s).
Writes x0143_ribbon_embedding_ledger.json next to itself.
"""

from __future__ import annotations

import json
import os
import sys
from collections import Counter
from fractions import Fraction as Fr
from math import gcd

import numpy as np

from mtft import combinatorial as C

N = 143
ARCHIVED_TREES = 3518081582959364640


def build_flags():
    units = [u for u in range(1, N) if gcd(u, N) == 1]
    valid = lambda c, d: not (c % 11 == 0 and d % 11 == 0) and \
                         not (c % 13 == 0 and d % 13 == 0)
    canon = lambda c, d: min(((c * u) % N, (d * u) % N) for u in units)
    P1 = sorted({canon(c, d) for c in range(N) for d in range(N)
                 if valid(c, d)})
    idx = {p: i for i, p in enumerate(P1)}
    perm = lambda f: [idx[canon(*f(*p))] for p in P1]
    sS = perm(lambda c, d: (d % N, (-c) % N))
    sT = perm(lambda c, d: (c, (c + d) % N))
    sR = perm(lambda c, d: (d % N, (d - c) % N))
    return P1, sS, sT, sR


def orbits(s):
    seen, out = set(), []
    for i in range(len(s)):
        if i in seen:
            continue
        o, j = [], i
        while j not in seen:
            seen.add(j)
            o.append(j)
            j = s[j]
        out.append(o)
    return out


def bareiss_det(M):
    M = [r[:] for r in M]
    n, prev, sg = len(M), 1, 1
    for k in range(n - 1):
        if M[k][k] == 0:
            for r in range(k + 1, n):
                if M[r][k]:
                    M[k], M[r] = M[r], M[k]
                    sg = -sg
                    break
            else:
                return 0
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                M[i][j] = (M[i][j] * M[k][k] - M[i][k] * M[k][j]) // prev
        prev = M[k][k]
    return sg * M[n - 1][n - 1]


def exact_charpoly(Mint):
    """Monic integer charpoly det(xI - M) by Bareiss interpolation."""
    n = len(Mint)
    pts = list(range(n + 1))
    vals = [bareiss_det([[(x if i == j else 0) - Mint[i][j]
                          for j in range(n)] for i in range(n)])
            for x in pts]
    work = [Fr(v) for v in vals]
    newton = [work[0]]
    for k in range(1, n + 1):
        work = [(work[i + 1] - work[i]) / Fr(k) for i in range(len(work) - 1)]
        newton.append(work[0])
    acc = [Fr(0)] * (n + 1)
    base = [Fr(1)] + [Fr(0)] * n
    for k in range(n + 1):
        for i in range(n + 1):
            acc[i] += newton[k] * base[i]
        nb = [Fr(0)] * (n + 1)
        for i in range(n):
            nb[i + 1] += base[i]
            nb[i] -= Fr(k) * base[i]
        base = nb
    assert all(c.denominator == 1 for c in acc)
    cp = [int(c) for c in acc]
    assert cp[n] == 1
    return cp


def poly_divmod_z(P, Q):
    """Exact division of integer polynomials (both monic)."""
    P = P[:]
    dq, out = len(Q) - 1, [0] * (len(P) - len(Q) + 1)
    for i in range(len(P) - 1, dq - 1, -1):
        c = P[i]
        if c:
            out[i - dq] = c
            for j, q in enumerate(Q):
                P[i - dq + j] -= c * q
    return out, P


def snf(Mrows):
    M = [r[:] for r in Mrows]
    Rn, Cn = len(M), len(M[0])
    inv, t = [], 0
    while t < min(Rn, Cn):
        pr = [(abs(M[i][j]), i, j) for i in range(t, Rn)
              for j in range(t, Cn) if M[i][j]]
        if not pr:
            break
        _, pi, pj = min(pr)
        M[t], M[pi] = M[pi], M[t]
        for r in M:
            r[t], r[pj] = r[pj], r[t]
        again = True
        while again:
            again = False
            for i in range(t + 1, Rn):
                if M[i][t]:
                    q = M[i][t] // M[t][t]
                    M[i] = [a - q * b for a, b in zip(M[i], M[t])]
                    if M[i][t]:
                        M[t], M[i] = M[i], M[t]
                        again = True
            for j in range(t + 1, Cn):
                if M[t][j]:
                    q = M[t][j] // M[t][t]
                    for r in range(Rn):
                        M[r][j] -= q * M[r][t]
                    if M[t][j]:
                        for r in range(Rn):
                            M[r][t], M[r][j] = M[r][j], M[r][t]
                        again = True
        inv.append(abs(M[t][t]))
        t += 1
    return inv


def main() -> int:
    ledger = {"study": "x0143_ribbon_embedding", "gates": {}}
    ok = True

    def gate(name, passed, **info):
        nonlocal ok
        ok &= bool(passed)
        ledger["gates"][name] = {"passed": bool(passed), **info}
        print(f"[{'PASS' if passed else 'FAIL'}] {name}  " +
              "  ".join(f"{k}={v}" for k, v in info.items()))

    P1, sS, sT, sR = build_flags()
    _units = [u for u in range(1, N) if gcd(u, N) == 1]
    canonP = lambda p: min(((p[0] * u) % N, (p[1] * u) % N) for u in _units)
    idxP = {p: i for i, p in enumerate(P1)}

    # G1 face permutation and fans
    face = [sR[sS[i]] for i in range(168)]
    fans = sorted(orbits(sT), key=len)
    gate("G1_face_perm_is_T", face == sT
         and [len(f) for f in fans] == [1, 11, 13, 143])
    cusp_of, width_of = {}, {}
    for k, o in enumerate(fans):
        width_of[k] = len(o)
        for f in o:
            cusp_of[f] = k

    tris = orbits(sR)
    tri_of = {}
    for ti, t in enumerate(tris):
        for f in t:
            tri_of[f] = ti

    # G2 genus
    gate("G2_genus", len(tris) - len(orbits(sS)) + len(fans) == 2 - 2 * 13,
         chi=len(tris) - 84 + 4)

    # G3 primal cusp multigraph
    W = np.zeros((4, 4), dtype=np.int64)
    loops = [0] * 4
    inf_edges = []
    for o in orbits(sS):
        a, b = cusp_of[o[0]], cusp_of[sS[o[0]]]
        if a == b:
            loops[a] += 1
            W[a, a] += 2
        else:
            W[a, b] += 1
            W[b, a] += 1
        if 0 in (a, b):
            inf_edges.append(sorted(P1[f] for f in o))
    gate("G3_primal_duality",
         W.sum(1).tolist() == [1, 11, 13, 143]
         and int(W[0, 3]) == 1 and int(W[1, 2]) == 1
         and int(W[1, 3]) == 10 and int(W[2, 3]) == 12
         and loops == [0, 0, 0, 60]
         and inf_edges == [[(0, 1), (1, 0)]],
         degrees=W.sum(1).tolist(), loops=loops,
         infinity_edge=str(inf_edges))

    # G4 width-gradient commutator + capacity
    comm = C.commutator_DL(W)
    routes = np.array_equal(comm, C.degree_gradient_matrix(W))
    P = sum(int(comm[i, j]) ** 2 for i in range(4) for j in range(i + 1, 4))
    Pf = int(comm[0, 1] * comm[2, 3] - comm[0, 2] * comm[1, 3]
             + comm[0, 3] * comm[1, 2])
    mu2 = (P + (P * P - 4 * Pf * Pf) ** 0.5) / 2
    cap_closed = 0.5 * mu2 ** 0.5
    cap_svd = 0.5 * max(np.linalg.svd(comm.astype(float),
                                      compute_uv=False))
    ent = sorted(abs(int(comm[i, j])) for i in range(4)
                 for j in range(i + 1, 4) if comm[i, j] != 0)
    # First run of this gate compared ALL six upper-triangle entries
    # (including the two zero pairs oo-w11, oo-w13) against the four
    # nonzero magnitudes and failed; the miss was in the gate, not the
    # mathematics.  Preserved per append-only correction discipline.
    gate("G4_primal_commutator", routes and ent == [2, 142, 1320, 1560]
         and P == 4196168 and Pf == 284
         and abs(cap_closed - cap_svd) < 1e-9,
         entries=ent, P=P, Pf=Pf, capacity=f"{cap_closed:.9f}")

    # G5 homology split 29 = 26 + 3
    E_orb = orbits(sS)
    eid, emin = {}, {}
    for k, o in enumerate(E_orb):
        for f in o:
            eid[f] = k
        emin[k] = min(o)
    B = [[0] * 4 for _ in range(len(E_orb))]
    for k, fc in enumerate(fans):
        for f in fc:
            B[eid[f]][k] += (1 if f == emin[eid[f]] else -1)
    rel = all(sum(r) == 0 for r in B)
    invf = snf(B)
    gate("G5_homology_26", rel and invf == [1, 1, 1]
         and 29 - len(invf) == 26 and 26 == 2 * 13,
         snf=invf, H1_rank=29 - len(invf))

    # G6 fan geometry and type census
    visits = [len({tri_of[f] for f in fc}) for fc in fans]
    types = Counter(tuple(sorted(width_of[cusp_of[f]] for f in t))
                    for t in tris)
    census = {str(k): v for k, v in sorted(types.items())}
    gate("G6_fans", visits == [1, 11, 13, 56]
         and types[(1, 143, 143)] == 1 and types[(143, 143, 143)] == 33
         and types[(13, 143, 143)] == 11 and types[(11, 143, 143)] == 9
         and types[(11, 13, 143)] == 2,
         fan_visits=visits, census=census)

    # G7 WL refinement -> 30-class equitable quotient, division certificate
    Ad = np.zeros((56, 56), dtype=np.int64)
    dual_loop_v = None
    for o in E_orb:
        a, b = tri_of[o[0]], tri_of[sS[o[0]]]
        if a != b:
            Ad[a, b] += 1
            Ad[b, a] += 1
        else:
            dual_loop_v = a
    col = {ti: tuple(sorted(width_of[cusp_of[f]] for f in tris[ti]))
           for ti in range(56)}
    while True:
        sig = {ti: (col[ti], tuple(sorted(col[int(u)]
               for u in np.nonzero(Ad[ti])[0]
               for _ in range(int(Ad[ti, u])))))
               for ti in range(56)}
        lab = {}
        for ti in range(56):
            if sig[ti] not in lab:
                lab[sig[ti]] = len(lab)
        nc = {ti: lab[sig[ti]] for ti in range(56)}
        stable = len(set(nc.values())) == len(set(col.values()))
        col = nc
        if stable:
            break
    kcl = len(set(col.values()))
    classes = [[ti for ti in range(56) if col[ti] == c] for c in range(kcl)]
    # equitability assert + quotient
    BQ = [[0] * kcl for _ in range(kcl)]
    equit = True
    for c, mem in enumerate(classes):
        rows = [[int(Ad[v, u]) for u in range(56)] for v in mem]
        for c2, mem2 in enumerate(classes):
            counts = {sum(r[u] for u in mem2) for r in rows}
            equit &= len(counts) == 1
            BQ[c][c2] = counts.pop()
    degs = Ad.sum(1)
    dclass = [int(degs[mem[0]]) for mem in classes]
    LQ = [[(dclass[i] if i == j else 0) - BQ[i][j]
           for j in range(kcl)] for i in range(kcl)]
    Ld = np.diag(degs) - Ad
    cpL = exact_charpoly([[int(Ld[i, j]) for j in range(56)]
                          for i in range(56)])
    cpQ = exact_charpoly(LQ)
    quo, rem = poly_divmod_z(cpL, cpQ)
    gate("G7_WL_quotient_spectrum", kcl == 30 and equit
         and all(r == 0 for r in rem),
         classes=kcl, division_remainder_zero=all(r == 0 for r in rem),
         class_sizes=sorted(len(m) for m in classes))
    ledger["quotient_charpoly_coeffs"] = cpQ

    # G8 third route to kappa
    kap_cp = -cpL[1] // 56
    L55 = [[int(Ld[i, j]) for j in range(1, 56)] for i in range(1, 56)]
    kap_bar = bareiss_det(L55)
    gate("G8_three_route_kappa", kap_cp == kap_bar == ARCHIVED_TREES,
         kappa=kap_cp)
    ledger["charpoly_L_coeffs"] = cpL

    # G9 primal trees + AG-D5 filing
    Wp = W.copy()
    np.fill_diagonal(Wp, 0)
    Lp = np.diag(Wp.sum(1)) - Wp
    ks = bareiss_det([[int(Lp[i, j]) for j in range(1, 4)]
                      for i in range(1, 4)])
    gate("G9_primal_trees", ks == 142 and ks == 120 + 10 + 12
         and gcd(ARCHIVED_TREES, ks) == 2 and ARCHIVED_TREES % 71 != 0,
         kappa_star=ks,
         ag_d5="142 = kappa* and 142 = h(0)-h(oo): proximity filed, "
               "no mechanism, dismissed")

    # G10 the real involution explains the WL profile (mechanism, not proximity)
    iota = {i: idxP[canonP(((-P1[i][0]) % N, P1[i][1]))] for i in range(168)}
    alg = (all(iota[iota[i]] == i for i in range(168))
           and all(iota[sS[i]] == sS[iota[i]] for i in range(168)))
    sTinv = [0] * 168
    for i, j in enumerate(sT):
        sTinv[j] = i
    alg &= all(iota[sT[iota[i]]] == sTinv[i] for i in range(168))
    psi = [sT[iota[i]] for i in range(168)]        # iota then (.)*T
    tri_ok = (all(psi[psi[i]] == i for i in range(168))
              and all(tri_of[psi[t[0]]] == tri_of[psi[f]]
                      for t in tris for f in t))
    tau = {ti: tri_of[psi[tris[ti][0]]] for ti in range(56)}
    autm = (all(tau[tau[a]] == a for a in tau)
            and all(int(Ad[tau[a], tau[b]]) == int(Ad[a, b])
                    for a in range(56) for b in range(56)))
    fixed = sorted(a for a in range(56) if tau[a] == a)
    tpairs = {tuple(sorted((a, tau[a]))) for a in range(56) if tau[a] != a}
    singles = sorted(m[0] for m in classes if len(m) == 1)
    wpairs = {tuple(sorted(m)) for m in classes if len(m) == 2}
    efix = sum(1 for o in E_orb
               if {psi[f] for f in o} == {o[0], sS[o[0]]})
    gate("G10_real_involution", alg and tri_ok and autm
         and len(fixed) == 4 and len(tpairs) == 26
         and singles == fixed and wpairs == tpairs
         and tau[dual_loop_v] == dual_loop_v and efix == 0,
         fixed_triangles=fixed, two_cycles=len(tpairs),
         WL_equals_tau_orbits=(singles == fixed and wpairs == tpairs),
         fixed_types=[str(tuple(sorted(width_of[cusp_of[f]]
                                       for f in tris[a]))) for a in fixed])

    ledger["all_passed"] = ok
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "x0143_ribbon_embedding_ledger.json")
    with open(out, "w") as f:
        json.dump(ledger, f, indent=2, default=str)
    print(f"\nledger -> {out}")
    print("ALL GATES PASS" if ok else "GATE FAILURE")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
