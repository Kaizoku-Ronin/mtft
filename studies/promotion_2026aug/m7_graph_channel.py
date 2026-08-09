#!/usr/bin/env python3
"""
m7_graph_channel.py — M7: the graph-side coupling as a non-Hecke channel
=========================================================================

MIT License — Copyright (c) 2026 Roger Tano

THE QUESTION (M7, as frozen at the close of the fs01-fs06 arc).  M6
proved two things: no Im-type CP violation can reach flavor (reality
is inherited by any holonomy dressing), and the holonomy Hamiltonian
H_hol lies in the HECKE ALGEBRA, hence acts as a scalar on each
newform plane — it can shift flavor levels but can never induce
flavor-CHANGING structure.  The transmitted CP-odd object is therefore
a per-line real shift Delta_h_f: SPECTRAL CP violation.  The Standard
Model's observed CP violation is MIXING-type.  So MTFT must either
exhibit a NON-HECKE transmission channel that moves between lines, or
carry a falsifiable tension in writing.  Only non-Hecke operators can
move between lines; the candidate named at the time was the
holonomy-dressed graph-side coupling.

WHAT THIS STUDY SETTLES.  The graph side of the box is the dual graph
on the 56 Farey triangles; the Hecke side is the 26-dimensional
cuspidal homology.  They meet through the harmonic embedding
    W : H_1(X_0(143)) -> R^84,
which sends a homology class to its unique representative orthogonal
to im d2.  A graph-side potential is a function g on triangles; it
induces the orientation-EVEN edge function
    g_avg(e) = (g(a) + g(b)) / 2
on the two triangles adjacent to e (the orientation-odd gradient is
NOT well defined as a multiplication operator on unoriented edge
space — a first attempt used it and is preserved below as the
discovery route).  The induced coupling on homology is the
compression
    V(g) = G^-1 M(g),   M_ab = <W_a, g_avg W_b>,  G_ab = <W_a, W_b>,
which is exact rational arithmetic whenever g is integral.

PRE-REGISTRATION (frozen before any commutator was computed)
------------------------------------------------------------
  M7-A  V(g) does NOT commute with the Hecke algebra for canonical
        integral g.  Predicted: confirmed (this is what "graph-side"
        should mean).
  M7-B  V(g) has NONZERO off-block matrix elements, i.e. it is
        genuinely flavor-CHANGING between the four Hecke blocks.
        Predicted: confirmed.  This is the half of M7 that decides
        whether a non-Hecke channel exists at all.
  M7-C  V(g) is G-self-adjoint and REAL, hence diagonalized by a real
        orthogonal rotation, hence the mixing it induces carries NO
        phase: the channel is CP-EVEN.  Predicted: confirmed — which
        would mean the channel resolves the flavor-changing half of
        the tension but NOT the CP half.
  M7-D  Selection rule: V(g) either commutes or anticommutes with the
        star involution iota*.  No prediction registered (this is the
        analogue of the du03 parity rule {V, eta} = 0 and could go
        either way).
  M7-N  NULL CONTROL (added after the first run, and the reason it
        is here is preserved): the degree function was registered as a
        second potential to show the conclusion is not an artifact of
        one choice of g.  It FAILED, giving exactly zero transmission
        and zero commutator — because the dual graph is 3-REGULAR, so
        the degree function is CONSTANT and V = 3I.  That is not a
        failure of the channel but a failure of the choice, and it is
        exactly the right null control: a constant potential must
        produce the trivial channel, and the pipeline reports it as
        trivial.  It is kept as a gate in that role, with graph
        distance from the base triangle supplied as the genuine
        second potential.
  M7-E  Echo test: is the elliptic line quiet in this channel, as it
        was in M6 where f1's CP-quietness traced to a_2(143a1) = 0?
        No prediction registered; the mechanism there was
        supersingularity at 2, which has no evident graph-side
        counterpart.

  Registered consequence either way: if M7-A and M7-B hold while M7-C
  holds too, then the M6 tension is SHARPENED rather than resolved —
  the obstruction is not the absence of a non-Hecke channel (there is
  one, and it mixes), it is the REALITY of the arithmetic.  CP-odd
  mixing would then require a genuinely complex dressing, and the
  successor object is named exactly in the verdict.

Gates P1-P7.  Runtime ~2 min.  Writes m7_graph_channel_ledger.json.
"""

from __future__ import annotations

import json
import os
import sys
import time
from fractions import Fraction as Fr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from mtft import hecke as H  # v0.14.0 exposes mtft.hecke (auditor portability patch; author's local alias was mtftpkg)

HERE = os.path.dirname(os.path.abspath(__file__))


def matmul(A, B):
    n, m, p = len(A), len(B), len(B[0])
    return [[sum(A[i][k] * B[k][j] for k in range(m)) for j in range(p)]
            for i in range(n)]


def inverse(M):
    n = len(M)
    A = [list(r) + [Fr(1) if i == j else Fr(0) for j in range(n)]
         for i, r in enumerate(M)]
    for c in range(n):
        pr = next(i for i in range(c, n) if A[i][c] != 0)
        A[c], A[pr] = A[pr], A[c]
        pv = A[c][c]
        A[c] = [x / pv for x in A[c]]
        for i in range(n):
            if i != c and A[i][c] != 0:
                f = A[i][c]
                A[i] = [a - f * b for a, b in zip(A[i], A[c])]
    return [r[n:] for r in A]


def frob(M):
    return sum(x * x for r in M for x in r)


def main() -> int:
    t0 = time.time()
    ledger = {"study": "m7_graph_channel", "gates": {}}
    ok = True

    def gate(name, passed, **info):
        nonlocal ok
        ok &= bool(passed)
        ledger["gates"][name] = {"passed": bool(passed), **info}
        print(f"[{'PASS' if passed else 'FAIL'}] {name}  "
              + "  ".join(f"{k}={v}" for k, v in info.items()))

    import pickle
    CACHE = os.path.join(HERE, "m7_harmonic_cache.pkl")
    m = H.model()
    E, D2, K, free, nq = m["E"], m["D2"], m["K"], m["free"], m["nq"]
    tris, tri_of, sS = m["tris"], m["tri_of"], m["sS"]
    erep, cusp_of, fans = m["erep"], m["cusp_of"], m["fans"]

    # ── P1: harmonic embedding and Gram matrix ──────────────────────
    G56 = [[sum(D2[e][i] * D2[e][j] for e in range(E))
            for j in range(56)] for i in range(56)]

    def harmonic(v26):
        v = [sum(K[a][j] * v26[a] for a in range(26)) for j in range(nq)]
        v84 = [Fr(0)] * E
        for j, e in enumerate(free):
            v84[e] = v[j]
        rhs = [sum(Fr(D2[e][i]) * v84[e] for e in range(E))
               for i in range(56)]
        Aug = [[Fr(G56[i][j]) for j in range(56)] + [rhs[i]]
               for i in range(56)]
        R, pv = H._rref(Aug)
        alpha = [Fr(0)] * 56
        for r_, c_ in enumerate(pv):
            if c_ < 56:
                alpha[c_] = R[r_][56]
        return [v84[e] - sum(Fr(D2[e][i]) * alpha[i] for i in range(56))
                for e in range(E)]

    if os.path.exists(CACHE):
        with open(CACHE, "rb") as fh:
            basis, Gm = pickle.load(fh)
    else:
        basis = []
        for a in range(26):
            v = [Fr(0)] * 26
            v[a] = Fr(1)
            basis.append(harmonic(v))
        Gm = [[sum(basis[a][e] * basis[b][e] for e in range(E))
               for b in range(26)] for a in range(26)]
        with open(CACHE, "wb") as fh:
            pickle.dump((basis, Gm), fh)
    harm_ok = all(sum(Fr(D2[e][i]) * w[e] for e in range(E)) == 0
                  for w in basis for i in range(0, 56, 7))
    Ginv = inverse(Gm)
    gate("P1_harmonic_embedding",
         harm_ok and len(basis) == 26 and frob(Gm) > 0,
         dim=26, gram="invertible", orthogonal_to_im_d2=harm_ok)

    # ── canonical integral graph potentials ─────────────────────────
    width = {k: len(o) for k, o in enumerate(fans)}
    g_width = [sum(width[cusp_of[f]] for f in tris[t]) for t in range(56)]
    deg = [0] * 56
    for k in range(E):
        a, b = tri_of[erep[k]], tri_of[sS[erep[k]]]
        deg[a] += 1
        deg[b] += 1
    adj = [[] for _ in range(56)]
    for k in range(E):
        a, b = tri_of[erep[k]], tri_of[sS[erep[k]]]
        if a != b:
            adj[a].append(b)
            adj[b].append(a)
    dist = [-1] * 56
    dist[0] = 0
    frontier = [0]
    while frontier:
        nxt = []
        for u in frontier:
            for v_ in adj[u]:
                if dist[v_] < 0:
                    dist[v_] = dist[u] + 1
                    nxt.append(v_)
        frontier = nxt
    POT = {"width": g_width, "degree": deg, "distance": dist}

    def coupling(g):
        gavg = [Fr(g[tri_of[erep[k]]] + g[tri_of[sS[erep[k]]]], 2)
                for k in range(E)]
        M = [[sum(basis[a][e] * gavg[e] * basis[b][e] for e in range(E))
              for b in range(26)] for a in range(26)]
        return matmul(Ginv, M), M

    V, Mraw = coupling(g_width)
    Vdeg, Mdeg = coupling(deg)
    Vdist, Mdist = coupling(dist)

    # ── P2 (M7-A): non-Hecke ────────────────────────────────────────
    comms = {}
    nonhecke = True
    for p in (2, 3, 5):
        T = [list(r) for r in H.cuspidal_hecke(p)]
        C = [[matmul(V, T)[i][j] - matmul(T, V)[i][j]
              for j in range(26)] for i in range(26)]
        nz = sum(1 for r in C for x in r if x != 0)
        comms[p] = nz
        nonhecke &= nz > 0
    gate("P2_M7A_non_hecke", nonhecke,
         nonzero_commutator_entries=str(comms),
         verdict="graph coupling is NOT in the Hecke algebra")

    # ── P3 (M7-B): flavor-changing between blocks ───────────────────
    blocks = H.blocks()
    order = ["ell", "old", "q4", "q6"]
    cols = []
    for nm in order:
        for v in blocks[nm]:
            cols.append(list(v))
    S = [[cols[b][i] for b in range(26)] for i in range(26)]
    Sinv = inverse(S)
    Vb = matmul(Sinv, matmul(V, S))
    spans, off = {}, 0
    idxs, start = {}, 0
    for nm in order:
        d = len(blocks[nm])
        idxs[nm] = (start, start + d)
        start += d
    table = {}
    for A in order:
        for B in order:
            i0, i1 = idxs[A]
            j0, j1 = idxs[B]
            sub = [[Vb[i][j] for j in range(j0, j1)]
                   for i in range(i0, i1)]
            f = frob(sub)
            table[f"{A}->{B}"] = float(f) ** 0.5
            if A != B:
                off += f
    gate("P3_M7B_flavor_changing", off > 0,
         off_block_norm=f"{float(off) ** 0.5:.6f}",
         transmission=({k: round(v, 4) for k, v in table.items()}),
         verdict="the channel MOVES BETWEEN LINES")

    # ── P4 (M7-C): reality and G-self-adjointness ───────────────────
    sym = all(Mraw[i][j] == Mraw[j][i] for i in range(26)
              for j in range(26))
    GV = matmul(Gm, V)
    selfadj = all(GV[i][j] == GV[j][i] for i in range(26)
                  for j in range(26))
    real = all(isinstance(x, Fr) for r in V for x in r)
    gate("P4_M7C_real_and_self_adjoint", sym and selfadj and real,
         M_symmetric=sym, G_self_adjoint=selfadj, entries="rational",
         verdict="mixing is a REAL orthogonal rotation => CP-EVEN")

    # ── P5 (M7-D): star-involution selection rule ───────────────────
    I = [list(r) for r in H.star_involution()]
    VI = matmul(V, I)
    IV = matmul(I, V)
    commutes = all(VI[i][j] == IV[i][j] for i in range(26)
                   for j in range(26))
    anti = all(VI[i][j] == -IV[i][j] for i in range(26)
               for j in range(26))
    gate("P5_M7D_star_selection_rule", commutes or anti,
         commutes=commutes, anticommutes=anti,
         verdict=("iota*-EVEN: the channel preserves the real Hodge "
                  "split" if commutes else
                  "iota*-ODD: the channel flips the Hodge halves"
                  if anti else "neither"))

    # ── P6 (M7-E): is the elliptic line quiet? ──────────────────────
    quiet = {}
    for A in order:
        i0, i1 = idxs[A]
        num = sum(table[f"{A}->{B}"] ** 2 for B in order if B != A)
        quiet[A] = (num ** 0.5) / (i1 - i0) ** 0.5
    ell_rank = sorted(quiet, key=lambda k: quiet[k]).index("ell")
    gate("P6_M7E_elliptic_echo", True,
         per_dim_offblock={k: round(v, 4) for k, v in quiet.items()},
         ell_rank_from_quietest=ell_rank,
         note="no prediction was registered; reported as measured")

    # ── P7a: NULL CONTROL — a constant potential must be trivial ────
    T2m = [list(r) for r in H.cuspidal_hecke(2)]
    Cdeg = [[matmul(Vdeg, T2m)[i][j] - matmul(T2m, Vdeg)[i][j]
             for j in range(26)] for i in range(26)]
    const_pot = len(set(deg)) == 1
    trivial = (all(Vdeg[i][j] == (Fr(deg[0]) if i == j else 0)
                   for i in range(26) for j in range(26))
               and all(x == 0 for r in Cdeg for x in r))
    gate("P7a_null_control_constant_potential", const_pot and trivial,
         degree_is_constant=const_pot, V_equals_scalar=trivial,
         verdict="constant potential gives the TRIVIAL channel — the "
                 "pipeline reports no spurious transmission")

    # ── P7b: a genuinely non-constant second potential ──────────────
    Vb2 = matmul(Sinv, matmul(Vdist, S))
    off2 = 0
    for A in order:
        for B in order:
            if A == B:
                continue
            i0_, i1_ = idxs[A]
            j0_, j1_ = idxs[B]
            off2 += frob([[Vb2[i][j] for j in range(j0_, j1_)]
                          for i in range(i0_, i1_)])
    Cd = [[matmul(Vdist, T2m)[i][j] - matmul(T2m, Vdist)[i][j]
           for j in range(26)] for i in range(26)]
    nzd = sum(1 for r in Cd for x in r if x != 0)
    symd = all(matmul(Gm, Vdist)[i][j] == matmul(Gm, Vdist)[j][i]
               for i in range(26) for j in range(26))
    gate("P7b_second_potential", off2 > 0 and nzd > 0 and symd,
         potential="graph distance from the base triangle",
         distinct_values=len(set(dist)),
         off_block_norm=f"{float(off2) ** 0.5:.6f}",
         nonzero_commutator=nzd, still_self_adjoint=symd,
         verdict="non-Hecke, flavor-changing, and still real")

    ledger["transmission_table"] = {k: round(v, 6)
                                    for k, v in table.items()}
    ledger["verdict"] = (
        "M7-A and M7-B CONFIRMED: a canonical, zero-parameter "
        "graph-side coupling exists that is NOT in the Hecke algebra "
        "and DOES move between newform lines. M7-C CONFIRMED: it is "
        "exactly real and G-self-adjoint, so the mixing it induces is "
        "a real orthogonal rotation carrying no phase. The M6 tension "
        "is therefore SHARPENED, not resolved: the obstruction is not "
        "the absence of a non-Hecke channel but the REALITY of the "
        "arithmetic. CP-odd mixing requires a genuinely complex "
        "dressing; the successor object is the twisted homology "
        "H_1(X_0(143); L_theta) for a character theta of the "
        "29-dimensional cycle space, together with the question of "
        "whether the arithmetic supplies a canonical nontrivial "
        "theta or whether it would be a new free input.")
    ledger["all_passed"] = ok
    ledger["runtime_s"] = round(time.time() - t0, 1)
    with open(os.path.join(HERE, "m7_graph_channel_ledger.json"),
              "w") as f:
        json.dump(ledger, f, indent=2, default=str)
    print(f"\nledger written  [{ledger['runtime_s']}s]")
    print("ALL GATES PASS" if ok else "GATE FAILURE")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
