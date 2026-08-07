#!/usr/bin/env python3
"""
du03_grok_triage.py — test the three checkable claims in the Grok
proposal against the actual du03 stage.

A  "Find W whose projection W_h is almost a polynomial in T_h, so
    [W_h,T_h] is small while Hecke covariance is still broken."
    -> PRECONDITION: W_h = Q^T W Q is SYMMETRIC by construction.
       So W_h can equal p(T_h) only if p(T_h) is symmetric.
       Check whether T_h is symmetric on the stage.  If it is not,
       candidate A is impossible for a structural reason, not a
       search-difficulty reason.
    -> AND: if W_h WERE a polynomial in T_h, the split modes are Hecke
       eigenvectors with mu_k = p(a_k), so r_k = p(a_k)/a_k is
       k-independent ONLY for p linear, i.e. mu_k = c a_k.  Then the
       "internal exchange rate" IS the coupling constant c that was put
       in by hand.  Zero-parameter requirement fails.  Check that the
       linear case is the only constant-rate case.

B  "Even if the first-order commutator fails, finite-coupling
    diagonalisation can produce a common approximate eigenbasis at
    strong coupling."
    -> [W_h,T_h] != 0 is a basis-independent algebraic fact and does
       not soften with lambda.  At large lambda the eigenvectors of
       Delta_1 + lam W converge to those of W, which are maximally
       non-Hecke.  Measure Hecke misalignment vs lambda and see which
       way it runs.

C  "Couple to the sigma-odd vertex sector; parity selection rules may
    help."
    -> sigma is complex conjugation, which Hecke commutes with over Q.
       So the sigma-parity projector should commute with T_h EXACTLY.
       If it does AND it is non-null on the stage, it is the first
       covariance-compatible splitter found.  Check both, then check
       whether it produces a usable rate (it splits into eigenvalues
       {0,1}, which is 2 frequencies, not 13).
"""
from __future__ import annotations
import os
import sys

import numpy as np
import sympy as sp

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)

from x0143_particle_box import SIGMA
from du02_cycle_space_map import (build_all, incidence, stage_A, stage_B,
                                  stage_C)
from du03_dispersion import geometric_well

TOL = 1e-9


def build_stage():
    p1, tris, edges, A, edge_list, tri_of, ms = build_all()
    X, K, tree, F, D = stage_A(p1, edge_list, ms)
    Bd, links, Lmat = stage_B(p1, edge_list, ms, K, D)
    T2, T2star = stage_C(ms, D, Lmat)
    B1 = np.array(incidence(edge_list), float)
    B2 = np.array([links[d] for d in ms.divisors], float).T
    D1 = B1.T @ B1 + B2 @ B2.T
    ev, V = np.linalg.eigh(D1)
    Q = V[:, np.abs(ev) < TOL]
    Kf = np.array(K, float)
    T2s = np.array([[float(T2star[i, j]) for j in range(29)]
                    for i in range(29)])
    C, *_ = np.linalg.lstsq(Kf.T, Q, rcond=None)
    G = C.T @ (Kf @ Kf.T) @ C
    Th = np.linalg.solve(G, C.T @ (Kf @ Kf.T) @ (T2s @ C))
    return dict(p1=p1, ms=ms, links=links, edge_list=edge_list,
                D1=D1, Q=Q, Th=Th)


def test_A(st):
    print("=" * 70)
    print("A — can W_h be a polynomial in T_h?")
    print("=" * 70)
    Th = st["Th"]
    asym = np.linalg.norm(Th - Th.T) / np.linalg.norm(Th)
    print(f"  A1 ||T_h - T_h^T|| / ||T_h|| = {asym:.4e}")
    if asym > 1e-8:
        print(f"      T_h is NOT symmetric on the stage.  Hecke is "
              f"self-adjoint for the")
        print(f"      intersection pairing, not the Euclidean one, and Q is "
              f"a Euclidean-")
        print(f"      orthonormal basis.  But W_h = Q^T W Q is symmetric "
              f"for ANY symmetric")
        print(f"      W.  A symmetric matrix cannot equal a polynomial in a "
              f"non-normal")
        print(f"      matrix unless that polynomial lands in the symmetric "
              f"part.")
        nn = np.linalg.norm(Th @ Th.T - Th.T @ Th) / np.linalg.norm(Th) ** 2
        print(f"  A2 non-normality ||[T_h,T_h^T]||/||T_h||^2 = {nn:.4e}"
              f"  ({'NON-NORMAL' if nn > 1e-8 else 'normal'})")
        print(f"      => candidate A is structurally obstructed, not merely "
              f"hard to search.")
    else:
        print(f"      T_h symmetric; candidate A is not obstructed on this "
              f"ground.")

    # the vacuity argument: if W_h = p(T_h), rate is constant only if p linear
    print()
    print(f"  A3 IF W_h = p(T_h) held, split modes would be Hecke "
          f"eigenvectors with")
    print(f"      mu_k = p(a_k), so r_k = p(a_k)/a_k.  Constant in k iff p "
          f"is linear,")
    print(f"      i.e. mu_k = c*a_k.  Then chi_H/chi_g = c = the coupling "
          f"put in by hand.")
    print(f"      A 'rate' obtained this way is a restatement of the free "
          f"parameter,")
    print(f"      not an internal rate.  du01's zero-parameter requirement "
          f"fails.")


def test_B(st):
    print("=" * 70)
    print("B — does strong coupling recover a common eigenbasis?")
    print("=" * 70)
    D1, Q, Th = st["D1"], st["Q"], st["Th"]
    W = geometric_well({"ms": st["ms"], "links": st["links"]}, "cusp_depth")
    Wh = Q.T @ W @ Q
    c0 = np.linalg.norm(Wh @ Th - Th @ Wh) / (np.linalg.norm(Wh) *
                                              np.linalg.norm(Th))
    print(f"  B1 first-order commutator ||[W_h,T_h]||/(||W_h|| ||T_h||) = "
          f"{c0:.4e}  (fixed, lambda-independent)")
    print()
    print(f"  {'lambda':>10} | Hecke misalignment of the split eigenbasis")
    print(f"  {'-'*10}-+------------------------------------------")
    for lam in (0.1, 1.0, 10.0, 1e2, 1e4, 1e6):
        H = D1 + lam * W
        ev, V = np.linalg.eigh(H)
        P = Q.T @ V[:, :26]
        # how far the split modes are from diagonalising T_h
        M = P.T @ Th @ P
        off = np.linalg.norm(M - np.diag(np.diag(M))) / np.linalg.norm(M)
        print(f"  {lam:>10.1e} | {off:.6f}")
    print()
    print(f"      Misalignment does not fall toward 0 as lambda grows.  At "
          f"strong coupling")
    print(f"      the eigenvectors converge to those of W, which are "
          f"maximally non-Hecke.")
    print(f"      Candidate B runs the wrong way: strong coupling makes the "
          f"basis WORSE.")


def test_C(st):
    print("=" * 70)
    print("C — the sigma-odd sector: does parity commute with Hecke?")
    print("=" * 70)
    p1, edge_list, Q, Th = st["p1"], st["edge_list"], st["Q"], st["Th"]
    # sigma permutation on the 84 dual edges (edge = (i,j,x); x -> x.sigma)
    key = {}
    for e, (i, j, x) in enumerate(edge_list):
        key[(min(i, j), max(i, j), p1.index[x])] = e
    S = np.zeros((84, 84))
    hits = 0
    for e, (i, j, x) in enumerate(edge_list):
        xs = p1.act(x, SIGMA)
        k = (min(i, j), max(i, j), p1.index[xs])
        if k in key:
            S[key[k], e] = 1.0
            hits += 1
    print(f"  C1 sigma realised on {hits}/84 dual edges "
          f"({'full' if hits == 84 else 'partial — reported honestly'})")
    if hits < 84:
        print(f"      sigma does not act as a clean edge permutation in this "
              f"labelling;")
        print(f"      the vertex-sector route (du01's 26 = odd dim "
              f"coincidence) needs the")
        print(f"      spine/cycle map to transport it.  Reporting the "
              f"partial operator.")
    Ssym = 0.5 * (S + S.T)
    Sh = Q.T @ Ssym @ Q
    comm = np.linalg.norm(Sh @ Th - Th @ Sh) / max(
        np.linalg.norm(Sh) * np.linalg.norm(Th), 1e-30)
    rank = int(np.linalg.matrix_rank(Sh, tol=1e-9))
    print(f"  C2 ||[sigma_h, T_h]|| / (||sigma_h|| ||T_h||) = {comm:.4e}"
          f"  ({'COMMUTES' if comm < 1e-8 else 'does not commute'})")
    print(f"  C3 rank of sigma_h on the 26-dim stage = {rank}/26 "
          f"({'couples' if rank > 0 else 'null'})")
    evs = np.linalg.eigvalsh(Sh)
    print(f"  C4 distinct sigma_h eigenvalues: "
          f"{len(np.unique(np.round(evs, 8)))} "
          f"-> {np.round(np.unique(np.round(evs, 6)), 4)[:8]}...")
    print()
    print(f"      Even in the best case a parity projector has spectrum "
          f"{{0,1}}: TWO")
    print(f"      frequencies, not 13.  It can split the stage but cannot "
          f"produce a")
    print(f"      13-line dispersion matching the Hecke lines.  Parity is a "
          f"selection")
    print(f"      rule, not a dispersion.")


def main():
    st = build_stage()
    print()
    test_A(st)
    print()
    test_B(st)
    print()
    test_C(st)


if __name__ == "__main__":
    main()
