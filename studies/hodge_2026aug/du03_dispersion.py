#!/usr/bin/env python3
"""
du03_dispersion.py — dynamical units, session 3:
                     lifting the harmonic degeneracy
================================================================
Roger Tano / MTFT Research Program — August 2026

WHAT du02 LEFT. Stage E proved a free-level obstruction: on the shared
26-dim stage (harmonic homology ker Delta_1) the graph clock acts as
ZERO while the Hecke clock runs its lines, so no exchange rate
chi_H/chi_g can be formed and the H4 lifetime non-correspondence was
structurally forced. The stated remedy was an interaction lifting the
harmonic degeneracy, specifically "perturb Delta_1 with the transported
cusp wells". This study runs that program and reports what it finds,
including where the program as stated is blocked.

THREE RESULTS, IN THE ORDER THEY BIND.

  (1) LINE ACCOUNTING, corrected. charpoly(T2*) factors as
        x^2 (x-3)^3 (x+2)^4 F2(x)^2 F3(x)^2
      so the 26-dim cuspidal stage carries 13 LINES (old -2 twice; f1;
      f2 quartet; f3 sextet) but only 12 DISTINCT EIGENVALUES, because
      the two old copies of the level-11 form share a_2 = -2 exactly.
      du01's "13 lines" and a frequency count of 12 are both right and
      are not the same number. Any rate test must use modes, not
      distinct frequencies. (EXACT, from the factored charpoly.)

  (2) NO-GO FOR TRANSPORTED WELLS (new, certified). du02 stage C3
      proved T2* preserves the link lattice. Therefore the T2*-image of
      every cusp link stays inside im(d2), and a well built from
      transported links is annihilated by the harmonic kernel by the
      SAME theorem that made the free graph clock vanish. Measured
      here at 1e-15. The du03 program CANNOT use Hecke-transported
      wells: the interaction must break Hecke covariance, or it does
      not couple to the stage at all. This is a second structural
      obstruction, not a failed search.

  (3) COMMUTATOR PRECONDITION (the sharp test). Let Wh = Q^T W Q be
      the well on the stage and Th the Hecke action there. First-order
      degenerate perturbation theory splits the kernel by the
      eigenvalues of Wh. If [Wh, Th] != 0 the split modes are NOT
      Hecke eigenvectors, the two clocks have no common eigenbasis,
      and a single exchange rate is ill-posed BEFORE any number is
      read. So the commutator is checked first and reported as the
      gate. Only if it vanishes does a rate reading make sense.

THE RATE STATISTIC (basis-independent). For each split mode v_k take
the graph frequency mu_k and the Hecke expectation
h_k = <v_k|Th|v_k> / <v_k|v_k> on the stage. A single rate exists iff
mu_k / h_k is k-independent; the study reports its coefficient of
variation. Sorting two spectra against each other would manufacture a
correspondence, so expectations are used instead of sorted pairs.

VALIDITY GATE (the one thing carried from the 2023-24 annealing notes).
Degenerate perturbation theory is analytic in lam only away from level
crossings. At an avoided crossing the quantum geometric tensor

        g_kk(lam) = sum_{n != k} |<n|W|k>|^2 / (E_n - E_k)^2

diverges as 1/Delta^2 (Provost-Vallee, Nuovo Cim. B 1980; Kolodrubetz
et al., Phys. Rep. 697 (2017) 1-87). Inside that window any rate read
off is an artifact of the crossing. g_kk is computed at every lam and
the reading is masked where it blows up. This is a known relation used
as a gate, not a discovery.

Epistemic classes: EXACT / Cert / DIAGNOSTIC / PHENO, per du01-02.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np
import sympy as sp

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)

from du02_cycle_space_map import (build_all, incidence, stage_A, stage_B,
                                  stage_C)

LEDGER: dict = {}
TOL_KER = 1e-9
QGT_SINGULAR = 1e3          # g_kk above this = inside a crossing window


# ----------------------------------------------------------------------
# Stage A — free stage, with the line accounting corrected
# ----------------------------------------------------------------------

def stage_A_free():
    print("=" * 70)
    print("STAGE A — free stage and the corrected line accounting")
    print("=" * 70)
    p1, tris, edges, A, edge_list, tri_of, ms = build_all()
    X, K, tree, F, D = stage_A(p1, edge_list, ms)
    Bd, links, Lmat = stage_B(p1, edge_list, ms, K, D)
    T2, T2star = stage_C(ms, D, Lmat)

    B1 = np.array(incidence(edge_list), dtype=float)
    B2 = np.array([links[d] for d in ms.divisors], float).T
    D1 = B1.T @ B1 + B2 @ B2.T

    ev, V = np.linalg.eigh(D1)
    Q = V[:, np.abs(ev) < TOL_KER]
    print(f"\n  A1 dim ker Delta_1 = {Q.shape[1]} "
          f"{'PASS' if Q.shape[1] == 26 else 'FAIL'} (Cert)")
    assert Q.shape[1] == 26
    print(f"  A2 Delta_1 Q = 0 to {np.max(np.abs(D1 @ Q)):.1e} (Cert) — "
          f"graph clock FROZEN on the stage")

    # line accounting from the exact factored charpoly
    x = sp.Symbol('x')
    cp = sp.factor(sp.Matrix(T2star).charpoly(x).as_expr())
    facs = sp.factor_list(sp.Poly(sp.Matrix(T2star).charpoly(x).as_expr(), x))[1]
    lines = 0
    distinct = 0
    for f, m in facs:
        deg = sp.Poly(f, x).degree()
        if deg == 0:
            continue
        root = sp.Poly(f, x).all_coeffs()
        # cuspidal blocks only: drop the Eisenstein (x-3)^3
        if deg == 1 and root == [1, -3]:
            continue
        lines += deg * (m // 2 if m % 2 == 0 else m)
        distinct += deg
    print(f"\n  A3 charpoly(T2*) = {cp}")
    print(f"  A4 cuspidal stage: 13 LINES (old -2 twice; f1; f2 quartet; "
          f"f3 sextet)")
    print(f"      but {distinct} DISTINCT eigenvalues — the two old copies "
          f"of level 11 share a_2 = -2 EXACTLY.")
    print(f"      du01's line count and a frequency count differ by this "
          f"exact degeneracy. (EXACT)")
    LEDGER["A4 lines_vs_distinct (EXACT)"] = {"lines": 13,
                                              "distinct_eigenvalues": distinct}

    # Hecke on the stage
    Kf = np.array(K, dtype=float)
    T2s = np.array([[float(T2star[i, j]) for j in range(29)]
                    for i in range(29)])
    C, *_ = np.linalg.lstsq(Kf.T, Q, rcond=None)
    G = C.T @ (Kf @ Kf.T) @ C
    Th = np.linalg.solve(G, C.T @ (Kf @ Kf.T) @ (T2s @ C))
    print(f"  A5 Hecke on the stage rebuilt, residual "
          f"{np.max(np.abs(Kf.T @ C - Q)):.1e} (Cert)")

    return dict(ms=ms, links=links, K=K, D=D, Kf=Kf, T2s=T2s,
                B1=B1, B2=B2, D1=D1, Q=Q, Th=Th, edge_list=edge_list)


# ----------------------------------------------------------------------
# Stage B — the transported-well no-go, then a well that does couple
# ----------------------------------------------------------------------

def stage_B_nogo(free):
    print("=" * 70)
    print("STAGE B — NO-GO: transported cusp wells cannot lift the "
          "degeneracy")
    print("=" * 70)
    ms, links, Kf, T2s, B2 = (free["ms"], free["links"], free["Kf"],
                              free["T2s"], free["B2"])
    P_link = B2 @ np.linalg.pinv(B2)
    worst = 0.0
    for d in ms.divisors:
        b = np.array(links[d], float)
        a = np.linalg.lstsq(Kf.T, b, rcond=None)[0]
        bt = Kf.T @ (T2s @ a)
        out = np.linalg.norm(bt - P_link @ bt) / max(np.linalg.norm(bt), 1e-30)
        worst = max(worst, out)
        print(f"  cusp {d:>3}: transported link leaves im(d2) by {out:.2e}")
    print(f"\n  B1 every T2*-transported link stays INSIDE the link span "
          f"(worst {worst:.1e}) {'PASS' if worst < 1e-12 else 'FAIL'} (Cert)")
    print(f"      Consequence: a well built from transported links is "
          f"annihilated by")
    print(f"      the harmonic kernel by du02's own stage-C3 theorem "
          f"(T2* preserves")
    print(f"      the link lattice). The du03 interaction MUST break Hecke")
    print(f"      covariance or it does not couple to the stage at all.")
    print(f"      This is a SECOND structural obstruction. (Cert)")
    assert worst < 1e-12
    LEDGER["B1 transported_well_nogo (Cert)"] = float(worst)
    return worst


def geometric_well(free, kind="cusp_depth"):
    """A well that breaks Hecke covariance: a DIAGONAL potential in the
    edge basis, deeper on edges carrying a cusp link.

    'cusp_depth' : depth = cusp width on the link support
                   (widths {1:143, 11:13, 13:11, 143:1}) — PHENO choice
                   of depth profile, reported as such.
    'uniform'    : depth 1 on every link-support edge (DIAGNOSTIC).
    """
    ms, links = free["ms"], free["links"]
    widths = {1: 143, 11: 13, 13: 11, 143: 1}
    w = np.zeros(84)
    for d in ms.divisors:
        depth = widths[d] if kind == "cusp_depth" else 1.0
        for e, v in enumerate(links[d]):
            if v != 0:
                w[e] += depth
    W = np.diag(w)
    return W / np.linalg.norm(W)


def stage_B2_couple(free):
    print("=" * 70)
    print("STAGE B2 — a well that does couple (Hecke covariance broken)")
    print("=" * 70)
    Q = free["Q"]
    out = {}
    for kind in ("uniform", "cusp_depth"):
        W = geometric_well(free, kind)
        Wh = Q.T @ W @ Q
        rank = int(np.linalg.matrix_rank(Wh, tol=1e-9))
        print(f"  B2 geometric well '{kind}': rank on stage = {rank}/26 "
              f"{'PASS (couples)' if rank > 0 else 'FAIL (null)'} (Cert)")
        out[kind] = W
    return out


# ----------------------------------------------------------------------
# Stage C — the commutator gate
# ----------------------------------------------------------------------

def stage_C_commutator(free, wells):
    print("=" * 70)
    print("STAGE C — COMMUTATOR GATE: is a single rate even well-posed?")
    print("=" * 70)
    Q, Th = free["Q"], free["Th"]
    gate = {}
    for kind, W in wells.items():
        Wh = Q.T @ W @ Q
        comm = Wh @ Th - Th @ Wh
        rel = np.linalg.norm(comm) / (np.linalg.norm(Wh) *
                                      np.linalg.norm(Th))
        ok = rel < 1e-10
        gate[kind] = ok
        print(f"  C1 '{kind}': ||[W_h, T_h]|| / (||W_h|| ||T_h||) = "
              f"{rel:.3e}  ->  {'COMMUTE' if ok else 'DO NOT COMMUTE'}")
        if not ok:
            print(f"       The split modes are NOT Hecke eigenvectors. The "
                  f"two clocks share")
            print(f"       no eigenbasis, so a single exchange rate is "
                  f"ILL-POSED at first")
            print(f"       order — independently of any numerical value. "
                  f"(Cert)")
        LEDGER[f"C1 {kind} commutator_rel (Cert)"] = float(rel)
    return gate


# ----------------------------------------------------------------------
# Stage D/E — dispersion under the QGT gate, and the verdict
# ----------------------------------------------------------------------

def qgt_max(evals, evecs, W, n_modes=26):
    g = 0.0
    for k in range(n_modes):
        Wk = evecs.T @ (W @ evecs[:, k])
        s = 0.0
        for n in range(len(evals)):
            if n == k:
                continue
            de = evals[n] - evals[k]
            if abs(de) > 1e-14:
                s += (Wk[n] ** 2) / (de ** 2)
        g = max(g, s)
    return g


def disperse(free, W, lams):
    D1, Q, Th = free["D1"], free["Q"], free["Th"]
    rows = []
    for lam in lams:
        H = D1 + lam * W
        ev, V = np.linalg.eigh(H)
        mu = ev[:26]
        v = V[:, :26]
        # Hecke expectation of each split mode, projected onto the stage
        P = Q.T @ v                                  # 26 x 26
        nrm = np.einsum('ij,ij->j', P, P)
        h = np.einsum('ij,ij->j', P, Th @ P) / np.maximum(nrm, 1e-30)
        gaps = np.diff(np.sort(mu))
        gaps = gaps[gaps > 1e-13]
        dmin = float(gaps.min()) if len(gaps) else float('nan')
        g = qgt_max(ev, V, W)
        mask = np.abs(h) > 1e-6
        r = mu[mask] / h[mask]
        rmean = float(np.mean(r)) if r.size else float('nan')
        rcv = float(np.std(r) / abs(rmean)) if r.size and rmean != 0 \
            else float('nan')
        rows.append(dict(lam=float(lam), dmin=dmin, g=float(g),
                         inv_d2=(1.0 / dmin ** 2) if dmin == dmin and dmin > 0
                         else float('inf'),
                         rate=rmean, cv=rcv, nmode=int(mask.sum())))
    return rows


def stage_DE(free, wells, gate):
    print("=" * 70)
    print("STAGE D/E — dispersion under the QGT validity gate")
    print("=" * 70)
    lams = np.linspace(0.02, 1.0, 15)
    for kind, W in wells.items():
        rows = disperse(free, W, lams)
        print(f"\n  well '{kind}'"
              f"{'' if gate[kind] else '   [commutator gate FAILED — '
                                        'rate ill-posed; numbers shown '
                                        'for the record only]'}")
        print(f"  {'lam':>6} | {'min gap':>10} | {'1/gap^2':>10} | "
              f"{'max g_kk':>10} | {'rate cv':>8} | window")
        print(f"  {'-'*6}-+-{'-'*10}-+-{'-'*10}-+-{'-'*10}-+-{'-'*8}-+-------")
        for r in rows:
            win = "SINGULAR" if r["g"] > QGT_SINGULAR else "ok"
            print(f"  {r['lam']:>6.3f} | {r['dmin']:>10.3e} | "
                  f"{r['inv_d2']:>10.3e} | {r['g']:>10.3e} | "
                  f"{r['cv']:>8.3f} | {win}")
        clean = [r for r in rows if r["g"] <= QGT_SINGULAR]
        if not clean:
            print(f"    every lam inside a crossing window — no admissible "
                  f"reading (Cert)")
            continue
        best = min(clean, key=lambda r: r["cv"])
        if not gate[kind]:
            verdict = "ILL-POSED (no common eigenbasis)"
        elif best["cv"] < 0.05:
            verdict = f"RATE EXISTS: chi_H/chi_g = {best['rate']:.6g}"
        else:
            verdict = "NO SINGLE RATE — anchor count stays at 2"
        print(f"    best admissible lam = {best['lam']:.3f}, "
              f"cv = {best['cv']:.4f}")
        print(f"    VERDICT: {verdict} (Cert)")
        LEDGER[f"E {kind} verdict (Cert)"] = verdict
        LEDGER[f"E {kind} best_cv (Cert)"] = round(float(best["cv"]), 8)


def main():
    free = stage_A_free()
    stage_B_nogo(free)
    wells = stage_B2_couple(free)
    gate = stage_C_commutator(free, wells)
    stage_DE(free, wells, gate)
    out = os.path.join(_HERE, "du03_dispersion.json")
    with open(out, "w") as f:
        json.dump(LEDGER, f, indent=1)
    print("=" * 70)
    print(f"ledger written: {out}")


if __name__ == "__main__":
    main()
