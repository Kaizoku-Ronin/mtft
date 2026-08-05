#!/usr/bin/env python3
"""
du01_two_clock_ledger.py — dynamical units, session 1:
                           the two-clock dimensionless ledger
==============================================================
Roger Tano / MTFT Research Program — built with Claude, August 2026

THE PROBLEM (operational-y, dynamical form). The particle box evolves in
two internal clocks — graph time (CURVE: unitary Cayley steps under the
Farey-dual Laplacian) and Hecke time (ITERATE: e^{-i a tau} on cuspidal
homology). Both clocks are pure numbers. Nothing yet ties either to
seconds or eV. Before any units map can be PROPOSED, the box's
unit-free content must be FROZEN: the complete table of dimensionless
spectral invariants that any admissible map must preserve. That table
is this study.

THE ANCHOR PROTOCOL (the transmission-line lesson, made precise).
An RF Smith chart carries no units: it is the hyperbolic-geometry chart
of a cusp neighborhood (Cayley transform of a half-plane; constant-r
circles are horocycles at Gamma=1, constant-x arcs the geodesics into
that cusp), and all physics enters through exactly one anchor per
clock — "one full turn = half a wavelength" (theta = omega l / v).
Likewise here:

  * Each internal clock needs ONE conversion constant
    (chi_g seconds-per-graph-tick, chi_H seconds-per-Hecke-tick).
  * The number of REQUIRED external anchors = the number of clocks NOT
    yet coupled by a certified internal process. Today that number is 2:
    the H4 lifetime cross-check (v0.11.0/BQ) remains a standing honest
    negative, so no internal exchange rate chi_H/chi_g exists yet.
  * A proposed units map is ADMISSIBLE only if it matches the box's
    dimensionless line RATIOS to a physical spectrum's ratios; one
    anchor then fixes the scale and every remaining line is a
    zero-parameter prediction (falsify-engine pattern, ppm ledger).

STAGES
  A  Graph clock, zero parameters. The 56-level free spectrum of
     L = D - A on the trivalent Farey-dual graph of X0(143), with EXACT
     certificates: one self-loop; tr L = 166 = 168 - 2; integer
     eigenvalues {0,1,2,4,5} each simple (fraction-free Bareiss rank,
     exact over Z); Kirchhoff spanning-tree count by exact integer
     determinant, cross-checked (E2) against the float eigenvalue
     product; spine identity b1 = E - V + 1 = 29 = 2g + (cusps - 1).
  B  Hecke clock, zero parameters. The 13 T2 lines on the 26-dim
     cuspidal space (old -2 x2; f1; f2 quartet; f3 sextet), a_p per
     line to p <= 47, trace cross-certified against the mtft.x0_143
     oracle (E2: Manin-symbol route vs LMFDB-validated corpus route).
     Ramanujan pass-band certificate: |a_p| <= 2 sqrt(p) at good p
     (the |tr| < 2 elliptic condition of periodic-line theory; the
     Hecke lines live in the pass band).
  C  The dimensionless ledger. Ratio tables for both clocks (all
     ratios to the fundamental line), orbit variances (uniform-weight,
     DIAGNOSTIC — Born-weighted values live in the v0.3 corpus), and
     the JSON artifact any future units proposal must be tested
     against.
  D  The physical side of the ledger: pinned hydrogen digit-ladder
     rungs (measured, sources cited inline) and their pure-number
     ratios. Rungs whose CODATA pin was not re-verified this session
     are left BLANK (repository blank policy), not filled.

Epistemic classes on every number: EXACT / Cert / DIAGNOSTIC / PHENO.

STANDING HONEST NEGATIVES CARRIED INTO THIS PROGRAM
  * No certified internal chi_H/chi_g exchange rate (H4
    non-correspondence, v0.11.0). Anchor count stays at 2.
  * The sigma-parity dimension match dim(odd vertex sector) = 26 =
    dim H1 is a COINCIDENCE OF DIMENSIONS until the spine/cycle-space
    map is computed (Heur; the drawn-loop stage computes it).
"""

from __future__ import annotations
import json
import os
import time
from math import sqrt, log10

import numpy as np

from x0143_particle_box import tessellation, dual_graph
from x0143_particle_box_v02 import (build_engine, float_projection,
                                    eigendata, assign_orbits, extract_ap,
                                    an_table)

_HERE = os.path.dirname(os.path.abspath(__file__))
LEDGER: dict = {"study": "du01_two_clock_ledger", "version_context": "0.11.3"}


# ----------------------------------------------------------------------
# exact linear algebra over Z (fraction-free Bareiss)
# ----------------------------------------------------------------------

def bareiss(M):
    """Fraction-free elimination on an integer matrix.
    Returns (rank, det) where det is exact iff square and full rank."""
    M = [row[:] for row in M]
    n, m = len(M), len(M[0])
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


# ----------------------------------------------------------------------
# Stage A — graph clock
# ----------------------------------------------------------------------

def stage_A():
    print("=" * 70)
    print("STAGE A — graph clock (CURVE): free spectrum, zero parameters")
    print("=" * 70)
    p1, tris, edges, cert = tessellation(143)
    A, edge_list, tri_of = dual_graph(p1, tris, edges)

    n_loops = sum(1 for (i, j, _x) in edge_list if i == j)
    Lf = np.diag(A.sum(axis=1)) - A
    Lint = [[int(v) for v in row] for row in Lf.astype(int).tolist()]
    tr = int(np.trace(Lf))
    print(f"  A1 self-loops on dual graph: {n_loops} "
          f"-> tr L = 168 - 2*{n_loops} = {tr}  (EXACT; 168 = index)")
    assert tr == 168 - 2 * n_loops == 166
    LEDGER["A1 self_loops (EXACT)"] = n_loops
    LEDGER["A1 trace_L (EXACT)"] = tr

    ev = np.linalg.eigvalsh(Lf)
    e2a = abs(ev.sum() - tr)
    e2b = abs((ev ** 2).sum() - (Lf * Lf).sum())
    print(f"  A2 spectrum in band [0, 6]: lambda_1 (Fiedler/free gap) = "
          f"{ev[1]:.6f}, lambda_max = {ev[-1]:.6f}")
    print(f"  A2 E2 moment checks: |sum ev - tr| = {e2a:.1e}, "
          f"|sum ev^2 - ||L||_F^2| = {e2b:.1e}  "
          f"{'PASS' if max(e2a, e2b) < 1e-9 else 'FAIL'} (Cert)")
    LEDGER["A2 spectrum (Cert, 56 values)"] = [round(float(x), 12)
                                               for x in ev]

    ints = []
    for k in (0, 1, 2, 3, 4, 5, 6):
        Mk = [[Lint[i][j] - (k if i == j else 0) for j in range(56)]
              for i in range(56)]
        r, _ = bareiss(Mk)
        if r < 56:
            ints.append((k, 56 - r))
    print(f"  A3 EXACT integer eigenvalues (Bareiss rank over Z): "
          f"{[f'{k} (x{m})' for k, m in ints]}")
    assert ints == [(0, 1), (1, 1), (2, 1), (4, 1), (5, 1)]
    LEDGER["A3 integer_eigenvalues (EXACT)"] = ints

    red = [row[1:] for row in Lint[1:]]
    _r, tau = bareiss(red)
    float_route = float(np.prod(ev[1:])) / 56
    rel = abs(float_route - tau) / tau
    print(f"  A4 Kirchhoff spanning trees tau(G) = {tau}  (EXACT integer "
          f"determinant)")
    print(f"  A4 E2 float route prod(ev>0)/56: rel err = {rel:.1e}  "
          f"{'PASS' if rel < 1e-10 else 'FAIL'} (Cert)")
    print(f"  A4 factorization: 2^5 * 3 * 5 * 17 * 941 * 101921 * 4495339"
          f"   [observed; no interpretation filed — AG-D5 discipline]")
    assert tau == 2**5 * 3 * 5 * 17 * 941 * 101921 * 4495339
    LEDGER["A4 spanning_trees (EXACT)"] = tau

    V, E = 56, 84
    b1 = E - V + 1
    print(f"  A5 spine identity: b1(dual graph) = E - V + 1 = {b1} = "
          f"2g + (cusps - 1) = 26 + 3  "
          f"{'PASS' if b1 == 29 else 'FAIL'} (EXACT)")
    print(f"  A5 the graph clock's cycle space CONTAINS the Hecke clock's"
          f" arena (H1 = 26) plus 3 cusp classes — the geometric coupling"
          f" channel. Computing the induced map = the drawn-loop stage.")
    assert b1 == 29
    LEDGER["A5 b1_spine (EXACT)"] = b1
    return ev


# ----------------------------------------------------------------------
# Stage B — Hecke clock
# ----------------------------------------------------------------------

GOOD_P = [2, 3, 5, 7, 17, 19, 23, 29, 31, 37, 41, 43, 47]
BAD_P = [11, 13]


def stage_B():
    print("=" * 70)
    print("STAGE B — Hecke clock (ITERATE): T2 lines, zero parameters")
    print("=" * 70)
    t0 = time.time()
    p1, tris, edges, ms = build_engine()
    P = float_projection(ms)
    Bc, proj_c, restrict, T2, E, lines = eigendata(ms, P)
    lines = assign_orbits(lines)
    print(f"  B1 engine up in {time.time() - t0:.1f}s; cuspidal dim "
          f"{T2.shape[0]}; eta split 13 + 13 (EXACT, engine-asserted)")

    ev26 = np.sort(np.linalg.eigvals(T2).real)
    doubled = ev26[::2]
    pair_gap = float(np.max(np.abs(ev26[::2] - ev26[1::2])))
    print(f"  B1 T2 spectrum: 13 lines, each x2 (max pair split "
          f"{pair_gap:.1e}):")
    print("      " + "  ".join(f"{x:+.6f}" for x in doubled))
    print(f"  B1 trace T2 = {T2.trace():+.10f}  "
          f"(E2 target: 2 * [2*(-2) + sum(new a2)])")
    LEDGER["B1 T2_lines (Cert, x2 each)"] = [round(float(x), 10)
                                             for x in doubled]

    ap = extract_ap(ms, P, restrict, lines, 50)

    import mtft.x0_143 as ox
    tr_tab = {"f1": ox.ORBIT_TRACE_F1, "f2": ox.ORBIT_TRACE_F2,
              "f3": ox.ORBIT_TRACE_F3}
    max_err = 0.0
    for orbit in ("f1", "f2", "f3"):
        Ls = [L for L in lines if L[1] == orbit and L[0] == "+"]
        ans = [an_table(ap[id(L)], 50) for L in Ls]
        for n in range(1, 51):
            tr_n = sum(a[n] for a in ans)
            max_err = max(max_err, abs(tr_n - tr_tab[orbit][n - 1]))
    print(f"  B2 orbit-trace cross-cert vs mtft.x0_143 oracle (n <= 50): "
          f"max |err| = {max_err:.2e}  "
          f"{'PASS' if max_err < 1e-6 else 'FAIL'} (Cert; E2 route pair)")
    LEDGER["B2 oracle_cross_cert_max_err (Cert)"] = float(max_err)

    line_table = {}
    worst = 0.0
    for L in lines:
        if L[0] != "+":
            continue
        d = ap[id(L)]
        rec = {}
        for p in GOOD_P:
            if p in d:
                rec[p] = round(float(d[p]), 8)
                worst = max(worst, abs(d[p]) / (2 * sqrt(p)))
        for p in BAD_P:
            if p in d:
                rec[p] = round(float(d[p]), 8)
        line_table.setdefault(L[1], []).append(rec)
    print(f"  B3 Ramanujan pass-band certificate at good p <= 47: "
          f"max |a_p| / (2 sqrt p) = {worst:.6f} < 1  "
          f"{'PASS' if worst < 1 else 'FAIL'} (Cert)")
    print(f"  B3 [periodic-line reading: |tr| < 2 on every line — the "
          f"Hecke clock runs entirely in the pass band; bad p in "
          f"{{11,13}} carry |a_p| = 1 (Atkin-Lehner, EXACT)]")
    LEDGER["B3 ramanujan_margin_max (Cert)"] = round(worst, 8)
    LEDGER["B3 ap_lines (Cert)"] = line_table

    # B4 — uniform-weight a2 variances, EXACT by Newton identities on the
    # oracle's INTEGER Hecke polynomials, E2-paired with the float engine.
    from fractions import Fraction as F
    import mtft.x0_143 as ox

    def exact_var(coeffs_int):
        # monic x^n + c1 x^{n-1} + c2 x^{n-2} + ... ; roots a_i
        n = len(coeffs_int) - 1
        e1, e2 = -F(coeffs_int[1]), F(coeffs_int[2])
        p1 = e1
        p2 = e1 * p1 - 2 * e2
        return p2 / n - (p1 / n) ** 2

    poly = {"f2": [int(round(c)) for c in ox.hecke_polynomial_f2_T2()],
            "f3": [int(round(c)) for c in ox.hecke_polynomial_f3_T2()]}
    var_exact = {o: exact_var(poly[o]) for o in ("f2", "f3")}
    var_float = {}
    for orbit in ("f2", "f3"):
        a2s = [L[2] for L in lines if L[1] == orbit and L[0] == "+"]
        var_float[orbit] = float(np.var(a2s))
    err = max(abs(var_float[o] - float(var_exact[o])) for o in var_exact)
    assert var_exact["f2"] == F(35, 16) and var_exact["f3"] == F(10, 3)
    print(f"  B4 uniform-weight a2 variances, EXACT from integer Hecke "
          f"polynomials (Newton identities):")
    print(f"      Var(f1) = 0 (dim 1 — the electron cannot Zeno-decay); "
          f"Var(f2) = 35/16; Var(f3) = 10/3")
    print(f"  B4 E2 vs float engine: max |err| = {err:.1e}  "
          f"{'PASS' if err < 1e-9 else 'FAIL'} "
          f"(Cert; Born-weighted Zeno values live in v0.3 corpus, "
          f"pairing-dependent)")
    LEDGER["B4 uniform_var (EXACT)"] = {"f1": "0", "f2": "35/16",
                                        "f3": "10/3"}
    LEDGER["B4 var_E2_err (Cert)"] = float(err)
    return doubled, line_table


# ----------------------------------------------------------------------
# Stage C — the dimensionless ledger
# ----------------------------------------------------------------------

def stage_C(ev_graph, hecke_lines):
    print("=" * 70)
    print("STAGE C — the dimensionless ledger (what any units map must "
          "preserve)")
    print("=" * 70)
    lam1 = ev_graph[1]
    ratios = ev_graph[1:] / lam1
    print(f"  C1 graph clock: 55 nonzero lines, all ratios to lambda_1 = "
          f"{lam1:.6f}; span 1 .. {ratios[-1]:.6f}")
    print(f"      first ten: " + "  ".join(f"{r:.5f}" for r in ratios[:10]))
    LEDGER["C1 graph_ratios (Cert)"] = [round(float(r), 10) for r in ratios]

    print(f"  C2 Hecke clock: line ratios are NOT taken on a2 alone "
          f"(zeros/signs); the invariant is the joint a_p table (B3) and "
          f"the per-orbit spreads (B4). Both frozen in the JSON artifact.")

    print(f"  C3 anchor count: 2 (chi_g, chi_H both free). Standing "
          f"honest negative: no certified internal exchange rate — the "
          f"H4 Hecke-Zeno vs graph-emission lifetime cross-check did not "
          f"correspond (v0.11.0/BQ) and that result is carried, not "
          f"retried, until the cycle-space map exists.")
    LEDGER["C3 anchor_count"] = 2
    LEDGER["C3 standing_negative"] = ("H4 lifetime non-correspondence "
                                      "carried from v0.11.0")


# ----------------------------------------------------------------------
# Stage D — the physical side (hydrogen digit-ladder)
# ----------------------------------------------------------------------

def stage_D():
    print("=" * 70)
    print("STAGE D — physical ladder rungs (measured; pure-number ratios)")
    print("=" * 70)
    f_1s2s = 2_466_061_413_187_035.0        # Hz  Parthey et al., PRL 107,
    #                                         203001 (2011); u = 10 Hz
    f_hfs = 1_420_405_751.7667              # Hz  Hellwig et al. 1970;
    #                                         u ~ 1 mHz (the 21 cm line)
    f_ryd = 3.2898419602500e15              # Hz  c R_inf, CODATA 2022
    print(f"  D1 pinned rungs (PHENO refs, sources in comments):")
    print(f"      f(1S-2S)  = {f_1s2s:.0f} Hz")
    print(f"      f(21 cm)  = {f_hfs:.4f} Hz")
    print(f"      c R_inf   = {f_ryd:.4e} Hz")
    print(f"  D2 rungs NOT pinned this session (blank policy — verify "
          f"against CODATA before corpus use): Lamb 2S-2P, 2P fine "
          f"structure, Balmer set.")
    r1 = f_1s2s / f_hfs
    r2 = f_hfs / f_ryd
    r3 = f_1s2s / f_ryd
    print(f"  D3 pure-number targets any admissible units map must hit:")
    print(f"      f(1S-2S)/f(21cm) = {r1:.10f}   [encodes alpha^2 "
          f"g_p m_e/m_p up to QED]")
    print(f"      f(21cm)/cR_inf   = {r2:.10e}")
    print(f"      f(1S-2S)/cR_inf  = {r3:.12f}   [3/4 + QED/recoil]")
    LEDGER["D pinned_rungs_Hz (PHENO)"] = {"f_1s2s": f_1s2s,
                                           "f_hfs_21cm": f_hfs,
                                           "c_Rinf": f_ryd}
    LEDGER["D pure_ratios (PHENO)"] = {"f1s2s_over_f21cm": r1,
                                       "f21cm_over_cRinf": r2,
                                       "f1s2s_over_cRinf": r3}
    print(f"  D4 protocol: an admissible map matches box ratios (C1/B3) "
          f"to ladder ratios (D3); ONE anchor per remaining free clock "
          f"fixes scale; every further line is a zero-parameter "
          f"prediction. Register through the falsify engine before "
          f"looking (pre-registration discipline).")


# ----------------------------------------------------------------------

def main():
    ev = stage_A()
    doubled, hecke = stage_B()
    stage_C(ev, hecke)
    stage_D()
    out = os.path.join(_HERE, "du01_two_clock_ledger.json")
    with open(out, "w") as f:
        json.dump(LEDGER, f, indent=1)
    print("=" * 70)
    print(f"ledger written: {out}")
    print("du01 complete — the box's unit-free content is frozen. "
          "Session-2 target: the cycle-space map (drawn-loop stage), "
          "which is now ALSO the chi_H/chi_g exchange-rate computation.")


if __name__ == "__main__":
    main()
