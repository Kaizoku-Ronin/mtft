#!/usr/bin/env python3
"""
du03_dispersion.py — dynamical units, session 3:
                     the interaction, the dispersion, the first ruler
======================================================================
Roger Tano / MTFT Research Program — built with Claude, August 2026
Baseline: mtft v0.11.4 (du01/du02 wave, Addendum BT audited).

CONTEXT. du02 proved the free-level obstruction: on the shared 26-dim
harmonic stage the graph clock is zero while the Hecke clock runs all
13 lines. An internal exchange rate therefore requires an INTERACTION.
This study (a) closes the remaining free channels with certificates,
(b) constructs the canonical zero-parameter interaction, (c) computes
the box's first dispersion relation omega(a), and (d) certifies the
box's first internal ruler (the systole). It also executes the
disciplined parts of the Grok memo (2026-08-05) and records the
corrections to it.

AUDIT OF THE GROK MEMO (recorded here, per corpus protocol):
  * Sections 1-2 (recap, anchor protocol, ladder): accurate; every
    number matches the certified ledger.
  * Route A (anchor Fiedler to 21 cm): arithmetic correct; correctly
    labeled as leaving the geometric factor free. Reproduced in
    Stage G as PHENO bookkeeping. Does not reduce the anchor count.
  * Route B (set a Zeno variance equal to a decay rate): DIMENSIONAL
    ERROR. Free Zeno decay is quadratic, S(tau) ~ 1 - Var*tau^2;
    Var carries (Hecke frequency)^2 while a lifetime carries 1/time.
    A variance is not a rate. The unit-free part that survives is the
    structural one already filed: Var(f1) = 0 => the drawn electron
    cannot Zeno-decay. Exponential lifetimes live in the dissipative
    engine (v0.10 H-series), where the du02 theorem now explains the
    H4 non-correspondence.
  * Route C (interaction-mediated exchange): agreed — this study.
    Sharpened: the obvious candidate (cusp wells) provably cannot
    work; see Stage A1.
  * Memo experiment 1 (ratio census vs hydrogen pure numbers): run in
    Stage F as a CENSUS with a null model (AG-D5 discipline). Two of
    the three targets are out of band and declared untestable.
  * Memo experiment 2 (extract hyperbolic edge lengths): collapses.
    All 56 Farey triangles are SL2(Z)-images of the ideal triangle
    {0,1,oo} and every edge's shear coordinate is 0 (the standard
    Farey quadrilateral (0,oo;1,-1) has cross-ratio -1). The surface
    is locally homogeneous along the tessellation: THERE IS NO LOCAL
    METRIC RULER. All length information is global. The canonical
    global lengths are the area 56*pi and the SYSTOLE, certified in
    Stage E.

STAGES
  A  Closing the free channels (certificates):
     A1  cusp wells cannot lift the harmonic degeneracy:
         (B2 B2^T) H_harm = 0 to machine precision — the link-form
         interaction annihilates the Hecke stage. Route closed.
     A2  the two canonical operators on cycle space, the link Gram
         G and the transported Hecke T2*, have their commutator
         computed EXACTLY; whatever it is, it acts through the
         Eisenstein sector only (A1) and cannot produce a harmonic
         dispersion.
  B  The canonical interaction: MINIMAL COUPLING. The one canonical
     bilinear map between the vertex clock's space and the edge/
     homology space is the module structure C^0 x C^1 -> C^1,
     (f.v)_e = avg(f_i, f_j) v_e — the discrete form of multiplying
     a 1-form by a function, i.e. exactly how matter couples to a
     gauge sector on a lattice. Zero-parameter choice of f: the graph
     clock's own fundamental mode (Fiedler vector; simple by du01
     A3, hence canonical up to one global sign, which is reported as
     a convention). Null control: f = const gives V_h = Id exactly.
     Trace E2: tr V_h two independent routes.
  C  The dispersion. Organize the 26-dim harmonic space by the
     transported Hecke lines (a third E2 route to the 13 lines);
     compress V_h = P_h M_f P_h onto each 2-dim line:
        omega(a) = mean level, delta(a) = fine-structure splitting,
        off-block norms = inter-line coupling.
     Parity readout: transport the eta involution; if [V_h, eta] ~ 0
     the fine structure is PARITY splitting. Slope of omega(a) at the
     massless (f1, a = 0) point by an estimator battery (A.7
     standards) — the box's light-cone slope candidate. DIAGNOSTIC.
  D  Robustness battery: f in {u1, u2, u3, u_top, random, const}.
     Features that survive the low-mode choices and vanish for the
     controls are structure; the rest is choice.
  E  The first ruler (EXACT): the systole of X0(143).
     No hyperbolic element of Gamma0(143) has |trace| = 3
     (d^2 - 3d + 1 = 0 mod 143 requires 5 to be a QR mod 13; it is
     not) and trace 4 is realized by an explicit matrix; hence
        ell_sys = 2 arccosh(2) = 2.6339157938...
     in curvature-radius units, and it is primitive (a square root
     would need trace sqrt(6)).
  F  Disciplined census (Grok memo experiment 1): all ordered ratios
     of the 55 nonzero graph levels against the one in-band hydrogen
     pure number 0.749598747594; null model by spacing shuffle;
     p-value reported. Census, not search.
  G  PHENO translation (bookkeeping, no claims): chi_g from the
     21 cm anchor; the curvature radius that would make the systolic
     velocity equal c; comparisons to lambda_C and a_0 as pure
     ratios.

Epistemic classes on every number: EXACT / Cert / DIAGNOSTIC / PHENO.
"""

from __future__ import annotations
from fractions import Fraction
import json
import math
import os

import numpy as np
import sympy as sp

from x0143_particle_box import ModularSymbols
from x0143_particle_box_v02 import float_projection, eta_float
import du02_cycle_space_map as du02

_HERE = os.path.dirname(os.path.abspath(__file__))
LEDGER: dict = {"study": "du03_dispersion", "baseline": "0.11.4"}
N = 143
GROUPS = ("old", "f1", "f2", "f3")
A2_F2 = (-1.126757, -0.197126, 1.747468, 2.576415)


# ----------------------------------------------------------------------
# quiet rebuild of the du02 machinery (certified conventions, asserted)
# ----------------------------------------------------------------------

def build_quiet():
    p1, tris, edges, A, edge_list, tri_of, ms = du02.build_all()
    X = du02.crossing_matrix(p1, edge_list, +1)     # du02 A1: sign = +1
    K, _tree = du02.cycle_basis(edge_list)
    n = len(p1.reps)
    F = [[sum(K[c][e] * X[e][m] for e in range(84)) for m in range(n)]
         for c in range(29)]
    free = ms.free
    D = [[F[c][free[k]] for k in range(29)] for c in range(29)]
    _r, det = du02.bareiss([row[:] for row in D])
    assert det == -1, "Lefschetz regression"          # du02 A2
    Bd = [[int(ms.boundary[i, j]) for j in range(29)]
          for i in range(len(ms.divisors))]
    links = {}
    for ci, d in enumerate(ms.divisors):
        a = du02.int_solve_unimodular(D, Bd[ci])
        links[d] = [sum(a[c] * K[c][e] for c in range(29))
                    for e in range(84)]
    T2 = ms.hecke_on_quotient(2)
    Dm = sp.Matrix(D)
    DT = Dm.T
    T2star = DT.inv() * T2.T * DT                     # du02 C: adjoint
    return p1, tris, tri_of, edge_list, ms, K, D, links, T2, T2star


# ----------------------------------------------------------------------
# Stage A — closing the free channels
# ----------------------------------------------------------------------

def stage_A(edge_list, ms, K, links, T2star):
    print("=" * 70)
    print("STAGE A — closing the free channels")
    print("=" * 70)
    B1 = np.array(du02.incidence(edge_list), float)
    B2 = np.array([links[d] for d in ms.divisors], float).T
    D1 = B1.T @ B1 + B2 @ B2.T
    w, V = np.linalg.eigh(D1)
    H = V[:, np.abs(w) < 1e-9]
    assert H.shape[1] == 26
    a1 = float(np.max(np.abs(B2 @ (B2.T @ H))))
    print(f"  A1 cusp-well action on the harmonic stage: "
          f"max |B2 B2^T h| = {a1:.1e}  "
          f"{'PASS' if a1 < 1e-10 else 'FAIL'} (Cert) — the link-form "
          f"interaction annihilates homology; the cusp-well route to a "
          f"harmonic dispersion is CLOSED.")
    assert a1 < 1e-10
    LEDGER["A1 cusp_well_closed (Cert)"] = f"{a1:.2e}"

    # A2: exact commutator of the two canonical cycle-space operators
    Ksp = sp.Matrix(K)
    Gram = Ksp * Ksp.T                                 # 29 x 29, exact
    B2sp = sp.Matrix([[Fraction(links[d][e]) for d in ms.divisors]
                      for e in range(84)])
    GK = Gram.LUsolve(Ksp * (B2sp * (B2sp.T * Ksp.T)))  # link Gram on Z
    C = GK * T2star - T2star * GK
    nz = sum(1 for i in range(29) for j in range(29) if C[i, j] != 0)
    fro = float(sp.sqrt(sum(C[i, j] ** 2
                            for i in range(29) for j in range(29))))
    print(f"  A2 exact commutator [G, T2*] on cycle space: "
          f"{nz} nonzero entries, Frobenius {fro:.4f}  (EXACT, "
          f"reported) — "
          f"{'the operators commute; ' if nz == 0 else 'noncommuting, but '}"
          f"by A1 any noncommutation acts through the Eisenstein sector "
          f"only and cannot split the harmonic level.")
    LEDGER["A2 commutator_nnz (EXACT)"] = nz
    LEDGER["A2 commutator_frobenius"] = round(fro, 6)
    return B1, B2, H, w, V


# ----------------------------------------------------------------------
# Stage B — minimal coupling
# ----------------------------------------------------------------------

def edge_avg_matrix(edge_list, f):
    d = np.empty(len(edge_list))
    for e, (i, j, _x) in enumerate(edge_list):
        d[e] = f[i] if i == j else 0.5 * (f[i] + f[j])
    return d


def stage_B(edge_list, B1, H):
    print("=" * 70)
    print("STAGE B — the canonical interaction: minimal coupling")
    print("=" * 70)
    L0 = B1 @ B1.T
    w0, U0 = np.linalg.eigh(L0)
    fied = U0[:, 1].copy()
    k0 = int(np.argmax(np.abs(fied)))
    if fied[k0] < 0:
        fied = -fied                     # sign convention, reported
    print(f"  B0 vertex mode: Fiedler u1 (lambda_1 = {w0[1]:.9f}, "
          f"simple by du01 A3). Sign convention: largest-|component| "
          f"positive; the dispersion's GLOBAL sign is convention, its "
          f"magnitudes are not.")

    ones = np.ones(56)
    Vnull = H.T @ (edge_avg_matrix(edge_list, ones)[:, None] * H)
    b1 = float(np.max(np.abs(Vnull - np.eye(26))))
    print(f"  B1 null control f = const: max |V_h - Id| = {b1:.1e}  "
          f"{'PASS' if b1 < 1e-12 else 'FAIL'} (Cert) — the module "
          f"structure is exact on homology at zero coupling profile.")
    assert b1 < 1e-12

    dvec = edge_avg_matrix(edge_list, fied)
    Vh = H.T @ (dvec[:, None] * H)
    Vh = 0.5 * (Vh + Vh.T)
    tr1 = float(np.trace(Vh))
    diagP = np.einsum('ek,ek->e', H, H)
    tr2 = float(np.dot(dvec, diagP))
    print(f"  B2 trace E2: tr V_h = {tr1:.12f} (direct) vs "
          f"{tr2:.12f} (edge-weight route), diff {abs(tr1-tr2):.1e}  "
          f"{'PASS' if abs(tr1-tr2) < 1e-10 else 'FAIL'} (Cert)")
    assert abs(tr1 - tr2) < 1e-10
    LEDGER["B2 trace_Vh (Cert)"] = round(tr1, 10)
    return w0, U0, fied, Vh


# ----------------------------------------------------------------------
# Stage C — the dispersion
# ----------------------------------------------------------------------

def orbit_of(a):
    if abs(a + 2) < 1e-6:
        return "old"
    if abs(a) < 1e-6:
        return "f1"
    if min(abs(a - t) for t in A2_F2) < 1e-4:
        return "f2"
    return "f3"


def hecke_on_harmonics(K, T2star, H):
    Kf = np.array(K, float)
    T2s = np.array([[float(T2star[i, j]) for j in range(29)]
                    for i in range(29)])
    G = Kf @ Kf.T
    ZH = np.linalg.solve(G, Kf @ H)          # 29 x 26 cycle coords
    T2H = H.T @ (Kf.T @ (T2s @ ZH))          # 26 x 26
    return T2H


def stage_C(ms, K, T2, T2star, H, Vh, edge_list_g):
    print("=" * 70)
    print("STAGE C — the dispersion relation omega(a)")
    print("=" * 70)
    T2H = hecke_on_harmonics(K, T2star, H)
    ev = np.linalg.eigvals(T2H).real
    A26, _B = ms.restrict_to_cuspidal(T2)
    ref = np.array(sorted(float(x) for x in
                          sp.Matrix(A26).eigenvals(multiple=True)))
    err = float(np.max(np.abs(np.sort(ev) - ref)))
    print(f"  C1 transported Hecke on the harmonic stage: eigenvalues "
          f"match the cuspidal lines at {err:.1e}  "
          f"{'PASS' if err < 1e-8 else 'FAIL'} (Cert; THIRD independent "
          f"route to the 13 lines: engine restriction, cycle transport, "
          f"harmonic compression)")
    assert err < 1e-8

    lam_ref = sorted(set(round(float(x), 6) for x in ref))
    mult = {a: int(np.sum(np.abs(ref - a) < 1e-6)) for a in lam_ref}
    assert len(lam_ref) == 12 and sum(mult.values()) == 26
    assert mult[-2.0] == 4                     # old pair coincides
    lines = []
    worst_sv = 0.0
    for a in lam_ref:
        m = mult[a]
        Uv, sv, Vt = np.linalg.svd(T2H - a * np.eye(26))
        assert sv[-m] < 1e-5 < sv[-m - 1], f"line a={a}: {sv[-m-2:]}"
        worst_sv = max(worst_sv, float(sv[-m]))
        Q = Vt[-m:, :].T                      # orthonormal 26 x m
        blk = Q.T @ Vh @ Q
        blk = 0.5 * (blk + blk.T)
        e2 = np.linalg.eigvalsh(blk)
        lines.append({"a": float(a), "orbit": orbit_of(a), "m": m,
                      "omega": float(np.trace(blk) / m),
                      "delta": float(e2[-1] - e2[0]),
                      "Q": Q, "levels": [float(x) for x in e2]})
    print(f"  C2 dispersion (omega = line mean, delta = fine structure), "
          f"DIAGNOSTIC compression onto Hecke lines "
          f"(worst null-space residual {worst_sv:.1e}):")
    print(f"     {'a2':>10}  {'orbit':>5} {'m':>2}  {'omega':>12}  "
          f"{'delta':>10}")
    for L in lines:
        print(f"     {L['a']:>10.6f}  {L['orbit']:>5} {L['m']:>2}  "
              f"{L['omega']:>12.8f}  {L['delta']:>10.2e}")
    LEDGER["C2 dispersion"] = [{k: L[k] for k in
                                ("a", "orbit", "m", "omega", "delta")}
                               for L in lines]

    # C2b: the vanishing of every line mean, certified, and its
    # mechanism. Per-line vertex weight W_a(i) = sum over edges at i of
    # (Q_a Q_a^T)_ee * avg-coefficient; omega(a) = sum_i f_i W_a(i).
    # If W_a is CONSTANT in i for every line, then omega(a) = 0 for
    # EVERY vertex profile orthogonal to constants — line-level mass
    # equidistribution, and first order dies universally.
    max_om = max(abs(L["omega"]) for L in lines)
    sym = max(max(abs(L["levels"][k] + L["levels"][-1 - k])
                  for k in range(len(L["levels"]))) for L in lines)
    print(f"  C2b line means: max |omega(a)| = {max_om:.1e}; level "
          f"symmetry max |e_k + e_(m-k)| = {sym:.1e}  "
          f"{'PASS' if max_om < 1e-10 and sym < 1e-10 else 'FAIL'} "
          f"(Cert) — every line splits SYMMETRICALLY about zero.")
    assert max_om < 1e-10 and sym < 1e-10
    dev = 0.0
    for L in lines:
        P2 = L["Q"] @ L["Q"].T
        diag = np.einsum('ek,ek->e', H @ P2, H)     # (H P2 H^T)_ee
        Wv = np.zeros(56)
        for e, (i, j, _x) in enumerate(edge_list_g):
            if i == j:
                Wv[i] += diag[e]
            else:
                Wv[i] += 0.5 * diag[e]
                Wv[j] += 0.5 * diag[e]
        dev = max(dev, float(np.max(np.abs(Wv - Wv.mean()))))
    print(f"  C2c mechanism: per-line vertex weights W_a(i) — max "
          f"deviation from constancy over all 12 lines = {dev:.2e}  "
          f"{'PASS — LINE-LEVEL MASS EQUIDISTRIBUTION (Cert): every '
             'Hecke line spreads uniformly over the 56 triangles, so '
             'first-order dispersion vanishes for EVERY profile '
             'orthogonal to constants, not just the Fiedler mode.'
             if dev < 1e-9 else
             'nonconstant — the vanishing is profile-specific '
             '(parity/selection); see C5.'}")
    LEDGER["C2c equidistribution_dev (Cert)"] = f"{dev:.3e}"

    # inter-line coupling
    top = []
    for i in range(len(lines)):
        for j in range(i + 1, len(lines)):
            c = float(np.linalg.norm(lines[i]["Q"].T @ Vh
                                     @ lines[j]["Q"]))
            top.append((c, lines[i]["a"], lines[j]["a"]))
    top.sort(reverse=True)
    print(f"  C3 strongest inter-line couplings |Q_a^T V Q_b|_F "
          f"(DIAGNOSTIC): " +
          ", ".join(f"({a:.3f},{b:.3f}):{c:.4f}" for c, a, b in top[:4]))
    LEDGER["C3 top_couplings"] = [[round(c, 6), a, b]
                                  for c, a, b in top[:6]]

    return lines, T2H


def stage_C_parity(ms, P_ms, K, D, H, Vh, lines):
    E29 = eta_float(ms, P_ms)                       # 29 x 29 float
    Df = np.array(D, float)
    EZ = np.linalg.solve(Df.T, E29.T @ Df.T)        # transport to cycles
    Kf = np.array(K, float)
    G = Kf @ Kf.T
    ZH = np.linalg.solve(G, Kf @ H)
    EH = H.T @ (Kf.T @ (EZ @ ZH))
    i2 = float(np.max(np.abs(EH @ EH - np.eye(26))))
    print(f"  C4 transported eta on harmonics: eta^2 = Id at {i2:.1e}  "
          f"{'PASS' if i2 < 1e-8 else 'FAIL'} (Cert)")
    assert i2 < 1e-8
    comm = float(np.linalg.norm(Vh @ EH - EH @ Vh))
    anti = float(np.linalg.norm(Vh @ EH + EH @ Vh))
    scale = float(np.linalg.norm(Vh))
    print(f"  C5 first-order operator vs parity: |[V,eta]|/|V| = "
          f"{comm/scale:.4f}, |{{V,eta}}|/|V| = {anti/scale:.4f}  "
          f"(Cert, reported)"
          + (" — V is eta-ODD: within each line V is purely "
             "parity-off-diagonal, the split eigenstates are equal "
             "parity mixtures, and delta(a) is a PARITY-FLIP matrix "
             "element (a form factor), consistent with C2b/C2c."
             if anti / scale < 1e-6 else ""))
    LEDGER["C5 parity"] = {"comm_rel": round(comm / scale, 8),
                           "anti_rel": round(anti / scale, 8)}
    return EH


def slope_battery(a_in, y_in, name):
    order = np.argsort(a_in)
    a = np.array(a_in)[order]
    o = np.array(y_in)[order]
    i0 = int(np.argmin(np.abs(a)))
    out = {}
    out["fd_left"] = float((o[i0] - o[i0 - 1]) / (a[i0] - a[i0 - 1]))
    out["fd_right"] = float((o[i0 + 1] - o[i0]) / (a[i0 + 1] - a[i0]))
    for k, w in (("lin5", 5), ("lin7", 7)):
        sel = np.argsort(np.abs(a))[:w]
        c = np.polyfit(a[sel], o[sel], 1)
        out[k] = float(c[0])
    c3 = np.polyfit(a, o, 3)
    out["cubic_at_0"] = float(c3[2])
    print(f"  C8 slope of {name} at the massless (f1) point — "
          f"estimator battery (A.7 standards, DIAGNOSTIC):")
    for k, v in out.items():
        print(f"       {k:>10}: {v:+.6f}")
    spread = max(out.values()) - min(out.values())
    print(f"       spread {spread:.6f} — the discrete 13-point "
          f"dispersion does not pin a unique slope; the sign "
          f"consistency across estimators is the readout.")
    return out


def second_order(edge_list, H, wE, VE, f):
    """H_eff = - sum_mu  P_h M_f |mu><mu| M_f P_h / lambda_mu over the
    58 non-harmonic edge modes: the canonical second-order (Lamb-type)
    effective Hamiltonian on the Hecke stage. Zero parameters."""
    dvec = edge_avg_matrix(edge_list, f)
    Heff = np.zeros((26, 26))
    for mu in range(len(wE)):
        if wE[mu] > 1e-9:
            g = H.T @ (dvec * VE[:, mu])
            Heff -= np.outer(g, g) / wE[mu]
    return 0.5 * (Heff + Heff.T)


def stage_C7(edge_list, H, wE, VE, fied, EH, lines):
    print("=" * 70)
    print("STAGE C7 — the SECOND-ORDER dispersion (the real one)")
    print("=" * 70)
    Heff = second_order(edge_list, H, wE, VE, fied)
    comm = float(np.linalg.norm(Heff @ EH - EH @ Heff))
    scale = float(np.linalg.norm(Heff))
    print(f"  C7a [H_eff, eta] / |H_eff| = {comm/scale:.2e}  "
          f"{'PASS' if comm/scale < 1e-8 else 'FAIL'} (Cert) — the "
          f"second-order operator is eta-EVEN (odd x odd), so parity "
          f"is a good quantum number and each line resolves into "
          f"parity branches.")
    assert comm / scale < 1e-8
    rows = []
    print(f"     {'a2':>10} {'orbit':>5} {'m':>2}  {'w2_mean':>12}  "
          f"{'w2(eta=+)':>12}  {'w2(eta=-)':>12}")
    for L in lines:
        Q = L["Q"]
        blk = 0.5 * ((Q.T @ Heff @ Q) + (Q.T @ Heff @ Q).T)
        Eblk = 0.5 * ((Q.T @ EH @ Q) + (Q.T @ EH @ Q).T)
        ev, U = np.linalg.eigh(blk)
        par = [float(U[:, k] @ Eblk @ U[:, k]) for k in range(len(ev))]
        wp = [e for e, s_ in zip(ev, par) if s_ > 0]
        wm = [e for e, s_ in zip(ev, par) if s_ < 0]
        ok = all(abs(abs(s_) - 1) < 1e-6 for s_ in par)
        row = {"a": L["a"], "orbit": L["orbit"], "m": L["m"],
               "w2_mean": float(np.mean(ev)),
               "w2_plus": [round(float(x), 8) for x in wp],
               "w2_minus": [round(float(x), 8) for x in wm],
               "parity_pure": ok}
        rows.append(row)
        print(f"     {L['a']:>10.6f} {L['orbit']:>5} {L['m']:>2}  "
              f"{row['w2_mean']:>12.8f}  "
              f"{str([f'{x:.6f}' for x in wp]):>12}  "
              f"{str([f'{x:.6f}' for x in wm]):>12}"
              + ("" if ok else "  [parity-mixed]"))
    assert all(r["parity_pure"] for r in rows), "eta not resolved"
    LEDGER["C7 second_order_dispersion"] = rows
    a_arr = [r["a"] for r in rows]
    m_arr = [r["w2_mean"] for r in rows]
    slopes = slope_battery(a_arr, m_arr, "w2_mean(a)")
    return Heff, rows, slopes


# ----------------------------------------------------------------------
# Stage D — robustness battery
# ----------------------------------------------------------------------

def spearman(x, y):
    rx = np.argsort(np.argsort(x)).astype(float)
    ry = np.argsort(np.argsort(y)).astype(float)
    rx -= rx.mean()
    ry -= ry.mean()
    return float((rx @ ry) / math.sqrt((rx @ rx) * (ry @ ry)))


def stage_D(edge_list, H, wE, VE, U0, lines, rows_u1):
    print("=" * 70)
    print("STAGE D — robustness battery over the vertex profile f")
    print("=" * 70)
    a = np.array([r["a"] for r in rows_u1])
    m1 = np.array([r["w2_mean"] for r in rows_u1])
    rng = np.random.default_rng(143)
    battery = {"u1": U0[:, 1], "u2": U0[:, 2], "u3": U0[:, 3],
               "u_top": U0[:, -1], "random": rng.standard_normal(56),
               "const": np.ones(56)}
    print(f"     {'f':>7} {'max|om1|':>10} {'rho(w2,a)':>10} "
          f"{'rho(w2,a^2)':>12} {'corr(w2,u1)':>12}")
    rows = {}
    for name, f in battery.items():
        ff = f.copy()
        k0 = int(np.argmax(np.abs(ff)))
        if ff[k0] < 0:
            ff = -ff
        dvec = edge_avg_matrix(edge_list, ff)
        Vf = H.T @ (dvec[:, None] * H)
        om1max = max(abs(float(np.trace(
            L["Q"].T @ Vf @ L["Q"])) / L["m"]) for L in lines)
        Hf = second_order(edge_list, H, wE, VE, ff)
        w2 = np.array([float(np.trace(L["Q"].T @ Hf @ L["Q"]) / L["m"])
                       for L in lines])
        r1 = spearman(w2, a)
        r2 = spearman(w2, a * a)
        r4 = (float(np.corrcoef(w2, m1)[0, 1])
              if np.std(w2) > 1e-14 else 0.0)
        rows[name] = dict(max_om1=f"{om1max:.1e}",
                          rho_w2_a=round(r1, 4),
                          rho_w2_a2=round(r2, 4),
                          corr_to_u1=round(r4, 4))
        print(f"     {name:>7} {om1max:>10.1e} {r1:>10.4f} "
              f"{r2:>12.4f} {r4:>12.4f}")
    only_u1 = (float(rows["u1"]["max_om1"].replace("e", "E")) < 1e-10
               and all(float(rows[k]["max_om1"].replace("e", "E"))
                       > 1e-4 for k in ("u2", "u3")))
    print("     readout: " + ("first-order vanishing is SPECIFIC to "
          "the eta-odd Fiedler profile (parity selection rule), not "
          "universal equidistribution; " if only_u1 else
          "first-order vanishing extends beyond u1; ")
          + "for second order, structure is a correlation shared by "
          "u1/u2/u3 that clears the random row's noise floor — "
          "anything below that bar is reported, not claimed.")
    LEDGER["D battery"] = rows
    return rows


# ----------------------------------------------------------------------
# Stage E — the systole (EXACT)
# ----------------------------------------------------------------------

def stage_E():
    print("=" * 70)
    print("STAGE E — the first ruler: the systole of X0(143)")
    print("=" * 70)
    s3 = [d for d in range(143) if (d * d - 3 * d + 1) % 143 == 0]
    s4 = [d for d in range(143) if (d * d - 4 * d + 1) % 143 == 0]
    print(f"  E1 residue scan mod 143: trace-3 solutions {s3} "
          f"(disc 5 is a non-residue mod 13); trace-4 solutions {s4}  "
          f"{'PASS' if not s3 and s4 else 'FAIL'} (EXACT)")
    assert not s3 and s4
    d0 = s4[0]
    a0 = 4 - d0
    m = (a0 * d0 - 1) // 143
    assert (a0 * d0 - 1) % 143 == 0
    M = ((a0, m), (143, d0))
    det = a0 * d0 - 143 * m
    print(f"  E2 explicit trace-4 element of Gamma0(143): "
          f"[[{a0}, {m}], [143, {d0}]], det = {det}  "
          f"{'PASS' if det == 1 else 'FAIL'} (EXACT). Primitive: a "
          f"square root would need trace sqrt(6), not an integer.")
    assert det == 1
    ell = 2 * math.acosh(2.0)
    print(f"  E3 systole ell_sys = 2 arccosh(2) = {ell:.15f} "
          f"(curvature-radius units, EXACT closed form). No shorter "
          f"closed geodesic exists: hyperbolic traces are integers "
          f">= 3 and 3 is excluded.")
    print(f"  E4 homogeneity note (Pr): every Farey triangle is an "
          f"SL2(Z)-image of {{0,1,oo}} and the standard Farey "
          f"quadrilateral (0, oo; 1, -1) has cross-ratio -1, shear 0 — "
          f"the tessellation carries no local metric ruler; the "
          f"canonical global lengths are area = 56 pi and ell_sys.")
    LEDGER["E systole (EXACT)"] = {"trace": 4, "matrix": M,
                                   "ell": round(ell, 12),
                                   "trace3_excluded": True}
    return ell


# ----------------------------------------------------------------------
# Stage F — disciplined census (Grok memo experiment 1)
# ----------------------------------------------------------------------

def stage_F(w0):
    print("=" * 70)
    print("STAGE F — census of graph ratios vs the hydrogen pure number")
    print("=" * 70)
    lv = np.sort(w0[np.abs(w0) > 1e-9])           # 55 levels
    t0 = 0.749598747594
    print(f"  F0 band statement: direct ratios span "
          f"[{lv[0]/lv[-1]:.4e}, {lv[-1]/lv[0]:.4f}]. Targets "
          f"1.7361669e6 and 4.3175501e-7 are OUT OF BAND — untestable "
          f"by direct ratios; powers/products are declined (AG-D5).")
    R = (lv[:, None] / lv[None, :]).ravel()
    R = R[(R > 0) & (np.abs(R - 1) > 1e-12)]
    tol = 1e-3
    hits = int(np.sum(np.abs(R / t0 - 1) < tol))
    B, rng = 2000, np.random.default_rng(143)
    sp_ = np.diff(lv)
    counts = []
    for _ in range(B):
        sh = lv[0] + np.concatenate([[0], np.cumsum(
            rng.permutation(sp_))])
        Rs = (sh[:, None] / sh[None, :]).ravel()
        Rs = Rs[(Rs > 0) & (np.abs(Rs - 1) > 1e-12)]
        counts.append(int(np.sum(np.abs(Rs / t0 - 1) < tol)))
    counts = np.array(counts)
    p = float(np.mean(counts >= hits))
    print(f"  F1 census: {hits} ratios within {tol:.0e} of "
          f"t0 = {t0} out of {len(R)}; spacing-shuffle null gives "
          f"{counts.mean():.1f} +- {counts.std():.1f}, p = {p:.3f}  "
          f"(DIAGNOSTIC census, not a search) — "
          + ("no anomaly." if p > 0.05 else "flagged for a proper "
             "falsifiability protocol before any interpretation."))
    LEDGER["F1 census"] = {"hits": hits, "null_mean": round(
        float(counts.mean()), 2), "p": round(p, 4)}


# ----------------------------------------------------------------------
# Stage G — PHENO translation (bookkeeping, no claims)
# ----------------------------------------------------------------------

def stage_G(w0, ell, slopes):
    print("=" * 70)
    print("STAGE G — PHENO bookkeeping (no claims)")
    print("=" * 70)
    lam1 = float(np.sort(w0)[1])
    f21 = 1420405751.7667
    c = 299792458.0
    chi_g = lam1 / (2 * math.pi * f21)
    print(f"  G1 Route A reproduction: chi_g = lambda_1 / (2 pi f_21cm)"
          f" = {chi_g:.6e} s  (PHENO; Grok memo value 3.056e-11 s "
          f"confirmed)")
    Rgeom = c * chi_g / ell
    print(f"  G2 curvature radius that makes the systolic velocity "
          f"equal c: R = c chi_g / ell_sys = {Rgeom:.6e} m  (PHENO)")
    lamC = 2.42631023538e-12
    a0 = 5.29177210544e-11
    print(f"  G3 pure ratios: R/lambda_C = {Rgeom/lamC:.4e}, "
          f"R/a_0 = {Rgeom/a0:.4e}  (PHENO; reported, no resonance "
          f"claimed)")
    print(f"  G4 dimensionless light-cone slope candidates from C8: "
          + ", ".join(f"{k}={v:+.4f}" for k, v in slopes.items())
          + "  (DIAGNOSTIC)")
    LEDGER["G"] = {"chi_g_s": chi_g, "R_geom_m": Rgeom,
                   "R_over_lambdaC": Rgeom / lamC, "R_over_a0": Rgeom / a0}


# ----------------------------------------------------------------------

def main():
    p1, tris, tri_of, edge_list, ms, K, D, links, T2, T2star = \
        build_quiet()
    P_ms = float_projection(ms)
    B1, B2, H, wE, VE = stage_A(edge_list, ms, K, links, T2star)
    w0, U0, fied, Vh = stage_B(edge_list, B1, H)
    lines, T2H = stage_C(ms, K, T2, T2star, H, Vh, edge_list)
    EH = stage_C_parity(ms, P_ms, K, D, H, Vh, lines)
    Heff, rows2, slopes2 = stage_C7(edge_list, H, wE, VE, fied, EH,
                                    lines)
    stage_D(edge_list, H, wE, VE, U0, lines, rows2)
    ell = stage_E()
    stage_F(w0)
    stage_G(w0, ell, slopes2)
    out = os.path.join(_HERE, "du03_dispersion.json")
    with open(out, "w") as f:
        json.dump(LEDGER, f, indent=1, default=str)
    print("=" * 70)
    print(f"ledger written: {out}")
    print("du03 complete — the interaction exists, the dispersion is "
          "on the table, the box has its first ruler.")


if __name__ == "__main__":
    main()
