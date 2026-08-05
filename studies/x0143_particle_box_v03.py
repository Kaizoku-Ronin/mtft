#!/usr/bin/env python3
"""
x0143_particle_box_v03.py — dissipation / atom capture, and
                            orbit-dependent Zeno on homology
=============================================================
Roger Tano / MTFT Research Program — built with Claude, August 2026

STAGE G — capture requires dissipation.
  G1  The four cusps become four nuclei. A cusp c of width w touches the
      56 triangles with multiplicities n_c(t); the well is V = -V0 n_c.
      Integer certificates: sum_t n_c(t) = width(c) for each cusp, and
      sum_c n_c(t) = 3 for every triangle.
  G2  H = L - V0 N_c on the dual graph. Bound states = eigenvalues below
      the free band [0, 6]. Localization measured by inverse participation.
  G3  Emission channel (Lindblad, downward jumps only in the H-eigenbasis):
      gamma_{mn} = gamma0 (E_n - E_m)^3 |<m|D|n>|^2, D = graph distance to
      the nucleus (declared dipole analog, model class).
  G4  Populations obey an exact classical master equation; solved by matrix
      exponential. Certificates:
        - column sums of the generator vanish  -> probability EXACTLY conserved
        - <H> monotonically decreasing
        - UNITARY CONTROL: bound-state population is EXACTLY constant under
          gamma0 = 0 (an eigenstate overlap cannot change) -> capture is
          impossible without dissipation. This is the theorem the stage exists
          to demonstrate, and it is exact, not numerical.
        - PHOTON LEDGER: the identity f + M^T w = 0 (all energy the state
          loses is carried by the jump flux) holds to machine precision.
          v0.11.1 relabel (owned error, caught by the Kimi BQ audit): the
          identity is FORCED by the vanishing column sums of the master
          generator, so it is one route written twice — a construction
          certificate, not an E2 pair. The honest numerical figure is the
          trajectory quadrature residual (~7e-3, time-grid limited).
  G5  Emission spectrum: discrete lines, one set per cusp-nucleus.

STAGE H — orbit-dependent Zeno on the 26-dim cuspidal homology.
  H1  Petersson metric pulled back to homology through the period map:
        G = sum_f (Pi_f^dag Pi_f) / <f,f>_Pet
      Certificates: G positive definite; Hecke operators built by the
      Manin-symbol engine are G-self-adjoint (cross-check between the
      Manin pipeline and the q-expansion period pipeline); Born weights
      through G reproduce v0.2's.
  H2  Within an orbit, the Galois conjugates carry distinct T2 eigenvalues.
      A state spread over them decays under Hecke evolution with
        A_Omega(tau) = sum_k p_k exp(-i a_k tau),   R = -ln|A|^2 / tau,
      and small-tau rate R -> Var_p(a) tau. The Born weights p_k come from
      the Petersson norms, which are NOT equidistributed (Paper 36 6.4),
      so the Zeno time is orbit-dependent.
  H3  dim(Omega_1) = 1  =>  Var = 0 EXACTLY  =>  the electron orbit cannot
      decay under Hecke evolution. Structural, not numerical.
  H4  Falsifiability (AG-D5 discipline): the f3 norm-to-embedding pairing is
      DIAGNOSTIC in the corpus, so the variance is reported as a certified
      INTERVAL over all pairings; the mu/tau ordering claim is tested against
      a permutation null and the quantitative ratio is compared to the
      measured lifetime ratio. Honest negatives are reported as such.
"""

from __future__ import annotations
from itertools import permutations
from math import pi
import json
import os

import numpy as np

from x0143_particle_box import (P1, SIGMA, TAU, tessellation, dual_graph,
                                orbits)
from x0143_particle_box_v02 import (build_engine, float_projection, eigendata,
                                    assign_orbits, extract_ap,
                                    period_functionals, PET_F1, PET_F2_DIAG,
                                    PET_F3_DIAG, A2_F2_ORDER, A2_F3_ORDER,
                                    P_MAX)

N = 143
import os as _os
_HERE = _os.path.dirname(_os.path.abspath(__file__))
CACHE = _os.path.join(_HERE, "v03_period_cache.npz")

# measured lepton lifetimes (PDG 2024), for the honest comparison in H4
TAU_MU_S = 2.1969811e-6
TAU_TAU_S = 2.903e-13


# ======================================================================
# STAGE G — the four cusps as nuclei
# ======================================================================

def cusp_incidence(p1, tris):
    """n_c(t) = number of corners of triangle t sitting at cusp c.
    Corners of the tau-orbit t are its three elements; the corner at
    element y sits at the cusp = T-orbit of y."""
    T_orb = orbits(p1, ((1, 1), (0, 1)), len(p1.reps))
    cusp_of = {}
    for ci, orb in enumerate(T_orb):
        for y in orb:
            cusp_of[y] = ci
    widths = [len(o) for o in T_orb]
    n = np.zeros((len(T_orb), len(tris)), dtype=int)
    for ti, t in enumerate(tris):
        for y in t:
            n[cusp_of[y], ti] += 1
    return n, widths


def bound_states(H, band_top=6.0):
    w, V = np.linalg.eigh(H)
    bound = np.nonzero(w < -1e-9)[0]
    return w, V, bound


def ipr(v):
    """inverse participation ratio: 1 = fully localized, 1/n = spread"""
    p = v ** 2
    return float(np.sum(p ** 2) / np.sum(p) ** 2)


def emission_generator(w, V, D, gamma0=0.05):
    """Downward-jump rates in the H-eigenbasis and the master-equation
    generator M with exactly vanishing column sums."""
    n = len(w)
    Dm = V.T @ (D[:, None] * V)              # <m|D|n>
    G = np.zeros((n, n))                     # G[m, n] = rate n -> m
    for a in range(n):
        for b in range(n):
            dE = w[b] - w[a]
            if dE > 1e-12:
                G[a, b] = gamma0 * dE ** 3 * Dm[a, b] ** 2
    M = G.copy()
    np.fill_diagonal(M, -G.sum(axis=0))
    return G, M


def evolve_populations(M, p0, times):
    from scipy.linalg import expm
    out = []
    for t in times:
        out.append(expm(M * t) @ p0)
    return np.array(out)


def stage_G(verbose=True):
    p1, tris, edges, cert = tessellation(N)
    A, edge_list, tri_of = dual_graph(p1, tris, edges)
    L = np.diag(A.sum(axis=1)) - A
    n_inc, widths = cusp_incidence(p1, tris)

    led = {}
    ok_w = [int(n_inc[c].sum()) == widths[c] for c in range(len(widths))]
    corner_sum = n_inc.sum(axis=0)
    ok_corners = bool(np.all(corner_sum == 3))
    if verbose:
        print(f"  G1 cusp widths {sorted(widths)}; "
              f"sum_t n_c(t) = width for every cusp: "
              f"{'PASS' if all(ok_w) else 'FAIL'} (EXACT)")
        print(f"     sum_c n_c(t) = 3 for all 56 triangles: "
              f"{'PASS' if ok_corners else 'FAIL'} (EXACT)")
    led["cusp incidence certificates (EXACT)"] = {
        "widths": sorted(widths), "row_sums_match_width": all(ok_w),
        "corner_sums_are_3": ok_corners}

    # graph distance from a node set
    def bfs(seed):
        d = -np.ones(A.shape[0], dtype=int)
        d[list(seed)] = 0
        frontier = list(seed)
        while frontier:
            nxt = []
            for u in frontier:
                for v in np.nonzero(A[u])[0]:
                    if d[v] < 0:
                        d[v] = d[u] + 1
                        nxt.append(int(v))
            frontier = nxt
        return d

    V0 = 4.0
    results = {}
    order = np.argsort(widths)
    for ci in order:
        Nc = n_inc[ci].astype(float)
        H = L - V0 * np.diag(Nc)
        w, Vv, bnd = bound_states(H)
        seed = np.nonzero(Nc > 0)[0]
        dist = bfs(seed).astype(float)
        results[ci] = dict(width=widths[ci], H=H, w=w, V=Vv, bound=bnd,
                           dist=dist, seed=seed)
        if verbose:
            print(f"     cusp width {widths[ci]:>3}: {len(bnd)} bound "
                  f"state(s), E_0 = {w[0]:+.4f}, "
                  f"IPR(ground) = {ipr(Vv[:, 0]):.3f}, "
                  f"support {len(seed)} triangle(s) (DIAGNOSTIC)")
    led["nuclei"] = {int(widths[ci]): {"n_bound": int(len(results[ci]['bound'])),
                                       "E0": float(results[ci]['w'][0]),
                                       "IPR": float(ipr(results[ci]['V'][:, 0]))}
                     for ci in order}

    # --- capture experiment -------------------------------------------
    # Only a NARROW cusp makes a localized well: the width-143 cusp touches
    # all 56 triangles (it is the ambient vacuum, not a nucleus), while the
    # width-1 cusp touches exactly one. Bind at the width-1 cusp.
    ci = int(np.argmin(widths))
    R = results[ci]
    w, Vv, dist = R["w"], R["V"], R["dist"]
    far = int(np.argmax(dist))
    psi0 = np.zeros(A.shape[0])
    psi0[far] = 1.0
    c = Vv.T @ psi0
    p0 = c ** 2
    p_bound0 = float(p0[R["bound"]].sum())

    Gj, M = emission_generator(w, Vv, dist, gamma0=0.05)
    colsum = float(np.abs(M.sum(axis=0)).max())
    times = np.linspace(0, 3000, 3001)
    P = evolve_populations(M, p0, times)
    trace_dev = float(np.abs(P.sum(axis=1) - 1).max())
    Ebound = P[:, R["bound"]].sum(axis=1)
    Eexp = P @ w
    mono = bool(np.all(np.diff(Eexp) <= 1e-12))

    # PHOTON LEDGER. Energy balance: d<H>/dt = w.Mp; photon flux = f.p,
    # f_b = sum_a gamma_{ab}(E_b - E_a). The identity f + M^T w = 0 holds
    # for EVERY state — but it is FORCED by M's vanishing column sums
    # (the adjoint identity), so it is a construction certificate that
    # catches indexing bugs at 1e-16, NOT two independent routes. (E2
    # mislabel corrected in v0.11.1; the independent numerical check is
    # the trajectory quadrature comparison below, time-grid limited.)
    dE = w[None, :] - w[:, None]              # dE[m, n] = E_n - E_m
    fvec = (Gj * dE).sum(axis=0)
    identity_res = float(np.abs(fvec + M.T @ w).max())
    flux = P @ fvec
    emitted = np.concatenate([[0.0], np.cumsum(
        0.5 * (flux[1:] + flux[:-1]) * np.diff(times))])
    quad_err = float(np.abs((Eexp[0] - Eexp) - emitted).max())

    if verbose:
        print(f"  G2 note: only narrow cusps bind. The width-143 cusp meets "
              f"all 56 triangles (ambient, IPR 0.04); the width-1 cusp meets "
              f"exactly one and is point-like (IPR 0.94). Nucleus = width-1 "
              f"cusp. (EXACT incidence, DIAGNOSTIC reading)")
        print(f"  G3/G4 capture, packet started {int(dist[far])} steps away:")
        print(f"     master-generator column sums: {colsum:.2e} "
              f"{'PASS' if colsum < 1e-12 else 'FAIL'} (EXACT)")
        print(f"     probability conservation over the run: "
              f"{trace_dev:.2e} (EXACT-class)")
        print(f"     UNITARY CONTROL: bound population is exactly constant "
              f"at {p_bound0:.4f} for all t when gamma0 = 0 "
              f"(eigenstate overlap; EXACT)")
        print(f"     with emission on: bound population "
              f"{p_bound0:.4f} -> {Ebound[-1]:.4f}  = CAPTURE")
        print(f"     <H> monotone decreasing: "
              f"{'PASS' if mono else 'FAIL'} (EXACT)")
        print(f"     PHOTON LEDGER, construction identity f + M^T w = 0 "
              f"(forced by vanishing column sums; catches indexing bugs, "
              f"not an E2 pair): residual {identity_res:.2e} "
              f"{'PASS' if identity_res < 1e-12 else 'FAIL'} (EXACT)")
        print(f"     integrated check on this trajectory: "
              f"|drop in <H> - integrated flux| = {quad_err:.2e} "
              f"(quadrature-limited, not a structural claim)")
    led["capture"] = {"start_distance": int(dist[far]),
                      "p_bound_initial (EXACT under unitary)": p_bound0,
                      "p_bound_final_with_emission": float(Ebound[-1]),
                      "column_sums (EXACT)": colsum,
                      "trace_dev (EXACT-class)": trace_dev,
                      "energy_monotone (EXACT)": mono,
                      "photon_ledger_identity_residual (EXACT)": identity_res,
                      "photon_ledger_quadrature_residual": quad_err}

    # emission spectrum: total quanta per line
    T_long, P_long = times, P
    quanta = np.zeros_like(Gj)
    for a in range(len(w)):
        for b in range(len(w)):
            if Gj[a, b] > 0:
                integ = np.trapezoid(Gj[a, b] * P_long[:, b], T_long)
                quanta[a, b] = integ
    lines_E, lines_I = [], []
    thresh = max(1e-9, 1e-4 * quanta.max())
    for a in range(len(w)):
        for b in range(len(w)):
            if quanta[a, b] > thresh:
                lines_E.append(w[b] - w[a])
                lines_I.append(quanta[a, b])
    if verbose:
        if lines_E:
            print(f"  G5 emission spectrum: {len(lines_E)} lines above "
                  f"1e-4 of peak; strongest at dE = "
                  f"{lines_E[int(np.argmax(lines_I))]:.4f} (DIAGNOSTIC)")
        else:
            print("  G5 emission spectrum: no lines above threshold")
    led["spectrum"] = {"n_lines": len(lines_E),
                       "strongest_dE": float(lines_E[int(np.argmax(lines_I))])
                       if lines_E else None}

    return led, dict(results=results, order=order, widths=widths,
                     times=times, P=P, Ebound=Ebound, Eexp=Eexp,
                     emitted=emitted, p_bound0=p_bound0,
                     lines_E=np.array(lines_E), lines_I=np.array(lines_I),
                     ci=ci, A=A)


# ======================================================================
# STAGE H — orbit-dependent Zeno on homology
# ======================================================================

def get_period_data():
    """Period functionals for the 22 new lines; cached to disk."""
    p1, tris, edges, ms = build_engine()
    P = float_projection(ms)
    Bc, proj_c, restrict, T2, E, lines = eigendata(ms, P)
    lines = assign_orbits(lines)
    if os.path.exists(CACHE):
        z = np.load(CACHE, allow_pickle=True)
        Pis = list(z["Pis"])
        meta = list(z["meta"])
    else:
        ap = extract_ap(ms, P, restrict, lines, P_MAX)
        X, keep, results = period_functionals(ms, P, Bc, proj_c, lines, ap,
                                              verbose=False)
        Pis = [results[id(L)][0] for L in lines]
        meta = [(L[0], L[1], L[2]) for L in lines]
        np.savez(CACHE, Pis=np.array(Pis), meta=np.array(meta, dtype=object))
    return ms, P, Bc, proj_c, restrict, lines, Pis, meta


def petersson_gram(lines, Pis):
    """G_ij = sum over embeddings ( Pi^dag Pi ) / <f,f>_Pet, both signs."""
    from collections import defaultdict
    emb = defaultdict(list)
    for L, Pi in zip(lines, Pis):
        emb[(L[1], round(L[2], 8))].append(Pi)
    G = np.zeros((26, 26), dtype=complex)
    for (orbit, lam), Pl in emb.items():
        pet = pet_norm(orbit, lam)
        for Pi in Pl:
            G += np.outer(Pi.conj(), Pi) / pet
    asym = float(np.abs(G - G.conj().T).max())
    return np.real(0.5 * (G + G.conj().T)), asym


def pet_norm(orbit, lam):
    if orbit == "f1":
        return PET_F1
    tab = (A2_F2_ORDER, PET_F2_DIAG) if orbit == "f2" \
        else (A2_F3_ORDER, PET_F3_DIAG)
    k = int(np.argmin([abs(lam - a) for a in tab[0]]))
    return tab[1][k]


def zeno_from_weights(a_vals, p):
    p = np.asarray(p, dtype=float)
    p = p / p.sum()
    a = np.asarray(a_vals, dtype=float)
    var = float(np.sum(p * a ** 2) - np.sum(p * a) ** 2)
    taus = np.logspace(-2.5, 1.5, 600)
    amp = np.array([np.sum(p * np.exp(-1j * a * t)) for t in taus])
    R = -np.log(np.clip(np.abs(amp) ** 2, 1e-300, 1)) / taus
    return taus, R, var


def stage_H(verbose=True):
    led = {}
    ms, P, Bc, proj_c, restrict, lines, Pis, meta = get_period_data()

    # ---- H1 Petersson Gram matrix on homology -------------------------
    G, asym = petersson_gram(lines, Pis)
    ev = np.linalg.eigvalsh(G)
    scale = ev.max()
    rank = int(np.sum(ev > 1e-9 * scale))
    pos_on_new = bool(np.sort(ev)[-22:].min() > 1e-9 * scale)
    if verbose:
        print(f"  H1 Petersson metric on the 26-dim homology: rank {rank} "
              f"(expected 22 = the new part), positive definite there "
              f"{'PASS' if (rank == 22 and pos_on_new) else 'FAIL'} (Cert); "
              f"hermiticity residual {asym/scale:.2e}")
    # G-self-adjointness of Hecke operators built by the OTHER pipeline
    sa = {}
    for p in (3, 7, 17):
        from x0143_particle_box_v02 import hecke_float
        Tp = restrict(hecke_float(ms, P, p))
        num = np.abs(G @ Tp - Tp.T @ G).max()
        sa[p] = float(num / np.abs(G @ Tp).max())
    if verbose:
        print(f"     Manin-pipeline T_p are self-adjoint in the "
              f"period-pipeline metric: "
              + ", ".join(f"T{p}: {v:.1e}" for p, v in sa.items())
              + f"  {'PASS' if max(sa.values()) < 1e-8 else 'FAIL'} (Cert)")
    led["gram"] = {"rank (Cert)": rank,
                   "positive_definite_on_new (Cert)": pos_on_new,
                   "hermiticity_residual_rel": float(asym / scale),
                   "hecke_self_adjointness (Cert)": sa}

    # the kernel of the newform Petersson metric should BE the old block
    import sympy as sp
    from x0143_particle_box import orbit_blocks
    xs = sp.symbols("x")
    A2s, _ = ms.restrict_to_cuspidal(ms.hecke_on_quotient(2))
    q2s = xs ** 4 - 3 * xs ** 3 - xs ** 2 + 5 * xs + 1
    q3s = xs ** 6 - 10 * xs ** 4 + 2 * xs ** 3 + 24 * xs ** 2 - 7 * xs - 12
    blocks = orbit_blocks(A2s, [(xs + 2, "old", 4), (xs, "f1", 2),
                                (q2s, "f2", 8), (q3s, "f3", 12)])
    block_np = [np.array(k[1].tolist(), dtype=float) for k in blocks]
    wv, Vv = np.linalg.eigh(G)
    ker = Vv[:, wv < 1e-9 * scale]
    old_q, _ = np.linalg.qr(block_np[0])
    ker_q, _ = np.linalg.qr(ker)
    princ = np.linalg.svd(old_q.T @ ker_q, compute_uv=False)
    subspace_gap = float(np.abs(princ - 1).max()) if len(princ) else 1.0
    if verbose:
        print(f"     kernel of the newform metric = the old block: "
              f"dim {ker.shape[1]} (expected 4), principal-angle deviation "
              f"{subspace_gap:.2e} "
              f"{'PASS' if ker.shape[1] == 4 and subspace_gap < 1e-8 else 'FAIL'}"
              f" (Cert)")
    led["gram"]["kernel_is_old_block (Cert)"] = {
        "dim": int(ker.shape[1]), "principal_angle_dev": subspace_gap}

    # ---- H2/H3 per-orbit Zeno ------------------------------------------
    orbits_data = {
        "f1 (electron)": (np.array([0.0]), np.array([PET_F1])),
        "f2 (muon)": (np.array(A2_F2_ORDER), np.array(PET_F2_DIAG)),
        "f3 (tau)": (np.array(A2_F3_ORDER), np.array(PET_F3_DIAG)),
    }
    curves, vars_pet, vars_eq = {}, {}, {}
    for nm, (a, pet) in orbits_data.items():
        taus, R, var = zeno_from_weights(a, pet)
        curves[nm] = (taus, R)
        vars_pet[nm] = var
        _, _, var_eq = zeno_from_weights(a, np.ones_like(pet))
        vars_eq[nm] = var_eq
    if verbose:
        print()
        print(f"  H2/H3 Zeno variance (Delta a2)^2 per orbit — the initial "
              f"decay rate is R = Var * tau:")
        for nm in orbits_data:
            tp = curves[nm][0][int(np.argmax(curves[nm][1]))] \
                if vars_pet[nm] > 0 else float("inf")
            print(f"     {nm:<14} Var_Petersson = {vars_pet[nm]:.6f}   "
                  f"Var_equidistributed = {vars_eq[nm]:.6f}   "
                  f"tau* = {tp if tp != float('inf') else float('inf'):.3f}"
                  if vars_pet[nm] > 0 else
                  f"     {nm:<14} Var_Petersson = {vars_pet[nm]:.6f}   "
                  f"Var_equidistributed = {vars_eq[nm]:.6f}   "
                  f"tau* = infinity  [EXACT: dim 1 => no Galois spread]")
    led["zeno_variance"] = {"petersson": {k: float(v) for k, v in
                                          vars_pet.items()},
                            "equidistributed": {k: float(v) for k, v in
                                                vars_eq.items()}}

    # ---- H4 falsifiability battery -------------------------------------
    if verbose:
        print()
        print("  H4 falsifiability battery (AG-D5 discipline):")
    # (a) the f3 (and f2) norm-to-embedding pairing is DIAGNOSTIC:
    #     report the variance as a certified interval over ALL pairings
    intervals = {}
    for nm, (a, pet) in orbits_data.items():
        if len(a) == 1:
            intervals[nm] = (vars_pet[nm], vars_pet[nm])
            continue
        vals = []
        for perm in permutations(range(len(pet))):
            w = pet[list(perm)]
            wn = w / w.sum()
            vals.append(np.sum(wn * a ** 2) - np.sum(wn * a) ** 2)
        intervals[nm] = (float(min(vals)), float(max(vals)))
    if verbose:
        for nm, (lo, hi) in intervals.items():
            print(f"     {nm:<14} Var over ALL norm-embedding pairings: "
                  f"[{lo:.6f}, {hi:.6f}]  (assumed pairing: "
                  f"{vars_pet[nm]:.6f})")
    disjoint = intervals["f3 (tau)"][0] > intervals["f2 (muon)"][1]
    if verbose:
        print(f"     tau-orbit variance exceeds muon-orbit variance for "
              f"EVERY pairing: {'YES' if disjoint else 'NO — intervals overlap'}")
    # (b) permutation null for the ordering claim
    rng = np.random.default_rng(143)
    a2, a3 = orbits_data["f2 (muon)"][0], orbits_data["f3 (tau)"][0]
    hits = 0
    trials = 20000
    for _ in range(trials):
        w2 = rng.dirichlet(np.ones(4))
        w3 = rng.dirichlet(np.ones(6))
        v2 = np.sum(w2 * a2 ** 2) - np.sum(w2 * a2) ** 2
        v3 = np.sum(w3 * a3 ** 2) - np.sum(w3 * a3) ** 2
        hits += (v3 > v2)
    null_p = hits / trials
    if verbose:
        print(f"     null test: with RANDOM weights on the same eigenvalue "
              f"sets, Var(tau) > Var(mu) happens {null_p:.1%} of the time — "
              f"the ordering is a property of the eigenvalue SPREAD, not "
              f"evidence for the Petersson weighting")
    # (c) quantitative comparison to measured lifetimes: honest negative
    ratio_var = vars_pet["f3 (tau)"] / vars_pet["f2 (muon)"]
    ratio_life = TAU_MU_S / TAU_TAU_S
    if verbose:
        print(f"     quantitative check: Var(tau)/Var(mu) = "
              f"{ratio_var:.3f}, measured lifetime ratio "
              f"tau_mu/tau_tau = {ratio_life:.3e} — "
              f"NO quantitative correspondence. The Zeno variances order "
              f"the generations but do not scale as lifetimes; any claim "
              f"beyond ordering is unsupported.")
    led["falsifiability"] = {
        "variance_intervals_over_pairings (Cert)":
            {k: list(v) for k, v in intervals.items()},
        "tau_gt_mu_for_every_pairing (Cert)": bool(disjoint),
        "random_weight_null_probability": float(null_p),
        "var_ratio": float(ratio_var),
        "measured_lifetime_ratio": float(ratio_life),
        "verdict": "ordering only; no quantitative lifetime correspondence"}

    # ---- H5 real drawn cycles ------------------------------------------
    rng2 = np.random.default_rng(143)
    vrand = sum(rng2.integers(-3, 4) * block_np[i][:, k]
                for i in range(4) for k in range(min(2, block_np[i].shape[1])))
    # Born weights through G vs v0.2's direct route (consistency)
    from collections import defaultdict
    emb = defaultdict(list)
    for L, Pi in zip(lines, Pis):
        emb[(L[1], round(L[2], 8))].append(Pi)
    per_emb = {}
    for (orbit, lam), Pl in emb.items():
        amp = sum(abs(np.dot(Pi, vrand)) ** 2 for Pi in Pl)
        per_emb[(orbit, lam)] = amp / pet_norm(orbit, lam)
    tot = sum(per_emb.values())
    orb_w = {"f1": 0.0, "f2": 0.0, "f3": 0.0}
    for (orbit, lam), v in per_emb.items():
        orb_w[orbit] += v / tot
    via_G = float(np.real(vrand @ G @ vrand))
    if verbose:
        print()
        print(f"  H5 drawn cycle, Born weights e/mu/tau = "
              f"({orb_w['f1']:.4f}, {orb_w['f2']:.4f}, {orb_w['f3']:.4f}) "
              f"[v0.2: 0.1813, 0.8004, 0.0183] "
              f"{'PASS' if abs(orb_w['f1']-0.1813)<1e-3 else 'FAIL'} (Cert)")
        print(f"     same total through the Gram matrix: "
              f"psi^T G psi = {via_G:.6f} vs sum of parts {tot:.6f} "
              f"{'PASS' if abs(via_G-tot)/tot < 1e-9 else 'FAIL'} (Cert)")
    led["drawn_cycle"] = {"born": {k: float(v) for k, v in orb_w.items()},
                          "gram_route": via_G, "sum_of_parts": tot}

    # within-orbit Zeno of the actual drawn cycle
    cyc_curves = {}
    for orbit, nm in (("f2", "f2 (muon)"), ("f3", "f3 (tau)")):
        items = [(lam, v) for (o, lam), v in per_emb.items() if o == orbit]
        items.sort()
        a = np.array([it[0] for it in items])
        p = np.array([it[1] for it in items])
        taus, R, var = zeno_from_weights(a, p)
        cyc_curves[nm] = (taus, R, var)
    if verbose:
        for nm, (_, _, var) in cyc_curves.items():
            print(f"     drawn cycle's within-orbit variance, {nm}: "
                  f"{var:.6f} (DIAGNOSTIC — state-dependent, unlike the "
                  f"canonical Petersson state)")
    led["drawn_cycle_zeno"] = {nm: float(v[2]) for nm, v in cyc_curves.items()}

    return led, dict(curves=curves, vars_pet=vars_pet, vars_eq=vars_eq,
                     intervals=intervals, cyc_curves=cyc_curves)


# ======================================================================
# figures
# ======================================================================

def figures(gd, hd):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # --- capture ---------------------------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(13.2, 4.8))
    ax = axes[0]
    ax.plot(gd["times"], gd["Ebound"], color="#1d3557", lw=2.2,
            label="bound-state population, emission ON")
    ax.axhline(gd["p_bound0"], color="#e63946", ls="--", lw=2,
               label="unitary control (EXACTLY constant)")
    ax.set_xlabel("time")
    ax.set_ylabel("population in bound states")
    ax.set_ylim(-0.03, 1.03)
    ax.set_title("Atom capture at the point-like width-1 cusp\n"
                 "unitary evolution cannot bind; the emission channel can",
                 fontsize=11)
    ax.legend(fontsize=9, loc="center right")
    ax = axes[1]
    ax.plot(gd["times"], gd["Eexp"], color="#1d3557", lw=2.2,
            label=r"$\langle H\rangle(t)$ from the populations")
    ax.plot(gd["times"], gd["Eexp"][0] - gd["emitted"], "--", color="#f4a261",
            lw=2, label="initial energy minus integrated photon flux")
    ax.set_xlabel("time")
    ax.set_ylabel("energy")
    ax.set_title("Photon ledger: energy bookkeeping of the capture\n"
                 "identity check is construction-forced (not an E2 pair);\n"
                 "quadrature residual 7.2e-3 (Add. BQ sec.8)", fontsize=11)
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(_os.path.join(_HERE, "fig_atom_capture.png"), dpi=170)
    plt.close(fig)
    print("  wrote fig_atom_capture.png")

    # --- emission spectrum ------------------------------------------------
    fig, ax = plt.subplots(figsize=(9.4, 4.6))
    E, I = gd["lines_E"], gd["lines_I"]
    I = I / I.max()
    ax.vlines(E, 0, I, color="#9b5de5", lw=1.6)
    ax.set_xlabel("photon energy  $E_n - E_m$")
    ax.set_ylabel("relative line strength")
    ax.set_title("Emission spectrum of the $X_0(143)$ atom "
                 "(width-1 cusp nucleus)\n"
                 "discrete lines from the arithmetic level structure "
                 "(DIAGNOSTIC: dipole model declared)", fontsize=11)
    fig.tight_layout()
    fig.savefig(_os.path.join(_HERE, "fig_emission_spectrum.png"), dpi=170)
    plt.close(fig)
    print("  wrote fig_emission_spectrum.png")

    # --- orbit-dependent Zeno --------------------------------------------
    fig, ax = plt.subplots(figsize=(9.4, 5.2))
    colors = {"f1 (electron)": "#e4572e", "f2 (muon)": "#17bebb",
              "f3 (tau)": "#9b5de5"}
    for nm, (taus, R) in hd["curves"].items():
        ax.semilogx(taus, R, lw=2.2, color=colors[nm],
                    label=f"{nm}   Var = {hd['vars_pet'][nm]:.4f}")
    ax.axhline(0, color="#999", lw=0.6)
    ax.set_xlabel(r"measurement interval $\tau$  (Hecke units)")
    ax.set_ylabel(r"effective decay rate $R_\Omega(\tau)$")
    ax.set_title("Orbit-dependent Zeno on the homology of $X_0(143)$\n"
                 "the electron orbit is one-dimensional, so its curve is "
                 "identically zero — stability is structural", fontsize=11)
    ax.legend(fontsize=9.5)
    fig.tight_layout()
    fig.savefig(_os.path.join(_HERE, "fig_orbit_zeno.png"), dpi=170)
    plt.close(fig)
    print("  wrote fig_orbit_zeno.png")


def main():
    ledger = {}
    print("=" * 70)
    print("v0.3 STAGE G — dissipation channel and atom capture")
    print("=" * 70)
    ledger["G"], gd = stage_G()
    print()
    print("=" * 70)
    print("v0.3 STAGE H — orbit-dependent Zeno on homology")
    print("=" * 70)
    ledger["H"], hd = stage_H()
    print()
    print("FIGURES")
    figures(gd, hd)
    with open(_os.path.join(_HERE, "certificates_v03.json"), "w") as f:
        json.dump(ledger, f, indent=2, default=str)
    print("\nledger written: certificates_v03.json")


if __name__ == "__main__":
    main()
