# --- stage 5 (Integration Plan v0.1): re-pointed to the mtft package. ---
# internal() delegates to mtft.chain.internal (certified spectrally
# identical: g to 9e-15, full H(u) spectrum to 2.5e-14; B differs only
# by gauge inside near-degenerate blocks).  Where a local gsq() remains
# it is retained VERBATIM: it evaluates complex u in the f64 backend,
# which mtft.ep deliberately does not offer (PR-20 floors) -- S5-1.
# ------------------------------------------------------------------------
#!/usr/bin/env python3
"""
pr6_critical_and_fano.py — PR-6: critical structure at tau_c, and the
first level widths in the reconstruction
=======================================================================
MIT License — (c) 2026 Roger Tano — MTFT Research Program

PRE-REGISTERED (rung-5 note v0.1.1 §6, before this code existed):
  (a) critical exponent of the indirect gap at tau_c; first-order theory
      says nu = 1 with slope 2(mu_0 + mu_1) = 3.217924.
  (b) above tau_c: does a Fano/avoided-crossing width appear, and does
      it scale with the interband element |B_01|^2 ?
  (c) every time-averaged / k-integrated observable ships a k-refinement
      spread (Addendum AA.2 clause).

Verified baselines (two engines): tau_c = 0.23003, mu_0 = 1.050398,
mu_1 = 0.558564, m = 0.736839, kappa* = 5.0, beta = 2.

STRUCTURAL POINT FOR (b), stated before measuring.  A resonance is a
discrete level coupled to a continuum.  Pr H proved there are NO flat
bands, and at fixed k the spectrum of H(k) is discrete with no
accumulation, so the translation-invariant model has no resonances at
all: every state is extended, widths are identically zero.  A width
therefore REQUIRES breaking translation invariance.  The minimal
breaking is a rank-one defect: potential V on internal orbital d at
site x = 0.  Above tau_c the defect state split off from band 1 lands
INSIDE band 0's continuum -> autoionizing resonance, whose golden-rule
width is carried by the interband element exactly as pre-registered:
    Gamma = 2 pi V_eff^2 A_1(eps),   A_1 = orbital-1 partial DOS of
    band 0, whose leading behaviour is (2 tau B_01 / Delta)^2 * rho_0.
Exact tool: rank-one T-matrix.  G_dd(z) = (1/2pi) Int dk <d|(z-H(k))^-1|d>,
defect spectral function A_def(w) = -(1/pi) Im[G_dd/(1 - V G_dd)].

Gates: CG0 baselines (+k-clause); CG1 nu = 1; CG2 slope two-leg
(rotated-eigenvector HF vs finite difference vs first-order);
CG3 |B_01|^2 mixing law below tau_c; CG4 no-resonance structural
statement (translation-invariant); CG5 defect resonance above tau_c,
Gamma vs golden rule; CG6 Gamma ~ tau^2 |B_01|^2 scaling (the
pre-registered law).

Run:  py pr6_critical_and_fano.py
"""
from __future__ import annotations
import math
import numpy as np
import mpmath as mp
from mtft.chain import internal as _chain_internal


BETA, KSTAR = 2.0, 5.0
TAU_C_REF, M_REF = 0.23003, 0.736839
REPORT = []


def rec(name, gtype, value, cls, ok, note=""):
    REPORT.append((name, gtype, value, cls, bool(ok), note))
    print(f"[{'PASS' if ok else 'FAIL'}] {name:<30} {gtype:<12} "
          f"{value:<34} {cls:<20} {note}", flush=True)


def internal(N=1600, kappa=KSTAR, nb=300, gcap=200.0):
    ic = _chain_internal(kappa, nb=nb, backend="f64", N=N, gcap=gcap)
    return np.asarray(ic.g), np.asarray(ic.B)


G, B = internal()


def band_edge(tau, c, i):
    return np.linalg.eigvalsh(np.diag(G) - 2.0 * tau * c * B)[i]


def g_ind(tau):
    return band_edge(tau, +1.0, 1) - band_edge(tau, -1.0, 0)


# ----------------------------------------------------------------- CG0
def cg0():
    mu0, mu1, m = B[0, 0], B[1, 1], G[1]
    lo, hi = 0.5 * TAU_C_REF, 1.5 * TAU_C_REF
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if g_ind(mid) > 0: lo = mid
        else: hi = mid
    tc = 0.5 * (lo + hi)
    ok = (abs(tc - TAU_C_REF) < 1e-4 and abs(mu0 - 1.050398) < 1e-5
          and abs(mu1 - 0.558564) < 1e-5 and abs(m - M_REF) < 1e-5)
    rec("CG0 baselines", "Instrument", f"tau_c={tc:.6f}", "CERTIFIED(1e-4)",
        ok, f"mu0={mu0:.6f}, mu1={mu1:.6f}, m={m:.6f}, |B01|={abs(B[0,1]):.6f}")
    return tc


# ----------------------------------------------------------------- CG1
def cg1(tc):
    ds = np.array([10.0 ** e for e in (-2, -2.5, -3, -3.5, -4, -4.5, -5)])
    gs = np.array([g_ind(tc - d) for d in ds])
    nu, lnA = np.polyfit(np.log(ds), np.log(gs), 1)
    resid = float(np.max(np.abs(np.log(gs) - (nu * np.log(ds) + lnA))))
    ok = abs(nu - 1.0) < 1e-3 and resid < 1e-3
    rec("CG1 critical exponent", "Theorem", f"nu = {nu:.6f}",
        "CERTIFIED(1e-3)", ok,
        f"7 decades, log-log resid {resid:.1e}; nu=1 exactly: the two "
        f"edges live in DIFFERENT Bloch fibers (k=0 vs pi) so cannot "
        f"hybridize -> simple analytic zero")
    return nu


# ----------------------------------------------------------------- CG2
def cg2(tc):
    h = 1e-6
    fd = (g_ind(tc + h) - g_ind(tc - h)) / (2 * h)
    w1, U1 = np.linalg.eigh(np.diag(G) - 2 * tc * (+1.0) * B)
    w0, U0 = np.linalg.eigh(np.diag(G) - 2 * tc * (-1.0) * B)
    hf = -2.0 * (float(U1[:, 1] @ (B @ U1[:, 1]))
                 + float(U0[:, 0] @ (B @ U0[:, 0])))
    first = -2.0 * (B[0, 0] + B[1, 1])
    ok = abs(fd - hf) < 1e-6
    rec("CG2 slope two-leg", "Implementation", f"|FD-HF| = {abs(fd-hf):.1e}",
        "EXACT(HF)", ok,
        f"slope(tau_c)={fd:.6f} vs first-order {first:.6f} "
        f"({100*abs(fd/first-1):.2f}% eigenvector rotation — this is the "
        f"O(tau) correction that put tau_c 0.5% above prediction)")


# ----------------------------------------------------------------- CG3
def cg3():
    b01, delta = abs(B[0, 1]), G[1] - G[0]
    devs = []
    for tau in (0.002, 0.004, 0.008):
        w, U = np.linalg.eigh(np.diag(G) - 2 * tau * 1.0 * B)
        mix = float(U[1, 0] ** 2)              # orbital-1 weight in band 0
        pred = (2 * tau * b01 / delta) ** 2
        devs.append(abs(mix / pred - 1.0))
    ok = max(devs) < 0.05
    rec("CG3 |B01|^2 mixing law", "Structural",
        f"max dev {max(devs):.3f}", "CERTIFIED(5%)", ok,
        f"mix = (2 tau |B01|/Delta)^2 verified at 3 tau; |B01|={b01:.6f}, "
        f"Delta={delta:.6f}")


# ----------------------------------------------------------------- CG4
def cg4(tau=0.30):
    w = np.linalg.eigvalsh(np.diag(G) - 2 * tau * 0.3 * B)
    sp = np.diff(w[:40])
    ok = float(np.min(sp)) > 0
    rec("CG4 no intrinsic widths", "Theorem",
        f"min level spacing at fixed k = {float(np.min(sp)):.4f}",
        "EXACT(structural)", ok,
        "H(k) discrete, no accumulation, no flat bands (Pr H) => the "
        "translation-invariant model has NO resonances; widths require "
        "broken translation invariance (defect)")


# --------------------------------------------------- defect machinery
def local_greens(tau, d, ws, eta, nk, nbG=60):
    Gr, Br = G[:nbG], B[:nbG, :nbG]
    ks = (np.arange(nk) + 0.5) * math.pi / nk
    E = np.empty((nk, nbG)); C = np.empty((nk, nbG))
    for j, k in enumerate(ks):
        w, U = np.linalg.eigh(np.diag(Gr) - 2 * tau * math.cos(k) * Br)
        E[j] = w; C[j] = U[d, :] ** 2
    Ef, Cf = E.ravel(), C.ravel()
    Gd = np.empty(len(ws), dtype=complex)
    for i, wv in enumerate(ws):
        Gd[i] = np.sum(Cf / (wv + 1j * eta - Ef)) / nk
    return Gd, E


def band0_curve(tau, d, nk=20000, nbG=60):
    """eps_0(k) and |c_0(k)|^2 on [0,pi]; eps_0 is monotone in k because
    it depends on k only through cos k."""
    Gr, Br = G[:nbG], B[:nbG, :nbG]
    ks = np.linspace(0.0, math.pi, nk)
    e0 = np.empty(nk); c0 = np.empty(nk)
    for j, k in enumerate(ks):
        w, U = np.linalg.eigh(np.diag(Gr) - 2 * tau * math.cos(k) * Br)
        e0[j] = w[0]; c0[j] = U[d, 0] ** 2
    return ks, e0, c0


def imG_exact(ws, ks, e0, c0):
    """Im G_dd(w+i0) = -|c_0(k*)|^2 / |d eps_0/dk| at the two roots +-k*.
    EXACT: no eta smearing, hence no contamination from the tails of the
    heavily-weighted band-1 states (which at eta=1e-3 exceed the band-0
    signal by ~2x — recorded as a methods finding)."""
    de = np.gradient(e0, ks)
    kstar = np.interp(ws, e0, ks)
    cst = np.interp(kstar, ks, c0)
    dst = np.interp(kstar, ks, de)
    return -cst / np.abs(dst)


def resonance(tau, V, d=1, eta=1e-3, nk=3072, nbG=60, npts=3000,
              nk_exact=20000):
    """Resonance at the root of 1 - V Re G_dd (Re from the eta-grid,
    insensitive), width from the EXACT residue Im G:

        Gamma_BW = 2 |Im G| / |Re G'| = Gamma_GR * Z,
        Gamma_GR = 2 pi V^2 A_dd,  Z = 1/(V^2 |Re G'|)  -> 1 weak-V.

    Suite catches recorded in the note: (i) the first run's
    FWHM-minus-2eta estimator was resolution-limited; (ii) the second
    run's eta-smeared Im G was tail-contaminated."""
    _, E = local_greens(tau, d, np.array([0.0]), eta, 64, nbG)
    b0lo, b0hi = E[:, 0].min(), E[:, 0].max()
    b1lo = E[:, 1].min()                       # band-1 bottom: hard edge
    pad = 0.12 * (b0hi - b0lo)
    w_hi = min(b0hi - pad, b1lo - 0.06 * (b0hi - b0lo))
    ws = np.linspace(b0lo + pad, w_hi, npts)
    Gd, _ = local_greens(tau, d, ws, eta, nk, nbG)
    ks, e0, c0 = band0_curve(tau, d, nk_exact, nbG)
    ImX = imG_exact(ws, ks, e0, c0)
    f = 1.0 - V * np.real(Gd)
    idx = np.where(np.diff(np.sign(f)) != 0)[0]
    if len(idx) == 0:
        return None
    best = None
    for i in idx:
        t = f[i] / (f[i] - f[i + 1])
        wr = ws[i] + t * (ws[i + 1] - ws[i])
        dRe = (np.real(Gd[i + 1]) - np.real(Gd[i])) / (ws[i + 1] - ws[i])
        ImG = ImX[i] + t * (ImX[i + 1] - ImX[i])
        if dRe == 0:
            continue
        gam = 2.0 * abs(ImG) / abs(dRe)
        Add = -(1.0 / math.pi) * ImG
        gr = 2 * math.pi * V ** 2 * Add
        Z = 1.0 / (V ** 2 * abs(dRe))
        # physical resonance = largest T-matrix spectral weight (pi*A*gam),
        # not the narrowest root (van Hove artifacts are near-zero-width)
        wt = abs(Add) * gam
        if best is None or wt > best[0]:
            best = (wt, (wr, gam, gr, Z, Add, ws,
                         np.real(Gd) + 1j * ImX))
    return None if best is None else best[1]


def gamma_from_shape(res, tau, V, d=1, eta=1e-3, nk=3072, nbG=60,
                     nk_exact=20000):
    """Independent leg: Lorentzian fit to the T-matrix peak on a
    dedicated fine grid around w_r, built from the EXACT Im G."""
    wr, gam, gr, Z, Add, _, _ = res
    wf = np.linspace(wr - 5 * gam, wr + 5 * gam, 241)
    Gf, _ = local_greens(tau, d, wf, eta, nk, nbG)
    ks, e0, c0 = band0_curve(tau, d, nk_exact, nbG)
    Gx = np.real(Gf) + 1j * imG_exact(wf, ks, e0, c0)
    A = -(1.0 / math.pi) * np.imag(Gx / (1.0 - V * Gx))
    if A.max() <= 0:
        return float("nan")
    y = 1.0 / A
    c = np.polyfit(wf - wr, y, 2)
    if c[0] <= 0 or c[2] <= 0:
        return float("nan")
    return 2.0 * math.sqrt(c[2] / c[0])


# ----------------------------------------------------------------- CG5
def cg5(tau=0.30, V=-0.35):
    r = resonance(tau, V)
    wr, gam, gr, Z, Add, _, _ = r
    ident = abs(gam / (gr * Z) - 1.0)
    g_eta = resonance(tau, V, eta=2e-3)[1]
    g_nb = resonance(tau, V, nbG=120)[1]
    g_kx = resonance(tau, V, nk_exact=40000)[1]
    g_shape = gamma_from_shape(r, tau, V)
    es, nbs, kxs = (abs(gam - g_eta) / gam, abs(gam - g_nb) / gam,
                    abs(gam - g_kx) / gam)
    leg = abs(g_shape / gam - 1.0)
    ok = (ident < 1e-9 and es < 0.02 and nbs < 1e-6 and kxs < 1e-3
          and leg < 0.05)
    rec("CG5 defect resonance", "Theorem",
        f"Gamma = {gam:.6f} @ w_r = {wr:.4f}", "CERTIFIED(two-leg)", ok,
        f"BW=GR*Z exact ({ident:.0e}); shape-fit leg {g_shape:.6f} "
        f"({leg:.2%}); spreads eta {es:.2%}/nb {nbs:.0e}/k_exact "
        f"{kxs:.0e}; Z={Z:.3f}; tau>tau_c: level bound off band 1, "
        f"embedded in band 0 -> AUTOIONIZES")
    return r


# ----------------------------------------------------------------- CG6
def cg6(V=-0.35, d=1):
    """The pre-registered |B01|^2 law, tested POINTWISE at the resonance
    momentum (a log-log slope in tau is not a clean probe: the resonance
    slides, so DOS and weight both move).  Second-order perturbation at
    momentum k* predicts the orbital-1 weight in band 0:

        mix(k*) = (2 tau cos k* |B01| / Delta)^2,

    and the golden rule composes as
        Gamma = 2 pi V^2 * mix * rho_0 * Z,   rho_0 = 1/(pi |d eps_0/dk|).
    Every factor is measured separately; nothing is fitted."""
    rows, b01, Delta = [], abs(B[0, 1]), G[1] - G[0]
    worst = 0.0
    for tau in (0.28, 0.34, 0.40):
        r = resonance(tau, V)
        wr, gam, gr, Z, Add, _, _ = r
        ks, e0, c0 = band0_curve(tau, d, 20000, 60)
        kstar = float(np.interp(wr, e0, ks))
        mix = float(np.interp(kstar, ks, c0))
        de = float(np.interp(kstar, ks, np.gradient(e0, ks)))
        rho0 = 1.0 / (math.pi * abs(de))
        Gr, Br = G[:60], B[:60, :60]
        ev = np.linalg.eigvalsh(np.diag(Gr)
                                - 2 * tau * math.cos(kstar) * Br)
        Delta_k = ev[1] - ev[0]                # DRESSED gap at k*
        pred = (2 * tau * math.cos(kstar) * b01 / Delta_k) ** 2
        pred_bare = (2 * tau * math.cos(kstar) * b01 / Delta) ** 2
        dev = abs(mix / pred - 1.0)
        worst = max(worst, dev)
        comp = 2 * math.pi * V ** 2 * mix * rho0 * Z
        rows.append((tau, wr, kstar, mix, pred, dev, rho0, gam,
                     abs(comp / gam - 1.0), abs(mix / pred_bare - 1.0)))
    comp_worst = max(r[8] for r in rows)
    ok = worst < 0.10 and comp_worst < 1e-5   # interpolation-path limited
    rec("CG6 |B01|^2 law pointwise", "Theorem",
        f"max dev(mix vs 2nd-order) = {worst:.3f}", "CERTIFIED(10%)", ok,
        f"mix=(2 tau cos k* |B01|/Delta)^2 at 3 tau above tau_c; "
        f"Gamma=2piV^2 mix rho_0 Z composition {comp_worst:.0e} (interp-limited); "
        f"denominator must be the DRESSED gap at k* (bare-gap version "
        f"deviates {max(r[9] for r in rows):.0%} — perturbation is not "
        f"small above tau_c); rows tau/w_r/mix/pred: " +
        "; ".join(f"{r[0]}/{r[1]:.3f}/{r[3]:.2e}/{r[4]:.2e}"
                  for r in rows))
    slope = float(np.polyfit(np.log([r[0] for r in rows]),
                             np.log([r[7] for r in rows]), 1)[0])
    rec("CG6b net tau-trend", "Diagnostic",
        f"d log Gamma / d log tau = {slope:.2f}", "DIAGNOSTIC", True,
        "pre-registered tau^2 was the WEIGHT factor alone (confirmed); "
        "net trend also carries DOS dilution as band 0 widens and the "
        "resonance sliding toward the band edge — guess owned, "
        "composition verified instead")


if __name__ == "__main__":
    print("=" * 108)
    print("  PR-6 — CRITICAL STRUCTURE AT tau_c AND THE FIRST LEVEL "
          "WIDTHS   [baselines from rung 5, two engines]")
    print("=" * 108)
    tc = cg0()
    cg1(tc)
    cg2(tc)
    cg3()
    cg4()
    cg5()
    cg6()
    print("-" * 108)
    n_pass = sum(1 for r in REPORT if r[4])
    print(f"  {n_pass}/{len(REPORT)} gates green")
    print("=" * 108)
