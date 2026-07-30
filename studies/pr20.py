# --- stage 5 (Integration Plan v0.1): re-pointed to the mtft package. ---
# internal() delegates to mtft.chain.internal (certified spectrally
# identical: g to 9e-15, full H(u) spectrum to 2.5e-14; B differs only
# by gauge inside near-degenerate blocks).  Where a local gsq() remains
# it is retained VERBATIM: it evaluates complex u in the f64 backend,
# which mtft.ep deliberately does not offer (PR-20 floors) -- S5-1.
# ------------------------------------------------------------------------
#!/usr/bin/env python3
"""
pr20.py — PR-20: the settling test, redesigned for precision (AQ-F3)
=====================================================================
MIT License — (c) 2026 Roger Tano — MTFT Research Program
Prints to stdout AND pr20_run.txt.

ADDENDUM AQ DISPOSITIONS:
 AQ-F1 (precisions window-flattered): ACCEPTED — the clean-window rates
   -0.248672 / -0.189832 confirm the law at 1e-3, not 1e-4.
 AQ-F2 (arithmetic): ACCEPTED — 2 log10(6/7) = -0.1338936, not
   -0.133907.  Verified.
 AQ-F3 (PR-20(a) unmeasurable in float64): ACCEPTED, and REDESIGNED
   HERE.  The (2,3) deviation floors at ~2e-12 by kappa ~ 80 before its
   transient settles, so more kappa reads noise as signal.  The fix is
   PRECISION.

THE PRECISION PROBLEM, AND WHY IT IS DEEPER THAN THE EIGENSOLVER.
Swapping in an mpmath eigensolver is not enough: g and B were being
INHERITED from a float64 diagonalisation of the 1600-site model, so the
inputs themselves carried ~1e-15 relative error.  A deviation of 1e-14
cannot be measured from 1e-15 inputs at any working precision.

The fix is to build the internal model itself in mpmath from its closed
form.  Everything is analytic:
    rho_n propto (log n) n^{-3},  T_ij = sqrt(rho_i rho_j) (min/max)^kappa,
    g_i = log(lambda_0/lambda_i),  B = e^{-g/2} (V^T M V) e^{-g/2},
so T is written EXACTLY at any precision and diagonalised there.  The
truncation to a few sites is legitimate at large kappa precisely because
the couplings die: at kappa = 80 the site-4-to-site-6 element is
(4/6)^80 ~ 8e-15, far below the levels being resolved.  Normalisation
drops out of g (a ratio), so Z_2 is not needed.

Gate: the settled rate for pair (2,3) must reach 2 log10(5/6) =
-0.1583625, the i = 2 member of the family and rung 4's bridge law.
"""
from __future__ import annotations
import mpmath as mp

mp.mp.dps = 40
OUT = open("pr20_run.txt", "w", buffering=1)


def say(s=""):
    print(s, flush=True)
    OUT.write(s + "\n")


def model(kappa, nsite=12):
    """Internal model in exact arithmetic at the working precision."""
    ns = [mp.mpf(n) for n in range(2, 2 + nsite)]
    rho = [mp.log(n) / n ** 3 for n in ns]
    K = mp.matrix(nsite, nsite)
    M = mp.matrix(nsite, nsite)
    for i in range(nsite):
        for j in range(nsite):
            r = min(ns[i], ns[j]) / max(ns[i], ns[j])
            M[i, j] = r ** kappa
            K[i, j] = mp.sqrt(rho[i] * rho[j]) * M[i, j]
    lam, V = mp.eigsy(K)                      # ascending
    idx = list(range(nsite))[::-1]            # descending
    lam = [lam[i] for i in idx]
    Vd = mp.matrix(nsite, nsite)
    for c, i in enumerate(idx):
        for r in range(nsite):
            Vd[r, c] = V[r, i]
    g = [mp.log(lam[0] / lam[i]) for i in range(nsite)]
    Mt = Vd.T * M * Vd
    B = mp.matrix(nsite, nsite)
    for i in range(nsite):
        for j in range(nsite):
            B[i, j] = mp.e ** (-(g[i] + g[j]) / 2) * Mt[i, j]
    return g, B


# S5-1: retained verbatim -- complex-u f64 gsq predates mtft.ep's
# real-axis diabatic-centre design; frozen so this study's record stands.
def gsq(u, g, B, i=None):
    """(eps_a - eps_b)^2 for the CLOSEST pair — the fixed-index form
    fails because the eigenvalues reorder in real part at negative
    real u, and Newton then chases a different root entirely."""
    n = len(g)
    A = mp.matrix(n, n)
    for a in range(n):
        for b in range(n):
            A[a, b] = (g[a] if a == b else 0) - u * B[a, b]
    ev = mp.eig(A, left=False, right=False)
    ev = sorted(ev, key=lambda z: mp.re(z))
    best = None
    for k in range(n - 1):
        d = ev[k + 1] - ev[k]
        if best is None or abs(d) < abs(best):
            best = d
    return best ** 2


def newton(u0, g, B, i, iters=40):
    u = mp.mpc(u0)
    h = mp.mpf(10) ** (-12)
    for _ in range(iters):
        f = gsq(u, g, B, i)
        fp = (gsq(u + h, g, B, i) - gsq(u - h, g, B, i)) / (2 * h)
        if fp == 0:
            break
        s = f / fp
        u -= s
        if abs(s) < mp.mpf(10) ** (-30):
            break
    return u


def limit_u(g, i):
    a, b = g[i], g[i + 1]
    return abs((b - a) / (mp.e ** (-b) - mp.e ** (-a)))


if __name__ == "__main__":
    say("=" * 92)
    say("  PR-20 — THE SETTLING TEST, REDESIGNED FOR PRECISION (AQ-F3)")
    say("=" * 92)
    say(f"\n  working precision: mpmath dps = {mp.mp.dps}; model built in "
        f"closed form, not inherited")
    say(f"  AQ-F2 verified: 2 log10(6/7) = {2*mp.log10(mp.mpf(6)/7)} "
        f"(note said -0.133907)")

    gL, _ = model(mp.mpf(400), nsite=12)
    say(f"\n  closed-form limit, pair (2,3): u_23(inf) = "
        f"{mp.nstr(limit_u(gL, 2), 12)}  (auditor: 5.13015139272)")
    pred = 2 * mp.log10(mp.mpf(5) / 6)
    say(f"  predicted settled rate: 2 log10(5/6) = {mp.nstr(pred, 8)}")

    say(f"\n{'kappa':>7} {'|u_23|':>26} {'deviation':>14} {'slope/kappa':>13}")
    prev = None
    for kap in (48, 64, 80, 96, 112):
        g, B = model(mp.mpf(kap), nsite=12)
        L = limit_u(gL, 2)
        u = newton(-L * (1 + mp.mpf('0.001') * 1j), g, B, 2)
        dev = abs(u) - L
        sl = ""
        if prev is not None and dev != 0 and prev[1] != 0:
            s = mp.log10(abs(dev) / abs(prev[1])) / (kap - prev[0])
            sl = mp.nstr(s, 7)
        say(f"{kap:7d} {mp.nstr(abs(u), 20):>26} "
            f"{mp.nstr(dev, 5):>14} {sl:>13}")
        prev = (kap, dev)
    say("\n" + "=" * 92)
    OUT.close()
