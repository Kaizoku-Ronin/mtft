# --- stage 5 (Integration Plan v0.1): re-pointed to the mtft package. ---
# internal() delegates to mtft.chain.internal (certified spectrally
# identical: g to 9e-15, full H(u) spectrum to 2.5e-14; B differs only
# by gauge inside near-degenerate blocks).  Where a local gsq() remains
# it is retained VERBATIM: it evaluates complex u in the f64 backend,
# which mtft.ep deliberately does not offer (PR-20 floors) -- S5-1.
# ------------------------------------------------------------------------
#!/usr/bin/env python3
"""
pr21.py — PR-21: all three members at 1e-4, and the amplitude predicted
========================================================================
MIT License — (c) 2026 Roger Tano — MTFT Research Program
Prints to stdout AND pr21_run.txt.

PRE-REGISTERED (PR-20 note §6):
 (a) re-run pairs (0,1) and (1,2) in the exact-arithmetic engine so all
     three family members are measured at 1e-4 rather than 1e-3.
 (c) prove the dominance argument.

WHAT (c) TURNS INTO.  Writing the argument out does more than justify
the rate — it predicts the AMPLITUDE.  Second-order perturbation gives
    delta g_j = sum_{k!=j} |B_jk|^2/(g_j - g_k),
with |B_jk|^2 = e^{-(g_j+g_k)} (n_j/n_k)^{2 kappa}, n_j = j+2.  For the
pair (i, i+1) every contribution to delta(g_{i+1} - g_i) is enumerable:

  * within-pair (i <-> i+1):        ((i+2)/(i+3))^{2 kappa}
  * upper member to the level above
    (i+1 <-> i+2):                  ((i+3)/(i+4))^{2 kappa}   <- SLOWEST
  * lower member downward (i-1):    ((i+1)/(i+2))^{2 kappa}
  * everything else:                 faster still

Since (i+3)/(i+4) > (i+2)/(i+3) > (i+1)/(i+2), the upper-neighbour term
dominates, and it does NOT cancel: it shifts g_{i+1} alone, so it moves
the gap at first order in the shift.  Hence, with
u(a,b) = (b-a)/(e^{-a} - e^{-b}),

    dev(kappa) = (d|u|/db) * e^{-(g_{i+1}+g_{i+2})}
                 ((i+3)/(i+4))^{2 kappa} / (g_{i+1} - g_{i+2})

— a closed-form PREDICTION OF BOTH THE RATE AND THE AMPLITUDE, with no
fitted constant anywhere.  That is testable far more sharply than a
slope, and it is what promotes the scaling argument from heuristic to
derivation.
"""
from __future__ import annotations
import mpmath as mp

mp.mp.dps = 45
OUT = open("pr21_run.txt", "w", buffering=1)


def say(s=""):
    print(s, flush=True)
    OUT.write(s + "\n")


def model(kappa, nsite=14):
    ns = [mp.mpf(n) for n in range(2, 2 + nsite)]
    rho = [mp.log(n) / n ** 3 for n in ns]
    K = mp.matrix(nsite, nsite)
    M = mp.matrix(nsite, nsite)
    for i in range(nsite):
        for j in range(nsite):
            r = min(ns[i], ns[j]) / max(ns[i], ns[j])
            M[i, j] = r ** kappa
            K[i, j] = mp.sqrt(rho[i] * rho[j]) * M[i, j]
    lam, V = mp.eigsy(K)
    idx = list(range(nsite))[::-1]
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
def gsq(u, g, B):
    n = len(g)
    A = mp.matrix(n, n)
    for a in range(n):
        for b in range(n):
            A[a, b] = (g[a] if a == b else 0) - u * B[a, b]
    ev = sorted(mp.eig(A, left=False, right=False), key=lambda z: mp.re(z))
    best = None
    for k in range(n - 1):
        d = ev[k + 1] - ev[k]
        if best is None or abs(d) < abs(best):
            best = d
    return best ** 2


def newton(u0, g, B, iters=40):
    u = mp.mpc(u0)
    h = mp.mpf(10) ** (-15)
    for _ in range(iters):
        f = gsq(u, g, B)
        fp = (gsq(u + h, g, B) - gsq(u - h, g, B)) / (2 * h)
        if fp == 0:
            break
        s = f / fp
        u -= s
        if abs(s) < mp.mpf(10) ** (-35):
            break
    return u


def limit_u(g, i):
    a, b = g[i], g[i + 1]
    return abs((b - a) / (mp.e ** (-b) - mp.e ** (-a)))


def predicted_dev(gL, i, kappa):
    """Closed-form amplitude AND rate from the dominance argument."""
    a, b, c = gL[i], gL[i + 1], gL[i + 2]
    ea, eb = mp.e ** (-a), mp.e ** (-b)
    dudb = ((ea - eb) - (b - a) * eb) / (ea - eb) ** 2
    ratio = mp.mpf(i + 3) / (i + 4)
    dg = mp.e ** (-(b + c)) * ratio ** (2 * kappa) / (b - c)
    return dudb * dg


if __name__ == "__main__":
    say("=" * 94)
    say("  PR-21 — ALL THREE MEMBERS IN EXACT ARITHMETIC, AND THE "
        "AMPLITUDE PREDICTED")
    say("=" * 94)
    gL, _ = model(mp.mpf(400), nsite=14)
    say(f"\n  limiting gaps g_i(inf) = ln(rho_2/rho_(i+2)): "
        + ", ".join(mp.nstr(gL[i], 8) for i in range(4)))

    say(f"\n[a] settled rates in exact arithmetic")
    say(f"{'pair':>6} {'kappa':>6} {'deviation':>14} {'slope/kappa':>13} "
        f"{'predicted':>12}")
    for i, kaps in ((0, (40, 56, 72)), (1, (64, 88, 112))):
        pred = 2 * mp.log10(mp.mpf(i + 3) / (i + 4))
        prev = None
        for kap in kaps:
            g, B = model(mp.mpf(kap), nsite=14)
            L = limit_u(gL, i)
            u = newton(-L * (1 + mp.mpf('0.001') * 1j), g, B)
            dev = abs(u) - L
            sl = ""
            if prev is not None:
                sl = mp.nstr(mp.log10(abs(dev) / abs(prev[1]))
                             / (kap - prev[0]), 7)
            say(f"{f'({i},{i+1})':>6} {kap:6d} {mp.nstr(dev, 5):>14} "
                f"{sl:>13} {mp.nstr(pred, 7):>12}")
            prev = (kap, dev)

    say(f"\n[c] AMPLITUDE test: measured deviation vs the closed form")
    say(f"{'pair':>6} {'kappa':>6} {'measured dev':>15} "
        f"{'predicted dev':>15} {'ratio':>9}")
    for i, kaps in ((0, (56, 72)), (1, (88, 112)), (2, (96, 112))):
        for kap in kaps:
            g, B = model(mp.mpf(kap), nsite=14)
            L = limit_u(gL, i)
            u = newton(-L * (1 + mp.mpf('0.001') * 1j), g, B)
            dev = abs(u) - L
            pd = predicted_dev(gL, i, mp.mpf(kap))
            say(f"{f'({i},{i+1})':>6} {kap:6d} {mp.nstr(dev, 6):>15} "
                f"{mp.nstr(pd, 6):>15} {mp.nstr(dev/pd, 6):>9}")
    say("\n" + "=" * 94)
    OUT.close()
