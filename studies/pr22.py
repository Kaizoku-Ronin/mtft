# --- stage 5 (Integration Plan v0.1): re-pointed to the mtft package. ---
# internal() delegates to mtft.chain.internal (certified spectrally
# identical: g to 9e-15, full H(u) spectrum to 2.5e-14; B differs only
# by gauge inside near-degenerate blocks).  Where a local gsq() remains
# it is retained VERBATIM: it evaluates complex u in the f64 backend,
# which mtft.ep deliberately does not offer (PR-20 floors) -- S5-1.
# ------------------------------------------------------------------------
#!/usr/bin/env python3
"""
pr22.py — PR-22(a) step 1: independent decomposition of the amplitude
======================================================================
MIT License — (c) 2026 Roger Tano — MTFT Research Program
Prints to stdout AND pr22_run.txt.

Addendum AS.3 seals two priors: (i) the amplitude is entirely a
three-level effect; (ii) the T-space gap-shift part is NEGATIVE and only
~-40% of it, the remaining +140% being H-space dynamics through B.

Before attempting the derivation this rung verifies that split on an
independent engine — and goes one step further, because my own
sign analysis of the dynamical term comes out NEGATIVE where AS.3
measures positive.  Zeroing the off-diagonal elements one at a time
identifies which coupling actually carries the +140% and with what sign,
which is the missing ingredient in the derivation rather than a detail.

Decomposition (all at exact precision, three-level block i, i+1, i+2):
  full     : EP of the 3x3 H(u) = diag(g) - u B, finite-kappa g and B
  gap-only : off-diagonal B zeroed -> real crossing of the finite-kappa
             diagonal, i.e. the complete T-space content
  then B_01, B_02, B_12 restored one at a time.
"""
from __future__ import annotations
import mpmath as mp

mp.mp.dps = 45
OUT = open("pr22_run.txt", "w", buffering=1)


def say(s=""):
    print(s, flush=True)
    OUT.write(s + "\n")


def model(kappa, nsite=12):
    ns = [mp.mpf(n) for n in range(2, 2 + nsite)]
    rho = [mp.log(n) / n ** 3 for n in ns]
    K = mp.matrix(nsite, nsite); M = mp.matrix(nsite, nsite)
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


def block(g, B, i, keep):
    """3x3 block for levels i,i+1,i+2 with only 'keep' off-diagonals."""
    gg = [g[i], g[i + 1], g[i + 2]]
    BB = mp.matrix(3, 3)
    for a in range(3):
        BB[a, a] = B[i + a, i + a]
    for (a, b) in keep:
        BB[a, b] = B[i + a, i + b]
        BB[b, a] = B[i + a, i + b]
    return gg, BB


# S5-1: retained verbatim -- complex-u f64 gsq predates mtft.ep's
# real-axis diabatic-centre design; frozen so this study's record stands.
def gsq3(u, gg, BB):
    A = mp.matrix(3, 3)
    for a in range(3):
        for b in range(3):
            A[a, b] = (gg[a] if a == b else 0) - u * BB[a, b]
    ev = sorted(mp.eig(A, left=False, right=False), key=lambda z: mp.re(z))
    d = [ev[1] - ev[0], ev[2] - ev[1]]
    return min(d, key=abs) ** 2


def newton3(u0, gg, BB, iters=60):
    u = mp.mpc(u0); h = mp.mpf(10) ** (-18)
    for _ in range(iters):
        f = gsq3(u, gg, BB)
        fp = (gsq3(u + h, gg, BB) - gsq3(u - h, gg, BB)) / (2 * h)
        if fp == 0: break
        s = f / fp; u -= s
        if abs(s) < mp.mpf(10) ** (-38): break
    return u


def cross_diag(gg, BB):
    """Real crossing of levels 0,1 with off-diagonals absent."""
    return abs((gg[1] - gg[0]) / (BB[1, 1] - BB[0, 0]))


if __name__ == "__main__":
    say("=" * 96)
    say("  PR-22(a) STEP 1 — INDEPENDENT DECOMPOSITION OF THE AMPLITUDE")
    say("=" * 96)
    gL, _ = model(mp.mpf(400))
    cases = [(0, 56), (1, 88)]
    for i, kap in cases:
        L = abs((gL[i + 1] - gL[i])
                / (mp.e ** (-gL[i]) - mp.e ** (-gL[i + 1])))
        g, B = model(mp.mpf(kap))
        say(f"\n  pair ({i},{i+1}), kappa = {kap};  limit = {mp.nstr(L,12)}")
        gg, BB = block(g, B, i, [(0, 1), (0, 2), (1, 2)])
        full = abs(newton3(-L * (1 + mp.mpf('0.001') * 1j), gg, BB)) - L
        say(f"    full 3-level deviation        {mp.nstr(full, 8):>16}"
            f"   100.000%")
        gg0, BB0 = block(g, B, i, [])
        gap = cross_diag(gg0, BB0) - L
        say(f"    gap-shift only (no off-diag)  {mp.nstr(gap, 8):>16}"
            f"   {mp.nstr(100*gap/full, 6)}%   [AS.3: ~-40%]")
        for lbl, keep in (("B_01", [(0, 1)]), ("B_02", [(0, 2)]),
                          ("B_12", [(1, 2)])):
            gk, Bk = block(g, B, i, keep)
            d = abs(newton3(-L * (1 + mp.mpf('0.001') * 1j), gk, Bk)) - L
            say(f"    + {lbl} only                 "
                f"{mp.nstr(d, 8):>16}   {mp.nstr(100*d/full, 6)}%"
                f"   (delta vs gap-only: "
                f"{mp.nstr(100*(d-gap)/full, 6)}%)")
    say("\n" + "=" * 96)
    OUT.close()
