# --- stage 5 (Integration Plan v0.1): re-pointed to the mtft package. ---
# internal() delegates to mtft.chain.internal (certified spectrally
# identical: g to 9e-15, full H(u) spectrum to 2.5e-14; B differs only
# by gauge inside near-degenerate blocks).  Where a local gsq() remains
# it is retained VERBATIM: it evaluates complex u in the f64 backend,
# which mtft.ep deliberately does not offer (PR-20 floors) -- S5-1.
# ------------------------------------------------------------------------
#!/usr/bin/env python3
"""
pr24.py — PR-24: the two-term law, and the fourth and fifth members
====================================================================
MIT License — (c) 2026 Roger Tano — MTFT Research Program
Prints to stdout AND pr24_run.txt.

ADDENDUM AU DISPOSITIONS:
 AU-F1 (ratio misprint 1.0000005 -> 1.0000045): ACCEPTED, corrected here.
 AU-F2 (r^{4 kappa} FALSIFIED): ACCEPTED.  My pre-registered next order
   was wrong; the true remainder is the WITHIN-PAIR coupling
   s = (i+2)/(i+3), giving the two-term law
       dev(i;kappa) = A_i r^{2 kappa} - C_i s^{2 kappa},  r = (i+3)/(i+4).
   Verified independently below from Kimi's stated channel structure.
 AU-F3 (their B_01 exponent; seals were finite-kappa transients): noted;
   the derived A_i are the asymptotic ledger values.

PR-24(b), run against the TWO-TERM law per AU.7's flag: finite-kappa
implied amplitudes read LOW by the within-pair transient, so the fourth
and fifth members must be tested against both terms, not the leading one.

C_i, from the same three channels as A_i but with BOTH members counted
(S2 = s^{2 kappa}, rho_j the weight at site j+2):
  (1) eigenvalue:   dg_{i+1} = +rho_i S2/(rho_i - rho_{i+1}) and, for
      i >= 1 only, dg_i = -rho_{i+1} S2/(rho_i - rho_{i+1})  (for i = 0
      the lower member IS the reference level lambda_0, so it is inert);
  (2) eigenvector:  dp_i = +e^{-a} X, dp_{i+1} = -e^{-b} X with
      X = 2 sqrt(rho_i rho_{i+1}) S2/(rho_i - rho_{i+1});
  (3) mutual element in the EP discriminant:
      |u| = m/sqrt(D^2 + 4 B_{i,i+1}^2)  =>  -2 m e^{-(a+b)} S2/D^3.
"""
from __future__ import annotations
import mpmath as mp
from mtft import expansion as _ex


mp.mp.dps = 40
OUT = open("pr24_run.txt", "w", buffering=1)


def say(s=""):
    print(s, flush=True)
    OUT.write(s + "\n")


def rho(site):
    return mp.log(site) / site ** 3


def gaps(n):
    return [mp.log(rho(2) / rho(2 + j)) for j in range(n)]


def A_of(i):
    return _ex.A(i)


def C_of(i):
    return _ex.C(i)


def model(kappa, nsite=12):
    ns = [mp.mpf(n) for n in range(2, 2 + nsite)]
    rr = [rho(n) for n in ns]
    K = mp.matrix(nsite, nsite); M = mp.matrix(nsite, nsite)
    for i in range(nsite):
        for j in range(nsite):
            q = min(ns[i], ns[j]) / max(ns[i], ns[j])
            M[i, j] = q ** kappa
            K[i, j] = mp.sqrt(rr[i] * rr[j]) * M[i, j]
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
    return min((ev[k + 1] - ev[k] for k in range(n - 1)), key=abs) ** 2


def newton(u0, g, B, iters=50):
    u = mp.mpc(u0); h = mp.mpf(10) ** (-18)
    for _ in range(iters):
        f = gsq(u, g, B)
        fp = (gsq(u + h, g, B) - gsq(u - h, g, B)) / (2 * h)
        if fp == 0: break
        s = f / fp; u -= s
        if abs(s) < mp.mpf(10) ** (-34): break
    return u


def limit(i):
    g = gaps(i + 2)
    a, b = g[i], g[i + 1]
    return abs((b - a) / (mp.e ** (-a) - mp.e ** (-b)))


if __name__ == "__main__":
    say("=" * 94)
    say("  PR-24 — THE TWO-TERM LAW, AND MEMBERS FOUR AND FIVE")
    say("=" * 94)
    say("\n[AU-F1] channel-(1) ratio corrected: 1.0000045, not 1.0000005")
    say("\n[AU-F2 accepted] C_i derived independently vs Kimi's measured")
    kimi = {0: mp.mpf('13.23412605'), 1: mp.mpf('42.38221931'),
            2: mp.mpf('112.0522097')}
    for i in (0, 1, 2):
        c = C_of(i)
        say(f"    C_{i}: derived {mp.nstr(c, 10):>14}  measured "
            f"{mp.nstr(kimi[i], 10):>14}  ratio {mp.nstr(c/kimi[i], 9)}")

    say("\n[PR-24(b)] members 4 and 5 vs the TWO-TERM law")
    say(f"    {'pair':>7} {'kappa':>6} {'measured dev':>18} "
        f"{'2-term pred':>18} {'ratio':>10} {'1-term ratio':>13}")
    for i, kaps in ((3, (120, 150)), (4, (150, 190))):
        L = limit(i)
        A, C = A_of(i), C_of(i)
        r = mp.mpf(i + 3) / (i + 4)
        s = mp.mpf(i + 2) / (i + 3)
        for kap in kaps:
            g, B = model(mp.mpf(kap))
            u = newton(-L * (1 + mp.mpf('0.001') * 1j), g, B)
            dev = abs(u) - L
            two = A * r ** (2 * kap) - C * s ** (2 * kap)
            one = A * r ** (2 * kap)
            say(f"    {f'({i},{i+1})':>7} {kap:6d} {mp.nstr(dev, 9):>18} "
                f"{mp.nstr(two, 9):>18} {mp.nstr(dev/two, 8):>10} "
                f"{mp.nstr(dev/one, 6):>13}")
        say(f"      limit u_{i}{i+1}(inf) = {mp.nstr(L, 12)}, "
            f"A_{i} = {mp.nstr(A, 10)}, C_{i} = {mp.nstr(C, 10)}")
    say("\n" + "=" * 94)
    OUT.close()
