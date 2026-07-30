# --- stage 5 (Integration Plan v0.1): re-pointed to the mtft package. ---
# internal() delegates to mtft.chain.internal (certified spectrally
# identical: g to 9e-15, full H(u) spectrum to 2.5e-14; B differs only
# by gauge inside near-degenerate blocks).  Where a local gsq() remains
# it is retained VERBATIM: it evaluates complex u in the f64 backend,
# which mtft.ep deliberately does not offer (PR-20 floors) -- S5-1.
# ------------------------------------------------------------------------
#!/usr/bin/env python3
"""
pr34.py — PR-34(b): the fifth order, and the directionality prediction
=======================================================================
MIT License — (c) 2026 Roger Tano — MTFT Research Program

PRE-REGISTERED BY BOTH ENGINES (BE.2): the fifth order is the (i, i+2)
coupling — upward from the pair's LOWER member — so the directionality
hypothesis (zeros only in the upward-from-UPPER class) predicts

    K5_i has NO zero on continuous i.

Partial test here: compute K5 at i = 2, 3, 4.  Same sign at all three is
consistent with the prediction; a sign change falsifies it outright.

WHY THIS EXTRACTION IS HARDER (BE.2).  P = (i,i+1) couples the crossing
members themselves, so level i cannot be decoupled and the object is a
genuine EP, not a diabatic crossing — Newton on the analytic square
replaces bisection on a level difference.  The count-3 structure
survives because p*r = t is formal in i:
    P = (i+2)/(i+3), R = (i+3)/(i+4), T = (i+2)/(i+4),
    monomials at t^{2 kappa}:  T^2,  P R T,  P^2 R^2.
Each is invariant under the sign-flip symmetries of the block
(conjugation by diag(+,-,+) flips P and R together), as it must be.
"""
from __future__ import annotations
import mpmath as mp

mp.mp.dps = 50
OUT = open("pr34_run.txt", "w", buffering=1)


def say(s):
    print(s, flush=True)
    OUT.write(s + "\n")


def rho(n):
    return mp.log(n) / n ** 3


r2 = rho(2)


def make(i):
    """block {i, i+1, i+2} — the crossing members plus the level above."""
    p0, p1, p2 = rho(i + 2), rho(i + 3), rho(i + 4)

    def block(P, R, T):
        t01 = mp.sqrt(p0 * p1) * P
        t12 = mp.sqrt(p1 * p2) * R
        t02 = mp.sqrt(p0 * p2) * T
        Tm = mp.matrix([[p0, t01, t02], [t01, p1, t12], [t02, t12, p2]])
        M = mp.matrix([[1, P, T], [P, 1, R], [T, R, 1]])
        lam, V = mp.eigsy(Tm)
        o = sorted(range(3), key=lambda j: -lam[j])
        lam = [lam[j] for j in o]
        Vd = mp.matrix(3, 3)
        for c, j in enumerate(o):
            for rw in range(3):
                Vd[rw, c] = V[rw, j]
        g = [mp.log(r2 / lam[j]) for j in range(3)]
        Mt = Vd.T * M * Vd
        B = mp.matrix(3, 3)
        for a in range(3):
            for b in range(3):
                B[a, b] = mp.e ** (-(g[a] + g[b]) / 2) * Mt[a, b]
        return g, B

    def gsq(u, g, B):
        H = mp.matrix(3, 3)
        for a in range(3):
            for b in range(3):
                H[a, b] = (g[a] if a == b else 0) - u * B[a, b]
        ev = sorted(mp.eig(H, left=False, right=False),
                    key=lambda z: mp.re(z))
        return min((ev[k + 1] - ev[k] for k in range(2)), key=abs) ** 2

    def uep(P, R, T, seed):
        g, B = block(P, R, T)
        u = mp.mpc(seed)
        h = mp.mpf(10) ** (-20)
        for _ in range(60):
            f = gsq(u, g, B)
            fp = (gsq(u + h, g, B) - gsq(u - h, g, B)) / (2 * h)
            if fp == 0:
                break
            s = f / fp
            u -= s
            if abs(s) < mp.mpf(10) ** (-40):
                break
        return abs(u)

    g0, _ = block(mp.mpf(0), mp.mpf(0), mp.mpf(0))
    L = abs((g0[1] - g0[0]) / (mp.e ** (-g0[0]) - mp.e ** (-g0[1])))
    seed = -L * (1 + mp.mpf('0.001') * 1j)
    return (lambda P, R, T: uep(P, R, T, seed) - L), L


def coeffs(dev):
    h = mp.mpf(10) ** (-5)
    a = (dev(0, 0, h) + dev(0, 0, -h)) / (2 * h ** 2)
    b = sum(sx * sy * sz * dev(sx * h, sy * h, sz * h)
            for sx in (1, -1) for sy in (1, -1) for sz in (1, -1)) \
        / (8 * h ** 3)
    hh = mp.mpf(10) ** (-4)
    d2 = lambda y: (dev(hh, y, 0) - 2 * dev(0, y, 0) + dev(-hh, y, 0)) \
        / hh ** 2
    c = (d2(hh) - 2 * d2(0) + d2(-hh)) / (4 * hh ** 2)
    return a, b, c


if __name__ == "__main__":
    say("=" * 86)
    say("  PR-34(b) — THE FIFTH ORDER: (i, i+2), UPWARD FROM THE LOWER "
        "MEMBER")
    say("=" * 86)
    say("\n  prediction on record (both engines): K5_i has NO zero")
    say(f"\n{'i':>3} {'limit u':>15} {'[T^2]':>15} {'[PRT]':>15} "
        f"{'[P^2R^2]':>15} {'K5_i':>14}")
    vals = {}
    for i in (2, 3, 4):
        dev, L = make(i)
        a, b, c = coeffs(dev)
        K5 = a + b + c
        vals[i] = K5
        say(f"{i:3d} {mp.nstr(L, 9):>15} {mp.nstr(a, 8):>15} "
            f"{mp.nstr(b, 8):>15} {mp.nstr(c, 8):>15} "
            f"{mp.nstr(K5, 8):>14}")
    sgn = [vals[i] > 0 for i in (2, 3, 4)]
    say(f"\n  signs: {['+' if s else '-' for s in sgn]}")
    say(f"  same sign at all three: {len(set(sgn)) == 1}"
        f"  ->  {'consistent with NO zero' if len(set(sgn)) == 1 else 'SIGN CHANGE: prediction FALSIFIED'}")
    say("\n" + "=" * 86)
    OUT.close()
