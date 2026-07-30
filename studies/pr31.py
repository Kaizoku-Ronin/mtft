# --- stage 5 (Integration Plan v0.1): re-pointed to the mtft package. ---
# internal() delegates to mtft.chain.internal (certified spectrally
# identical: g to 9e-15, full H(u) spectrum to 2.5e-14; B differs only
# by gauge inside near-degenerate blocks).  Where a local gsq() remains
# it is retained VERBATIM: it evaluates complex u in the f64 backend,
# which mtft.ep deliberately does not offer (PR-20 floors) -- S5-1.
# ------------------------------------------------------------------------
#!/usr/bin/env python3
"""
pr31.py — PR-31(b): K_i for general i (scope i >= 2)
=====================================================
MIT License — (c) 2026 Roger Tano — MTFT Research Program

BB dispositions: K_4 = -54.78257 adopted as the ledger value (direct,
Richardson-extrapolated); the closure's -54.7822 demoted to "kappa-fit,
+-5e-4" because slowly separated contaminant classes (0.9604^kappa,
0.9795^kappa) bias any attainable kappa-window at exactly the 1e-4 scale.
X_rep -> -729.3137.

PR-31(b).  Same construction, pair (i, i+1), scope i >= 2 (AZ's
completeness bound): decouple level i, keep only the (i+1,i+2)=R,
(i+2,i+3)=S, (i+1,i+3)=Q couplings, and monomial-extract the crossing.
Every stencil validated separately (PR-30 sec.3), and the R^2 S^2
stencil carries the corrected /4.

Control: i = 3 must reproduce -54.7826.
"""
from __future__ import annotations
import mpmath as mp

mp.mp.dps = 50
OUT = open("pr31_run.txt", "w", buffering=1)


def say(s):
    print(s, flush=True)
    OUT.write(s + "\n")


def rho(n):
    return mp.log(n) / n ** 3


r2 = rho(2)


def make(i):
    """levels i (decoupled) and the block {i+1, i+2, i+3}."""
    pl, p1, p2, p3 = rho(i + 2), rho(i + 3), rho(i + 4), rho(i + 5)
    gl = mp.log(r2 / pl)
    Bl = mp.e ** (-gl)

    def block(R, S, Q):
        t12 = mp.sqrt(p1 * p2) * R
        t23 = mp.sqrt(p2 * p3) * S
        t13 = mp.sqrt(p1 * p3) * Q
        T = mp.matrix([[p1, t12, t13], [t12, p2, t23], [t13, t23, p3]])
        M = mp.matrix([[1, R, Q], [R, 1, S], [Q, S, 1]])
        lam, V = mp.eigsy(T)
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

    def ucross(R, S, Q):
        g, B = block(R, S, Q)

        def f(u):
            H = mp.matrix(3, 3)
            for a in range(3):
                for b in range(3):
                    H[a, b] = (g[a] if a == b else 0) - u * B[a, b]
            return (gl - u * Bl) - min(mp.eigsy(H, eigvals_only=True))

        lo, hi = mp.mpf(-60), mp.mpf(-1)
        flo = f(lo)
        for _ in range(220):
            mid = (lo + hi) / 2
            if f(mid) * flo > 0:
                lo = mid
            else:
                hi = mid
            if hi - lo < mp.mpf(10) ** (-45):
                break
        return -(lo + hi) / 2

    L = ucross(mp.mpf(0), mp.mpf(0), mp.mpf(0))
    return (lambda R, S, Q: ucross(R, S, Q) - L), L


def coeffs(dev):
    h = mp.mpf(10) ** (-5)
    a = (dev(0, 0, h) + dev(0, 0, -h)) / (2 * h ** 2)
    b = sum(sx * sy * sz * dev(sx * h, sy * h, sz * h)
            for sx in (1, -1) for sy in (1, -1) for sz in (1, -1)) \
        / (8 * h ** 3)
    hh = mp.mpf(10) ** (-4)
    d2 = lambda y: (dev(hh, y, 0) - 2 * dev(0, y, 0) + dev(-hh, y, 0)) \
        / hh ** 2
    c = (d2(hh) - 2 * d2(0) + d2(-hh)) / (4 * hh ** 2)   # /4, per BA
    return a, b, c


if __name__ == "__main__":
    say("=" * 88)
    say("  PR-31(b) — K_i FOR GENERAL i  (scope i >= 2)")
    say("=" * 88)
    say(f"\n  ledger adopted from BB: K_4 = -54.78257, X_rep = -729.3137")
    say(f"\n{'i':>3} {'limit u':>16} {'[Q^2]':>16} {'[RSQ]':>16} "
        f"{'[R^2S^2]':>16} {'K_i':>15}")
    for i in (2, 3, 4):
        dev, L = make(i)
        a, b, c = coeffs(dev)
        K = a + b + c
        tag = "  <- control, target -54.7826" if i == 3 else ""
        say(f"{i:3d} {mp.nstr(L, 10):>16} {mp.nstr(a, 9):>16} "
            f"{mp.nstr(b, 9):>16} {mp.nstr(c, 9):>16} "
            f"{mp.nstr(K, 9):>15}{tag}")
    say("\n" + "=" * 88)
    OUT.close()
