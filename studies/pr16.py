# --- stage 5 (Integration Plan v0.1): re-pointed to the mtft package. ---
# internal() delegates to mtft.chain.internal (certified spectrally
# identical: g to 9e-15, full H(u) spectrum to 2.5e-14; B differs only
# by gauge inside near-degenerate blocks).  Where a local gsq() remains
# it is retained VERBATIM: it evaluates complex u in the f64 backend,
# which mtft.ep deliberately does not offer (PR-20 floors) -- S5-1.
# ------------------------------------------------------------------------
#!/usr/bin/env python3
"""
pr16.py — PR-16: the occupancy check, and a four-parameter radius
==================================================================
MIT License — (c) 2026 Roger Tano — MTFT Research Program
Prints to stdout AND pr16_run.txt.

PRE-REGISTERED (PR-15 note §5):
 (a) add the OCCUPANCY CHECK so "not found" can never masquerade as
     "not there"; rerun kappa = 3, 4, 6.  Addendum AM-F1 supplies the
     natural engine: the LEVEL-RESOLVED PAIR WINDING — track eps_0 and
     eps_1 around |u| = r by continuation and wind their difference;
     each order-2 branch point of the 0-1 connection contributes a
     half-winding, so (level-0/1 EPs inside) = 2 x winding.  No Newton,
     no basins, and it certifies OCCUPANCY and LEVEL IDENTITY in one
     step.  Gate 1 is an INDEPENDENT re-derivation of AM-F1's numbers.
 (b) a sharper contour-free estimator.  The singularity is a conjugate
     pair, so model the Cauchy coefficients directly:
         c_k ~ A k^{-3/2} R^{-k} cos(k theta + phi),
     four parameters against 26 coefficients.  For fixed (R, theta) the
     model is LINEAR in (A cos phi, A sin phi), so the fit reduces to a
     2-D scan with an exact inner solve.  If it returns 1.358 it is a
     genuine second route to the radius.
 (c) the PR-6 gate-string amendment (documentation, in the note).
"""
from __future__ import annotations
import math
import numpy as np
import mpmath as mp
from mtft.chain import internal as _chain_internal


BETA = 2.0
OUT = open("pr16_run.txt", "w", buffering=1)


def say(s=""):
    print(s, flush=True)
    OUT.write(s + "\n")


def internal(N=1600, kappa=5.0, nb=30, gcap=200.0):
    ic = _chain_internal(kappa, nb=nb, backend="f64", N=N, gcap=gcap)
    return np.asarray(ic.g), np.asarray(ic.B)


def pair_winding(r, g, B, npts=600):
    """Winding of eps_1 - eps_0 around |u| = r, sheets tracked by
    continuation.  Returns (level-0/1 EP count inside) = 2 * winding."""
    th = 2 * math.pi * np.arange(npts + 1) / npts
    ev = np.linalg.eigvalsh(np.diag(g) - r * B)
    a, b = complex(ev[0]), complex(ev[1])
    tot = 0.0
    qprev = b - a
    for t in th[1:]:
        u = r * np.exp(1j * t)
        cur = np.linalg.eigvals(np.diag(g) - u * B)
        ia = int(np.argmin(np.abs(cur - a)))
        a = cur[ia]
        cur2 = np.delete(cur, ia)
        b = cur2[int(np.argmin(np.abs(cur2 - b)))]
        q = b - a
        tot += np.angle(q / qprev)
        qprev = q
    return 2.0 * tot / (2 * math.pi)


def R_certified(kap, lo=0.30, hi=2.60, nb=30, npts=600):
    """Bisect on the pair winding: occupancy 0 below R, 2 above."""
    if pair_winding(lo, *internal(kappa=kap, nb=nb), npts=npts) > 1.0:
        return float('nan')
    g, B = internal(kappa=kap, nb=nb)
    if pair_winding(hi, g, B, npts) < 1.0:
        return float('nan')
    for _ in range(34):
        mid = 0.5 * (lo + hi)
        if pair_winding(mid, g, B, npts) < 1.0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def cauchy_coeffs(kap=5.0, rho=1.0, npts=512, kmax=40, nb=30):
    g, B = internal(kappa=kap, nb=nb)
    th = 2 * math.pi * np.arange(npts) / npts
    vals = np.empty(npts, dtype=complex)
    cur = complex(np.linalg.eigvalsh(np.diag(g) - rho * B)[1])
    for j, t in enumerate(th):
        u = rho * np.exp(1j * t)
        ev = np.linalg.eigvals(np.diag(g) - u * B)
        cur = ev[int(np.argmin(np.abs(ev - cur)))]
        vals[j] = cur
    return np.array([np.sum(vals * np.exp(-1j * k * th)) / npts / rho ** k
                     for k in range(kmax + 1)]).real


def fit_pair_model(cs, ks):
    """c_k = k^{-3/2} R^{-k} [a cos(k th) - b sin(k th)]; scan (R, th),
    solve (a, b) exactly by least squares at each node."""
    y = cs[ks]
    best = None
    for R in np.arange(1.20, 1.65, 0.0005):
        base = ks ** -1.5 * R ** (-ks.astype(float))
        for thh in np.arange(2.60, 3.00, 0.002):
            A = np.column_stack([base * np.cos(ks * thh),
                                 -base * np.sin(ks * thh)])
            sol, res, *_ = np.linalg.lstsq(A, y, rcond=None)
            r2 = float(np.sum((A @ sol - y) ** 2))
            if best is None or r2 < best[0]:
                best = (r2, R, thh, sol)
    return best


if __name__ == "__main__":
    say("=" * 92)
    say("  PR-16 — OCCUPANCY CHECK (pair winding) AND A FOUR-PARAMETER "
        "RADIUS")
    say("=" * 92)
    say("\n[a] R(kappa) by level-resolved pair winding — independent "
        "re-derivation of AM-F1")
    say(f"{'kappa':>6} {'R (this engine)':>18} {'AM-F1 (Kimi)':>16} "
        f"{'deviation':>12}")
    kimi = {3.0: 1.0242466134, 4.0: 1.2568164063, 5.0: 1.3578911871,
            6.0: 1.3995852406, 8.0: 1.4244407112}
    mine = {}
    for kap in (3.0, 4.0, 5.0, 6.0, 8.0):
        R = R_certified(kap)
        mine[kap] = R
        say(f"{kap:6.1f} {R:18.7f} {kimi[kap]:16.7f} "
            f"{abs(R - kimi[kap]):12.2e}")
    say(f"  monotone increasing in kappa: "
        f"{all(mine[a] < mine[b] for a, b in zip([3.,4.,5.,6.], [4.,5.,6.,8.]))}")
    say("  nb-stability at kappa=4 (nb 30 vs 40): "
        f"{abs(R_certified(4.0, nb=40) - mine[4.0]):.2e}")

    say("\n[b] four-parameter coefficient fit (conjugate-pair model)")
    cs = cauchy_coeffs()
    ks = np.arange(15, 41)
    r2, R, thh, sol = fit_pair_model(cs, ks)
    say(f"  best fit over k=15..40: R = {R:.4f}, theta = {thh:.4f} rad")
    say(f"  certified: R = 1.357891, theta = 2.7995 rad")
    say(f"  deviations: R {abs(R-1.357891):.4f} ({100*abs(R-1.357891)/1.357891:.2f}%), "
        f"theta {abs(thh-2.7995):.4f}")
    say(f"  residual sum sq = {r2:.3e}; amplitude (a,b) = "
        f"({sol[0]:.4e}, {sol[1]:.4e})")
    say(f"  PR-15 ratio estimator gave 1.449 +- 0.082 — the four-parameter "
        f"model is the sharper route as pre-registered: "
        f"{'YES' if abs(R-1.357891) < 0.082 else 'NO'}")
    say("\n" + "=" * 92)
    OUT.close()
