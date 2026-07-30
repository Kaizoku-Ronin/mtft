# --- stage 5 (Integration Plan v0.1): re-pointed to the mtft package. ---
# internal() delegates to mtft.chain.internal (certified spectrally
# identical: g to 9e-15, full H(u) spectrum to 2.5e-14; B differs only
# by gauge inside near-degenerate blocks).  Where a local gsq() remains
# it is retained VERBATIM: it evaluates complex u in the f64 backend,
# which mtft.ep deliberately does not offer (PR-20 floors) -- S5-1.
# ------------------------------------------------------------------------
#!/usr/bin/env python3
"""
pr17.py — PR-17: the tracker fix, and a closed form for R(infinity)
====================================================================
MIT License — (c) 2026 Roger Tano — MTFT Research Program
Prints to stdout AND pr17_run.txt.

PRE-REGISTERED (PR-16 note §6):
 (a) fix the kappa = 8 tracker dropout (two-sheet subspace tracking);
 (b) does R(kappa) approach a finite limit, and is it related to rung
     4's strong-coupling constant m_inf = log(27 ln2/(8 ln3))?
 (c) extend the four-parameter fit to kappa != 5 with k-window analysis
     (Addendum AN-F1).

THE PREDICTION FOR (b), DERIVED BEFORE MEASURING.  Rung 4's Pr G:
as kappa -> infinity, M -> I, so T -> diag(rho-hat) and the dressed
hopping B = e^{-h/2} M e^{-h/2} -> diag(e^{-g}).  Then

    H(u) = diag(g_i - u e^{-g_i})     — DIAGONAL,

so the levels are LINEAR in u and levels 0,1 meet at a genuine REAL
crossing (permitted: the limit is diagonal, hence maximally symmetric):

    u_cross = -(g_1 - g_0)/(e^{-g_1} - e^{-g_0}) = -m/(1 - e^{-m}).

Hence the closed form, in terms of the modular anchor alone:

    R(infinity) = m_inf/(1 - e^{-m_inf}),   e^{-m_inf} = 8 ln3/(27 ln2)
    R(infinity)/m_inf = 27 ln2/(27 ln2 - 8 ln3)   [dimensionless]

This is a bridge: the ANALYTIC structure (radius of convergence) is
fixed by the MODULAR constant (rung-4 strong-coupling gap).  Falsifiable
— R(kappa) must climb to it and stop.
"""
from __future__ import annotations
import math
import numpy as np
import mpmath as mp
from mtft.chain import internal as _chain_internal


BETA = 2.0
OUT = open("pr17_run.txt", "w", buffering=1)


def say(s=""):
    print(s, flush=True)
    OUT.write(s + "\n")


def internal(N=1600, kappa=5.0, nb=30, gcap=200.0):
    ic = _chain_internal(kappa, nb=nb, backend="f64", N=N, gcap=gcap)
    return np.asarray(ic.g), np.asarray(ic.B)


def pair_winding(r, g, B, npts=900):
    """(a) tracker fix: match against a LINEAR EXTRAPOLATION of the two
    tracked sheets rather than their last value, so fast-moving levels
    and third-level crowding no longer steal the assignment."""
    th = 2 * math.pi * np.arange(npts + 1) / npts
    ev = np.linalg.eigvalsh(np.diag(g) - r * B)
    a, b = complex(ev[0]), complex(ev[1])
    ap, bp = a, b
    tot, qprev = 0.0, b - a
    for t in th[1:]:
        u = r * np.exp(1j * t)
        cur = np.linalg.eigvals(np.diag(g) - u * B)
        pa, pb = 2 * a - ap, 2 * b - bp          # linear prediction
        ia = int(np.argmin(np.abs(cur - pa)))
        rest = np.delete(cur, ia)
        ib = int(np.argmin(np.abs(rest - pb)))
        ap, bp = a, b
        a, b = cur[ia], rest[ib]
        q = b - a
        tot += np.angle(q / qprev)
        qprev = q
    return 2.0 * tot / (2 * math.pi)


def R_of(kap, lo=0.30, hi=2.60, nb=30, npts=900, iters=40):
    g, B = internal(kappa=kap, nb=nb)
    if pair_winding(lo, g, B, npts) > 1.0 or pair_winding(hi, g, B, npts) < 1.0:
        return float('nan')
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        if pair_winding(mid, g, B, npts) < 1.0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


if __name__ == "__main__":
    say("=" * 92)
    say("  PR-17 — TRACKER FIX, AND A CLOSED FORM FOR R(infinity)")
    say("=" * 92)

    with mp.workdps(30):
        m_inf = float(mp.log(27 * mp.log(2) / (8 * mp.log(3))))
        emi = float(8 * mp.log(3) / (27 * mp.log(2)))
    R_inf = m_inf / (1 - emi)
    say(f"\n[b] prediction, derived before measuring:")
    say(f"    m_inf                = {m_inf:.7f}   (rung 4, Pr G)")
    say(f"    e^-m_inf = 8ln3/27ln2= {emi:.7f}")
    say(f"    R(inf) = m/(1-e^-m)  = {R_inf:.7f}")
    say(f"    R(inf)/m_inf         = {1/(1-emi):.7f} "
        f"= 27ln2/(27ln2-8ln3)")

    say(f"\n[a]+[b] R(kappa) with the fixed tracker")
    say(f"{'kappa':>7} {'R':>13} {'R_inf - R':>13} {'m(k)/(1-e^-m(k))':>18}")
    prev = None
    for kap in (5.0, 8.0, 12.0, 16.0, 24.0):
        g, B = internal(kappa=kap)
        R = R_of(kap)
        mk = g[1]
        loc = mk / (1 - math.exp(-mk))
        say(f"{kap:7.1f} {R:13.6f} {R_inf - R:13.6f} {loc:18.6f}")
        prev = R
    say(f"\n  ledger values for reference: R(5)=1.3578912, R(8)=1.4244407")
    say("\n" + "=" * 92)
    OUT.close()
