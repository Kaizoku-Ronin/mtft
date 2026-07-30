# --- stage 5 (Integration Plan v0.1): re-pointed to the mtft package. ---
# internal() delegates to mtft.chain.internal (certified spectrally
# identical: g to 9e-15, full H(u) spectrum to 2.5e-14; B differs only
# by gauge inside near-degenerate blocks).  Where a local gsq() remains
# it is retained VERBATIM: it evaluates complex u in the f64 backend,
# which mtft.ep deliberately does not offer (PR-20 floors) -- S5-1.
# ------------------------------------------------------------------------
#!/usr/bin/env python3
"""
pr19.py — PR-19: the correction rate is pairwise, and in closed form
=====================================================================
MIT License — (c) 2026 Roger Tano — MTFT Research Program
Prints to stdout AND pr19_run.txt.

PRE-REGISTERED (PR-18 note §6, Addendum AP.4): is the finite-kappa
correction rate GLOBAL or PAIRWISE?  AP measures -0.249 for pair (0,1)
(settled) and -0.172/-0.190 for pair (1,2) (still transient), and
recommends extending both to kappa = 64 before declaring.

INSTRUMENT UPGRADE (the reason this rung is possible here).  PR-18 §1
declined to measure overshoots because winding-bisection alone floors at
~2e-4.  Fixed: the winding CERTIFIES occupancy and level identity, then
NEWTON on the analytic square g(u) = (eps_a - eps_b)^2 refines the same
root to machine precision.  Certification and precision from different
instruments, which is the correct division of labour.

THE PREDICTION, derived before measuring.  In the diagonal limit the
levels are the sites n = 2, 3, 4, ... and the Mellin coupling between
sites n and n+1 is (n/(n+1))^kappa.  The crossing of levels (i, i+1)
sits at u = -(g_{i+1}-g_i)/(e^{-g_{i+1}} - e^{-g_i}), so its finite-
kappa correction is driven by whatever perturbs those two gaps SLOWEST.
Second-order repulsion inside the pair decays as ((i+2)/(i+3))^{2kappa};
repulsion of the UPPER member against the level above it decays as
((i+3)/(i+4))^{2kappa} — which is SLOWER, hence dominant.  Therefore

    rate(i, i+1) = 2 log10( (i+3)/(i+4) )   per unit kappa,

    pair (0,1): 2 log10(3/4) = -0.2498775   [AP measures -0.249]
    pair (1,2): 2 log10(4/5) = -0.1938200   [AP measures -0.190]
    pair (2,3): 2 log10(5/6) = -0.1583625   [= RUNG 4's (5/6)^{2kappa}]

The third line is the punchline: rung 4's bridge-contraction law is the
i = 2 member of this family.  PAIRWISE, not global, with a closed form
per pair — and PR-18(a)'s withdrawn guess was the right law applied to
the wrong pair.
"""
from __future__ import annotations
import math
import numpy as np
import mpmath as mp
from mtft.chain import internal as _chain_internal


BETA = 2.0
OUT = open("pr19_run.txt", "w", buffering=1)


def say(s=""):
    print(s, flush=True)
    OUT.write(s + "\n")


def internal(N=1600, kappa=5.0, nb=30, gcap=200.0):
    ic = _chain_internal(kappa, nb=nb, backend="f64", N=N, gcap=gcap)
    return np.asarray(ic.g), np.asarray(ic.B)


# S5-1: retained verbatim -- complex-u f64 gsq predates mtft.ep's
# real-axis diabatic-centre design; frozen so this study's record stands.
def gsq(u, g, B):
    ev = np.linalg.eigvals(np.diag(g) - u * B)
    ev = ev[np.argsort(ev.real)]
    d = ev[1:] - ev[:-1]
    return d[int(np.argmin(np.abs(d)))] ** 2


def newton(u0, g, B, iters=90, h=1e-8):
    u = complex(u0)
    for _ in range(iters):
        f = gsq(u, g, B)
        fp = (gsq(u + h, g, B) - gsq(u - h, g, B)) / (2 * h)
        if fp == 0:
            break
        s = f / fp
        u -= s
        if abs(s) < 1e-15:
            break
    return u


def levels_of(u, g, B, steps=300):
    ev = np.linalg.eigvals(np.diag(g) - u * B)
    ev = ev[np.argsort(ev.real)]
    d = ev[1:] - ev[:-1]
    i = int(np.argmin(np.abs(d)))
    pair = np.array([ev[i], ev[i + 1]])
    for t in np.linspace(1.0, 0.0, steps)[1:]:
        cur = np.linalg.eigvals(np.diag(g) - u * t * B)
        new, used = [], set()
        for p in pair:
            j = int(np.argmin([abs(p - c) if k not in used else 9e9
                               for k, c in enumerate(cur)]))
            used.add(j); new.append(cur[j])
        pair = np.array(new)
    return sorted(int(np.argmin(np.abs(g - p.real))) for p in pair)


def limit_u(g, i):
    a, b = g[i], g[i + 1]
    return abs((b - a) / (math.exp(-b) - math.exp(-a)))


if __name__ == "__main__":
    say("=" * 96)
    say("  PR-19 — THE CORRECTION RATE IS PAIRWISE, WITH A CLOSED FORM")
    say("=" * 96)
    say("\n  predicted rate(i,i+1) = 2 log10((i+3)/(i+4)) per unit kappa:")
    for i in (0, 1, 2):
        say(f"    pair ({i},{i+1}) <-> sites ({i+2},{i+3}):  "
            f"{2*math.log10((i+3)/(i+4)):+.7f}"
            + ("   [= rung 4's (5/6)^2k]" if i == 2 else ""))

    gL, _ = internal(kappa=200.0, nb=30)
    lim = {i: limit_u(gL, i) for i in (0, 1, 2)}
    say(f"\n  closed-form limits (kappa=200 gaps): "
        + ", ".join(f"u_{i}{i+1}={lim[i]:.7f}" for i in (0, 1, 2)))

    say(f"\n{'pair':>6} {'kappa':>6} {'|u| (Newton)':>16} "
        f"{'|u|-limit':>13} {'slope/kappa':>12}")
    rates = {}
    for i in (0, 1, 2):
        prev = None
        for kap in (32.0, 48.0, 64.0):
            g, B = internal(kappa=kap)
            L = limit_u(g, i)
            u = newton(-L * (1 + 0.002j), g, B)
            ok = abs(gsq(u, g, B)) < 1e-13 and levels_of(u, g, B) == [i, i + 1]
            dev = abs(u) - lim[i]
            sl = ""
            if prev is not None and dev != 0 and prev[1] != 0:
                s = math.log10(abs(dev) / abs(prev[1])) / (kap - prev[0])
                sl = f"{s:+.6f}"
                rates[(i, kap)] = s
            say(f"{f'({i},{i+1})':>6} {kap:6.1f} {abs(u):16.10f} "
                f"{dev:+13.3e} {sl:>12}" + ("" if ok else "  [CHECK]"))
            prev = (kap, dev)

    say(f"\n  measured vs predicted rate (kappa 48->64):")
    for i in (0, 1, 2):
        p = 2 * math.log10((i + 3) / (i + 4))
        m = rates.get((i, 64.0))
        if m is not None:
            say(f"    pair ({i},{i+1}): measured {m:+.6f}  predicted "
                f"{p:+.6f}  dev {abs(m-p):.6f}")
    say("\n" + "=" * 96)
    OUT.close()
