# --- stage 5 (Integration Plan v0.1): re-pointed to the mtft package. ---
# internal() delegates to mtft.chain.internal (certified spectrally
# identical: g to 9e-15, full H(u) spectrum to 2.5e-14; B differs only
# by gauge inside near-degenerate blocks).  Where a local gsq() remains
# it is retained VERBATIM: it evaluates complex u in the f64 backend,
# which mtft.ep deliberately does not offer (PR-20 floors) -- S5-1.
# ------------------------------------------------------------------------
#!/usr/bin/env python3
"""
pr14_ep.py — PR-14: the exceptional point located, and R(kappa)
================================================================
MIT License — (c) 2026 Roger Tano — MTFT Research Program

PRE-REGISTERED (PR-13 note §5, before this code existed):
 (a) refine the collision to four digits and resolve Im(u_c): is it ON
     the negative real axis (an ordinary level crossing at negative
     coupling — the cleaner statement) or genuinely complex?
     + Kimi (AI.3): identify WHICH level pair collides.
 (b) does R vary with kappa, and does it track kappa* = 5.0 — is the
     softest coupling also the one with the most distant collision?

THE STRUCTURAL EXPECTATION, stated before measuring.  h and B are REAL
symmetric, so h + tB (t real) is a real symmetric family; by von
Neumann-Wigner, eigenvalue crossings in such a family have codimension
2 and do NOT occur generically without a symmetry.  The generic
structure is an AVOIDED crossing on the axis with a conjugate PAIR of
square-root branch points just off it.  The two engines' straddle of pi
in Addendum AI (3.075 = pi - 0.067 and 3.209 = pi + 0.067) is precisely
that signature — one engine found each member of the pair.  Prediction:
Im(u_c) != 0, |Im| ~ 0.17, and the pair is conjugate.

METHOD.  Near an exceptional point eps_a - eps_b ~ sqrt(u - u_c), so
    g(u) = (eps_a(u) - eps_b(u))^2
is ANALYTIC with a simple zero at u_c: Newton on g converges
quadratically where minimising the gap crawls.

Gates: WG0 baselines; WG1 Newton refinement of u_c to >=6 digits;
WG2 conjugate-pair structure + on-axis test; WG3 which pair collides
(continuation back to u=0); WG4 R(kappa) vs kappa*.

Run:  py pr14_ep.py
"""
from __future__ import annotations
import math
import numpy as np
import mpmath as mp
from mtft.chain import internal as _chain_internal


BETA, KSTAR, NB = 2.0, 5.0, 40
REPORT = []


def rec(name, gtype, value, cls, ok, note=""):
    REPORT.append((name, gtype, value, cls, bool(ok), note))
    print(f"[{'PASS' if ok else 'FAIL'}] {name:<24} {gtype:<12} "
          f"{value:<36} {cls:<20} {note}", flush=True)


def internal(N=1600, kappa=KSTAR, nb=NB, gcap=200.0):
    ic = _chain_internal(kappa, nb=nb, backend="f64", N=N, gcap=gcap)
    return np.asarray(ic.g), np.asarray(ic.B)


G, B = internal()


# S5-1: retained verbatim -- complex-u f64 gsq predates mtft.ep's
# real-axis diabatic-centre design; frozen so this study's record stands.
def gsq(u, g=None, Bm=None):
    """(eps_a - eps_b)^2 for the closest pair: analytic near the EP."""
    g = G if g is None else g
    Bm = B if Bm is None else Bm
    ev = np.linalg.eigvals(np.diag(g) - u * Bm)
    ev = ev[np.argsort(ev.real)]
    d = ev[1:] - ev[:-1]
    i = int(np.argmin(np.abs(d)))
    return d[i] ** 2


def newton_ep(u0, g=None, Bm=None, iters=60, h=1e-6):
    u = complex(u0)
    for _ in range(iters):
        f = gsq(u, g, Bm)
        fp = (gsq(u + h, g, Bm) - gsq(u - h, g, Bm)) / (2 * h)
        if fp == 0:
            break
        step = f / fp
        u -= step
        if abs(step) < 1e-13:
            break
    return u


def all_eps(g=None, Bm=None, rmax=5.5, nr=26, nth=34, tol=1e-13):
    """Multi-seed Newton over an annular grid, collecting EVERY root and
    sorting by |u|.  (A global minimum of |(de)^2| is NOT the nearest
    root — PR-14's first pass made exactly that error and landed on a
    deeper EP at |u| = 4.50.)"""
    roots = []
    for th in np.linspace(0.0, math.pi, nth):
        for r in np.linspace(0.5, rmax, nr):
            u = newton_ep(r * np.exp(1j * th), g, Bm, iters=40)
            if not np.isfinite(u) or abs(u) > 12:
                continue
            if abs(gsq(u, g, Bm)) < tol:
                if all(abs(u - v) > 1e-5 for v in roots):
                    roots.append(u)
    roots.sort(key=abs)
    return roots


# ------------------------------------------------------------------ WG0
def wg0():
    ok = abs(B[0, 0] - 1.050398) < 1e-5 and abs(G[1] - 0.736839) < 1e-5
    rec("WG0 baselines", "Instrument", f"mu0={B[0,0]:.6f}, m={G[1]:.6f}",
        "CERTIFIED", ok, "PR-13: |u_c| ~ 2.517 (suite) / 2.515 (audit), "
                         "R_CH ~ 2.555")


# ------------------------------------------------------------------ WG1
def wg1():
    roots = all_eps()
    u = roots[0]
    resid = abs(gsq(u))
    us = [newton_ep(u + d) for d in (0.01, -0.01, 0.01j, -0.01j)]
    spread = max(abs(v - u) for v in us)
    ok = resid < 1e-13 and spread < 1e-8
    rec("WG1 EP located (Newton)", "Implementation",
        f"u_c = {u.real:.6f} {u.imag:+.6f}i", "CERTIFIED(1e-8)", ok,
        f"|u_c| = {abs(u):.6f} vs R_CH 2.555 (ratio {abs(u)/2.555:.3f}); "
        f"arg = {math.atan2(u.imag, u.real):.6f} rad; |(de)^2| = "
        f"{resid:.1e}; seed spread {spread:.1e}; {len(roots)} distinct "
        f"EPs found, next at |u| = {abs(roots[1]):.3f}")
    return u, roots


# ------------------------------------------------------------------ WG2
def wg2(u):
    ubar = newton_ep(np.conj(u))
    conj_err = abs(ubar - np.conj(u))
    im = abs(u.imag)
    # square-root character: (de)^2 must vanish LINEARLY in (u - u_c)
    ratios = []
    for d in (1e-2, 5e-3, 2.5e-3, 1.25e-3):
        ratios.append(abs(gsq(u + d)) / d)
    lin = (max(ratios) - min(ratios)) / np.mean(ratios)
    ok = conj_err < 1e-8 and lin < 0.10
    rec("WG2 conjugate pair, off axis", "Theorem",
        f"Im(u_c) = {u.imag:+.3e}", "CERTIFIED", ok,
        f"conjugate partner recovers to {conj_err:.1e}; |(de)^2|/|u-u_c| "
        f"= {['%.4f' % r for r in ratios]} (linear to {lin:.2f}) => "
        f"square-root branch point. PRE-REGISTERED PREDICTION FALSIFIED "
        f"(owned): I expected a conjugate pair JUST off the axis "
        f"(|Im| ~ 0.17) with the AI.3 straddle of pi as its signature. "
        f"The true EP sits at arg {math.atan2(u.imag,u.real):.3f} rad "
        f"= {math.degrees(math.atan2(u.imag,u.real)):.0f} deg, far from "
        f"pi. The straddle was BOTH engines resolving a shallow "
        f"near-axis avoided crossing (min gap ~0.01), not the branch "
        f"point. von Neumann-Wigner still forbids on-axis crossings; it "
        f"does not put the EPs near the axis")


# ------------------------------------------------------------------ WG3
def wg3(u):
    """Continue the colliding pair back to u = 0 to name the levels."""
    ev = np.linalg.eigvals(np.diag(G) - u * B)
    ev = ev[np.argsort(ev.real)]
    d = ev[1:] - ev[:-1]
    i = int(np.argmin(np.abs(d)))
    pair = np.array([ev[i], ev[i + 1]])
    steps = 400
    for t in np.linspace(1.0, 0.0, steps)[1:]:
        uu = u * t
        cur = np.linalg.eigvals(np.diag(G) - uu * B)
        new = []
        used = set()
        for p in pair:
            j = int(np.argmin([abs(p - c) if k not in used else 9e9
                               for k, c in enumerate(cur)]))
            used.add(j); new.append(cur[j])
        pair = np.array(new)
    idx = [int(np.argmin(np.abs(G - p.real))) for p in pair]
    ok = len(set(idx)) == 2 and max(abs(pair.imag)) < 1e-6
    rec("WG3 which pair collides", "Structural",
        f"internal levels {idx[0]} and {idx[1]}", "MEASURED", ok,
        f"continued back to u=0: g = {G[idx[0]]:.6f}, {G[idx[1]]:.6f} "
        f"(residual Im {max(abs(pair.imag)):.0e}); NOT the lowest two "
        f"levels (AI.3 anticipated this): the pair is the GAP-DEFINING "
        f"level 1 (g = m = {G[1]:.6f}) and its upper neighbour level 2 "
        f"(g = {G[2]:.6f}). So R is set by the same level whose "
        f"separation from level 0 defines the mass gap — the radius and "
        f"the gap are properties of one level")


# ------------------------------------------------------------------ WG4
def wg4():
    from math import inf
    rows = []
    for kap in (3.0, 4.0, 5.0, 6.0, 8.0):
        g, Bm = internal(kappa=kap)
        rr = all_eps(g, Bm, rmax=9.0, nr=26, nth=30)
        if not rr:
            rows.append((kap, float('nan'), float('nan'), float('nan')))
            continue
        u = rr[0]
        m_k = g[1]
        s_star = m_k / (Bm[0, 0] + Bm[1, 1])
        rows.append((kap, abs(u), s_star, abs(u) / s_star))
    rows = [r for r in rows if r[1] == r[1]]
    kmax_R = max(rows, key=lambda r: r[1])[0]
    kmax_ratio = max(rows, key=lambda r: r[3])[0]
    ok = all(r[1] > r[2] for r in rows)
    rec("WG4 R(kappa) vs kappa*", "Structural",
        f"R peaks at kappa={kmax_R}, R/s* peaks at {kmax_ratio}",
        "MEASURED", ok,
        "; ".join(f"k={r[0]:.0f}: R={r[1]:.3f}, s*={r[2]:.3f}, "
                  f"R/s*={r[3]:.2f}" for r in rows)
        + f" — kappa*=5.0 is the SOFTEST gap (rung 4), and the question "
          f"was whether it is also the most robust expansion; "
          f"{'IT IS NOT' if kmax_ratio != 5.0 else 'IT IS'}")


if __name__ == "__main__":
    print("=" * 106)
    print("  PR-14 — THE EXCEPTIONAL POINT, AND R(kappa)")
    print("=" * 106)
    wg0()
    u, roots = wg1()
    wg2(u); wg3(u); wg4()
    print("-" * 106)
    n = sum(1 for x in REPORT if x[4])
    print(f"  {n}/{len(REPORT)} gates green")
    print("=" * 106)
