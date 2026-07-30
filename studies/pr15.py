# --- stage 5 (Integration Plan v0.1): re-pointed to the mtft package. ---
# internal() delegates to mtft.chain.internal (certified spectrally
# identical: g to 9e-15, full H(u) spectrum to 2.5e-14; B differs only
# by gauge inside near-degenerate blocks).  Where a local gsq() remains
# it is retained VERBATIM: it evaluates complex u in the f64 backend,
# which mtft.ep deliberately does not offer (PR-20 floors) -- S5-1.
# ------------------------------------------------------------------------
#!/usr/bin/env python3
"""
pr15.py — PR-15: R(kappa) by certified census, and a contour-free radius
========================================================================
MIT License — (c) 2026 Roger Tano — MTFT Research Program
Prints to stdout AND pr15_run.txt (Addendum AL procedural note).

CARRIED ITEMS (PR-14 v0.1.2 §6, all pre-registered):
 (1) R(kappa) via the staircase localizer inside the nb-certified scope
     r <~ 4.8, WITH LEVEL IDENTIFICATION: only level-0/1 collisions
     count toward the radius (standing rule 8).
 (2) Cauchy-Hadamard extension: an independent, CONTOUR-FREE estimate of
     R at kappa = 5, to test the certified 1.35789.  Coefficients are
     taken by Cauchy integral on a circle of radius rho < R (spectrally
     accurate, no polynomial fitting), NOT by least-squares fits.
     Note the singularity is a CONJUGATE PAIR at arg = +-2.799 rad, so
     c_k oscillates; estimators must average over the oscillation.
 (3) audit-trail re-examination for unpinned-number "confirmations"
     (standing rule 9) — tabulated in the note.

Method note for (1): the staircase brackets moduli without Newton
basins, then Newton is seeded INSIDE a bracket whose occupancy the
census already fixed — search is legitimate once the census says how
many roots are there.
"""
from __future__ import annotations
import math
import numpy as np
import mpmath as mp
from mtft.chain import internal as _chain_internal


BETA = 2.0
OUT = open("pr15_run.txt", "w", buffering=1)


def say(s=""):
    print(s, flush=True)
    OUT.write(s + "\n")


def internal(N=1600, kappa=5.0, nb=30, gcap=200.0):
    ic = _chain_internal(kappa, nb=nb, backend="f64", N=N, gcap=gcap)
    return np.asarray(ic.g), np.asarray(ic.B)


def count(r, g, B, npts=400):
    th = 2 * math.pi * np.arange(npts) / npts
    tot = 0j
    for t in th:
        u = r * math.cos(t) + 1j * r * math.sin(t)
        ev, V = np.linalg.eig(np.diag(g) - u * B)
        nrm = np.einsum('ij,ij->j', V, V)
        dep = -np.einsum('ij,jk,ki->i', V.T, B, V) / nrm
        d = ev[:, None] - ev[None, :]
        dd = dep[:, None] - dep[None, :]
        np.fill_diagonal(d, 1.0); np.fill_diagonal(dd, 0.0)
        tot += np.sum(dd / d) * 1j * u
    return (tot * (2 * math.pi / npts) / (2j * math.pi)).real


# S5-1: retained verbatim -- complex-u f64 gsq predates mtft.ep's
# real-axis diabatic-centre design; frozen so this study's record stands.
def gsq(u, g, B):
    ev = np.linalg.eigvals(np.diag(g) - u * B)
    ev = ev[np.argsort(ev.real)]
    d = ev[1:] - ev[:-1]
    return d[int(np.argmin(np.abs(d)))] ** 2


def newton(u0, g, B, iters=70, h=1e-7):
    u = complex(u0)
    for _ in range(iters):
        f = gsq(u, g, B)
        fp = (gsq(u + h, g, B) - gsq(u - h, g, B)) / (2 * h)
        if fp == 0:
            break
        s = f / fp
        u -= s
        if abs(s) < 1e-14:
            break
    return u


def levels_of(u, g, B, steps=400):
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


def R_of_kappa(kap):
    g, B = internal(kappa=kap)
    rs = np.arange(0.4, 3.01, 0.10)
    prev, brackets = None, []
    for r in rs:
        c = count(r, g, B)
        if prev is not None and round(c) > round(prev):
            brackets.append((r - 0.10, r, round(c) - round(prev)))
        prev = c
    for lo, hi, jump in brackets:
        roots = []
        for th in np.linspace(0, 2 * math.pi, 40, endpoint=False):
            for rr in (lo + 0.25 * (hi - lo), 0.5 * (lo + hi),
                       lo + 0.75 * (hi - lo)):
                u = newton(rr * np.exp(1j * th), g, B)
                if (lo - 0.05 < abs(u) < hi + 0.05
                        and abs(gsq(u, g, B)) < 1e-12
                        and all(abs(u - v) > 1e-6 for v in roots)):
                    roots.append(u)
        for u in sorted(roots, key=abs):
            lv = levels_of(u, g, B)
            if lv == [0, 1]:
                return abs(u), lv, brackets
    return float('nan'), None, brackets


def taylor_contour(kap=5.0, rho=1.0, npts=512, kmax=40):
    """c_k by Cauchy integral on |u| = rho < R: spectrally accurate,
    no fitting.  eps_1 tracked by continuation around the circle."""
    g, B = internal(kappa=kap)
    th = 2 * math.pi * np.arange(npts) / npts
    vals = np.empty(npts, dtype=complex)
    ev0 = np.linalg.eigvals(np.diag(g) - rho * B)
    cur = np.sort(ev0.real)[1] + 0j
    for j, t in enumerate(th):
        u = rho * np.exp(1j * t)
        ev = np.linalg.eigvals(np.diag(g) - u * B)
        cur = ev[int(np.argmin(np.abs(ev - cur)))]
        vals[j] = cur
    cs = []
    for k in range(kmax + 1):
        cs.append(np.sum(vals * np.exp(-1j * k * th)) / npts / rho ** k)
    return np.array(cs)


if __name__ == "__main__":
    say("=" * 96)
    say("  PR-15 — R(kappa) BY CERTIFIED CENSUS, AND A CONTOUR-FREE RADIUS")
    say("=" * 96)

    say("\n[1] R(kappa): staircase brackets -> Newton inside a bracket "
        "of known occupancy -> level ID")
    say(f"{'kappa':>6} {'R (levels 0-1)':>16} {'first brackets':>34}")
    rows = []
    for kap in (3.0, 4.0, 5.0, 6.0, 8.0):
        R, lv, br = R_of_kappa(kap)
        rows.append((kap, R))
        bs = "; ".join(f"({a:.2f},{b:.2f})+{j}" for a, b, j in br[:3])
        say(f"{kap:6.1f} {R:16.5f}   {bs}")
    say(f"  kappa=5 check vs certified 1.35789: "
        f"{abs(rows[2][1] - 1.35789):.5f}")

    say("\n[2] contour-free radius at kappa = 5 (Cauchy coefficients, "
        "rho = 1.0)")
    cs = taylor_contour()
    say(f"  c_2 = {cs[2].real:+.6f} (RSPT S_1 = -0.005790 expected)")
    say(f"  c_3 = {cs[3].real:+.6f} (RSPT T_1 = -0.003084 expected)")
    ests = []
    for m in (8, 12, 16):
        for k in (16, 20, 24):
            if k + m <= 40 and abs(cs[k + m]) > 0:
                ests.append((abs(cs[k]) / abs(cs[k + m])) ** (1.0 / m))
    ests = np.array(ests)
    say(f"  |c_k/c_(k+m)|^(1/m) over k=16,20,24 x m=8,12,16: "
        f"mean {ests.mean():.5f}, sd {ests.std():.5f}, "
        f"range [{ests.min():.4f}, {ests.max():.4f}]")
    say(f"  certified census value: 1.35789  ->  deviation of mean "
        f"{abs(ests.mean() - 1.35789):.5f}")
    say("\n" + "=" * 96)
    OUT.close()
