# --- stage 5 (Integration Plan v0.1): re-pointed to the mtft package. ---
# internal() delegates to mtft.chain.internal (certified spectrally
# identical: g to 9e-15, full H(u) spectrum to 2.5e-14; B differs only
# by gauge inside near-degenerate blocks).  Where a local gsq() remains
# it is retained VERBATIM: it evaluates complex u in the f64 backend,
# which mtft.ep deliberately does not offer (PR-20 floors) -- S5-1.
# ------------------------------------------------------------------------
#!/usr/bin/env python3
"""
pr10_origins.py — PR-10: one origin for two constants, and the level test
=========================================================================
MIT License — (c) 2026 Roger Tano — MTFT Research Program
Engines: numpy/mpmath + PARI/GP 2.15.4 (mftraceform) + mtft 0.9.1.

PRE-REGISTERED (PR-9 note §5, before this code existed):
 (a) are the two geometry-free constants — 1.00459 (Pr Y, tau direction)
     and 0.9485 (Pr Z, tau_2 direction) — produced by ONE perturbative
     calculation?  If it produces only one, the mechanism is incomplete.
 (b) does the equidistribution of Pr AA sharpen with level (N = 187, 209)?

THE ONE CALCULATION.  Second-order perturbation of H = h - u B in u:
    eps_i(u) = g_i - u B_ii + u^2 S_i,   S_i = sum_{j!=i} |B_ij|^2/(g_i-g_j).
Both edges are eps_1(u_+) and eps_0(u_-), and with the second shell
u_+- = +-tau x_max + 3p tau_2, so with  D = S_1 - S_0,  Sg = S_1 + S_0:

    gap = m - s(mu_0+mu_1) + s^2 D + c(mu_0-mu_1) + 2 s c Sg + c^2 D,
    s = tau x_max,  c = 3p tau_2.

  c = 0  =>  ROT1 = s*(mu_0+mu_1)/m  with s* the root of the quadratic
             (depends on D ALONE);
  d/dc   =>  ROT2 = [1 + 2 s* Sg/(mu_0-mu_1)] / [1 - 2 s* D/(mu_0+mu_1)]
             (depends on D AND Sg).
Two constants, one pair of sums.  Falsifiable: compute S_0 and S_1 from
the internal data and both numbers must come out.

THE LEVEL TEST, with its pre-registration corrected in the open.  A
z-score is ALREADY normalized by 1/sqrt(d); under equidistribution it
stays O(1) while the RAW deviations shrink.  The PR-9 wording ("z-scores
shrinking like 1/sqrt(dim)") was wrong as written; the corrected test is
run here and the error is owned.

Gates: RG0 baselines; RG1 (PR-10a) both constants from one calculation;
RG2 cross-engine (PARI trace forms vs the corpus at level 143);
RG3 (PR-10b) three levels, corrected test; RG4 the rigidity finding.

Run:  py pr10_origins.py     (needs gp on PATH)
"""
from __future__ import annotations
import math, subprocess, json
import numpy as np
import mpmath as mp
from mtft.chain import internal as _chain_internal


BETA, KSTAR, NB = 2.0, 5.0, 60
ROT1_MEAS, ROT2_MEAS = 1.004590, 0.9485
REPORT = []


def rec(name, gtype, value, cls, ok, note=""):
    REPORT.append((name, gtype, value, cls, bool(ok), note))
    print(f"[{'PASS' if ok else 'FAIL'}] {name:<24} {gtype:<12} "
          f"{value:<34} {cls:<20} {note}", flush=True)


def internal(N=1600, kappa=KSTAR, nb=NB, gcap=200.0):
    ic = _chain_internal(kappa, nb=nb, backend="f64", N=N, gcap=gcap)
    return np.asarray(ic.g), np.asarray(ic.B)


G, B = internal()
MU0, MU1, MGAP = B[0, 0], B[1, 1], G[1]


def second_order(g, Bm, nb):
    S0 = sum(Bm[0, j] ** 2 / (g[0] - g[j]) for j in range(1, nb))
    S1 = (Bm[1, 0] ** 2 / (g[1] - g[0])
          + sum(Bm[1, j] ** 2 / (g[1] - g[j]) for j in range(2, nb)))
    return float(S0), float(S1)


def rots(S0, S1, mu0, mu1, m):
    D, Sg = S1 - S0, S1 + S0
    a, b = mu0 + mu1, mu0 - mu1
    disc = a ** 2 - 4 * m * D
    s = (a - math.sqrt(disc)) / (2 * D) if abs(D) > 1e-14 else m / a
    rot1 = s * a / m
    rot2 = (1 + 2 * s * Sg / b) / (1 - 2 * s * D / a)
    return s, rot1, rot2, D, Sg


# ------------------------------------------------------------------ RG0
def rg0():
    ok = abs(MU0 - 1.050398) < 1e-5 and abs(MGAP - 0.736839) < 1e-5
    rec("RG0 baselines", "Instrument", f"mu0={MU0:.6f}, m={MGAP:.6f}",
        "CERTIFIED", ok, f"targets: ROT1={ROT1_MEAS}, ROT2={ROT2_MEAS}")


# ------------------------------------------------------------------ RG1
def rg1():
    conv = {}
    for nb in (30, 45, 60):
        g, Bm = internal(nb=nb)
        S0, S1 = second_order(g, Bm, nb)
        conv[nb] = rots(S0, S1, Bm[0, 0], Bm[1, 1], g[1])
    s, r1, r2, D, Sg = conv[60]
    drift = max(abs(conv[60][1] - conv[45][1]), abs(conv[60][2] - conv[45][2]))
    e1 = abs(r1 / ROT1_MEAS - 1)
    e2 = abs(r2 / ROT2_MEAS - 1)
    ok = e1 < 0.01 and e2 < 0.03 and drift < 1e-3
    rec("RG1 one origin, two constants", "Theorem",
        f"ROT1 {r1:.5f} / ROT2 {r2:.5f}", "CERTIFIED(1%,3%)", ok,
        f"measured {ROT1_MEAS} / {ROT2_MEAS} (dev {e1:.4f} / {e2:.4f}); "
        f"from ONE pair S0={S0_:.6f}, S1={S1_:.6f} -> D={D:.6f}, "
        f"Sg={Sg:.6f}; s*={s:.5f}; nb-drift {drift:.1e}")


# ------------------------------------------------------------------ RG2
def pari_traces(levels=(143, 187, 209)):
    src = "default(parisize, 800000000);\n{\nNs = %s;\nfor(i = 1, %d,\n  NN = Ns[i];\n  tf = mftraceform([NN, 2], 0);\n  dd = mfdim([NN, 2], 0);\n  cc = mfcoefs(tf, 32);\n  print(NN, \" \", dd, \" \", cc[3], \" \", cc[5], \" \", cc[9], \" \", cc[17], \" \", cc[33]);\n);\n}\nquit;\n" % (list(levels), len(levels))
    open("tf_run.gp", "w").write(src)
    out = subprocess.run(["gp", "-q", "tf_run.gp"], capture_output=True,
                         text=True, timeout=1800).stdout
    res = {}
    for line in out.strip().splitlines():
        f = line.split()
        if len(f) == 7 and f[0].isdigit():
            res[int(f[0])] = (int(f[1]), [int(v) for v in f[2:]])
    return res


def rg2(pari):
    import mtft
    tt = mtft.TRACE_TOTALS_50
    corpus = (tt[0], [tt[1], tt[3], tt[7], tt[15], tt[31]])
    ok = pari[143] == corpus
    rec("RG2 cross-engine (corpus)", "Identity",
        f"PARI {pari[143][1]} vs corpus {corpus[1]}", "EXACT", ok,
        f"dim {pari[143][0]} = {corpus[0]}; the shipped TRACE_TOTALS_50 "
        f"reproduced independently by PARI mftraceform")


# ------------------------------------------------------------------ RG3
def power_sums(d, T):
    p1 = T[0]
    p2 = T[1] + 2 * d
    p3 = T[2] + 4 * p1
    p4 = T[3] + 6 * p2 - 4 * d
    p5 = T[4] + 8 * p3 - 12 * p1
    return [p1, p2, p3, p4, p5]


def kesten_m(p=2, kmax=10, n=900):
    t, w = np.polynomial.legendre.leggauss(n)
    th = 0.5 * math.pi * t
    xm, q = 2 * math.sqrt(p), p + 1
    x = xm * np.sin(th)
    dens = q * xm ** 2 * np.cos(th) ** 2 / (2 * math.pi * (q ** 2 - x ** 2))
    wt = w * 0.5 * math.pi * dens
    return {k: float(np.sum(wt * x ** k)) for k in range(1, kmax + 1)}


def rg3(pari):
    mk = kesten_m()
    rows, allz = [], []
    for N in (143, 187, 209):
        d, T = pari[N]
        ps = power_sums(d, T)
        z = []
        for k in range(1, 6):
            sd = math.sqrt(max(mk[2 * k] - mk[k] ** 2, 0.0) / d)
            z.append((ps[k - 1] / d - mk[k]) / sd)
        rows.append((N, d, ps, z, max(abs(v) for v in z),
                     abs(ps[1] / d - mk[2])))
        allz += z

    ok = all(r[4] < 2.0 for r in rows)
    rec("RG3 level test (corrected)", "Theorem",
        f"max|z| = " + "/".join(f"{r[4]:.2f}" for r in rows),
        "CERTIFIED(2 sigma)", ok,
        f"N=143/187/209, dim {[r[1] for r in rows]}; power sums "
        f"{[r[2] for r in rows]}; PRE-REGISTRATION CORRECTED IN THE OPEN: "
        f"a z-score is already 1/sqrt(d)-normalized, so under "
        f"equidistribution it stays O(1) while raw deviations shrink — "
        f"the PR-9 wording was wrong and is owned")
    return allz


# ------------------------------------------------------------------ RG4
def rg4(allz):
    rms = math.sqrt(sum(v * v for v in allz) / len(allz))
    mx = max(abs(v) for v in allz)
    # probability that 15 iid N(0,1) all fall inside |z| < mx
    from math import erf, sqrt
    pin = erf(mx / sqrt(2))
    prob = pin ** len(allz)
    ok = rms < 0.6
    rec("RG4 arithmetic rigidity", "Structural",
        f"RMS z = {rms:.3f} over {len(allz)} moments", "MEASURED", ok,
        f"iid model predicts 1.000; max|z| = {mx:.2f}, and P(all "
        f"{len(allz)} inside) = {prob:.1e} under iid — the arithmetic "
        f"spectrum is ~{1/rms:.1f}x MORE regular than random, consistent "
        f"with trace-formula error terms (elliptic/hyperbolic, bounded) "
        f"rather than random fluctuation. Kimi AE.5-F2 anticipated the "
        f"iid caveat; this quantifies it")


if __name__ == "__main__":
    print("=" * 108)
    print("  PR-10 — ONE ORIGIN FOR TWO CONSTANTS, AND THE LEVEL TEST")
    print("=" * 108)
    rg0()
    g_, B_ = internal()
    S0_, S1_ = second_order(g_, B_, NB)
    rg1()
    pari = pari_traces()
    rg2(pari)
    z = rg3(pari)
    rg4(z)
    print("-" * 108)
    n = sum(1 for x in REPORT if x[4])
    print(f"  {n}/{len(REPORT)} gates green")
    print("=" * 108)
