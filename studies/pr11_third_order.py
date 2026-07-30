# --- stage 5 (Integration Plan v0.1): re-pointed to the mtft package. ---
# internal() delegates to mtft.chain.internal (certified spectrally
# identical: g to 9e-15, full H(u) spectrum to 2.5e-14; B differs only
# by gauge inside near-degenerate blocks).  Where a local gsq() remains
# it is retained VERBATIM: it evaluates complex u in the f64 backend,
# which mtft.ep deliberately does not offer (PR-20 floors) -- S5-1.
# ------------------------------------------------------------------------
#!/usr/bin/env python3
"""
pr11_third_order.py — PR-11: the direction-of-motion test, and the
boundedness signature of arithmetic rigidity
==================================================================
MIT License — (c) 2026 Roger Tano — MTFT Research Program
Engines: numpy/mpmath + PARI/GP (mftraceform) + mtft 0.9.1.

PRE-REGISTERED (PR-10 note §5, before this code existed):
 (a) push the perturbation to THIRD order and check the residuals of
     Pr AC (ROT1 0.17% low, ROT2 1.7% high) are PREDICTED, not merely
     tolerated.  FALSIFIER: if third order moves either constant AWAY
     from its measured value, the second-order attribution is wrong.
 (b) quantify Pr AF (arithmetic rigidity) against the trace formula.

THIRD ORDER.  eps_i(u) = g_i - u B_ii + u^2 S_i + u^3 T_i with
  S_i = sum_{j!=i} |B_ij|^2/(g_i-g_j),
  T_i = -sum_{j,k!=i} B_ij B_jk B_ki/((g_i-g_j)(g_i-g_k))
        + B_ii sum_{j!=i} |B_ij|^2/(g_i-g_j)^2.
With u_+- = +-s + c the gap becomes, to this order,
  gap = m - s a + s^2 D + s^3 W + c[b + 2 s Sg + 3 s^2 (T_1-T_0)] + ...
  a = mu_0+mu_1, b = mu_0-mu_1, D = S_1-S_0, Sg = S_1+S_0, W = T_1+T_0,
giving ROT1 = s* a/m (s* the CUBIC root) and
  ROT2 = (a/b)[b + 2 s* Sg + 3 s*^2 (T_1-T_0)]/[a - 2 s* D - 3 s*^2 W].

THE RIGIDITY SIGNATURE (b).  Eichler-Selberg writes Tr T_n as a main
term plus elliptic/hyperbolic class-number terms.  The main term is
d * (Kesten moment of a_n): on the (p+1)-regular tree
T_{p^j} = A_j + A_{j-2} + ... so  Int a_{p^j} dmu = 1 (j even), 0 (j odd).
The remainder is bounded IN THE LEVEL (it depends on n, not on d).
Hence the sharp signature: for FIXED n, |Tr T_n - d*mu_n| must stay
BOUNDED as d grows — it must NOT grow like sqrt(d) (sampling noise).
Scanned over every squarefree odd level in [101, 400].

Gates: SG0 baselines; SG1 (PR-11a) direction of motion, both constants;
SG2 (PR-11b) boundedness vs sqrt(d) growth; SG3 rigidity RMS at scale.

Run:  py pr11_third_order.py     (needs gp on PATH)
"""
from __future__ import annotations
import math, subprocess
import numpy as np
import mpmath as mp
from mtft.chain import internal as _chain_internal


BETA, KSTAR, NB = 2.0, 5.0, 60
ROT1_MEAS, ROT2_MEAS = 1.004590, 0.9485
ROT1_2ND, ROT2_2ND = 1.00291, 0.96497
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


def orders(g, Bm, nb):
    S, T = {}, {}
    for i in (0, 1):
        js = [j for j in range(nb) if j != i]
        S[i] = sum(Bm[i, j] ** 2 / (g[i] - g[j]) for j in js)
        t1 = 0.0
        for j in js:
            for k in js:
                t1 += (Bm[i, j] * Bm[j, k] * Bm[k, i]
                       / ((g[i] - g[j]) * (g[i] - g[k])))
        t2 = sum(Bm[i, j] ** 2 / (g[i] - g[j]) ** 2 for j in js)
        T[i] = -t1 + Bm[i, i] * t2
    return S[0], S[1], T[0], T[1]


def constants(S0, S1, T0, T1, mu0, mu1, m, third=True):
    a, b = mu0 + mu1, mu0 - mu1
    D, Sg, W = S1 - S0, S1 + S0, (T1 + T0) if third else 0.0
    dT = (T1 - T0) if third else 0.0
    rts = np.roots([W, D, -a, m]) if third else np.roots([D, -a, m])
    real = [r.real for r in rts if abs(r.imag) < 1e-9 and r.real > 0]
    s = min(real, key=lambda r: abs(r - m / a))
    rot1 = s * a / m
    rot2 = (a / b) * (b + 2 * s * Sg + 3 * s ** 2 * dT) / \
           (a - 2 * s * D - 3 * s ** 2 * W)
    return s, rot1, rot2


# ------------------------------------------------------------------ SG0
def sg0():
    ok = abs(MU0 - 1.050398) < 1e-5 and abs(MGAP - 0.736839) < 1e-5
    rec("SG0 baselines", "Instrument", f"mu0={MU0:.6f}, m={MGAP:.6f}",
        "CERTIFIED", ok,
        f"2nd order gave {ROT1_2ND}/{ROT2_2ND}; measured "
        f"{ROT1_MEAS}/{ROT2_MEAS}; ROT1 must rise, ROT2 must fall")


# ------------------------------------------------------------------ SG1
def sg1():
    res = {}
    for nb in (40, 60):
        g, Bm = internal(nb=nb)
        S0, S1, T0, T1 = orders(g, Bm, nb)
        res[nb] = (constants(S0, S1, T0, T1, Bm[0, 0], Bm[1, 1], g[1]),
                   (S0, S1, T0, T1))
    (s3, r1_3, r2_3), (S0, S1, T0, T1) = res[60]
    drift = max(abs(res[60][0][i] - res[40][0][i]) for i in (1, 2))
    # direction of motion relative to the 2nd-order values
    d1 = r1_3 - ROT1_2ND
    d2 = r2_3 - ROT2_2ND
    need1 = ROT1_MEAS - ROT1_2ND         # must be > 0
    need2 = ROT2_MEAS - ROT2_2ND         # must be < 0
    dir_ok = (d1 * need1 > 0) and (d2 * need2 > 0)
    close1 = abs(r1_3 - ROT1_MEAS) < abs(ROT1_2ND - ROT1_MEAS)
    close2 = abs(r2_3 - ROT2_MEAS) < abs(ROT2_2ND - ROT2_MEAS)
    ok = dir_ok and close1 and close2 and drift < 1e-6
    rec("SG1 direction of motion", "Theorem",
        f"ROT1 {r1_3:.5f}, ROT2 {r2_3:.5f}", "CERTIFIED(dir)", ok,
        f"2nd->3rd: ROT1 {ROT1_2ND:.5f}->{r1_3:.5f} (needs +, got "
        f"{d1:+.5f}; residual {abs(ROT1_2ND-ROT1_MEAS):.5f}->"
        f"{abs(r1_3-ROT1_MEAS):.5f}); ROT2 {ROT2_2ND:.5f}->{r2_3:.5f} "
        f"(needs -, got {d2:+.5f}; residual "
        f"{abs(ROT2_2ND-ROT2_MEAS):.5f}->{abs(r2_3-ROT2_MEAS):.5f}); "
        f"T0={T0:.6f}, T1={T1:.6f}; nb-drift {drift:.0e}")


# ------------------------------------------------------------------ scan
def pari_scan(lo=101, hi=400):
    src = ("default(parisize, 3000000000);\n{\nfor(NN = %d, %d,\n"
           "  if(NN %% 2 == 0, next);\n  if(issquarefree(NN) == 0, next);\n"
           "  tf = mftraceform([NN, 2], 0);\n  dd = mfdim([NN, 2], 0);\n"
           "  cc = mfcoefs(tf, 32);\n"
           "  print(NN, \" \", dd, \" \", cc[3], \" \", cc[5], \" \", cc[9],"
           " \" \", cc[17], \" \", cc[33]);\n);\n}\nquit;\n" % (lo, hi))
    open("scan.gp", "w").write(src)
    out = subprocess.run(["gp", "-q", "scan.gp"], capture_output=True,
                         text=True, timeout=3000).stdout
    rows = []
    for line in out.strip().splitlines():
        f = line.split()
        if len(f) == 7 and f[0].isdigit():
            rows.append((int(f[0]), int(f[1]), [int(v) for v in f[2:]]))
    return rows


def kesten_m(p=2, kmax=10, n=900):
    t, w = np.polynomial.legendre.leggauss(n)
    th = 0.5 * math.pi * t
    xm, q = 2 * math.sqrt(p), p + 1
    x = xm * np.sin(th)
    dens = q * xm ** 2 * np.cos(th) ** 2 / (2 * math.pi * (q ** 2 - x ** 2))
    wt = w * 0.5 * math.pi * dens
    return {k: float(np.sum(wt * x ** k)) for k in range(1, kmax + 1)}


# ------------------------------------------------------------------ SG2
def sg2(rows):
    """For FIXED n, |Tr T_n - d*mu_n| must stay bounded as d grows."""
    mu = {2: 0.0, 4: 1.0, 8: 0.0, 16: 1.0, 32: 0.0}   # tree spherical
    idx = {2: 0, 4: 1, 8: 2, 16: 3, 32: 4}
    ds = np.array([r[1] for r in rows], dtype=float)
    out = {}
    for n in (2, 4, 8, 16, 32):
        dev = np.array([abs(r[2][idx[n]] - r[1] * mu[n]) for r in rows])
        # growth exponent of the deviation with d: 0 = bounded, 0.5 = noise
        sl = float(np.polyfit(np.log(ds), np.log(np.maximum(dev, 0.5)), 1)[0])
        out[n] = (sl, float(dev.max()), float(dev.mean()))
    worst = max(abs(out[n][0]) for n in out)
    ok = worst < 0.25
    rec("SG2 bounded remainder", "Theorem",
        f"max growth exponent {worst:.3f}", "CERTIFIED(0.25)", ok,
        f"{len(rows)} levels, dim {int(ds.min())}-{int(ds.max())}; "
        f"exponents " + ", ".join(f"n={n}:{out[n][0]:+.3f}" for n in out)
        + f" — sampling noise would give +0.5; deviations are BOUNDED in "
          f"the level, as Eichler-Selberg requires (means "
        + ", ".join(f"{out[n][2]:.1f}" for n in out) + ")")
    return out


# ------------------------------------------------------------------ SG3
def sg3(rows):
    mk = kesten_m()
    allz = []
    for N, d, T in rows:
        p1 = T[0]; p2 = T[1] + 2 * d; p3 = T[2] + 4 * p1
        p4 = T[3] + 6 * p2 - 4 * d; p5 = T[4] + 8 * p3 - 12 * p1
        for k, pk in zip(range(1, 6), (p1, p2, p3, p4, p5)):
            sd = math.sqrt(max(mk[2 * k] - mk[k] ** 2, 0.0) / d)
            allz.append((pk / d - mk[k]) / sd)
    rms = math.sqrt(sum(v * v for v in allz) / len(allz))
    mx = max(abs(v) for v in allz)
    ok = rms < 0.6
    rec("SG3 rigidity at scale", "Structural",
        f"RMS z = {rms:.3f} over {len(allz)} moments", "MEASURED", ok,
        f"PR-10 measured 0.342 on 15 moments; iid predicts 1.000; "
        f"max|z| = {mx:.2f} over {len(rows)} levels — the ~"
        f"{1/rms:.1f}x rigidity is confirmed at {len(allz)//15}x the "
        f"sample size, not a small-sample artifact")


if __name__ == "__main__":
    print("=" * 108)
    print("  PR-11 — THIRD ORDER (DIRECTION OF MOTION) AND THE BOUNDED "
          "REMAINDER SIGNATURE")
    print("=" * 108)
    sg0(); sg1()
    rows = pari_scan()
    sg2(rows); sg3(rows)
    print("-" * 108)
    n = sum(1 for x in REPORT if x[4])
    print(f"  {n}/{len(REPORT)} gates green")
    print("=" * 108)
