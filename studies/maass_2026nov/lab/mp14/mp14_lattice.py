#!/usr/bin/env python3
"""MP14 Route B (test T-G10, check B6): the Eisenstein series of Gamma_0(143) at s = 2, summed over the group.

No Fourier expansion, no Bessel function and no formula from level 1.  The series at the cusp infinity is the sum over
the cosets Gamma_inf \\ Gamma_0(143), one coset for each bottom row (C, D) with gcd(C, D) = 1 and 143 | C, up to sign:

    E_inf(z, 2) = y^2 + sum_{C >= 1, 143 | C} sum_{gcd(C, D) = 1} y^2 / |C z + D|^4.

The other cusps come from the Atkin-Lehner matrices, E_a(z) = E_inf(W_Q z) with W_Q(inf) = a:
W_143 = [[0, -1], [143, 0]] (a = 0), W_13 = [[78, 7], [143, 13]] (a = 6/11, the cusp 1/11),
W_11 = [[66, 5], [143, 11]] (a = 6/13, the cusp 1/13).  Composing the coset sum with W_Q sorts the rows by gcd(C, 143)
(derivation in the Result of MP14):

    E_0 = 143^-2 sum_{gcd(C, 143) = 1},   E_1/11 = 13^-2 sum_{gcd(C, 143) = 11},   E_1/13 = 11^-2 sum_{gcd(C, 143) = 13},

each of the same inner sum over D, and the level-1 series is y^2 plus the sum over every C.  Both forms are computed here:
the four series at z by the sorted rows, and E_inf at W_Q z, so that the cusp labels are checked from their definition.

Method for the rows: gcd(C, D) = 1 is removed by Moebius inversion over e | C, and each whole row
sum_{D in Z} |M z + D|^-4 = S(Mx, My) is done in closed form from the classical sum_m 1/((u + m)^2 + Y^2) = (pi/Y) g,
g = sinh(2 pi Y)/(cosh(2 pi Y) - cos(2 pi u)), so that S(u, Y) = pi g/(2 Y^3) - pi g_Y/(2 Y^2).  Collected by M = C/e,

    sum over the rows C of one class of sum_{gcd(C, D) = 1} |C z + D|^-4 = sum_{M >= 1} w(M) S(Mx, My),

where w(M) = sum of mu(e) e^-4 over squarefree e with eM in the class: K times one factor for each of p = 11, 13, with
K = prod over primes p other than 11, 13 of (1 - p^-4) = 1/(zeta(4)(1 - 11^-4)(1 - 13^-4)).  M runs to 143 * 7000;
beyond that S(Mx, My) = pi/(2 (My)^3) to within exp(-2 pi M y), and that tail is summed with Hurwitz zeta functions.

Checks of the machinery itself: S against brute summation over m, and E_inf at one point against the plainest coset sum
(rows C = 143k, k <= 40, |D| <= 10^5, gcd tested one pair at a time).

Usage: python3 mp14_lattice.py mp14_numbers.json   (reads Route A's values and adds the key "route_B" to the file)
"""
import json
import math
import sys
import time
from fractions import Fraction

import mpmath as mp
import numpy as np

T0 = time.time()
N = 143
XT = 7000
X = N * XT
ZETA4 = math.pi ** 4 / 90
K = 1 / (ZETA4 * (1 - 11.0 ** -4) * (1 - 13.0 ** -4))
CHECKS = []


def check(tag, text, ok, detail=""):
    CHECKS.append({"tag": tag, "text": text, "ok": bool(ok), "detail": detail})
    print(f"{'ok  ' if ok else 'FAIL'} {tag:4s} {text}" + (f"   [{detail}]" if detail else ""))


def S2(u, Y):
    """sum over m in Z of 1/((u + m)^2 + Y^2)^2 in closed form: with t = 2 pi Y and c = cos(2 pi u),
    S = pi/(2 Y^3) (g - t dg/dt), g = sinh t/(cosh t - c), and g - t dg/dt = (A + c B)/(cosh t - c)^2 with
    A = sinh t cosh t - t, B = t cosh t - sinh t.  A and B by their Taylor series for t < 1/2 (no cancellation),
    directly up to t = 20, and in q = exp(-t) beyond (no overflow); cosh t - c = 2 sinh^2(t/2) + 2 sin^2(pi u)."""
    scalar = np.ndim(u) == 0 and np.ndim(Y) == 0
    u, Y = np.broadcast_arrays(np.atleast_1d(np.asarray(u, dtype=float)), np.atleast_1d(np.asarray(Y, dtype=float)))
    f = u - np.round(u)
    t = 2 * np.pi * Y
    c = np.cos(2 * np.pi * f)
    out = np.empty(t.shape)
    lo = t <= 20
    tl, cl = t[lo], c[lo]
    A, B = np.zeros(tl.shape), np.zeros(tl.shape)
    sm = tl < 0.5
    ts = tl[sm]
    termA, termB, As, Bs = ts.copy(), ts.copy(), np.zeros(ts.shape), np.zeros(ts.shape)
    for k in range(1, 15):                           # 4^k t^(2k+1)/(2k+1)!  and  2k t^(2k+1)/(2k+1)!
        termA = termA * 4 * ts * ts / ((2 * k) * (2 * k + 1))
        termB = termB * ts * ts / ((2 * k) * (2 * k + 1))
        As += termA
        Bs += 2 * k * termB
    A[sm], B[sm] = As, Bs
    tm = tl[~sm]
    A[~sm] = np.sinh(tm) * np.cosh(tm) - tm
    B[~sm] = tm * np.cosh(tm) - np.sinh(tm)
    den = 2 * np.sinh(tl / 2) ** 2 + 2 * np.sin(np.pi * f[lo]) ** 2
    out[lo] = np.pi / (2 * Y[lo] ** 3) * (A + cl * B) / den ** 2
    hi = ~lo
    th, ch = t[hi], c[hi]
    q = np.exp(-th)
    dq = 1 + q * q - 2 * q * ch
    g = (1 - q * q) / dq
    dg = 2 * q * (2 * q - (1 + q * q) * ch) / dq ** 2
    out[hi] = np.pi / (2 * Y[hi] ** 3) * (g - th * dg)
    return float(out[0]) if scalar else out


# 1. the closed form against brute summation over m
worst = 0.0
for u, Y in ((0.3, 0.9), (-0.1, 0.05), (0.5, 5e-4), (0.001, 5e-4), (0.27, 3.0)):
    m = np.arange(-2_000_000, 2_000_001, dtype=float)
    brute = math.fsum(1 / ((u + m) ** 2 + Y * Y) ** 2)
    worst = max(worst, abs(float(S2(u, Y)) - brute) / brute)
check("S", "the closed form of sum_m ((u + m)^2 + Y^2)^-2 agrees with brute summation over |m| <= 2e6 at five (u, Y)", worst < 1e-12, f"{worst:.1e}")

# 2. the weights w(M) by class of gcd(C, 143), and their periodic tail
M = np.arange(1, X + 1, dtype=np.int64)
d11, d13 = (M % 11 == 0), (M % 13 == 0)
f11 = np.where(d11, 1 - 11.0 ** -4, -(11.0 ** -4))
f13 = np.where(d13, 1 - 13.0 ** -4, -(13.0 ** -4))
W = {"inf": K * f11 * f13, "0": K * (~d11) * (~d13), "1/11": K * f11 * (~d13), "1/13": K * f13 * (~d11),
     "level1": np.full(X, 1 / ZETA4)}
PREF = {"inf": 1.0, "0": 143.0 ** -2, "1/11": 13.0 ** -2, "1/13": 11.0 ** -2, "level1": 1.0}
part = max(abs(W["inf"] + W["0"] + W["1/11"] + W["1/13"] - 1 / ZETA4))
check("W", "the four classes of rows partition the rows of level 1: their weights add up to 1/zeta(4) for every M", part < 1e-15, f"{part:.1e}")
mp.mp.dps = 30
TAIL = {a: float(sum(mp.mpf(float(W[a][j - 1])) * mp.zeta(3, XT + mp.mpf(j) / N) for j in range(1, N + 1)) / mp.mpf(N) ** 3)
        for a in W}


def frac_parts(x):
    """M x mod 1 for M = 1 .. X: exact, from numerator and denominator, when x is a fraction; else in floating point.
    First run (finished 12:21:57): the label check L missed its tolerance 1e-10 with 2.1e-10 at z = 0.3 + 0.9i for the
    cusp 1/13 (2.7e-11 for 1/11), because M x was rounded in floating point at W_Q z, which lies within 6e-4 of the cusps
    6/13 and 6/11, where the Moebius rows cancel each other by a factor of about 600.  A fault of the apparatus, recorded in
    MP14; the fractional parts are now exact at rational points."""
    if isinstance(x, Fraction):
        p, q = x.numerator % x.denominator, x.denominator
        return (M * p % q) / q
    return M * x


def series(x, y):
    """The five series at s = 2 at z = x + iy: the four cusps of Gamma_0(143) and level 1, by the sorted rows."""
    yf = float(y)
    S = S2(frac_parts(x), M * yf)
    out = {}
    for a in W:
        body = math.fsum(W[a] * S) + math.pi / (2 * yf ** 3) * TAIL[a]
        out[a] = (yf * yf if a in ("inf", "level1") else 0.0) + PREF[a] * yf * yf * body
    return out


def mob(mat, x, y):
    """W z for z = x + iy, exactly when x and y are fractions."""
    a, b, c, d = mat
    nr, ni, dr, di = a * x + b, a * y, c * x + d, c * y
    den = dr * dr + di * di
    return (nr * dr + ni * di) / den, (ni * dr - nr * di) / den


WQ = {"0": (0, -1, 143, 0), "1/11": (78, 7, 143, 13), "1/13": (66, 5, 143, 11)}

# 3. the plainest coset sum, at one point, against the method
z = complex(0.3, 0.9)
plain = z.imag ** 2
Dall = np.arange(-100_000, 100_001, dtype=np.int64)
for k in range(1, 41):
    C = N * k
    D = Dall[np.gcd(Dall, C) == 1]
    plain += math.fsum(z.imag ** 2 / ((C * z.real + D) ** 2 + (C * z.imag) ** 2) ** 2)
tail40 = math.pi / (2 * z.imag) * sum(1.0 / (N * k) ** 3 for k in range(41, 10 ** 6))   # rows beyond k = 40, bound
viaS = series(Fraction(3, 10), Fraction(9, 10))["inf"]
dplain = abs(plain - viaS)
check("P", "at z = 0.3 + 0.9i the plainest coset sum (C = 143k, k <= 40, |D| <= 1e5) agrees with the method to within the rows it leaves out",
      dplain < tail40, f"difference {dplain:.2e}, rows left out at most {tail40:.2e}")

# 4. Route A's values (mp14_tg10_verify.py, 40 digits) against Route B
path = sys.argv[1] if len(sys.argv) > 1 else "mp14_numbers.json"
data = json.load(open(path))
routeA = data["numbers"]["B"]["s2"]
rows, b6, extra, labels = [], 0.0, 0.0, 0.0
for rec in routeA:
    zx, zy = (Fraction(str(v)) for v in rec["z"])          # the decimals Route A used, exactly
    Bv = series(zx, zy)
    A = dict(rec["E"])
    A["level1"] = rec["E1"]
    rel = {a: abs(Bv[a] - A[a]) / abs(A[a]) for a in A}
    b6 = max(b6, rel["inf"])
    extra = max(extra, max(rel[a] for a in ("0", "1/11", "1/13", "level1")))
    viaW = {a: series(*mob(WQ[a], zx, zy))["inf"] for a in WQ}
    lab = {a: abs(viaW[a] - Bv[a]) / abs(Bv[a]) for a in WQ}
    labels = max(labels, max(lab.values()))
    rows.append({"z": rec["z"], "B": Bv, "A": A, "rel": rel, "E_inf_at_WQz": viaW, "rel_label": lab})
check("B6", "Route B sums E_inf(z, 2) over the cosets of Gamma_0(143) and agrees with the formula from level 1 at both points (1e-8)", b6 < 1e-8, f"{b6:.1e}")
check("B6+", "not predicted: E_0, E_1/11, E_1/13 and the level-1 series agree as well (1e-8)", extra < 1e-8, f"{extra:.1e}")
check("L", "not predicted: E_inf(W_Q z) = E_a(z) for Q = 143, 13, 11, so the cusp labels 0, 1/11, 1/13 hold by their definition (1e-10)",
      labels < 1e-10, f"{labels:.1e}")

# 5. the Pell line by Route B (B3 again, at s = 2 only)
r0 = 1 / math.sqrt(143)
pell = []
for th in (0.4, 1.1, 2.3):
    z = r0 * complex(math.cos(th), math.sin(th))
    E = series(z.real, z.imag)
    sc = max(abs(E[a]) for a in ("inf", "0", "1/11", "1/13"))
    pell.append(max(abs(E["inf"] - E["0"]), abs(E["1/11"] - E["1/13"])) / sc)
    E2 = series(1.3 * z.real, 1.3 * z.imag)
    pell.append(-abs(E2["inf"] - E2["0"]) / max(abs(E2[a]) for a in ("inf", "0", "1/11", "1/13")))
on = max(p for p in pell if p >= 0)
off = min(-p for p in pell if p < 0)
check("B3'", "not predicted: on the Pell line E_inf = E_0 and E_1/11 = E_1/13 at s = 2 by Route B as well; off it they differ", on < 1e-11 and off > 1e-3,
      f"on the line {on:.1e}, off it at least {off:.2e}")

npass = sum(c["ok"] for c in CHECKS)
print(f"{npass} of {len(CHECKS)} checks passed in {time.time() - T0:.0f} s")
for r in rows:
    print(f"  z = {r['z']}: " + ", ".join(f"{a} {r['B'][a]:.15g} (A {r['A'][a]:.15g})" for a in ("inf", "0", "1/11", "1/13", "level1")))
data["route_B"] = {"checks": CHECKS, "M_max": X, "points": rows, "pell_line": on, "off_line": off, "S_closed_form": worst,
                   "plain_coset_sum": {"value": plain, "difference": dplain, "bound": tail40}}
json.dump(data, open(path, "w"), indent=1)
