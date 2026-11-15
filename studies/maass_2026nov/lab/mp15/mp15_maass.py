#!/usr/bin/env python3
"""MP15, H24's test: the lowest Maass cusp forms of Gamma_0(143), sorted by parity and Atkin-Lehner signs (Hejhal's method).

Predictions S1, S2, N0 to N4 and the forecasts F1 to F5 are in MP15_H24_Plan.md, saved at 13:31:34 before this file existed.

The group is G = Gamma_0*(N) = Gamma_0(N) with all its Atkin-Lehner matrices: the integer matrices [[Q a, b], [N c, Q d]] of
determinant Q for Q | N.  A Maass form with Atkin-Lehner signs (eps_p) satisfies f(M z) = eps_Q f(z) for M of determinant Q,
eps_Q the product of eps_p over p | Q.  A point is carried into the Ford domain of G by applying, again and again, the element
that raises it most; the character is multiplied at each step.  The form is f = sum_{n >= 1} c_n sqrt(y) K_iR(2 pi n y) trig(2 pi n x),
trig = cos (even) or sin (odd), c_1 = 1.  On a horocycle y = Y below the domain, the coefficients satisfy
    c_n sqrt(Y) K_iR(2 pi n Y) = (2/Qs) sum_m chi_m f(z_m*) trig(2 pi n x_m),
which is a linear system for c_2 .. c_M.  An eigenvalue is a value of R at which the solutions from two horocycles agree.

Stages (each can be run alone):  python3 mp15_maass.py checks     V0 the K-Bessel routine, V1 the reduction, V2 SL(2, Z)
                                 python3 mp15_maass.py scan R0 R1 dR out.npz
                                 python3 mp15_maass.py roots scan1.npz [scan2.npz ...] [sectors=0,1,2,3] [Rmax=1.6] [out=roots.json]
The run for MP15: checks; scans 0.01..1.255 and 1.26..2.5 in steps of 0.005, and 0.0025..0.4725 and 0.475..0.95 in steps of 0.0025;
roots with Rmax=1.6, even and odd sectors apart; mp15_extra.py dips for the three the sign-change filter missed.  mp15_numbers.json
gathers the results.
"""
import itertools
import json
import math
import sys
import time

import mpmath as mp
import numpy as np
from scipy.optimize import brentq

T0 = time.time()

# ------------------------------------------------------------------ K_iR(x), real R, by the trapezoid rule
S_GRID = np.arange(0.0, 11.0 + 1e-9, 0.2)
W_S = np.full(S_GRID.size, 0.2)
W_S[0] = 0.1


def Kir(R, x):
    """K_{iR}(x) = int_0^inf exp(-x cosh t) cos(R t) dt with t = s/sqrt(1 + x), the trapezoid rule in s (spectrally accurate:
    the integrand is even and analytic).  Entries with x > 60 are set to 0 (below e^-60 of the rest)."""
    x = np.asarray(x, dtype=float)
    out = np.zeros(x.shape)
    keep = x <= 60.0
    xs = x[keep]
    sc = 1.0 / np.sqrt(1.0 + xs)
    t = sc[:, None] * S_GRID[None, :]
    out[keep] = (np.exp(-xs[:, None] * np.cosh(t)) * np.cos(R * t)) @ W_S * sc
    return out


def Kir_mp(R, x):
    return float(mp.besselk(1j * R, x).real)


# ------------------------------------------------------------------ the group Gamma_0*(N) and the reduction
def egcd(a, b):
    if b == 0:
        return (a, 1, 0) if a >= 0 else (-a, -1, 0)
    g, u, v = egcd(b, a % b)
    return g, v, u - (a // b) * v


def primes_of(N):
    return [p for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37) if N % p == 0]


def reduce_point(x, y, N, use_al=True):
    """Carry z = x + iy into the Ford domain of Gamma_0*(N) (of Gamma_0(N) if use_al is False).
    Returns x*, y*, and for each prime p | N the parity of the number of times an element with p | Q was applied:
    f(z) = prod eps_p^(parity_p) f(z*)."""
    ps = primes_of(N)
    Qs = [1]
    if use_al:
        Qs = sorted({math.prod(s) for k in range(len(ps) + 1) for s in itertools.combinations(ps, k)})
    par = [0] * len(ps)
    for _ in range(10000):
        x -= math.floor(x + 0.5)
        best, arg = 1.0 - 1e-13, None
        for Q in Qs:
            Nq = N // Q
            cmax = int(math.sqrt(Q) / (N * y)) + 1
            for c in range(1, cmax + 1):
                d0 = -N * c * x / Q
                fl = math.floor(d0)
                for d in range(fl - 1, fl + 3):
                    if math.gcd(Q * d, Nq * c) != 1:
                        continue
                    v = ((N * c * x + Q * d) ** 2 + (N * c * y) ** 2) / Q
                    if v < best:
                        best, arg = v, (Q, c, d)
        if arg is None:
            return x, y, tuple(par)
        Q, c, d = arg
        g, u, w = egcd(Q * d, (N // Q) * c)            # u Q d + w (N/Q) c = 1
        a, b = u, -w                                     # Q a d - (N/Q) b c = 1, so det [[Qa, b], [Nc, Qd]] = Q
        den = (N * c * x + Q * d) ** 2 + (N * c * y) ** 2
        xn = ((Q * a * x + b) * (N * c * x + Q * d) + Q * a * N * c * y * y) / den
        y = Q * y / den
        x = xn
        for i, p in enumerate(ps):
            if Q % p == 0:
                par[i] ^= 1
    raise RuntimeError("no convergence")


# ------------------------------------------------------------------ Eisenstein series at s = 2 by sorted rows (MP14 Route B)
def S2(u, Y):
    u = np.asarray(u, dtype=float)
    Y = np.asarray(Y, dtype=float)
    f = u - np.round(u)
    t = 2 * np.pi * Y
    c = np.cos(2 * np.pi * f)
    out = np.empty(t.shape)
    lo = t <= 20
    tl, cl = t[lo], c[lo]
    A, B = np.zeros(tl.shape), np.zeros(tl.shape)
    sm = tl < 0.5
    ts = tl[sm]
    tA, tB, As, Bs = ts.copy(), ts.copy(), np.zeros(ts.shape), np.zeros(ts.shape)
    for k in range(1, 15):
        tA = tA * 4 * ts * ts / ((2 * k) * (2 * k + 1))
        tB = tB * ts * ts / ((2 * k) * (2 * k + 1))
        As += tA
        Bs += 2 * k * tB
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
    out[hi] = np.pi / (2 * Y[hi] ** 3) * ((1 - q * q) / dq - th * 2 * q * (2 * q - (1 + q * q) * ch) / dq ** 2)
    return out


class Eis143:
    """E_inf, E_0, E_1/11, E_1/13 of Gamma_0(143) at s = 2, as in mp14_lattice.py."""
    def __init__(self, XT=3000):
        N = 143
        self.M = np.arange(1, N * XT + 1, dtype=np.int64)
        d11, d13 = self.M % 11 == 0, self.M % 13 == 0
        z4 = math.pi ** 4 / 90
        K = 1 / (z4 * (1 - 11.0 ** -4) * (1 - 13.0 ** -4))
        f11 = np.where(d11, 1 - 11.0 ** -4, -(11.0 ** -4))
        f13 = np.where(d13, 1 - 13.0 ** -4, -(13.0 ** -4))
        self.W = {"inf": K * f11 * f13, "0": K * (~d11) * (~d13), "1/11": K * f11 * (~d13), "1/13": K * f13 * (~d11)}
        self.P = {"inf": 1.0, "0": 143.0 ** -2, "1/11": 13.0 ** -2, "1/13": 11.0 ** -2}
        mp.mp.dps = 30
        self.T = {a: float(sum(mp.mpf(float(self.W[a][j - 1])) * mp.zeta(3, XT + mp.mpf(j) / N) for j in range(1, N + 1)) / mp.mpf(N) ** 3)
                  for a in self.W}

    def __call__(self, x, y):
        S = S2(self.M * x, self.M * y)
        out = {}
        for a in self.W:
            out[a] = (y * y if a == "inf" else 0.0) + self.P[a] * y * y * (math.fsum(self.W[a] * S) + math.pi / (2 * y ** 3) * self.T[a])
        return out


# ------------------------------------------------------------------ Hejhal's system
class Horocycle:
    def __init__(self, N, Y, Qs, use_al=True):
        self.N, self.Y, self.Qs = N, Y, Qs
        m = np.arange(1, Qs + 1)
        self.xm = (m - 0.5) / (2 * Qs)
        red = [reduce_point(float(x), Y, N, use_al) for x in self.xm]
        self.xs = np.array([r[0] for r in red])
        self.ys = np.array([r[1] for r in red])
        npr = len(primes_of(N))
        self.par = np.array([list(r[2]) for r in red], dtype=int).reshape(Qs, npr) if npr else np.zeros((Qs, 0), dtype=int)
        assert self.ys.min() > Y * 1.0001, (self.ys.min(), Y)

    def chi(self, eps):
        if len(eps) == 0:
            return np.ones(self.Qs)
        e = np.array(eps, dtype=float)
        return np.prod(np.where(self.par == 1, e[None, :], 1.0), axis=1)


def solve_coeffs(R, hc, M, parity, chis, Kcache=None):
    """For one horocycle: the coefficient vectors c (c_1 = 1, length M) for each character in chis, and the residual of the
    dropped first equation."""
    l = np.arange(1, M + 1)
    if Kcache is None:
        Kst = Kir(R, 2 * np.pi * np.outer(hc.ys, l))                 # Qs x M at the reduced points
        Kdiag = Kir(R, 2 * np.pi * l * hc.Y)
    else:
        Kst, Kdiag = Kcache
    trig = np.cos if parity == 0 else np.sin
    P = np.sqrt(hc.ys)[:, None] * Kst * trig(2 * np.pi * np.outer(hc.xs, l))
    Tn = trig(2 * np.pi * np.outer(l, hc.xm)) * (2.0 / hc.Qs)        # M x Qs
    out = []
    for chi in chis:
        V = (Tn * chi[None, :]) @ P
        V[np.arange(M), np.arange(M)] -= np.sqrt(hc.Y) * Kdiag
        A = V[1:, 1:]
        sc = np.abs(A).max(axis=0)
        sc[sc == 0] = 1.0
        cr = np.linalg.solve(A / sc[None, :], -V[1:, 0]) / sc
        c = np.concatenate([[1.0], cr])
        res = V[0, :] @ c / max(abs(V[0, 0]), 1e-300)
        out.append((c, res))
    return out


SECTORS = [(par, e11, e13) for par in (0, 1) for e11 in (1, -1) for e13 in (1, -1)]


def sector_name(s):
    par, e11, e13 = s
    return f"{'even' if par == 0 else 'odd'} e11={e11:+d} e13={e13:+d}"


def setup143(M=None):
    """The two horocycles for N = 143 and the truncation, chosen from the lowest point of the domain met by the samples."""
    N = 143
    probe = [reduce_point(float(x), 0.002, N)[1] for x in np.linspace(0.0001, 0.4999, 4000)]
    Y0 = min(probe)
    Y1, Y2 = 0.96 * Y0, 0.88 * Y0
    if M is None:
        M = int(math.ceil(40.0 / (2 * math.pi * Y0)))
    Qs = int(math.ceil(1.08 * M / 0.88)) + 12
    return N, Y0, (Horocycle(N, Y1, Qs), Horocycle(N, Y2, Qs)), M


def scan(R0, R1, dR, out):
    N, Y0, hcs, M = setup143()
    chis = [[hc.chi((e11, e13)) for (_, e11, e13) in SECTORS[:4]] for hc in hcs]
    Rs = np.arange(R0, R1 + 1e-12, dR)
    nchk = 8
    D = np.zeros((len(Rs), len(SECTORS), nchk))
    RES = np.zeros((len(Rs), len(SECTORS), 2))
    l = np.arange(1, M + 1)
    t = time.time()
    for i, R in enumerate(Rs):
        sols = []
        for h, hc in enumerate(hcs):
            Kc = (Kir(R, 2 * np.pi * np.outer(hc.ys, l)), Kir(R, 2 * np.pi * l * hc.Y))
            row = []
            for parity in (0, 1):
                row += solve_coeffs(R, hc, M, parity, chis[h], Kc)
            sols.append(row)
        for s in range(len(SECTORS)):
            c1, r1 = sols[0][s]
            c2, r2 = sols[1][s]
            D[i, s, :] = (c1 - c2)[1:nchk + 1]
            RES[i, s, :] = (r1, r2)
        if i % 50 == 0:
            print(f"  R = {R:.4f}  ({i + 1}/{len(Rs)}, {time.time() - t:.0f} s)", flush=True)
    np.savez(out, Rs=Rs, D=D, RES=RES, Y0=Y0, Y=[hc.Y for hc in hcs], M=M, Qs=hcs[0].Qs)
    print(f"scan {R0}..{R1} step {dR}: Y0 = {Y0:.6f}, M = {M}, Qs = {hcs[0].Qs}, {time.time() - T0:.0f} s")


# ------------------------------------------------------------------ checks V0, V1, V2
def checks():
    res = {}
    # V0: the K-Bessel routine against mpmath
    worst = 0.0
    for R in (0.05, 0.37, 1.0, 2.2, 3.0):
        for x in (0.02, 0.1, 0.7, 2.5, 9.0, 30.0):
            a, b = float(Kir(R, np.array([x]))[0]), Kir_mp(R, x)
            worst = max(worst, abs(a - b) / max(abs(b), 1e-300))
    print(f"{'ok  ' if worst < 1e-10 else 'FAIL'} V0  K_iR(x) by the trapezoid rule against mpmath, R <= 3, 0.02 <= x <= 30   [{worst:.1e}]")
    res["V0"] = worst
    # V1: the reduction and the characters, on the four Eisenstein combinations at s = 2
    E = Eis143()
    N = 143
    rng = np.random.default_rng(1)
    worst, low = 0.0, 1.0
    for x in rng.uniform(-0.5, 0.5, 24):
        y = float(rng.uniform(0.004, 0.03))
        xs, ys, par = reduce_point(float(x), y, N)
        low = min(low, ys)
        Ez, Es = E(float(x), y), E(xs, ys)
        for e11 in (1, -1):
            for e13 in (1, -1):
                F = lambda v: v["inf"] + e11 * e13 * v["0"] + e13 * v["1/11"] + e11 * v["1/13"]
                chi = (e11 ** par[0]) * (e13 ** par[1])
                worst = max(worst, abs(F(Ez) - chi * F(Es)) / max(abs(F(Ez)), abs(F(Es))))
    print(f"{'ok  ' if worst < 1e-8 else 'FAIL'} V1  the reduction into the domain of Gamma_0*(143) and its characters: F_eps(z) = chi F_eps(z*) for 24 points and 4 characters   [{worst:.1e}; lowest reduced height {low:.4f}]")
    res["V1"] = worst
    # V2: SL(2, Z), with K from mpmath
    def g1(R, parity):
        cs = []
        for Y in (0.72, 0.62):
            hc = Horocycle(1, Y, 40, use_al=False)
            l = np.arange(1, 26)
            Kst = np.vectorize(lambda v: Kir_mp(R, v))(2 * np.pi * np.outer(hc.ys, l))
            Kd = np.array([Kir_mp(R, 2 * np.pi * k * Y) for k in l])
            cs.append(solve_coeffs(R, hc, 25, parity, [hc.chi(())], (Kst, Kd))[0][0])
        return cs[0][1] - cs[1][1]
    mp.mp.dps = 20
    r_odd = brentq(lambda R: g1(R, 1), 9.50, 9.56, xtol=1e-11)
    r_even = brentq(lambda R: g1(R, 0), 13.75, 13.81, xtol=1e-11)
    e1, e2 = abs(r_odd - 9.5336952613535575), abs(r_even - 13.779751351890738)
    print(f"{'ok  ' if max(e1, e2) < 1e-6 else 'FAIL'} V2  SL(2, Z): the solver finds R = {r_odd:.10f} (odd) and {r_even:.10f} (even)   [errors {e1:.1e}, {e2:.1e}]")
    res["V2"] = {"odd": r_odd, "even": r_even, "err": [e1, e2]}
    N, Y0, hcs, M = setup143()
    print(f"     setup at level 143: lowest point of the domain met {Y0:.6f}; horocycles {hcs[0].Y:.6f}, {hcs[1].Y:.6f}; M = {M}, Qs = {hcs[0].Qs}")
    res["setup"] = {"Y0": Y0, "Y": [hcs[0].Y, hcs[1].Y], "M": M, "Qs": hcs[0].Qs}
    json.dump(res, open("mp15_checks.json", "w"), indent=1)
    print(f"checks done in {time.time() - T0:.0f} s")


def roots(files, out="mp15_numbers.json", sectors=None, Rmax=99.0):
    """Brackets from the scan where c_2 of the two horocycles changes sign and the difference is small at the bracket
    (a zero, not a pole), refined by Brent's method; each root checked by the agreement of c_2 .. c_20 between the
    horocycles, by Hecke relations that were not imposed, and by |c_11| sqrt(11), |c_13| sqrt(13) (newform or oldform)."""
    Z = [np.load(f) for f in files]
    Rs = np.concatenate([z["Rs"] for z in Z])
    D = np.concatenate([z["D"] for z in Z])
    order = np.argsort(Rs)
    Rs, D = Rs[order], D[order]
    N, Y0, hcs, M = setup143()
    l = np.arange(1, M + 1)
    cache = {}

    def coeffs(R, s):
        par, e11, e13 = SECTORS[s]
        key = (round(R, 13), s)
        if key not in cache:
            cs = []
            for hc in hcs:
                Kc = (Kir(R, 2 * np.pi * np.outer(hc.ys, l)), Kir(R, 2 * np.pi * l * hc.Y))
                cs.append(solve_coeffs(R, hc, M, par, [hc.chi((e11, e13))], Kc)[0][0])
            cache[key] = cs
        return cache[key]

    g = lambda R, s: coeffs(R, s)[0][1] - coeffs(R, s)[1][1]
    found = {sector_name(SECTORS[s]): [] for s in range(len(SECTORS))}
    rejected = 0
    for s in (range(len(SECTORS)) if sectors is None else sectors):
        d = D[:, s, 0]
        for i in range(2, len(Rs) - 3):
            if np.sign(d[i]) == np.sign(d[i + 1]) or Rs[i] > Rmax:
                continue
            if max(abs(d[i]), abs(d[i + 1])) >= min(abs(d[i - 2]), abs(d[i + 3])):
                rejected += 1                                     # |difference| grows into the bracket: a pole
                continue
            a, b = Rs[i], Rs[i + 1]
            try:
                R = brentq(lambda r: g(r, s), a, b, xtol=1e-12, rtol=1e-14, maxiter=60)
            except ValueError:
                print(f"  {sector_name(SECTORS[s])}: bracket [{a:.3f}, {b:.3f}] has no sign change on recomputation", flush=True)
                rejected += 1
                continue
            c1, c2 = coeffs(R, s)
            agree = float(np.max(np.abs(c1[1:20] - c2[1:20])) / max(1.0, np.max(np.abs(c1[1:20]))))
            c = c1
            hecke = max(abs(c[1] * c[2] - c[5]), abs(c[1] * c[4] - c[9]), abs(c[3] - (c[1] ** 2 - 1)),
                        abs(c[2] * c[4] - c[14]), abs(c[8] - (c[2] ** 2 - 1)))
            if agree > 1e-6 or hecke > 1e-4:
                print(f"  {sector_name(SECTORS[s])}: bracket [{a:.3f}, {b:.3f}] -> R = {R:.6f} fails the checks (agree {agree:.1e}, hecke {hecke:.1e}): a pole", flush=True)
                rejected += 1
                continue
            new = abs(abs(c[10]) * math.sqrt(11) - 1) < 1e-4 and abs(abs(c[12]) * math.sqrt(13) - 1) < 1e-4
            found[sector_name(SECTORS[s])].append({"R": R, "lambda": 0.25 + R * R, "new": bool(new), "agree": agree, "hecke": float(hecke),
                                          "c2": float(c[1]), "c3": float(c[2]), "c11_sqrt11": float(c[10] * math.sqrt(11)),
                                          "c13_sqrt13": float(c[12] * math.sqrt(13))})
            print(f"  {sector_name(SECTORS[s])}: R = {R:.10f}  lambda = {0.25 + R * R:.8f}  {'new' if new else 'OLD'}  agree {agree:.1e}  hecke {hecke:.1e}  "
                  f"c11 sqrt11 = {c[10] * math.sqrt(11):+.6f}  c13 sqrt13 = {c[12] * math.sqrt(13):+.6f}", flush=True)
    json.dump({"setup": {"Y0": Y0, "Y": [hc.Y for hc in hcs], "M": M, "Qs": hcs[0].Qs, "R_range": [float(Rs[0]), float(Rs[-1])],
                         "step": float(Rs[1] - Rs[0])}, "found": found, "rejected_brackets": rejected},
              open(out, "w"), indent=1)
    print(f"roots: {sum(len(v) for v in found.values())} eigenvalues, {rejected} brackets rejected, {time.time() - T0:.0f} s")


if __name__ == "__main__":
    if sys.argv[1] == "checks":
        checks()
    elif sys.argv[1] == "scan":
        scan(float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]), sys.argv[5])
    elif sys.argv[1] == "roots":
        args = sys.argv[2:]
        opt = {a.split("=")[0]: a.split("=")[1] for a in args if "=" in a}
        files = [a for a in args if "=" not in a]
        roots(files, out=opt.get("out", "mp15_numbers.json"), sectors=[int(v) for v in opt["sectors"].split(",")] if "sectors" in opt else None,
              Rmax=float(opt.get("Rmax", 99)))
