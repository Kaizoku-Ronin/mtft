#!/usr/bin/env python3
"""MP19, the even sectors of X0(N): the Selberg trace formula with the continuum.  Plan: ../MP19_Even_Trace_Plan.md (saved at
13:39:36 New York time on 13 November 2026, before this file existed).

For squarefree N with n = 2^omega(N) Atkin-Lehner elements (and n cusps) and a character eps of them, the even and odd cusp
forms of the sector, the constant function (trivial eps only) and the continuum are the spectrum of the one-cusp group
Gamma_0*(N) with the character eps.  Its trace formula is the standard one; subtracting the odd sum of MP16 and MP17 gives

    Even_eps(h) = Full_eps(h) - Odd_eps(h) - C_eps(h) - delta_eps h(i/2),

    Full_eps = (1/n) sum_Q eps_Q T_{W_Q},
      T_1     = (vol/4 pi) int h(r) r tanh(pi r) dr + sum_t w_t g(l_t) + n [ h(0)/2 - g(0) log 2 - Psi(h) ],
      T_{W_Q} = E_Q (1/4) int h(r) sech(pi r) dr + sum_s w_s g(l_s)                       (Q != 1: no cusp term),
      Psi(h)  = (1/2 pi) int h(r) psi(1 + i r) dr,
    Odd_eps  = (1/2n) sum_Q eps_Q (T_{W_Q} - T_{R W_Q})                                    (MP16, MP17),
    C_eps    = -(1/4 pi) int h(r) (phi_eps'/phi_eps)(1/2 + i r) dr,   phi_eps = phi_1 prod_p f_p^(eps_p),
             = C_1 + sum_{p | N} log p [ g(0) + sum_{k >= 1} (-eps_p)^k p^(-k/2) g(k log p) ],
    C_1      = (1/2 pi) int h(r) psi(1/2 + i r) dr - g(0) log pi - 2 sum_n (Lambda(n)/n) g(2 log n) + h(0)/2.

Every geometric weight comes from MP16's and MP17's caches (lengths to 23.5).  The elliptic term of a class of order m is
sum_{l=1}^{m-1} (4 m sin(pi l/m))^-1 int h(r) cosh(pi (1 - 2l/m) r)/cosh(pi r) dr (MP16 measured the m = 2 constant).

Usage: python3 -I mp19_even.py checks                       V9, V10, V11                         -> mp19_checks.json
       python3 -I mp19_even.py level1                       V12: SL(2, Z)'s whole formula cancels -> mp19_level1.json
       python3 -I mp19_even.py level11                      P1: Hejhal at level 11 and the even formula -> mp19_level11.json
       python3 -I mp19_even.py first                        P2, P3 on MP15's even list             -> mp19_first.json
       python3 -I mp19_even.py scan R0 R1 dR out.npz        the pole-free detector, 8 untwisted sectors, horocycle 0.96 Y0
       python3 -I mp19_even.py find out.json part scan1.npz ...   dips confirmed on two horocycles (part: 0, 1 or all)
       python3 -I mp19_even.py fine which R0 R1 dR s1,s2 out.json   the full system on one horocycle on a fine grid (a dip test)
       python3 -I mp19_even.py final lists.json out.json [R_top]   P3, P4, P6 on the final lists (R_top 1.8 unless given)
Needs numpy, scipy, mpmath, cypari2 and the scripts and caches of MP15 to MP18 in ../mp15 ... ../mp18."""
import importlib.util
import json
import math
import os
import sys
import time

import mpmath as mp
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq, minimize_scalar
from scipy.special import digamma

HERE = os.path.dirname(os.path.abspath(__file__))
T0 = time.time()


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


t17 = _load("mp17_trace", os.path.join(HERE, "..", "mp17", "mp17_trace.py"))
t16 = t17.t16
EULER = t17.EULER
LMAX = 23.5
MP16 = os.path.join(HERE, "..", "mp16")


def dump(obj, name):
    json.dump(obj, open(os.path.join(HERE, name), "w"), indent=1)


# ------------------------------------------------------------------------------------------------- test functions
def testfun(r0, sig):
    h, g, gp = t17.testfun(r0, sig)
    hi2 = 2 * math.exp((0.25 - r0 * r0) / (2 * sig * sig)) * math.cos(r0 / (2 * sig * sig))        # h(i/2), real
    return h, g, gp, hi2


TESTS = ((0.0, 0.35), (0.0, 0.45), (0.445, 0.35), (0.9, 0.35))


def tau(h, Rtop):
    return quad(lambda r: 3.5 * r * h(r), Rtop, 60, limit=200)[0] + 1e-6


# ------------------------------------------------------------------------------------------------- cusp and continuum terms
def int_r_tanh(h):
    return 2 * quad(lambda r: h(r) * r * math.tanh(math.pi * r), 0, 40, limit=400, epsabs=1e-14, epsrel=1e-13)[0]


def psi_quad(h, a):
    """(1/2 pi) int h(r) psi(a + i r) dr over the real line, a = 1 or 1/2."""
    return quad(lambda r: h(r) * digamma(a + 1j * r).real, 0, 40, limit=400, epsabs=1e-14, epsrel=1e-13)[0] / math.pi


def psi1_from_g(g):
    """Second route for Psi(h): -gamma g(0) + int_0^inf (g(0) - g(u))/(e^u - 1) du."""
    g0 = float(g(0.0))
    f = lambda u: (g0 - float(g(u))) / math.expm1(u) if u > 1e-12 else 0.0
    return -EULER * g0 + quad(f, 0, 80, limit=400, epsabs=1e-14, epsrel=1e-13)[0]


_GL = None


def gl_nodes():
    """Gauss-Legendre nodes on [0, 14] (28 panels of 20) with theta_1'(r) = (phi_1'/phi_1)(1/2 + i r) at each, by mpmath."""
    global _GL
    if _GL is not None:
        return _GL
    cache = os.path.join(HERE, "theta1.npz")
    if os.path.exists(cache):
        z = np.load(cache)
        _GL = (z["x"], z["w"], z["th"])
        return _GL
    x0, w0 = np.polynomial.legendre.leggauss(20)
    edges = np.linspace(0.0, 14.0, 29)
    x = np.concatenate([(b - a) / 2 * x0 + (a + b) / 2 for a, b in zip(edges[:-1], edges[1:])])
    w = np.concatenate([(b - a) / 2 * w0 for a, b in zip(edges[:-1], edges[1:])])
    mp.mp.dps = 30
    th = np.array([theta1p(float(r)) for r in x])
    np.savez(cache, x=x, w=w, th=th)
    _GL = (x, w, th)
    return _GL


def theta1p(r):
    """theta_1'(r) = (phi_1'/phi_1)(1/2 + i r) = -2 Re psi(1/2 + i r) + 2 log pi - 4 Re (zeta'/zeta)(1 + 2 i r); real on the line."""
    s = mp.mpc(1, 2 * r)
    return float(-2 * mp.re(mp.digamma(mp.mpc(0.5, r))) + 2 * mp.log(mp.pi) - 4 * mp.re(mp.zeta(s, derivative=1) / mp.zeta(s)))


def C1_quad(h):
    x, w, th = gl_nodes()
    return -float(np.sum(w * h(x) * th)) / (2 * math.pi)


_PRIMES = None


def lambda_sum(g, umax=LMAX):
    """sum_n (Lambda(n)/n) g(2 log n) over prime powers n with 2 log n <= umax."""
    global _PRIMES
    nmax = int(math.exp(umax / 2)) + 1
    if _PRIMES is None:
        sieve = np.ones(nmax + 1, dtype=bool)
        sieve[:2] = False
        for i in range(2, int(nmax ** 0.5) + 1):
            if sieve[i]:
                sieve[i * i::i] = False
        _PRIMES = np.nonzero(sieve)[0]
    tot = 0.0
    for p in _PRIMES:
        lp, q = math.log(p), int(p)
        while q <= nmax:
            tot += lp / q * float(g(2 * math.log(q)))
            q *= int(p)
    return tot


def C1_closed(h, g):
    return psi_quad(h, 0.5) - float(g(0.0)) * math.log(math.pi) - 2 * lambda_sum(g) + float(h(0.0)) / 2


def level_quad(h, p, e):
    """-(1/4 pi) int h theta_p' dr with theta_p'(r) = -2 L (p + e sqrt(p) cos(r L)) / (1 + p + 2 e sqrt(p) cos(r L)), e = eps_p."""
    L, sp = math.log(p), math.sqrt(p)
    f = lambda r: h(r) * (-2 * L) * (p + e * sp * math.cos(r * L)) / (1 + p + 2 * e * sp * math.cos(r * L))
    return -quad(f, 0, 40, limit=800, epsabs=1e-14, epsrel=1e-13)[0] / (2 * math.pi)


def level_series(g, p, e):
    L = math.log(p)
    tot, k = float(g(0.0)), 1
    while k * L <= 60:
        tot += (-e) ** k * p ** (-k / 2) * float(g(k * L))
        k += 1
    return L * tot


def checks():
    out = {}
    for r0, sig in TESTS + ((1.0, 0.35),):
        h, g, gp, hi2 = testfun(r0, sig)
        a, b = C1_quad(h), C1_closed(h, g)
        row = {"V9 C1 by quadrature": a, "V9 C1 closed form": b, "V9 difference": a - b}
        worst = 0.0
        for p in (11, 13):
            for e in (1, -1):
                q_, s_ = level_quad(h, p, e), level_series(g, p, e)
                row[f"V10 p={p} eps={e:+d}"] = {"quadrature": q_, "series": s_, "difference": q_ - s_}
                worst = max(worst, abs(q_ - s_))
        pa, pb = psi_quad(h, 1.0), psi1_from_g(g)
        row.update({"V10 worst difference": worst, "V11 Psi by quadrature": pa, "V11 Psi from g": pb, "V11 difference": pa - pb, "h(i/2)": hi2})
        out[f"r0={r0} sigma={sig}"] = row
        print(f"(r0 {r0}, sigma {sig}): V9 {a:+.12f} vs {b:+.12f} ({a - b:+.1e})   V10 worst {worst:.1e}   V11 {pa:+.12f} vs {pb:+.12f} ({pa - pb:+.1e})", flush=True)
    out["pass"] = all(abs(v["V9 difference"]) < 1e-9 and v["V10 worst difference"] < 1e-10 and abs(v["V11 difference"]) < 1e-10 for v in out.values())
    dump(out, "mp19_checks.json")
    print("V9, V10, V11:", "pass" if out["pass"] else "FAIL", f"({time.time() - T0:.0f} s)")


def C_eps(h, g, eps):
    """The continuum of the sector: eps is a dict prime -> sign."""
    return C1_quad(h) + sum(level_series(g, p, e) for p, e in eps.items())


# ------------------------------------------------------------------------------------------------- level 1: V12
def level1():
    cache = os.path.join(HERE, "weights_1_1.npz")
    if os.path.exists(cache):
        hyp = np.load(cache)["hyp"]
    else:
        rows, t = [], 3
        while 2 * math.acosh(t / 2) <= LMAX:
            W, D, c0 = t16.W_of(t * t - 4, 1, 1)
            rows.append((t, 2 * math.acosh(t / 2), 2 * W))
            if t % 20000 == 0:
                print(f"  t = {t} ({time.time() - T0:.0f} s)", flush=True)
            t += 1
        hyp = np.array(rows)
        np.savez(cache, lmax=LMAX, hyp=hyp)
    c_ell, _ = t16.constants()
    rows = []
    for r0, sig in ((0.0, 0.35), (0.0, 0.45), (1.0, 0.35)):
        h, g, gp, hi2 = testfun(r0, sig)
        p = {"identity": (math.pi / 3) / (4 * math.pi) * int_r_tanh(h),
             "hyperbolic": float(np.sum(g(hyp[:, 1]) * hyp[:, 2])),
             "elliptic, order 2": 0.5 * c_ell * t16.int_h_sech(h),
             "elliptic, order 3": 2 * quad(lambda r: h(r) * math.cosh(math.pi * r / 3) / math.cosh(math.pi * r), 0, 40, limit=400,
                                           epsabs=1e-14, epsrel=1e-13)[0] / (3 * math.sqrt(3)),
             "cusp: h(0)/2": float(h(0.0)) / 2, "cusp: -g(0) log 2": -float(g(0.0)) * math.log(2), "cusp: -Psi": -psi_quad(h, 1.0),
             "continuum: -C1": -C1_quad(h), "constant: -h(i/2)": -hi2}
        tot = sum(p.values())
        rows.append({"r0": r0, "sigma": sig, "parts": p, "sum (cusp forms)": tot, "pass": abs(tot) < 1e-7})
        print(f"level 1 (r0 {r0}, sigma {sig}): sum = {tot:+.3e}   " + "  ".join(f"{k} {v:+.6f}" for k, v in p.items()), flush=True)
    dump({"hyperbolic classes (by t)": int(hyp.shape[0]), "rows": rows, "pass": all(r["pass"] for r in rows)}, "mp19_level1.json")
    print("V12:", "pass" if all(r["pass"] for r in rows) else "FAIL", f"({time.time() - T0:.0f} s)")


# ------------------------------------------------------------------------------------------------- the even sector at level N
class Level:
    def __init__(self, N):
        self.N = N
        self.ps = t16.primes_of(N)
        self.n = 2 ** len(self.ps)
        self.B = t17.OddTrace(N, LMAX)
        self.c_ell, c_refl = t16.constants()
        self.Qs = sorted(q for q in range(2, N + 1) if N % q == 0)
        self.cos = {Q: t16.Coset(N, Q, LMAX, self.c_ell, c_refl, cache=os.path.join(MP16, f"weights_{N}_{Q}.npz")) for Q in self.Qs}

    def epsQ(self, eps, Q):
        return math.prod(e for p, e in eps.items() if Q % p == 0)

    def pieces(self, r0, sig, eps, L=None):
        """Every piece of Even_eps(h); they sum to the even cusp forms of the sector.  eps: dict prime -> sign."""
        h, g, gp, hi2 = testfun(r0, sig)
        n = self.n
        bp = self.B.parts(h, g, gp, L)
        out = {"area": bp["identity"] / n, "geodesics": bp["geodesics"] / n, "glides": -bp["glides"] / n,
               "cusps: h(0)/2 - g(0) log 2 - Psi": float(h(0.0)) / 2 - float(g(0.0)) * math.log(2) - psi_quad(h, 1.0),
               "cusp mirrors (MP17's C, sign reversed)": -bp["cusps and cusp mirrors"] / n}
        for Q in self.Qs:
            c = self.cos[Q]
            if L is not None:
                saved = (c.hyp, c.gli)
                c.hyp, c.gli = [x for x in c.hyp if x[1] <= L], [x for x in c.gli if x[1] <= L]
            p = c.parts(h, g)
            if L is not None:
                c.hyp, c.gli = saved
            e = self.epsQ(eps, Q)
            out[f"W{Q} fixed points"] = e * p["elliptic"] / (2 * n)
            out[f"W{Q} mirrors"] = -e * p["mirror"] / (2 * n)
            out[f"W{Q} geodesics"] = e * p["hyperbolic"] / (2 * n)
            out[f"W{Q} glides"] = -e * p["glide"] / (2 * n)
        out["continuum: -C1"] = -C1_quad(h)
        for p_, e in eps.items():
            out[f"continuum: level factor at {p_}"] = -level_series(g, p_, e)
        out["constant: -h(i/2)"] = -hi2 if all(e == 1 for e in eps.values()) else 0.0
        return out, h

    def even(self, r0, sig, eps, L=None):
        p, h = self.pieces(r0, sig, eps, L)
        return sum(p.values()), p, h


def sector_name(par, eps):
    return ("even" if par == 0 else "odd") + " (" + ",".join("+" if e > 0 else "-" for e in eps) + ")"


# ------------------------------------------------------------------------------------------------- level 11: P1
def _m15():
    return _load("mp15_maass", os.path.join(HERE, "..", "mp15", "mp15_maass.py"))


def vmat(m, R, hc, fr, parity, chi, Kc):
    Kst, Kdiag = Kc
    M = fr.size
    trig = np.cos if parity == 0 else np.sin
    P = np.sqrt(hc.ys)[:, None] * Kst * trig(2 * np.pi * np.outer(hc.xs, fr))
    Tn = trig(2 * np.pi * np.outer(fr, hc.xm)) * (2.0 / hc.Qs)
    V = (Tn * chi[None, :]) @ P
    V[np.arange(M), np.arange(M)] -= np.sqrt(hc.Y) * Kdiag
    sc = np.abs(V).max(axis=0)
    sc[sc == 0] = 1.0
    return V / sc[None, :]


def kc(m, R, hc, fr):
    return m.Kir(R, 2 * np.pi * np.outer(hc.ys, fr)), m.Kir(R, 2 * np.pi * fr * hc.Y)


def smin_sys(m, R, hc, fr, parity, chi, Kc=None):
    sv = np.linalg.svd(vmat(m, R, hc, fr, parity, chi, Kc if Kc is not None else kc(m, R, hc, fr)), compute_uv=False)
    return sv[-1] / sv[0]


def tests_untw(c1, c2, ps):
    """MP15's rules: agreement of c_2..c_20 between the horocycles; c2 c3 = c6, c2 c5 = c10, c4 = c2^2 - 1, c3 c5 = c15,
    c9 = c3^2 - 1 (none imposed); c_p sqrt(p) for p | N."""
    agree = float(np.max(np.abs(c1[1:20] - c2[1:20])) / max(1.0, np.max(np.abs(c1[1:20]))))
    c = c1
    hecke = float(max(abs(c[1] * c[2] - c[5]), abs(c[1] * c[4] - c[9]), abs(c[3] - (c[1] ** 2 - 1)), abs(c[2] * c[4] - c[14]),
                      abs(c[8] - (c[2] ** 2 - 1))))
    return agree, hecke, {p: float(c[p - 1] * math.sqrt(p)) for p in ps}


def level11():
    m = _m15()
    N = 11
    probe = [m.reduce_point(float(x), 0.01, N)[1] for x in np.linspace(0.0001, 0.4999, 4000)]
    Y0 = min(probe)
    M = int(math.ceil(40.0 / (2 * math.pi * Y0)))
    Qs = int(math.ceil(1.08 * M / 0.88)) + 12
    hcs = {"A": (m.Horocycle(N, 0.96 * Y0, Qs), m.Horocycle(N, 0.88 * Y0, Qs)), "B": (m.Horocycle(N, 0.93 * Y0, Qs), m.Horocycle(N, 0.85 * Y0, Qs))}
    fr = np.arange(1, M + 1, dtype=float)
    print(f"level 11: lowest point met {Y0:.6f}, M = {M}, Qs = {Qs}", flush=True)
    Rs = np.arange(0.05, 3.0 + 1e-9, 0.005)
    found, unconfirmed = [], []

    def coeffs(R, pair, par, e):
        return [m.solve_coeffs(R, hc, M, par, [hc.chi((e,))], kc(m, R, hc, fr))[0][0] for hc in hcs[pair]]

    for par in (0, 1):
        for e in (1, -1):
            hc = hcs["A"][0]
            chi = hc.chi((e,))
            S = np.log(np.array([smin_sys(m, float(R), hc, fr, par, chi) for R in Rs]))
            dips = [i for i in range(1, len(Rs) - 1) if S[i] < S[i - 1] and S[i] <= S[i + 1]]
            for i in dips:
                res = minimize_scalar(lambda r: math.log(smin_sys(m, r, hc, fr, par, chi)), bounds=(Rs[i - 1], Rs[i + 1]), method="bounded",
                                      options={"xatol": 1e-10})
                if res.fun > -13:
                    continue
                r1, got = float(res.x), None
                for pair in ("A", "B"):
                    for d in (2e-4, 6e-4):
                        ca, cb = coeffs(r1 - d, pair, par, e), coeffs(r1 + d, pair, par, e)
                        ks = [k for k in (1, 2, 3, 4) if np.sign(ca[0][k] - ca[1][k]) != np.sign(cb[0][k] - cb[1][k])]
                        for k in ks:
                            F = lambda r, k=k: (lambda c: c[0][k] - c[1][k])(coeffs(r, pair, par, e))
                            try:
                                Rr = brentq(F, r1 - d, r1 + d, xtol=1e-12, rtol=1e-14)
                            except ValueError:
                                continue
                            c = coeffs(Rr, pair, par, e)
                            agree, hecke, cp = tests_untw(c[0], c[1], (11,))
                            if agree <= 1e-6 and hecke <= 1e-4:
                                got = (pair, Rr, agree, hecke, cp)
                                break
                        if got:
                            break
                    if got:
                        break
                if got is None:
                    other = hcs["A"][1]
                    r2 = minimize_scalar(lambda r: math.log(smin_sys(m, r, other, fr, par, other.chi((e,)))), bounds=(r1 - 0.003, r1 + 0.003),
                                         method="bounded", options={"xatol": 1e-8})
                    unconfirmed.append({"sector": sector_name(par, (e,)), "dip": r1, "log smin": float(res.fun),
                                        "on 0.88 Y0": {"R at minimum": float(r2.x), "log smin": float(r2.fun)}})
                    print(f"  {sector_name(par, (e,))}: dip at {r1:.6f} ({res.fun:.1f}) not confirmed; on 0.88 Y0: {r2.x:.6f} ({r2.fun:.1f})", flush=True)
                    continue
                pair, Rr, agree, hecke, cp = got
                found.append({"sector": sector_name(par, (e,)), "parity": "even" if par == 0 else "odd", "eps11": e, "R": Rr, "agree": agree, "hecke": hecke,
                              "c11_sqrt11": cp[11], "new": bool(abs(abs(cp[11]) - 1) < 1e-4), "confirmed on pair": pair})
                print(f"  {sector_name(par, (e,))}: R = {Rr:.10f}  agree {agree:.1e}  hecke {hecke:.1e}  c11 sqrt11 {cp[11]:+.6f}", flush=True)
    lev = Level(11)
    rows = []
    for r0, sig in ((0.0, 0.35), (1.0, 0.35), (0.0, 0.45), (2.0, 0.35)):
        for e in (1, -1):
            ev, p, h = lev.even(r0, sig, {11: e})
            lst = sum(float(h(x["R"])) for x in found if x["parity"] == "even" and x["eps11"] == e)
            S_odd, _, _ = odd_sum(lev, r0, sig, {11: e})
            lst_odd = sum(float(h(x["R"])) for x in found if x["parity"] == "odd" and x["eps11"] == e)
            rows.append({"r0": r0, "sigma": sig, "eps11": e, "even: geometry": ev, "even: Hejhal list": lst, "even: difference": lst - ev,
                         "odd: geometry (MP17)": S_odd, "odd: Hejhal list": lst_odd, "odd: difference": lst_odd - S_odd, "parts": p,
                         "registered": (r0, sig) in ((0.0, 0.35), (1.0, 0.35)), "pass": abs(lst - ev) <= 1e-6})
            print(f"level 11 (r0 {r0}, sigma {sig}) eps11 {e:+d}: even geometry {ev:+.3e}  list {lst:+.3e}  diff {lst - ev:+.2e}   |  odd geometry {S_odd:+.3e}  list {lst_odd:+.3e}  diff {lst_odd - S_odd:+.2e}", flush=True)
    dump({"Y0": Y0, "M": M, "Qs": Qs, "R range": [0.05, 3.0], "forms": found, "unconfirmed dips": unconfirmed, "rows": rows,
          "pass": all(r["pass"] for r in rows if r["registered"])}, "mp19_level11.json")
    print("P1:", "pass" if all(r["pass"] for r in rows if r["registered"]) else "FAIL", f"({time.time() - T0:.0f} s)")


def odd_sum(lev, r0, sig, eps, L=None):
    """MP16's and MP17's odd sum of the sector, from geometry."""
    h, g, gp, _ = testfun(r0, sig)
    bp = lev.B.parts(h, g, gp, L)
    tot = sum(bp.values()) / lev.n
    for Q in lev.Qs:
        tot += lev.epsQ(eps, Q) * sum(lev.cos[Q].parts(h, g).values()) / (2 * lev.n)
    return tot, bp, h


# ------------------------------------------------------------------------------------------------- level 143: P2, P3
SECT = ((1, 1), (1, -1), (-1, 1), (-1, -1))
MP15_EVEN = {(1, 1): [0.8671117067, 0.8992552254, 1.3097603252], (1, -1): [0.9926470124, 1.5760971285], (-1, 1): [1.4142825437],
             (-1, -1): [0.4449803610, 1.0820932551, 1.1450130025, 1.5755000386]}


def solve_R1(S_geo, others, r0, sig, lo=0.05, hi=0.9):
    h, _, _, _ = testfun(r0, sig)
    target = S_geo - sum(float(h(R)) for R in others)
    return brentq(lambda R: float(h(R)) - target, lo, hi, xtol=1e-13)


def first():
    lev = Level(143)
    rows = []
    for r0, sig in TESTS:
        for s in SECT:
            ev, p, h = lev.even(r0, sig, {11: s[0], 13: s[1]})
            lst = sum(float(h(R)) for R in MP15_EVEN[s])
            t_ = tau(h, 1.6)
            rows.append({"r0": r0, "sigma": sig, "sector": s, "geometry": ev, "MP15 list": lst, "list - geometry": lst - ev, "tau": t_,
                         "pass": abs(lst - ev) <= t_, "parts": p})
            print(f"(r0 {r0}, sigma {sig}) even {s}: geometry {ev:+.9f}  list {lst:+.9f}  diff {lst - ev:+.2e}  tau {t_:.1e}  "
                  f"{'PASS' if rows[-1]['pass'] else 'MISS'}", flush=True)
    ev, p, h = lev.even(0.0, 0.35, {11: -1, 13: -1})
    R1 = solve_R1(ev, MP15_EVEN[(-1, -1)][1:], 0.0, 0.35)
    print(f"P3 on MP15's list: R1 = {R1:.9f}, from 0.4449803610 by {R1 - 0.4449803610:+.2e}")
    dump({"rows": rows, "P3 on MP15's list": {"R1": R1, "difference": R1 - 0.4449803610, "pass (2e-3)": abs(R1 - 0.4449803610) <= 2e-3}}, "mp19_first.json")
    print(f"first: {time.time() - T0:.0f} s")


# ------------------------------------------------------------------------------------------------- the pole-free search at 143: P5
MD = 300
SUB = [(par, e) for par in (0, 1) for e in SECT]


def _tw():
    return _load("mp18_twisted", os.path.join(HERE, "..", "mp18", "mp18_twisted.py"))


def scan(R0, R1, dR, out):
    tw = _tw()
    m = tw.m
    Y0 = tw.setup(0.96, 0.88)[0]
    hc = tw.Horo(0.96 * Y0, int(math.ceil(1.08 * MD / 0.88)) + 12)
    fr = np.arange(1, MD + 1, dtype=float)
    chis = [hc.chi(e, twist=False) for e in SECT]
    Rs = np.arange(R0, R1 + 1e-12, dR)
    S = np.zeros((len(Rs), 8))
    for i, R in enumerate(Rs):
        K = kc(m, float(R), hc, fr)
        for s, (par, e) in enumerate(SUB):
            S[i, s] = smin_sys(m, float(R), hc, fr, par, chis[SECT.index(e)], K)
        if i % 50 == 0:
            print(f"  R = {R:.4f} ({i + 1}/{len(Rs)}, {time.time() - T0:.0f} s)", flush=True)
    np.savez(out, Rs=Rs, S=S, Y0=Y0, MD=MD)
    print(f"scan {R0}..{R1}: {time.time() - T0:.0f} s")


def find(out, part, files):
    tw = _tw()
    m = tw.m
    Z = [np.load(f) for f in files]
    Rs = np.concatenate([z["Rs"] for z in Z])
    S = np.concatenate([z["S"] for z in Z])
    Rs, first_ = np.unique(np.round(Rs, 9), return_index=True)
    S = S[first_]
    L = np.log(S)
    Y0, M, hcsA = tw.setup(0.96, 0.88)
    Qs = hcsA[0].Qs
    pairs = {"A": hcsA, "B": (tw.Horo(0.93 * Y0, Qs), tw.Horo(0.85 * Y0, Qs)), "C": (tw.Horo(0.90 * Y0, Qs), tw.Horo(0.82 * Y0, Qs))}
    hcd = tw.Horo(0.96 * Y0, int(math.ceil(1.08 * MD / 0.88)) + 12)
    frd = np.arange(1, MD + 1, dtype=float)
    fr = np.arange(1, M + 1, dtype=float)
    found, unconfirmed = [], []
    for s, (par, e) in enumerate(SUB):
        if part in ("0", "1") and s % 2 != int(part):
            continue
        chid = hcd.chi(e, twist=False)
        dips = [i for i in range(1, len(Rs) - 1) if L[i, s] < L[i - 1, s] and L[i, s] <= L[i + 1, s]]
        for i in dips:
            res = minimize_scalar(lambda r: math.log(smin_sys(m, r, hcd, frd, par, chid)), bounds=(Rs[i - 1], Rs[i + 1]), method="bounded",
                                  options={"xatol": 1e-10})
            if res.fun > -13:
                continue
            r1, got = float(res.x), None
            for name, hcs in pairs.items():
                cf = lambda r: tw.coeffs(r, hcs, fr, par, e, twist=False)
                for d in (2e-4, 6e-4):
                    ca, cb = cf(r1 - d), cf(r1 + d)
                    ks = [k for k in (1, 2, 3, 4) if np.sign(ca[0][k] - ca[1][k]) != np.sign(cb[0][k] - cb[1][k])]
                    for k in ks:
                        try:
                            Rr = brentq(lambda r, k=k: (lambda c: c[0][k] - c[1][k])(cf(r)), r1 - d, r1 + d, xtol=1e-12, rtol=1e-14)
                        except ValueError:
                            continue
                        c = cf(Rr)
                        agree, hecke, cp = tests_untw(c[0], c[1], (11, 13))
                        if agree <= 1e-6 and hecke <= 1e-4:
                            got = (name, Rr, agree, hecke, cp, c[0])
                            break
                    if got:
                        break
                if got:
                    break
            if got is None:
                row = {"subsector": sector_name(par, e), "dip": r1, "log smin (detector)": float(res.fun)}
                for nm, hc in (("0.96 Y0", hcsA[0]), ("0.88 Y0", hcsA[1])):                # the full system on both horocycles of pair A
                    chi = hc.chi(e, twist=False)
                    r2 = minimize_scalar(lambda r: math.log(smin_sys(m, r, hc, fr, par, chi)), bounds=(r1 - 0.003, r1 + 0.003), method="bounded",
                                         options={"xatol": 1e-8})
                    row[nm] = {"R at minimum": float(r2.x), "log smin": float(r2.fun)}
                unconfirmed.append(row)
                print(f"  {sector_name(par, e)}: dip at {r1:.6f} ({res.fun:.1f}) not confirmed; full system: 0.96 Y0 {row['0.96 Y0']['R at minimum']:.6f} "
                      f"({row['0.96 Y0']['log smin']:.1f}), 0.88 Y0 {row['0.88 Y0']['R at minimum']:.6f} ({row['0.88 Y0']['log smin']:.1f})", flush=True)
                continue
            name, Rr, agree, hecke, cp, c1 = got
            new = abs(abs(cp[11]) - 1) < 1e-4 and abs(abs(cp[13]) - 1) < 1e-4
            found.append({"subsector": sector_name(par, e), "parity": "even" if par == 0 else "odd", "eps11": e[0], "eps13": e[1], "R": Rr,
                          "lambda": 0.25 + Rr * Rr, "new": bool(new), "agree": agree, "hecke": hecke, "c2": float(c1[1]), "c3": float(c1[2]),
                          "c11_sqrt11": cp[11], "c13_sqrt13": cp[13], "confirmed on pair": name, "detector dip": r1, "log smin at dip": float(res.fun)})
            print(f"  {sector_name(par, e)}: R = {Rr:.10f}  {'new' if new else 'OLD'}  pair {name}  agree {agree:.1e}  hecke {hecke:.1e}  "
                  f"c11 {cp[11]:+.4f}  c13 {cp[13]:+.4f}", flush=True)
        print(f"  [{sector_name(par, e)} done, {time.time() - T0:.0f} s]", flush=True)
    json.dump({"found": sorted(found, key=lambda x: x["R"]), "unconfirmed dips": unconfirmed, "seconds": time.time() - T0}, open(out, "w"), indent=1)
    print(f"find: {len(found)} forms, {len(unconfirmed)} unconfirmed dips, {time.time() - T0:.0f} s")


def fine(which, R0, R1, dR, subs, out):
    """log(sigma_min / sigma_max) of the full system on one horocycle of pair A (which = 0: 0.96 Y0, 1: 0.88 Y0) on a fine grid,
    for the listed subsectors (0 to 7: even ++, +-, -+, --, then odd): a form dips at the same R on both horocycles."""
    tw = _tw()
    m = tw.m
    Y0, M, hcs = tw.setup(0.96, 0.88)
    fr = np.arange(1, M + 1, dtype=float)
    hc = hcs[which]
    Rs = np.arange(R0, R1 + 1e-9, dR)
    res = {}
    for s in subs:
        par, e = SUB[s]
        chi = hc.chi(e, twist=False)
        res[sector_name(par, e)] = [float(math.log(smin_sys(m, float(R), hc, fr, par, chi))) for R in Rs]
    json.dump({"horocycle": ["0.96 Y0", "0.88 Y0"][which], "full-system coefficients": M, "Rs": [float(r) for r in Rs], "log smin": res}, open(out, "w"))
    for k, v in res.items():
        i = int(np.argmin(v))
        print(f"  {['0.96 Y0', '0.88 Y0'][which]}  {k}: minimum {v[i]:.2f} at R = {Rs[i]:.4f}; median {np.median(v):.2f}; maximum {max(v):.2f}", flush=True)


# ------------------------------------------------------------------------------------------------- the final lists: P3, P4, P6
def final(listfile, out, Rtop=1.8):
    F = json.load(open(listfile))["found"]
    lst = {(par, e): sorted(x["R"] for x in F if x["parity"] == ("even" if par == 0 else "odd") and (x["eps11"], x["eps13"]) == e) for par, e in SUB}
    lev = Level(143)
    rows = []
    for r0, sig in TESTS + ((0.424, 0.35),):
        for s in SECT:
            eps = {11: s[0], 13: s[1]}
            ev, p, h = lev.even(r0, sig, eps)
            od, _, _ = odd_sum(lev, r0, sig, eps)
            t_ = tau(h, Rtop)
            le, lo = sum(float(h(R)) for R in lst[(0, s)]), sum(float(h(R)) for R in lst[(1, s)])
            rows.append({"r0": r0, "sigma": sig, "sector": s, "tau": t_,
                         "even": {"geometry": ev, "list": le, "list - geometry": le - ev, "pass": abs(le - ev) <= t_},
                         "odd": {"geometry": od, "list": lo, "list - geometry": lo - od, "pass": abs(lo - od) <= t_}, "even parts": p})
            print(f"(r0 {r0}, sigma {sig}) {s}: even geometry {ev:+.9f} list {le:+.9f} diff {le - ev:+.2e} | odd geometry {od:+.9f} list {lo:+.9f} "
                  f"diff {lo - od:+.2e} | tau {t_:.1e}  {'PASS' if rows[-1]['even']['pass'] and rows[-1]['odd']['pass'] else 'MISS'}", flush=True)
    ev, p, h = lev.even(0.0, 0.35, {11: -1, 13: -1})
    R1 = solve_R1(ev, lst[(0, (-1, -1))][1:], 0.0, 0.35)
    print(f"P3 on the final list: R1 = {R1:.10f}, from 0.4449803610 by {R1 - 0.4449803610:+.2e}")
    cut = []
    for Lc in (3, 4, 5, 6, 7, 8, 10, 12, 14, 17, 20, 23.5):
        evL, _, _ = lev.even(0.0, 0.35, {11: -1, 13: -1}, L=Lc)
        try:
            cut.append((Lc, solve_R1(evL, lst[(0, (-1, -1))][1:], 0.0, 0.35)))
        except ValueError:
            cut.append((Lc, None))
    print("R1 of even (-,-) with the geodesic sums cut at L:", cut)
    json.dump({"R_top": Rtop, "lists": {sector_name(par, e): v for (par, e), v in lst.items()}, "rows": rows,
               "P3 on the final list": {"R1": R1, "difference": R1 - 0.4449803610}, "R1 against the cut L": cut}, open(out, "w"), indent=1)
    print(f"final: {time.time() - T0:.0f} s")


if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd == "checks":
        checks()
    elif cmd == "level1":
        level1()
    elif cmd == "level11":
        level11()
    elif cmd == "first":
        first()
    elif cmd == "scan":
        scan(float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]), sys.argv[5])
    elif cmd == "find":
        find(sys.argv[2], sys.argv[3], sys.argv[4:])
    elif cmd == "fine":
        fine(int(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5]), [int(v) for v in sys.argv[6].split(",")], sys.argv[7])
    elif cmd == "final":
        final(sys.argv[2], sys.argv[3], float(sys.argv[4]) if len(sys.argv) > 4 else 1.8)
