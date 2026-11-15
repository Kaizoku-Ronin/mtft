#!/usr/bin/env python3
"""MP17: the untwisted odd trace B(h) = sum over odd Maass cusp forms of h(r_j) on X0(N), N squarefree and Gamma0(N) without
elliptic points (N = 143 and N = 11 here), and with MP16's twists the trace of each odd sector.
Plan: ../MP17_Odd_Trace_Plan.md (saved 18:41:25, 12 Nov 2026; Amendment 1 at 18:42).

B(h) = (1/2)(T_1 - T_R); the continuous spectrum cancels exactly (all Eisenstein series are even).
  identity:     (vol/8 pi) int h(r) r tanh(pi r) dr
  geodesics:    (1/2) sum_t w_t g(l_t),   2 cosh(l_t/2) = t >= 3,  w_t = 2 W_1(t^2 - 4)   [MP16's Q = 1 weights]
  glides:      -(1/2) sum_t w'_t g(l'_t), 2 sinh(l'_t/2) = t >= 1, w'_t = 2 W_1(t^2 + 4)
  cusps and the cusp-to-cusp mirrors of R (every reflection in R Gamma0(N) has discriminant 4):
      C = (h_c/2)[g(0)(log 2 + gamma_E) + 4K - 2J] - (g(0)/4) sum_sigma beta_sigma
      K = int_0^inf k(s^2) log s ds = -(1/4) int_0^inf g'(r) log(sinh(r/2)/2) dr   (the plan's 1/2 was a slip; V7 caught it)
      J = int_0^inf k(v^2) log(1 + v^2) dv = -int_0^inf g'(r) log cosh(r/4) dr
      beta_sigma = log(w1 w2) + 2 log|p2 q1 - p1 q2| for the mirror's cusps p1/q1, p2/q2 of widths w1, w2;
      sum_sigma beta_sigma = (1/2) sum over the 2 h_c mirror ends.
Sector traces: S_chi = B/4 + e11 A11 + e13 A13 + e143 A143, A_Q = (T_{W_Q} - T_{RW_Q})/8 (MP16).

Usage: python3 -I mp17_trace.py v7 | ends | level11 | level143 | locate
Needs mp16_trace.py one folder over (../mp16), numpy, scipy, mpmath, cypari2."""
import importlib.util
import json
import math
import os
import sys
import time
from fractions import Fraction

import mpmath as mp
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("t16", os.path.join(HERE, "..", "mp16", "mp16_trace.py"))
t16 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(t16)
pari = t16.pari
T0 = time.time()
EULER = 0.57721566490153286061


def dump(obj, name):
    with open(os.path.join(HERE, name), "w") as fh:
        json.dump(obj, fh, indent=1)


# ------------------------------------------------------------------------------------------------- K and J
def KJ_from_g(gp):
    """K and J from g'(r) (a callable), by the Abel-transform formulas."""
    K = -0.25 * quad(lambda r: gp(r) * math.log(math.sinh(r / 2) / 2), 0, 60, limit=400, epsabs=1e-14, epsrel=1e-12)[0]   # 1/4: see the plan's Amendment 2
    J = -quad(lambda r: gp(r) * math.log(math.cosh(r / 4)), 0, 60, limit=400, epsabs=1e-14, epsrel=1e-12)[0]
    return K, J


def v7():
    """K and J for k(u) = exp(-a u): directly from k, and from g(r) = 2 sqrt(pi/a) exp(-a sinh^2(r/2))."""
    out = {}
    for a in (1.0, 3.0):
        Kd = quad(lambda s: math.exp(-a * s * s) * math.log(s), 0, np.inf, limit=400, epsabs=1e-14)[0]
        Jd = quad(lambda v: math.exp(-a * v * v) * math.log1p(v * v), 0, np.inf, limit=400, epsabs=1e-14)[0]
        gp = lambda r: -2 * math.sqrt(math.pi / a) * a * math.sinh(r / 2) * math.cosh(r / 2) * math.exp(-a * math.sinh(r / 2) ** 2)
        Kg, Jg = KJ_from_g(gp)
        out[f"a={a}"] = {"K direct": Kd, "K from g": Kg, "J direct": Jd, "J from g": Jg, "K closed form (a=1)": -(math.sqrt(math.pi) / 4) * (EULER + 2 * math.log(2)) if a == 1 else None}
        print(f"k = exp(-{a}u): K {Kd:.15f} / {Kg:.15f}   J {Jd:.15f} / {Jg:.15f}", flush=True)
    dump(out, "mp17_v7.json")


# ------------------------------------------------------------------------------------------------- cusps and mirror ends
def cusp_data(N):
    """For each cusp class 1/d (d | N; d = N is the class of infinity): width w = N/d and A in SL(2,Z) with first column (1, d)
    and N | 2 d v, so that A diag(-1,1) A^-1 lies in R Gamma0(N) and the cusp's mirrors are x = 0 and x = 1/2 locally."""
    ps = t16.primes_of(N)
    divs = sorted({math.prod(s) for k in range(len(ps) + 1) for s in __import__("itertools").combinations(ps, k)})
    out = []
    for d in divs:
        w = N // d
        p, q = 1, d
        u = 0
        while True:                                   # v = 1 + u q, need N | 2 q v
            v = 1 + u * q
            if (2 * q * v) % N == 0:
                break
            u += 1
        out.append({"cusp": f"1/{d}" if d < N else "1/N ~ inf", "p": p, "q": q, "u": u, "v": v, "w": w})
    return out


def reduce_frac(a, b):
    if b < 0:
        a, b = -a, -b
    g = math.gcd(a, b)
    return (a // g, b // g) if g else (1, 0)


def mirror_ends(N):
    ends, checks = [], []
    for c in cusp_data(N):
        p, q, u, v, w = c["p"], c["q"], c["u"], c["v"], c["w"]
        A = np.array([[p, u], [q, v]], dtype=object)
        Ainv = np.array([[v, -u], [-q, p]], dtype=object)
        for j, L in ((0, np.array([[-1, 0], [0, 1]], dtype=object)), (1, np.array([[-1, w], [0, 1]], dtype=object))):
            M = A.dot(L).dot(Ainv)                    # the reflection, global coordinates (acts on conj z)
            det = M[0, 0] * M[1, 1] - M[0, 1] * M[1, 0]
            fix = (M[0, 0] * p + M[0, 1] * q, M[1, 0] * p + M[1, 1] * q)
            checks.append({"cusp": c["cusp"], "mirror": "x = 0" if j == 0 else "x = 1/2", "matrix": [[int(M[0, 0]), int(M[0, 1])], [int(M[1, 0]), int(M[1, 1])]],
                           "det": int(det), "trace": int(M[0, 0] + M[1, 1]), "N | lower-left": int(M[1, 0]) % N == 0,
                           "fixes the cusp": fix[0] * q == fix[1] * p})
            # the other end: A applied to the real point w*j/2 ... in local coordinates the mirror x = j/2 ends at x = j/2,
            # i.e. at sigma_a(j/2) = A(w j / 2)
            num, den = (u, v) if j == 0 else (p * w + 2 * u, q * w + 2 * v)
            p2, q2 = reduce_frac(num, den)
            w2 = N // math.gcd(q2, N) if q2 else 1
            cval = abs(p2 * q - p * q2)
            beta = math.log(w * w2) + 2 * math.log(cval)
            ends.append({"cusp": c["cusp"], "width": w, "mirror": "x = 0" if j == 0 else "x = 1/2", "other end": f"{p2}/{q2}", "other width": w2,
                         "|p2 q1 - p1 q2|": cval, "beta": beta})
    return ends, checks


def ends():
    out = {}
    for N in (143, 11):
        e, c = mirror_ends(N)
        out[str(N)] = {"ends": e, "matrix checks": c, "sum beta over classes": 0.5 * sum(x["beta"] for x in e)}
        print(f"level {N}:", flush=True)
        for x in e:
            print(f"   {x}", flush=True)
        print(f"   all matrix checks pass: {all(ch['det'] == -1 and ch['trace'] == 0 and ch['N | lower-left'] and ch['fixes the cusp'] for ch in c)}")
    dump(out, "mp17_ends.json")


# ------------------------------------------------------------------------------------------------- the odd trace
class OddTrace:
    def __init__(self, N, lmax):
        self.N = N
        ps = t16.primes_of(N)
        self.index = N * math.prod(1 + Fraction(1, p) for p in ps)
        self.vol = math.pi / 3 * float(self.index)
        e, _ = mirror_ends(N)
        self.hc = len(e) // 2
        self.sum_beta = 0.5 * sum(x["beta"] for x in e)
        cache = os.path.join(HERE, f"weights_{N}_1.npz")
        if os.path.exists(cache):
            z = np.load(cache)
            assert float(z["lmax"]) == lmax
            self.hyp, self.gli = z["hyp"], z["gli"]
        else:
            hyp, gli = t16.geodesic_weights(N, 1, lmax)
            self.hyp = np.array([r[:3] for r in hyp], dtype=float)
            self.gli = np.array([r[:3] for r in gli], dtype=float)
            np.savez(cache, lmax=lmax, hyp=self.hyp, gli=self.gli)

    def parts(self, h, g, gp, L=None):
        ident = self.vol / (8 * math.pi) * 2 * quad(lambda r: h(r) * r * math.tanh(math.pi * r), 0, 40, limit=400, epsabs=1e-14, epsrel=1e-13)[0]
        hy, gl = self.hyp, self.gli
        if L is not None:
            hy, gl = hy[hy[:, 1] <= L], gl[gl[:, 1] <= L]
        geo = 0.5 * float(np.sum(g(hy[:, 1]) * hy[:, 2]))
        glide = -0.5 * float(np.sum(g(gl[:, 1]) * gl[:, 2]))
        K, J = KJ_from_g(gp)
        g0 = float(g(0.0))
        cusp = (self.hc / 2) * (g0 * (math.log(2) + EULER) + 4 * K - 2 * J) - g0 / 4 * self.sum_beta
        return {"identity": ident, "geodesics": geo, "glides": glide, "cusps and cusp mirrors": cusp}


def testfun(r0, sig):
    h, g = t16.testfun(r0, sig)
    c = 2 * sig / math.sqrt(2 * math.pi)
    gp = lambda u: c * (-r0 * math.sin(r0 * u) - sig * sig * u * math.cos(r0 * u)) * math.exp(-sig * sig * u * u / 2)
    return h, g, gp


def level11():
    B = OddTrace(11, 23.5)
    rows = []
    for sig in (0.35, 0.45):
        for r0 in (0.0, 0.3, 0.6, 1.0):
            h, g, gp = testfun(r0, sig)
            p = B.parts(h, g, gp)
            rows.append({"sigma": sig, "r0": r0, "B": sum(p.values()), "parts": p})
            print(f"level 11 sigma {sig} r0 {r0}: B = {sum(p.values()):+.3e}   " + "  ".join(f"{k} {v:+.6f}" for k, v in p.items()), flush=True)
    dump({"index": float(B.index), "cusps": B.hc, "sum beta": B.sum_beta, "geodesic classes (by t)": len(B.hyp), "glide classes (by t)": len(B.gli), "rows": rows}, "mp17_level11.json")


# ------------------------------------------------------------------------------------------------- level 143: the four odd sectors
SECT = ((1, 1), (1, -1), (-1, 1), (-1, -1))


def odd_list():
    spec_ = t16.odd_spectrum()
    for R, a, b in ((1.3069477649, -1, 1), (1.3671856490, -1, 1), (1.3427626324, 1, 1), (1.5826864635, 1, -1)):
        spec_.setdefault((a, b), []).append(R)
    return {k: sorted(v) for k, v in spec_.items()}


def twists(h, g, L=None):
    c_ell, c_refl = t16.constants()
    out = {}
    for Q in (11, 13, 143):
        c = t16.Coset(143, Q, 23.5, c_ell, c_refl, cache=os.path.join(HERE, "..", "mp16", f"weights_143_{Q}.npz"))
        if L is not None:
            saved = (c.hyp, c.gli)
            c.hyp, c.gli = [x for x in c.hyp if x[1] <= L], [x for x in c.gli if x[1] <= L]
        out[Q] = c.parts(h, g)
        if L is not None:
            c.hyp, c.gli = saved
    return out


def sector_traces(B143, r0, sig, L=None):
    h, g, gp = testfun(r0, sig)
    bp = B143.parts(h, g, gp, L)
    Bv = sum(bp.values())
    tw = twists(h, g, L)
    A = {Q: sum(p.values()) / 8 for Q, p in tw.items()}
    eps = {11: lambda a, b: a, 13: lambda a, b: b, 143: lambda a, b: a * b}
    S = {s: Bv / 4 + sum(eps[Q](*s) * A[Q] for Q in A) for s in SECT}
    return S, bp, tw, h


def level143():
    B = OddTrace(143, 23.5)
    lst = odd_list()
    rows = []
    for sig, r0 in ((0.35, 0.0), (0.45, 0.0), (0.35, 0.374), (0.35, 0.15)):
        S, bp, tw, h = sector_traces(B, r0, sig)
        tau = quad(lambda r: 3.5 * r * h(r), 1.65, 60, limit=200)[0]
        for s in SECT:
            spec_side = sum(float(h(R)) for R in lst[s])
            rows.append({"sigma": sig, "r0": r0, "sector": s, "geometric": S[s], "spectral (list)": spec_side, "diff": spec_side - S[s], "tau_s": tau,
                         "pass": abs(spec_side - S[s]) <= tau + 1e-6})
            print(f"sigma {sig} r0 {r0:5.3f} sector {s}: geometric {S[s]:+.9f}  list {spec_side:+.9f}  diff {spec_side - S[s]:+.2e}  tau_s {tau:.1e}  {'PASS' if rows[-1]['pass'] else 'FAIL'}", flush=True)
        rows.append({"sigma": sig, "r0": r0, "B parts": bp, "B": sum(bp.values())})
    dump({"index": float(B.index), "vol": B.vol, "cusps": B.hc, "sum beta": B.sum_beta, "geodesic classes (by t)": len(B.hyp), "glide classes (by t)": len(B.gli),
          "odd list": {str(k): v for k, v in lst.items()}, "rows": rows}, "mp17_level143.json")


# ------------------------------------------------------------------------------------------------- P3 and P4: the place of the lowest (-,-) form
def solve_R1(S_geo, others, r0, sig, lo=0.05, hi=0.8):
    h, _ = t16.testfun(r0, sig)
    target = S_geo - sum(float(h(R)) for R in others)
    f = lambda R: float(h(R)) - target
    try:
        return brentq(f, lo, min(hi, r0) if r0 > lo + 0.01 else hi, xtol=1e-14) if r0 > lo + 0.01 else brentq(f, lo, hi, xtol=1e-14)
    except ValueError:
        return None


def locate():
    B = OddTrace(143, 23.5)
    lst = odd_list()
    others = [R for R in lst[(-1, -1)] if R > 0.5]
    out = {"others used (P3a)": others, "P3a": {}, "P3b": {}, "P4": [], "largest terms at r0 = 0.374": None}
    for r0 in (0.0, 0.15, 0.374):
        S, bp, tw, h = sector_traces(B, r0, 0.35)
        R1 = solve_R1(S[(-1, -1)], others, r0, 0.35, hi=0.8) if r0 == 0 else None
        if r0 > 0:                                    # below the centre h increases with R; take the root on the side of the known value
            hh, _ = t16.testfun(r0, 0.35)
            target = S[(-1, -1)] - sum(float(hh(R)) for R in others)
            f = lambda R: float(hh(R)) - target
            roots = []
            for a, b in ((0.05, r0), (r0, 0.9)):
                try:
                    roots.append(brentq(f, a, b, xtol=1e-14))
                except ValueError:
                    pass
            R1 = min(roots, key=lambda x: abs(x - 0.3738518157)) if roots else None
            out["P3a"][f"r0={r0} all roots"] = roots
        out["P3a"][f"r0={r0}"] = {"R1": R1, "minus Hejhal": (R1 - 0.3738518157) if R1 is not None else None, "S_geo": S[(-1, -1)]}
        print(f"P3a r0 {r0}: R1 = {R1}  (Hejhal 0.3738518157)", flush=True)
    # P3b: no spectral input.  Every (-,-) form above R1 is replaced by the sector's smooth part (identity, cusps and cusp
    # mirrors, fixed points, Pell mirrors) applied to h times a smooth step at Rcut; R1 then solves S_geo = h(R1) + that.
    sig, r0 = 0.35, 0.0
    S, bp, tw, hh = sector_traces(B, r0, sig)
    rg = np.linspace(0, 6, 6001)
    ug = np.linspace(0, 60, 3001)
    c_ell, c_refl = t16.constants()
    for Rcut in (0.6, 0.8):
        H = np.asarray(hh(rg)) * 0.5 * (1 + __import__("scipy.special", fromlist=["erf"]).erf((rg - Rcut) / 0.05))
        gH = np.array([np.trapezoid(H * np.cos(rg * u), rg) for u in ug]) / math.pi
        gpH = np.array([-np.trapezoid(H * rg * np.sin(rg * u), rg) for u in ug]) / math.pi
        with np.errstate(divide="ignore"):
            lk = np.where(ug > 0, np.log(np.sinh(np.maximum(ug, 1e-300) / 2) / 2), 0.0)
        KH = -0.25 * np.trapezoid(gpH[1:] * lk[1:], ug[1:])
        JH = -np.trapezoid(gpH * np.log(np.cosh(ug / 4)), ug)
        g0 = gH[0]
        ident = B.vol / (8 * math.pi) * 2 * np.trapezoid(H * rg * np.tanh(math.pi * rg), rg)
        cusp = (B.hc / 2) * (g0 * (math.log(2) + EULER) + 4 * KH - 2 * JH) - g0 / 4 * B.sum_beta
        sech = 2 * np.trapezoid(H / np.cosh(math.pi * rg), rg)
        twist_smooth = (-1) * 0 + (-1) * 4 * c_ell * sech / 8 + (+1) * (20 * c_ell * sech - 12.705252722366623 * c_refl * g0) / 8
        tail = (ident + cusp) / 4 + twist_smooth
        R1b = solve_R1(S[(-1, -1)] - tail, [], r0, sig, hi=Rcut)
        out["P3b"][f"smooth part above {Rcut}"] = {"tail": float(tail), "identity/4": float(ident / 4), "cusps/4": float(cusp / 4), "twisted smooth": float(twist_smooth),
                                                   "R1": R1b, "minus Hejhal": (R1b - 0.3738518157) if R1b else None}
        print(f"P3b (smooth part above {Rcut}): tail {tail:.5f}  R1 = {R1b}", flush=True)
    # P4: cut every geodesic and glide sum (untwisted and twisted) at length L
    for L in (3, 4, 5, 6, 7, 8, 10, 12, 14, 17, 20, 23.5):
        S, _, _, _ = sector_traces(B, 0.0, 0.35, L)
        R1 = solve_R1(S[(-1, -1)], others, 0.0, 0.35)
        out["P4"].append({"L": L, "S_geo": S[(-1, -1)], "R1": R1, "minus Hejhal": (R1 - 0.3738518157) if R1 is not None else None})
        print(f"P4 L {L:5.1f}: S = {S[(-1, -1)]:+.9f}  R1 = {R1}", flush=True)
    # the largest single terms of S_(-,-) at r0 = 0.374
    h, g, gp = testfun(0.374, 0.35)
    S, bp, tw, _ = sector_traces(B, 0.374, 0.35)
    terms = [(f"untwisted {k}", v / 4) for k, v in bp.items() if k in ("identity", "cusps and cusp mirrors")]
    terms += [(f"geodesic t={int(t)} (l={l:.4f})", 0.5 * float(g(l)) * w / 4) for t, l, w in B.hyp[:12]]
    terms += [(f"glide t={int(t)} (l={l:.4f})", -0.5 * float(g(l)) * w / 4) for t, l, w in B.gli[:12]]
    eps = {11: -1, 13: -1, 143: 1}
    c_ell, c_refl = t16.constants()
    for Q, p in tw.items():
        terms += [(f"W{Q} {k}", eps[Q] * v / 8) for k, v in p.items() if k in ("elliptic", "mirror")]
        c = t16.Coset(143, Q, 23.5, c_ell, c_refl, cache=os.path.join(HERE, "..", "mp16", f"weights_143_{Q}.npz"))
        terms += [(f"W{Q} geodesic s={int(x[0])} (l={x[1]:.4f})", eps[Q] * float(g(x[1])) * x[2] / 8) for x in c.hyp[:8]]
        terms += [(f"W{Q} glide s={int(x[0])} (l={x[1]:.4f})", -eps[Q] * float(g(x[1])) * x[2] / 8) for x in c.gli[:8]]
    total = S[(-1, -1)]
    terms.sort(key=lambda x: -abs(x[1]))
    out["largest terms at r0 = 0.374"] = {"total": total, "terms": terms[:14],
                                          "largest single geodesic or glide": max((abs(v) for k, v in terms if ("geodesic" in k or "glide" in k)), default=0) / abs(total)}
    print(json.dumps(out["largest terms at r0 = 0.374"], indent=1), flush=True)
    dump(out, "mp17_locate.json")


if __name__ == "__main__":
    {"v7": v7, "ends": ends, "level11": level11, "level143": level143, "locate": locate}[sys.argv[1]]()
    print(f"done in {time.time() - T0:.0f} s")
