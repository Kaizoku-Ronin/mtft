#!/usr/bin/env python3
"""MP21, T-G7: the shear of the flux tube at its core line.  Plan: ../MP21_TG7_Shear_Plan.md (written 14 November 2026, 08:27
New York time, before this file existed).

    python3 -I mp21_tg7.py        ->  mp21_numbers.json next to this file

The closed form under test:  d(iota)/du at the core line = nu [3pq(p+q) - 2(p^2 + pq + q^2)] / (12 pq (p+q)),  u = 1 - h/h_max.
E1, E2 exact (sympy); P1 the profile integral at 30 digits (mpmath); P2 magnetic lines integrated in space on the field
B = Im F_pq itself, at t = 0 and t = 0.8; P3 is mp21_core.js (the page's own profile function).  A failed check prints MISS."""
import json
import math
import pathlib

import mpmath as mp
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

HERE = pathlib.Path(__file__).resolve().parent
NPASS = NFAIL = 0
OUT = {}


def ok(name, cond, detail=""):
    global NPASS, NFAIL
    cond = bool(cond)
    NPASS += cond
    NFAIL += not cond
    print(("  ok   " if cond else "  MISS ") + name + (("   [" + detail + "]") if detail else ""), flush=True)
    return cond


nu_of = lambda p, q: math.sqrt(2 * p * q * (p + q))
hmax_of = lambda p, q: (p / (p + q)) ** (p / 2) * (q / (p + q)) ** (q / 2)
shear = lambda p, q: nu_of(p, q) * (3 * p * q * (p + q) - 2 * (p * p + p * q + q * q)) / (12 * p * q * (p + q))

# ---------------------------------------------------------------- E: the derivation, exactly
print("E the derivation")
eta, Th, x = sp.symbols("eta Theta x", real=True)
P, Q = sp.symbols("p q", positive=True)
# MP09's line field on the label sphere, and the reduced pair (x = ln cos^2 eta, Theta) with xi_1 as the time
d_eta = -sp.sin(Th) * sp.sin(eta) * sp.cos(eta)
d_xi1 = sp.cos(Th) * sp.sin(eta) ** 2
d_Th = sp.cos(Th) * (P * sp.sin(eta) ** 2 - Q * sp.cos(eta) ** 2)
dx = sp.simplify(sp.diff(sp.log(sp.cos(eta) ** 2), eta) * d_eta / d_xi1)                   # dx/dxi_1
dT = sp.simplify(d_Th / d_xi1)                                                             # dTheta/dxi_1
lam = P * x + Q * sp.log(1 - sp.exp(x))
Ham = lam + 2 * sp.log(sp.cos(Th))                                                         # 2 ln H
sub = {eta: sp.acos(sp.exp(x / 2))}
e1a = sp.simplify((dx - 2 * sp.tan(Th)).subs(sub)) == 0
e1b = sp.simplify(sp.simplify(dT.subs(sub)) - (P - Q * sp.exp(x) / (1 - sp.exp(x)))) == 0
e1c = sp.simplify(-sp.diff(Ham, Th) - 2 * sp.tan(Th)) == 0 and sp.simplify(sp.diff(Ham, x) - (P - Q * sp.exp(x) / (1 - sp.exp(x)))) == 0
ok("E1 dx/dxi_1 = 2 tan(Theta), dTheta/dxi_1 = p - q cos^2/sin^2, and these are Hamilton's equations for 2 ln H (EXACT, symbolic p, q)", e1a and e1b and e1c)
x0 = sp.log(P / (P + Q))
l2, l3, l4 = [sp.simplify(sp.diff(lam, x, k).subs(x, x0)) for k in (2, 3, 4)]
a, b, c = sp.simplify(-l2 / 2), sp.simplify(-l3 / 6), sp.simplify(-l4 / 24)
n = P + Q
e2a = (sp.simplify(a - P * n / (2 * Q)) == 0 and sp.simplify(b - P * n * (2 * P + Q) / (6 * Q ** 2)) == 0
       and sp.simplify(c - P * n * (n ** 2 + 4 * P * n + P ** 2) / (24 * Q ** 3)) == 0)
beta4 = sp.Rational(1, 6)                                                                   # -2 ln cos(Theta) = Theta^2 + Theta^4/6 + ...
ser = sp.series(-2 * sp.log(sp.cos(Th)), Th, 0, 6).removeO()
e2b = sp.simplify(ser - (Th ** 2 + Th ** 4 / 6)) == 0
# area A = (pi eps / sqrt(a)) [1 + eps K], K = 15 b^2/(32 a^3) - 3c/(8a^2) - 3 beta/8;  P = dA/d eps;  iota = 2 pi q / P;  eps = 2u
K = sp.Rational(15, 32) * b ** 2 / a ** 3 - sp.Rational(3, 8) * c / a ** 2 - sp.Rational(3, 8) * beta4
nu = sp.sqrt(2 * P * Q * n)
slope_sym = sp.simplify(2 * sp.pi * Q / (sp.pi / sp.sqrt(a)) * (-2 * K) * 2)               # d iota/du = (2 pi q sqrt(a)/pi) (-2K) (d eps/du = 2)
closed = nu * (3 * P * Q * n - 2 * (P ** 2 + P * Q + Q ** 2)) / (12 * P * Q * n)
e2c = sp.simplify(slope_sym - closed) == 0 and sp.simplify(2 * Q * sp.sqrt(a) - nu) == 0
ok("E2 the coefficients a, b, c, the Theta^4/6 term, and the bracket of the closed form from the area's second-order term (EXACT, symbolic p, q)", e2a and e2b and e2c)
# the area's second-order coefficient itself, checked exactly on a polynomial case by the same polar expansion, order by order
rho, phi, eps = sp.symbols("rho phi epsilon", positive=True)
aa, bb, cc, be = sp.symbols("a b c beta", positive=True)
r1 = -bb / aa ** sp.Rational(3, 2) * sp.cos(phi) ** 3
r2 = sp.Rational(3, 2) * (bb ** 2 / aa ** 3) * sp.cos(phi) ** 6 - (cc / aa ** 2) * sp.cos(phi) ** 4 - be * sp.sin(phi) ** 4
rho2 = 1 + sp.sqrt(eps) * r1 + eps * r2                                                     # rho^2 on the level curve, X = rho cos, Y = rho sin, scaled
lhs = rho2 + sp.sqrt(eps) * (bb / aa ** sp.Rational(3, 2)) * rho2 ** sp.Rational(3, 2) * sp.cos(phi) ** 3 + eps * rho2 ** 2 * ((cc / aa ** 2) * sp.cos(phi) ** 4 + be * sp.sin(phi) ** 4)
res = sp.series(lhs - 1, eps, 0, sp.Rational(3, 2)).removeO()
area2 = sp.integrate(r2, (phi, 0, 2 * sp.pi)) / 2
e2d = sp.simplify(res) == 0 and sp.simplify(area2 / sp.pi - (sp.Rational(15, 32) * bb ** 2 / aa ** 3 - sp.Rational(3, 8) * cc / aa ** 2 - sp.Rational(3, 8) * be)) == 0
ok("E2' the level curve to order eps and its area: A = (pi eps/sqrt(a)) [1 + eps (15 b^2/(32 a^3) - 3c/(8 a^2) - 3 beta/8)] (EXACT)", e2d)

# ---------------------------------------------------------------- P1: the profile integral at 30 digits
print("P1 the profile integral")
mp.mp.dps = 34


def iota_mp(p, q, u):
    """MP09's integral: iota = 2 pi q / int dx / sqrt(e^{lambda(x) - 2 ln h} - 1) between the two roots, h = h_max (1 - u)"""
    p, q = mp.mpf(p), mp.mpf(q)
    lw = lambda xx: p * xx + q * mp.log1p(-mp.exp(xx))
    xc = mp.log(p / (p + q))
    l2h = lw(xc) + 2 * mp.log1p(-mp.mpf(u))
    w = mp.sqrt(-2 * mp.log1p(-mp.mpf(u)) / (p * (p + q) / (2 * q)))                        # the harmonic half-width, as a first guess
    hi = min(xc + 10 * w, mp.mpf("-1e-30"))                                                  # lambda rises to x_c and falls after it: one root each side
    x1 = mp.findroot(lambda xx: lw(xx) - l2h, (xc - 10 * w, xc), solver="illinois")
    x2 = mp.findroot(lambda xx: lw(xx) - l2h, (xc, hi), solver="illinois")
    assert x1 < xc < x2 and abs(lw(x1) - l2h) < mp.mpf("1e-28") and abs(lw(x2) - l2h) < mp.mpf("1e-28")
    mid, half = (x1 + x2) / 2, (x2 - x1) / 2
    dmax = lw(xc) - l2h
    slope = lambda xe: abs(p - q * mp.exp(xe) / (1 - mp.exp(xe)))                            # |d lambda/dx| at an end

    def f(th):
        d = lw(mid + half * mp.sin(th)) - l2h
        if d < mp.mpf("1e-22") * dmax:                                                       # the last sliver by its tangent, free of cancellation
            d = slope(x2 if th > 0 else x1) * half * 2 * mp.sin((mp.pi / 2 - abs(th)) / 2) ** 2
        return half * mp.cos(th) / mp.sqrt(mp.expm1(d))
    return 2 * mp.pi * q / mp.quad(f, mp.linspace(-mp.pi / 2, mp.pi / 2, 9))


def richardson3(f1, f2, f4):
    """f(u), f(2u), f(4u) with f = f0 + c1 u + c2 u^2 + ...: f0 to O(u^3)"""
    return (8 * f1 - 6 * f2 + f4) / 3


P1 = {}
worst = 0
for p, q in ((1, 2), (2, 3), (3, 2), (2, 5), (1, 8), (3, 5), (1, 49)):
    nu_ = mp.sqrt(2 * p * q * (p + q))
    s = [(iota_mp(p, q, uu) - nu_) / mp.mpf(uu) for uu in ("1e-3", "2e-3", "4e-3")]
    ext = richardson3(*s)
    S = nu_ * (3 * p * q * (p + q) - 2 * (p * p + p * q + q * q)) / (12 * p * q * (p + q))
    P1["%d,%d" % (p, q)] = {"closed form": float(S), "extrapolated": float(ext), "difference / nu": float((ext - S) / nu_), "at u = 0.001, 0.002, 0.004": [float(v) for v in s]}
    worst = max(worst, abs(float((ext - S) / nu_)))
hop = max(abs(iota_mp(1, 1, uu) - 2) for uu in ("1e-3", "2e-3", "4e-3"))
OUT["P1"] = {"pairs": P1, "worst difference / nu": worst, "hopfion |iota - 2|": float(hop)}
ok("P1 PREDICTION (iota(u) - nu)/u at u = 1e-3, 2e-3, 4e-3, extrapolated, equals the closed form to 1e-8 of nu for seven pairs; the hopfion's iota - 2 is zero to 1e-20",
   worst < 1e-8 and hop < mp.mpf("1e-20"), "worst %.1e of nu; hopfion %.1e" % (worst, float(hop)))
# finer levels, to see the extrapolation settle (recorded, not a registered test)
fine = {}
for p, q in ((2, 3), (1, 8)):
    nu_ = mp.sqrt(2 * p * q * (p + q))
    s = [(iota_mp(p, q, uu) - nu_) / mp.mpf(uu) for uu in ("1e-5", "2e-5", "4e-5")]
    fine["%d,%d" % (p, q)] = float(richardson3(*s) - mp.mpf(shear(p, q)))
OUT["P1"]["extrapolated from u = 1e-5, 2e-5, 4e-5, minus the closed form"] = fine
# the second-order coefficient, measured (recorded for F1)
sec = {}
for p, q in ((1, 2), (2, 3), (3, 2), (2, 5), (1, 8)):
    nu_ = mp.sqrt(2 * p * q * (p + q))
    S = nu_ * (3 * p * q * (p + q) - 2 * (p * p + p * q + q * q)) / (12 * p * q * (p + q))
    g = [((iota_mp(p, q, uu) - nu_) / mp.mpf(uu) - S) / mp.mpf(uu) for uu in ("1e-3", "2e-3", "4e-3")]
    sec["%d,%d" % (p, q)] = float(richardson3(*g))
OUT["P1"]["second-order coefficient c2 in iota = nu + S u + c2 u^2 (measured)"] = sec

# ---------------------------------------------------------------- P2: lines in space
print("P2 lines in space")


def ab(tv, X):
    X = np.asarray(X, float)
    return X[..., 0] - 1j * X[..., 1], tv - 1j - X[..., 2]


def dd(tv, X):
    X = np.asarray(X, float)
    return np.sum(X * X, -1) - (tv - 1j) ** 2


def labels(tv, X):
    a_, b_ = ab(tv, X)
    d = dd(tv, X)
    return 1 - 2j * b_ / d, 2 * a_ / d


def nvec(tv, X):
    a_, b_ = ab(tv, X)
    return np.stack([4 * (b_ * b_ - a_ * a_), -4j * (a_ * a_ + b_ * b_), 8 * a_ * b_], -1)


def Fpq(tv, X, p=1, q=1):
    al, be_ = labels(tv, X)
    return (p * q * al ** (p - 1) * be_ ** (q - 1))[..., None] * nvec(tv, X) / dd(tv, X)[..., None] ** 3


def chart(tv, al, be_):
    w = 1 / (1 - al)
    tau = (tv - w.imag) / w.real
    xy = be_ * (1 + 1j * tau) * w
    return np.stack([xy.real, -xy.imag, tv - tau], -1)


def bun(tv, X, p, q):
    B = Fpq(tv, X, p, q).imag
    return B / np.linalg.norm(B, axis=-1)[..., None]


def roots(p, q, h):
    lw = lambda xx: p * xx + q * np.log1p(-np.exp(xx))
    l2h = 2 * math.log(h)
    x1 = brentq(lambda xx: lw(xx) - l2h, -700.0, math.log(p / (p + q)), xtol=1e-15)
    v2 = brentq(lambda v: p * math.log1p(-math.exp(v)) + q * v - l2h, -700.0, math.log(q / (p + q)), xtol=1e-15)
    return x1, math.log1p(-math.exp(v2))


def iota_space(p, q, tv, u, nret=6):
    """a magnetic line of the torus H = h_max (1 - u) at the time tv, integrated in space; iota from the advance of phi = -xi_2/p over nret returns to Theta = 0"""
    h = hmax_of(p, q) * (1 - u)
    x1, x2 = roots(p, q, h)
    xs = math.exp(x2)
    X0 = chart(tv, np.array(math.sqrt(xs) + 0j), np.array(math.sqrt(1 - xs) + 0j))
    f = lambda s, X: bun(tv, X, p, q) * (1 + X @ X) / 2

    def ev(s, X):
        al, be_ = labels(tv, X)
        return float((al ** p * be_ ** q).imag)
    ev.direction = -1 if ev(0, X0 + 1e-7 * f(0, X0)) < 0 else 1
    got, t0, Xc, ts, ph = [], 0.0, X0, [], []
    while len(got) < nret and t0 < 400.0:
        so = solve_ivp(f, (t0, t0 + 2.0), Xc, method="DOP853", rtol=1e-12, atol=1e-14, events=ev, dense_output=True)
        got += [(te, ye) for te, ye in zip(so.t_events[0], so.y_events[0]) if te > 1e-6]
        tt = np.linspace(so.t[0], so.t[-1], 1200, endpoint=False)
        ts.append(tt)
        ph.append(np.angle(labels(tv, so.sol(tt).T)[1]))
        t0, Xc = so.t[-1], so.y[:, -1]
    ts = np.concatenate(ts)
    x2u = np.unwrap(np.concatenate(ph))
    te, ye = got[nret - 1]
    w0 = float(np.angle(labels(tv, X0)[1]))
    w1 = float(np.angle(labels(tv, ye)[1]))
    w1 += 2 * np.pi * round((np.interp(te, ts, x2u) - w1) / (2 * np.pi))
    dphi = -(w1 - w0) / p
    Hs = [float((labels(tv, ye_)[0] ** p * labels(tv, ye_)[1] ** q).real) for _, ye_ in got[:nret]]
    return 2 * np.pi * nret / abs(dphi), max(abs(v / h - 1) for v in Hs)


P2 = {}
worst2, worst_t, hop2, drift = 0.0, 0.0, 0.0, 0.0
US = (0.004, 0.008, 0.016)
for p, q in ((2, 3), (3, 2), (2, 5), (1, 8)):
    row = {}
    for tv in (0.0, 0.8):
        io = []
        for u in US:
            v, dr = iota_space(p, q, tv, u)
            io.append(v)
            drift = max(drift, dr)
        s = [(io[k] - nu_of(p, q)) / US[k] for k in range(3)]
        ext = richardson3(*s)
        row["t = %.1f" % tv] = {"iota at u = 0.004, 0.008, 0.016": io, "extrapolated slope": ext, "difference / nu": (ext - shear(p, q)) / nu_of(p, q),
                                "profile integral at the same levels": [float(iota_mp(p, q, u)) for u in US]}
        worst2 = max(worst2, abs(ext - shear(p, q)) / nu_of(p, q))
    worst_t = max(worst_t, abs(row["t = 0.0"]["extrapolated slope"] - row["t = 0.8"]["extrapolated slope"]) / nu_of(p, q))
    row["closed form"] = shear(p, q)
    P2["%d,%d" % (p, q)] = row
    print("    (%d,%d): closed form %.7f; in space %.7f (t = 0), %.7f (t = 0.8)" % (p, q, shear(p, q), row["t = 0.0"]["extrapolated slope"], row["t = 0.8"]["extrapolated slope"]), flush=True)
for tv in (0.0, 0.8):
    for u in US:
        v, dr = iota_space(1, 1, tv, u)
        hop2 = max(hop2, abs(v - 2))
OUT["P2"] = {"pairs": P2, "worst difference / nu": worst2, "worst difference between the two times / nu": worst_t, "hopfion |iota - 2|": hop2, "worst drift of H at a return": drift}
ok("P2 PREDICTION the slope from lines integrated in space at t = 0 and t = 0.8 (u = 0.004, 0.008, 0.016, extrapolated) is the closed form to 1e-5 of nu for (2,3), (3,2), (2,5), (1,8); hopfion |iota - 2| < 1e-8",
   worst2 < 1e-5 and worst_t < 1e-5 and hop2 < 1e-8, "worst %.1e of nu; between the times %.1e; hopfion %.1e; drift of H %.1e" % (worst2, worst_t, hop2, drift))

# ---------------------------------------------------------------- F1: the next term
print("F1 the size of the next term")
F1 = {}
for p, q in ((2, 3), (3, 2), (2, 5), (1, 8)):
    u = 0.016
    rem = float(iota_mp(p, q, u)) - nu_of(p, q) - shear(p, q) * u
    F1["%d,%d" % (p, q)] = {"iota - nu - S u at u = 0.016": rem, "3 S u^2": 3 * shear(p, q) * u * u, "ratio to S u^2": rem / (shear(p, q) * u * u)}
OUT["F1"] = F1
f1 = all(abs(v["iota - nu - S u at u = 0.016"]) < v["3 S u^2"] for v in F1.values())
print("  %s F1 FORECAST |iota - nu - S u| < 3 S u^2 at u = 0.016 for the four knotted pairs   [ratios to S u^2: %s]" % ("held  " if f1 else "missed", ", ".join("%.3f" % v["ratio to S u^2"] for v in F1.values())))
OUT["F1 held"] = f1
OUT["table"] = {"%d,%d" % (p, q): {"nu": nu_of(p, q), "d iota/du": shear(p, q), "h_max": hmax_of(p, q), "d iota/dh": -shear(p, q) / hmax_of(p, q)} for p, q in ((1, 1), (1, 2), (2, 3), (3, 2), (2, 5), (1, 8))}
OUT["passed"], OUT["failed"] = NPASS, NFAIL
(HERE / "mp21_numbers.json").write_text(json.dumps(OUT, indent=1))
print("%d passed, %d failed" % (NPASS, NFAIL))
