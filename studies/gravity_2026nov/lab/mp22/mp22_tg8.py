#!/usr/bin/env python3
"""MP22, T-G8: the gravitational hopfion's lines of force after t = 0.  Plan: ../MP22_TG8_Hopfion_Plan.md (written 14 November
2026, 08:47 New York time, and filed, before this file existed).

    python3 -I mp22_tg8.py [E P0 P1 P2 P3 P4 P5 P6 P7]   ->  mp22_numbers.json next to this file (sections named, or all)
                                                         (P7 was added by Amendment 1 of the plan, written before P7 was run)

The claim under test: at the time t the lines of force (the lines of Im G and of Re G, G = sqrt(d) F_11) are the orbits, on the
level H = 1, of

    H = |P|^2 + |s^2 - sigma^2|^2 / t^2,   s = (1 - it(1 - alpha))^(1/2),   sigma^2 = 1 - it,   P = conj(beta) or i conj(beta),

and that Hamiltonian is integrable: with zeta = x + iy = s conj(sigma_hat), p = P conj(sigma_hat), m = |sigma|^2,

    I_2 = L^2 + 4 m p_x^2 + (4m/t^2) x^2 (r^2 - 2m)

is a second invariant.  Sections P0 to P6 work in space, on the field itself.  A failed check prints MISS and is counted."""
import json
import math
import pathlib
import sys
import time

import mpmath as mp
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp

HERE = pathlib.Path(__file__).resolve().parent
T0 = time.time()
NPASS = NFAIL = 0
OUTF = HERE / "mp22_numbers.json"
OUT = json.loads(OUTF.read_text()) if OUTF.exists() else {}
X0H = np.array([0.6, -0.3, 0.5])                      # MP10's and MP11's point


def ok(name, cond, detail=""):
    global NPASS, NFAIL
    cond = bool(cond)
    NPASS += cond
    NFAIL += not cond
    print(("  ok   " if cond else "  MISS ") + name + (("   [" + detail + "]") if detail else ""), flush=True)
    return cond


# ---------------------------------------------------------------- the field in space (MP08's conventions)
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


def Gf(tv, X):
    """G = sqrt(d) F_11, F_11 = n / d^3; the principal root, smooth for t >= 0"""
    d = dd(tv, X)
    return (np.sqrt(d + 0j) / d ** 3)[..., None] * nvec(tv, X)


def unit(V):
    return V / np.linalg.norm(V, axis=-1)[..., None]


def field(tv, kind):
    return (lambda X: unit(Gf(tv, X).imag)) if kind == "B" else (lambda X: unit(Gf(tv, X).real))


def chart(tv, al, be):
    w = 1 / (1 - al)
    tau = (tv - w.imag) / w.real
    xy = be * (1 + 1j * tau) * w
    return np.stack([xy.real, -xy.imag, tv - tau], -1)


# ---------------------------------------------------------------- the model
class Model:
    def __init__(self, tv):
        self.t = tv
        self.sig = np.sqrt(1 - 1j * tv)
        self.sh = self.sig / abs(self.sig)
        self.m = abs(self.sig) ** 2
        self.f = 2 * math.sqrt(self.m)
        self.eps = tv / (2 * math.sqrt(1 + tv * tv))
        self.g = 4 / self.eps ** 2
        self.Kmin = 1 / tv ** 2
        self.c0 = 0.5 - 0.25 * math.sqrt(1 + 4 * self.eps ** 2 / 3)
        self.Kmax = self.g * self.c0 * (self.c0 - 0.75) ** 2 - self.c0
        self.I2_axial, self.I2_cross = -4 * self.m * self.Kmin, -4 * self.m * self.Kmax
        self.I2_range = self.I2_axial - self.I2_cross

    def sP(self, al, be, kind):
        s = np.sqrt(1 - 1j * self.t * (1 - al))
        return s, (np.conj(be) if kind == "B" else 1j * np.conj(be))

    def zp(self, al, be, kind):
        s, P = self.sP(al, be, kind)
        return s * np.conj(self.sh), P * np.conj(self.sh)

    def H(self, z, p):
        return abs(p) ** 2 + abs(z * z - self.m) ** 2 / self.t ** 2

    def I2(self, z, p):
        x, y, px, py = z.real, z.imag, p.real, p.imag
        return (x * py - y * px) ** 2 + 4 * self.m * px ** 2 + (4 * self.m / self.t ** 2) * x * x * (x * x + y * y - 2 * self.m)

    def I2_space(self, X, kind):
        return self.I2(*self.zp(*labels(self.t, X), kind))

    def elliptic(self, z):
        """cosh(mu) and cos(nu) of zeta = f cosh(mu + i nu)"""
        a, b = abs(z + self.f), abs(z - self.f)
        return (a + b) / (2 * self.f), (a - b) / (2 * self.f)

    def labels_of(self, z, p, kind):
        """the label point of a phase point (zeta, p)"""
        s, P = z * self.sh, p * self.sh
        al = (s * s - self.sig ** 2) / (1j * self.t)
        be = np.conj(P) if kind == "B" else np.conj(-1j * P)
        return al, be


def periods(K, eps, dps=25):
    """T_mu, T_nu and the roots of Pi(X) = K + X - (4/eps^2) X (X - 3/4)^2, at mp precision"""
    mp.mp.dps = dps
    g = 4 / mp.mpf(eps) ** 2
    K = mp.mpf(K)
    r = sorted(mp.re(v) for v in mp.polyroots([-g, 3 * g / 2, 1 - 9 * g / 16, K], maxsteps=400, extraprec=400))
    X1, X2, X3 = r

    def fn(th):
        X = (X1 + X2) / 2 + (X2 - X1) / 2 * mp.sin(th)
        return 1 / mp.sqrt(g * X * (1 - X) * (X3 - X))

    def fm(th):
        X = (1 + X3) / 2 + (X3 - 1) / 2 * mp.sin(th)
        return 1 / mp.sqrt(g * X * (X - X1) * (X - X2))
    Tn = mp.quad(fn, mp.linspace(-mp.pi / 2, mp.pi / 2, 5))
    Tm = 2 * mp.quad(fm, mp.linspace(-mp.pi / 2, mp.pi / 2, 5))
    return Tm, Tn, (X1, X2, X3)


def W_axial(eps, dps=25):
    mp.mp.dps = dps
    e = mp.mpf(eps)
    Tn = (e / 2) * mp.quad(lambda ph: 1 / ((mp.mpf(3) / 4 - e / 2 * mp.sin(ph)) * mp.sqrt(mp.mpf(1) / 4 + e / 2 * mp.sin(ph))), mp.linspace(-mp.pi / 2, mp.pi / 2, 5))
    Tm = (4 * mp.pi * e / 3) / mp.sqrt(1 - 4 * e * e / 9)
    return Tn / Tm


# ---------------------------------------------------------------- lines in space, at unit speed on the label sphere
def rhs_line(tv, kind, mod=None, tau=False, h=1e-5):
    u = field(tv, kind)

    def f(l, Y):
        X = Y[:3]
        d = u(X)
        al, be = labels(tv, np.stack([X + h * d, X - h * d, X]))
        dal, dbe = (al[0] - al[1]) / (2 * h), (be[0] - be[1]) / (2 * h)
        sp_ = math.sqrt(abs(dal) ** 2 + abs(dbe) ** 2)                 # label length per unit of length in space
        if not tau:
            return d / sp_
        s, P = mod.sP(al[2], be[2], kind)
        ds = (1j * tv / (2 * s)) * dal / sp_
        dP = (np.conj(dbe) if kind == "B" else 1j * np.conj(dbe)) / sp_
        gV = 4 * (s * s - mod.sig ** 2) * np.conj(s) / tv ** 2
        dtau = math.sqrt(abs(ds) ** 2 + abs(dP) ** 2) / math.sqrt(4 * abs(P) ** 2 + abs(gV) ** 2)
        z = s * np.conj(mod.sh)
        return np.concatenate([d / sp_, [2 * dtau / abs(z * z - 4 * mod.m)]])
    return f


def run_line(tv, kind, X0, length, n_eval=None, **kw):
    te = None if n_eval is None else np.linspace(0, length, n_eval)
    return solve_ivp(rhs_line(tv, kind), (0, length), X0, method="DOP853", rtol=1e-12, atol=1e-13, t_eval=te, **kw)


SECT = set(sys.argv[1:]) or {"E", "P0", "P1", "P2", "P3", "P4", "P5", "P6", "P7"}
rng = np.random.default_rng(143)
STARTS = [X0H] + [rng.normal(0, 0.8, 3) for _ in range(4)]

# ================================================================ E: the derivation, exactly
if "E" in SECT:
    print("E the derivation")
    t, A = sp.symbols("t A")
    al = A + sp.I * t / 4 * (A - 1) ** 2
    ok("E1 alpha = A + (it/4)(A - 1)^2 gives (d alpha/dA)^2 = 1 - it(1 - alpha), so dA = lambda_t^(1/2) d alpha (EXACT)",
       sp.expand(sp.diff(al, A) ** 2 - (1 - sp.I * t * (1 - al))) == 0 and sp.diff(al, A).subs(A, 1) == 1)
    x, y, px, py, m, tt = sp.symbols("x y p_x p_y m t", real=True)
    r2 = x * x + y * y
    V = (r2 ** 2 - 2 * m * (x * x - y * y) + m * m) / tt ** 2
    Hs = px ** 2 + py ** 2 + V
    L = x * py - y * px
    I2s = L ** 2 + 4 * m * px ** 2 + (4 * m / tt ** 2) * x * x * (r2 - 2 * m)
    pb = lambda F, G: sum(sp.diff(F, q) * sp.diff(G, p) - sp.diff(F, p) * sp.diff(G, q) for q, p in ((x, px), (y, py)))
    zc = x + sp.I * y
    ok("E2 {H, I_2} = 0 for H = p^2 + [r^4 - 2m(x^2 - y^2) + m^2]/t^2, which is |zeta^2 - m|^2/t^2 (EXACT, symbolic t and m)",
       sp.expand(pb(Hs, I2s)) == 0 and sp.expand(V * tt ** 2 - sp.expand((zc ** 2 - m) * (sp.conjugate(zc) ** 2 - m))) == 0)
    # the second-order normal form of H = |pi|^2 + |u|^2 |1 + eps u/2|^2 at the 1:1 resonance
    a, ab_, b, bb, e = sp.symbols("a abar b bbar epsilon")
    br = lambda F, G: -2 * sp.I * (sp.diff(F, a) * sp.diff(G, ab_) - sp.diff(F, ab_) * sp.diff(G, a) + sp.diff(F, b) * sp.diff(G, bb) - sp.diff(F, bb) * sp.diff(G, b))
    u1, u2 = (a + ab_) / 2, (b + bb) / 2
    V1 = sp.expand(u1 * (u1 ** 2 + u2 ** 2))                      # the cubic term, without its epsilon
    V2 = sp.expand((u1 ** 2 + u2 ** 2) ** 2 / 4)                  # the quartic term, without its epsilon^2
    wt = lambda mono: sp.Poly(mono, a, ab_, b, bb).monoms()[0]    # a, b turn as e^{-2i tau}: the frequency of a monomial is 2 (holomorphic degree - antiholomorphic degree)

    def harmonic(expr, k):
        return sum(term for term in sp.Add.make_args(sp.expand(expr)) if (lambda mm: mm[0] + mm[2] - mm[1] - mm[3])(wt(term)) == k)
    avg2 = harmonic(V2, 0)
    K2 = avg2
    for k in (1, 3):                                              # frequencies 2k: K_2 = <V_2> - i sum_{omega > 0} {V^(omega), V^(-omega)} / omega
        K2 += -sp.I * br(harmonic(V1, k), harmonic(V1, -k)) / (2 * k)
    K2 = sp.expand(K2)
    na, nb = a * ab_, b * bb
    target = sp.expand((-12 * na ** 2 - sp.Rational(32, 3) * na * nb + sp.Rational(4, 3) * nb ** 2) / 32)
    e3a = sp.expand(K2 - target) == 0
    Na, Nb = sp.symbols("n_a n_b", positive=True)
    K2n = (-12 * Na ** 2 - sp.Rational(32, 3) * Na * Nb + sp.Rational(4, 3) * Nb ** 2) / 32
    e3b = sp.simplify(2 * (sp.diff(K2n, Nb) - sp.diff(K2n, Na)) - sp.Rational(5, 6) * (Na + Nb)) == 0
    # the 1-dof check of the formula for K_2: H = p^2 + x^2 + eps x^3 has K_2 = -(15/32) n^2
    c1 = sp.expand(u1 ** 3)
    k1 = sp.expand(sum(-sp.I * br(harmonic(c1, k), harmonic(c1, -k)) / (2 * k) for k in (1, 3)))
    e3c = sp.expand(k1 + sp.Rational(15, 32) * na ** 2) == 0
    ok("E3 the normal form: K_2 = (1/32)[-12 n_a^2 - (32/3) n_a n_b + (4/3) n_b^2], no phase-coupling term; omega_b - omega_a = (5/6) eps^2 N; the cubic oscillator's -15/32 as a check of the method (EXACT)",
       e3a and e3b and e3c)
    es = sp.symbols("e", positive=True)
    ph = sp.symbols("phi", real=True)
    yv = es / 2 * sp.sin(ph)
    ser = sp.series(1 / ((sp.Rational(3, 4) - yv) * sp.sqrt(sp.Rational(1, 4) + yv)), es, 0, 3).removeO()
    Tn = es / 2 * sp.integrate(sp.expand(ser), (ph, -sp.pi / 2, sp.pi / 2))
    Wser = sp.series(Tn / (4 * sp.pi * es / 3) * sp.sqrt(1 - 4 * es ** 2 / 9), es, 0, 3).removeO()
    ok("E4 on the axial line W = 1 + (5/12) eps^2 + O(eps^4), from the two periods (EXACT)", sp.simplify(Wser - 1 - sp.Rational(5, 12) * es ** 2) == 0)

# ================================================================ P0: the reduction, in space
if "P0" in SECT:
    print("P0 the reduction, in space")
    worst = 0.0
    worstH = 0.0
    r0 = np.random.default_rng(22)
    for tv in (0.1, 0.4, 0.8, 1.5):
        mod = Model(tv)
        for kind in "BE":
            u = field(tv, kind)
            for X in r0.normal(0, 1.5, (60, 3)):
                h = 1e-6
                d = u(X)
                al, be = labels(tv, np.stack([X + h * d, X - h * d, X]))
                s, P = mod.sP(al, be, kind)
                v = np.array([(s[0] - s[1]).real, (s[0] - s[1]).imag, (P[0] - P[1]).real, (P[0] - P[1]).imag])
                gV = 4 * (s[2] ** 2 - mod.sig ** 2) * np.conj(s[2]) / tv ** 2
                w = np.array([2 * P[2].real, 2 * P[2].imag, -gV.real, -gV.imag])
                cs = abs(v @ w) / (np.linalg.norm(v) * np.linalg.norm(w))
                worst = max(worst, math.sqrt(max(0.0, 1 - cs * cs)))
                worstH = max(worstH, abs(mod.H(*mod.zp(al[2], be[2], kind)) - 1))
    OUT["P0"] = {"worst angle": worst, "worst |H - 1|": worstH}
    ok("P0 PREDICTION a step along the field moves (s, P) parallel to (2P, -grad V): 60 points, four times, both types (1e-6)", worst < 1e-6, "worst angle %.1e; |H - 1| %.1e" % (worst, worstH))

# ================================================================ P1: the invariant
if "P1" in SECT:
    print("P1 the invariant")
    P1 = {}
    worst = 0.0
    for tv in (0.1, 0.2, 0.4, 0.8, 1.5):
        mod = Model(tv)
        for kind in "BE":
            for k, X0 in enumerate(STARTS):
                so = run_line(tv, kind, X0, 20 * 2 * math.pi, 2400)
                I = mod.I2_space(so.y.T, kind)
                dr = float((I.max() - I.min()) / mod.I2_range)
                Hd = float(np.abs(mod.H(*mod.zp(*labels(tv, so.y.T), kind)) - 1).max())
                P1["t=%.1f %s line %d" % (tv, kind, k)] = {"drift / range": dr, "place in the range": float((I.mean() - mod.I2_cross) / mod.I2_range),
                                                         "status": int(so.status), "max |X|": float(np.linalg.norm(so.y, axis=0).max()), "|H - 1|": Hd}
                worst = max(worst, dr if so.status == 0 else 1.0)
        print("    t = %.1f done, worst so far %.1e (%.0f s)" % (tv, worst, time.time() - T0), flush=True)
    OUT["P1"] = {"lines": P1, "worst drift / range": worst}
    ok("P1 PREDICTION I_2 is constant along lines of force integrated in space: five lines of each type at t = 0.1, 0.2, 0.4, 0.8, 1.5, twenty circuits (1e-8 of its range)", worst < 1e-8, "worst %.1e" % worst)

# ================================================================ P2: the two closed lines
if "P2" in SECT:
    print("P2 the two closed lines")
    P2 = {}
    w_tan = w_close = w_hyp = 0.0
    links, inside = [], True
    edge = 1.0

    def gauss(C1, C2):
        """the linking integral of two closed polygons, by the midpoint rule"""
        d1, d2 = np.roll(C1, -1, 0) - C1, np.roll(C2, -1, 0) - C2
        R = (C1 + d1 / 2)[:, None, :] - (C2 + d2 / 2)[None, :, :]
        tri = np.einsum("ijk,ijk->ij", R, np.cross(d1[:, None, :], d2[None, :, :]))
        return float(np.sum(tri / np.linalg.norm(R, axis=2) ** 3) / (4 * math.pi))

    def first_closure(so, X0):
        """the upward crossing of the plane through the start, within the first circuit and a half, that comes back nearest to it"""
        cand = [(l, Y) for l, Y in zip(so.t_events[0], so.y_events[0]) if 1 < l < 9.5]
        return min(cand, key=lambda v: np.linalg.norm(v[1] - X0))

    for tv in (0.1, 0.4, 0.8, 1.5):
        mod = Model(tv)
        for kind in "BE":
            u = field(tv, kind)
            # the axial line: the great circle (-i sh^2 cos phi, conj(sh) sin phi), with i conj(sh) sin phi for the electric type
            circ = lambda ph: chart(tv, -1j * mod.sh ** 2 * np.cos(ph) + 0j, (np.conj(mod.sh) if kind == "B" else 1j * np.conj(mod.sh)) * np.sin(ph) + 0j)
            ph = np.linspace(0, 2 * math.pi, 400, endpoint=False) + 0.1
            Tg = unit(circ(ph + 1e-6) - circ(ph - 1e-6))
            tan = float(np.sqrt(np.maximum(0, 1 - np.sum(Tg * u(circ(ph)), -1) ** 2)).max())
            Iax = float(np.abs(mod.I2_space(circ(ph), kind) - mod.I2_axial).max() / mod.I2_range)
            Xa = circ(np.array(0.3))
            n0 = u(Xa)
            ev = lambda l, Y: (Y - Xa) @ n0
            ev.direction = 1
            so = solve_ivp(rhs_line(tv, kind), (0, 12), Xa, method="DOP853", rtol=1e-12, atol=1e-13, events=ev)
            ret = [first_closure(so, Xa)]
            close_ax = float(np.linalg.norm(ret[0][1] - Xa))
            # the cross line: from x_0 = 2 sqrt(m c_0) on the axis, the momentum across it
            x0 = 2 * math.sqrt(mod.m * mod.c0)
            py = math.sqrt(1 - (x0 * x0 - mod.m) ** 2 / tv ** 2)
            Xc = chart(tv, *[np.array(v) for v in mod.labels_of(x0 + 0j, 1j * py, kind)])
            n0c = u(Xc)
            evc = lambda l, Y: (Y - Xc) @ n0c
            evc.direction = 1
            soc = solve_ivp(rhs_line(tv, kind), (0, 12), Xc, method="DOP853", rtol=1e-12, atol=1e-13, events=evc, dense_output=True)
            retc = [first_closure(soc, Xc)]
            close_cr = float(np.linalg.norm(retc[0][1] - Xc))
            Cc = soc.sol(np.linspace(0, retc[0][0], 900, endpoint=False)).T
            zc_, pc_ = mod.zp(*labels(tv, Cc), kind)
            cosn = mod.elliptic(zc_)[1]
            hyp = float(np.abs(cosn ** 2 - mod.c0).max())
            Icr = float(np.abs(mod.I2(zc_, pc_) - mod.I2_cross).max() / mod.I2_range)
            lk = gauss(circ(np.linspace(0, 2 * math.pi, 900, endpoint=False)), Cc)
            # 20,000 random points of the label sphere
            rr = np.random.default_rng(int(1000 * tv) + (kind == "E"))
            v = rr.normal(size=(20000, 4))
            v /= np.linalg.norm(v, axis=1)[:, None]
            Ir = mod.I2(*mod.zp(v[:, 0] + 1j * v[:, 1], v[:, 2] + 1j * v[:, 3], kind))
            lo, hi = float((Ir.min() - mod.I2_cross) / mod.I2_range), float((Ir.max() - mod.I2_cross) / mod.I2_range)
            inside = inside and lo > 0 and hi < 1
            edge = min(edge, lo, 1 - hi)
            P2["t=%.1f %s" % (tv, kind)] = {"axial: tangent off the field": tan, "axial: I_2 off its value / range": Iax, "axial: closure": close_ax, "axial: label length": float(ret[0][0]),
                                          "cross: closure": close_cr, "cross: cos^2 nu off c_0": hyp, "cross: I_2 off its value / range": Icr, "cross: label length": float(retc[0][0]),
                                          "linking number (Gauss)": lk, "random points, lowest and highest place in the range": [lo, hi]}
            w_tan, w_close, w_hyp = max(w_tan, tan), max(w_close, close_ax, close_cr), max(w_hyp, hyp, Iax, Icr)
            links.append(lk)
    OUT["P2"] = {"cases": P2, "worst tangent": w_tan, "worst closure": w_close, "worst off the hyperbola or off the extreme value": w_hyp, "linking numbers": links, "closest a random point came to an end of the range": edge}
    ok("P2 PREDICTION the axial great circle and the cross line are closed lines of force (tangent 1e-7, closure 1e-8), linked once; I_2 on 20,000 random points lies strictly between its values on them",
       w_tan < 1e-7 and w_close < 1e-8 and all(abs(abs(v) - 1) < 0.02 for v in links) and inside,
       "tangent %.1e; closure %.1e; off the hyperbola %.1e; linking %s; nearest to an end %.1e of the range" % (w_tan, w_close, w_hyp, [round(v, 4) for v in links], edge))

# ================================================================ P3: the periods
if "P3" in SECT:
    print("P3 the periods")
    P3 = {}
    worst = 0.0
    for tv in (0.1, 0.4, 0.8):
        mod = Model(tv)
        for k, X0 in enumerate(STARTS[:3]):
            Kv = -mod.I2_space(X0, "B") / (4 * mod.m)
            Tm, Tn, (X1, X2, X3) = periods(Kv, mod.eps)
            cmid = math.sqrt(float((X1 + X2) / 2))

            def emu(l, Y):
                return float(mod.zp(*labels(tv, Y[:3]), "B")[0].imag)

            def enu(l, Y):
                return float(mod.elliptic(mod.zp(*labels(tv, Y[:3]), "B")[0])[1] - cmid)
            emu.direction = enu.direction = 1
            so = solve_ivp(rhs_line(tv, "B", mod, tau=True), (0, 16 * 2 * math.pi), np.append(X0, 0.0), method="DOP853", rtol=1e-12, atol=1e-13, events=(emu, enu))
            dm, dn = np.diff(so.y_events[0][:, 3]), np.diff(so.y_events[1][:, 3])
            tm, tn = float(Tm) / mod.f, float(Tn) / mod.f
            row = {"K": float(Kv), "place in the range": float((Kv - mod.Kmin) / (mod.Kmax - mod.Kmin)), "roots": [float(X1), float(X2), float(X3)],
                   "T_mu / f (integral)": tm, "measured mean": float(dm.mean()), "measured spread": float(dm.max() - dm.min()), "periods measured": int(dm.size),
                   "T_nu / f (integral)": tn, "measured mean ": float(dn.mean()), "measured spread ": float(dn.max() - dn.min()), "W = T_nu / T_mu": float(Tn / Tm), "W measured": float(dn.mean() / dm.mean())}
            P3["t=%.1f line %d" % (tv, k)] = row
            worst = max(worst, abs(dm.mean() / tm - 1), abs(dn.mean() / tn - 1), (dm.max() - dm.min()) / tm, (dn.max() - dn.min()) / tn)
            print("    t = %.1f line %d: W = %.9f, measured %.9f; periods constant to %.1e, %.1e" % (tv, k, row["W = T_nu / T_mu"], row["W measured"], row["measured spread"] / tm, row["measured spread "] / tn), flush=True)
    OUT["P3"] = {"lines": P3, "worst relative difference": worst}
    ok("P3 PREDICTION the separated time between crossings of y = 0 is constant along a line and equals T_mu/f; between passages of nu through its middle, T_nu/f (1e-7; three lines at t = 0.1, 0.4, 0.8)", worst < 1e-7, "worst %.1e" % worst)

# ================================================================ P4: the leading law
if "P4" in SECT:
    print("P4 the leading law")
    P4 = {}
    mod = Model(0.05)
    rat = []
    for fr in (0.05, 0.25, 0.5, 0.75, 0.95):
        Kv = mp.mpf(mod.Kmin) + fr * (mp.mpf(mod.Kmax) - mp.mpf(mod.Kmin))
        Tm, Tn, _ = periods(Kv, mod.eps, 40)
        rat.append(float((Tn / Tm - 1) / mp.mpf(mod.eps) ** 2))
    wax = float(W_axial(Model(0.4).eps))
    table = {}
    for tv in (0.1, 0.2, 0.4, 0.8, 1.5, 3.0):
        md = Model(tv)
        row = {"eps": md.eps, "leading law (5/12) eps^2": 5 / 12 * md.eps ** 2, "axial W - 1": float(W_axial(md.eps) - 1)}
        for fr in (0.25, 0.5, 0.75, 0.999):
            Tm, Tn, _ = periods(mp.mpf(md.Kmin) + fr * (mp.mpf(md.Kmax) - mp.mpf(md.Kmin)), md.eps, 30)
            row["W - 1 at %.3f of the range" % fr] = float(Tn / Tm - 1)
        table["t=%.1f" % tv] = row
    OUT["P4"] = {"(W - 1)/eps^2 at t = 0.05, at 0.05, 0.25, 0.5, 0.75, 0.95 of the range of K": rat, "5/12": 5 / 12, "axial W at t = 0.4": wax, "table": table}
    ok("P4 PREDICTION (W - 1)/eps^2 is within 2 % of 5/12 on five tori at t = 0.05, and W on the axial line at t = 0.4 is 1.0156 +- 0.0002",
       all(abs(v / (5 / 12) - 1) < 0.02 for v in rat) and abs(wax - 1.0156) < 0.0002, "ratios %s; axial W(0.4) = %.6f" % ([round(v, 5) for v in rat], wax))

# ================================================================ P5: MP10's numbers re-read
if "P5" in SECT:
    print("P5 MP10's numbers re-read")

    def frame(u, Xg):
        n0 = u(Xg)
        e1 = unit(np.cross(n0, [0.3, 0.5, 0.8]))
        return n0, e1, np.cross(n0, e1)

    def returns10(u, Xc, n0, e1, e2, y, Lmax, smin=1.0):
        """MP10's return: of the upward crossings of the plane through Xc within the arc Lmax, the one closest to Xc; also which crossing it was and its arc length"""
        ev = lambda s, X: (X - Xc) @ n0
        ev.direction = 1
        so = solve_ivp(lambda s, X: u(X), (0, Lmax), Xc + y[0] * e1 + y[1] * e2, method="DOP853", rtol=1e-12, atol=1e-14, events=ev)
        cand = [(s, Xr) for s, Xr in zip(so.t_events[0], so.y_events[0]) if s > smin]
        if not cand:
            raise RuntimeError("no return within %.0f" % Lmax)
        k = min(range(len(cand)), key=lambda i: np.linalg.norm(cand[i][1] - Xc))
        s, Xr = cand[k]
        return np.array([(Xr - Xc) @ e1, (Xr - Xc) @ e2]), k + 1, float(s), len(cand)

    P5 = {}
    # (a) the defect over t^2 as t -> 0
    al0, be0 = labels(0.0, X0H)
    na, nb = al0.imag ** 2 + be0.real ** 2, al0.real ** 2 + be0.imag ** 2
    lim = 5 * math.pi / 24 * (1 + X0H @ X0H) / 2 * math.sqrt(na * nb)
    dq = {}
    for tv in (0.0125, 0.00625, 0.003125):
        u = field(tv, "B")
        n0, e1, e2 = frame(u, X0H)
        yr, k, s, nc = returns10(u, X0H, n0, e1, e2, np.zeros(2), 40.0)
        dq[tv] = float(np.linalg.norm(yr)) / tv ** 2
    ext = (8 * dq[0.003125] - 6 * dq[0.00625] + dq[0.0125]) / 3
    mp10 = {0.0125: 3.06751280175973e-05 / 0.0125 ** 2, 0.025: 0.00011918828836359899 / 0.025 ** 2, 0.05: 0.00044713515623104674 / 0.05 ** 2}
    P5["a"] = {"shares |a|^2, |b|^2": [float(na), float(nb)], "the limit, from the labels": float(lim), "defect / t^2": {str(k_): v for k_, v in dq.items()}, "extrapolated": float(ext),
               "MP10's recorded values extrapolated": (8 * mp10[0.0125] - 6 * mp10[0.025] + mp10[0.05]) / 3}
    oka = abs(ext - 0.2017) < 0.0010
    print("    (a) defect/t^2 = %s -> %.6f; the derivation's limit %.6f" % ([round(v, 6) for v in dq.values()], ext, lim), flush=True)
    # (b) the true winding at t = 0.4
    tv = 0.4
    mod = Model(tv)
    Kv = -mod.I2_space(X0H, "B") / (4 * mod.m)
    Tm, Tn, (X1, X2, X3) = periods(Kv, mod.eps)
    W = float(Tn / Tm)
    cmid = math.sqrt(float((X1 + X2) / 2))

    def enu(l, Y):
        return float(mod.elliptic(mod.zp(*labels(tv, Y), "B")[0])[1] - cmid)
    enu.direction = 1
    NRET = 200
    pts, l0, Xc = [], 0.0, X0H
    while len(pts) < NRET + 1:
        so = solve_ivp(rhs_line(tv, "B"), (l0, l0 + 40 * 2 * math.pi), Xc, method="DOP853", rtol=1e-12, atol=1e-13, events=enu)
        pts += list(so.y_events[0])
        l0, Xc = so.t[-1], so.y[:, -1]
    pts = np.array(pts[:NRET + 1])
    z, p = mod.zp(*labels(tv, pts), "B")
    ch, cn = mod.elliptic(z)
    mu = np.sign(z.imag) * np.arccosh(np.maximum(ch, 1.0))
    sn = np.sqrt(1 - cn ** 2)
    pmu = p.real * mod.f * np.sinh(mu) * cn + p.imag * mod.f * np.cosh(mu) * sn
    PiX = lambda X: Kv + X - mod.g * X * (X - 0.75) ** 2
    ang = np.unwrap(np.arctan2(pmu / (mod.f * math.sqrt(PiX(1.0))), mu / math.acosh(math.sqrt(float(X3)))))
    turns = abs(ang[-1] - ang[0]) / (2 * math.pi)
    Idr = float(np.ptp(mod.I2(z, p)) / mod.I2_range)
    P5["b"] = {"K": float(Kv), "place in the range": float((Kv - mod.Kmin) / (mod.Kmax - mod.Kmin)), "W": W, "W - 1": W - 1, "one turn in this many circuits": 1 / (W - 1),
               "turns counted over 200 returns": float(turns), "200 (W - 1)": 200 * (W - 1), "drift of I_2 over the 200 returns / range": Idr}
    okb = abs(turns - 200 * (W - 1)) < 0.05 and 0.0125 < W - 1 < 0.0175
    print("    (b) W - 1 = %.6f (one turn in %.1f circuits); turns in 200 returns %.4f against %.4f" % (W - 1, 1 / (W - 1), turns, 200 * (W - 1)), flush=True)
    # (c) MP10's sixty closest-returns
    u = field(tv, "B")
    n0, e1, e2 = frame(u, X0H)
    ys, ks, ss, ncs = [np.zeros(2)], [], [], []
    for _ in range(60):
        yr, k, s, nc = returns10(u, X0H, n0, e1, e2, ys[-1], 60.0)
        ys.append(yr)
        ks.append(k)
        ss.append(s)
        ncs.append(nc)
    ys = np.array(ys)
    c = ys - ys.mean(0)
    tm10 = float(np.diff(np.unwrap(np.arctan2(c[:, 1], c[:, 0]))).sum() / (2 * math.pi))
    circuits = int(sum(ks))
    P5["c"] = {"circuits taken by the sixty returns": circuits, "which crossing was taken, counts": {str(v): int(ks.count(v)) for v in sorted(set(ks))}, "arc length in space": float(sum(ss)),
               "arc length of the first return": float(ss[0]), "turns about the mean, as MP10 counted": tm10, "MP10's recorded value": 2.037620105774404,
               "circuits times (W - 1)": circuits * (W - 1), "upward crossings within the arc of 60, least and most": [int(min(ncs)), int(max(ncs))]}
    okc = 115 <= circuits <= 165
    print("    (c) the sixty returns took %d circuits (%s); turns as MP10 counted %.4f; circuits x (W - 1) = %.3f" % (circuits, P5["c"]["which crossing was taken, counts"], tm10, circuits * (W - 1)), flush=True)
    OUT["P5"] = P5
    ok("P5a PREDICTION defect/t^2 for MP10's line extrapolates to 0.2017 +- 0.0010", oka, "%.6f" % ext)
    ok("P5b PREDICTION at t = 0.4 the returns turn 200 (W - 1) in 200 circuits (+- 0.05), with W - 1 between 0.0125 and 0.0175", okb, "W - 1 = %.6f; turns %.4f against %.4f" % (W - 1, turns, 200 * (W - 1)))
    ok("P5c PREDICTION MP10's sixty closest-returns took between 115 and 165 circuits", okc, "%d circuits" % circuits)

# ================================================================ P6: the lost line
if "P6" in SECT:
    print("P6 the lost line at t = 0.8")
    tv = 0.8
    mod = Model(tv)
    so = run_line(tv, "B", X0H, 100 * 2 * math.pi, 60000)
    I = mod.I2_space(so.y.T, "B")
    dr = float(np.ptp(I) / mod.I2_range)
    dist = np.linalg.norm(so.y.T - X0H, axis=1)
    far = float(np.linalg.norm(so.y, axis=0).max())
    back = float(dist[so.t > 3].min())
    Kv = -I[0] / (4 * mod.m)
    Tm, Tn, roots = periods(Kv, mod.eps)
    # how near the torus comes to the label of infinity: alpha = 1 is s = 1, zeta = conj(sigma_hat)
    zinf = np.conj(mod.sh)
    chi, cni = mod.elliptic(zinf)
    OUT["P6"] = {"drift of I_2 / range": dr, "status": int(so.status), "label length": float(so.t[-1]), "farthest from the origin": far, "closest return to the start": back,
                 "K": float(Kv), "place in the range": float((Kv - mod.Kmin) / (mod.Kmax - mod.Kmin)), "W": float(Tn / Tm), "roots X1, X2, X3": [float(v) for v in roots],
                 "the label of infinity: cos^2 nu, cosh^2 mu": [float(cni ** 2), float(chi ** 2)], "K of the torus through infinity, place in the range": float((-mod.I2(zinf, 0j) / (4 * mod.m) - mod.Kmin) / (mod.Kmax - mod.Kmin))}
    ok("P6 PREDICTION at t = 0.8 the line through (0.6, -0.3, 0.5) keeps I_2 (1e-8) over 100 circuits and comes back within 0.3 of its start", dr < 1e-8 and back < 0.3 and so.status == 0,
       "drift %.1e; closest return %.4f; farthest %.1f" % (dr, back, far))
    print("  %s P6 FORECAST it reaches more than 10 from the origin   [%.2f]" % ("held  " if far > 10 else "missed", far))
    OUT["P6"]["forecast held"] = bool(far > 10)

# ================================================================ P7: a torus of closed lines (Amendment 1)
if "P7" in SECT:
    print("P7 a torus of closed lines")
    from scipy.optimize import brentq
    tv = 0.4
    mod = Model(tv)
    Wof = lambda K: float((lambda r: r[1] / r[0])(periods(K, mod.eps)))
    span = mod.Kmax - mod.Kmin
    K65 = brentq(lambda K: Wof(K) - 65 / 64, mod.Kmin + 1e-6 * span, mod.Kmax - 1e-6 * span, xtol=1e-14)
    K66 = brentq(lambda K: Wof(K) - 66 / 65, mod.Kmin + 1e-6 * span, mod.Kmax - 1e-6 * span, xtol=1e-14)

    def start_on(K):
        """a point of the torus K on the axis y = 0, at the middle of the swing of cos^2 nu, as a point of space (magnetic type)"""
        _, _, (X1, X2, X3) = periods(K, mod.eps)
        c = float((X1 + X2) / 2)
        Pi = lambda X: K + X - mod.g * X * (X - 0.75) ** 2
        sn = math.sqrt(1 - c)
        z = mod.f * math.sqrt(c) + 0j
        p = -math.sqrt(-Pi(c)) / sn + 1j * math.sqrt(Pi(1.0)) / sn
        assert abs(mod.H(z, p) - 1) < 1e-12 and abs(-mod.I2(z, p) / (4 * mod.m) - K) < 1e-11
        return chart(tv, *[np.array(v) for v in mod.labels_of(z, p, "B")])

    def closest_return(K, ncirc=70):
        X0 = start_on(K)
        u = field(tv, "B")
        n0 = u(X0)
        ev = lambda l, Y: (Y - X0) @ n0
        ev.direction = 1
        emu = lambda l, Y: float(mod.zp(*labels(tv, Y), "B")[0].imag)
        emu.direction = 1
        so = solve_ivp(rhs_line(tv, "B"), (0, ncirc * 2 * math.pi), X0, method="DOP853", rtol=1e-12, atol=1e-13, events=(ev, emu))
        cand = [(l, float(np.linalg.norm(Y - X0))) for l, Y in zip(so.t_events[0], so.y_events[0]) if l > 1]
        k = min(range(len(cand)), key=lambda i: cand[i][1])
        return {"closest return": cand[k][1], "it is plane crossing number": k + 1, "label length": cand[k][0], "crossings of y = 0 before it": int(np.sum(so.t_events[1] < cand[k][0] + 1e-6)),
                "second closest": sorted(v[1] for v in cand)[1], "W": Wof(K), "place in the range": (K - mod.Kmin) / span}
    r65 = closest_return(K65)
    rmid = closest_return((K65 + K66) / 2)
    OUT["P7"] = {"K with W = 65/64": K65, "K with W = 66/65": K66, "on the torus W = 65/64": r65, "halfway to the torus W = 66/65": rmid}
    ok("P7 PREDICTION (Amendment 1) at t = 0.4 a line on the torus W = 65/64 closes after 64 circuits (1e-7), with 65 crossings of y = 0; a line halfway to W = 66/65 comes no nearer than 1e-3 in 70 circuits",
       r65["closest return"] < 1e-7 and r65["it is plane crossing number"] == 64 and r65["crossings of y = 0 before it"] == 65 and rmid["closest return"] > 1e-3,
       "closure %.1e at plane crossing %d, %d crossings of y = 0; the neighbour's nearest return %.4f" % (r65["closest return"], r65["it is plane crossing number"], r65["crossings of y = 0 before it"], rmid["closest return"]))

if SECT & {"E", "P0", "P1", "P2", "P3", "P4", "P5", "P6", "P7"}:
    OUT.setdefault("checks", {})["run of sections %s" % " ".join(sorted(SECT))] = {"passed": NPASS, "failed": NFAIL, "seconds": time.time() - T0}
    OUTF.write_text(json.dumps(OUT, indent=1))
    print("%d passed, %d failed (%.0f s)" % (NPASS, NFAIL, time.time() - T0))
