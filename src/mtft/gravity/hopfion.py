r"""mtft.gravity.hopfion — the squared hopfion's lines of force after t = 0, as an integrable system (MODEL; v1.1.0, 2026-11-16).

The lines of force of the squared hopfion W = n ⊗ n/d⁵ are the lines of Im G (magnetic type, "B") and of Re G
(electric type, "E"), G = √d · F₁₁: the light's own lines turned by the local duality angle arg(d)/2, which is
Swearngin et al.'s result (arXiv:1302.1431).  At t = 0 they are the hopfion's circles.  The lab's record MP22
(14 November 2026; derived there, tested in space, unreviewed) reduces them at any later time to a particle in a plane:

    s = (1 − it(1 − α))^{1/2},   σ² = 1 − it,   P = β̄ (magnetic) or iβ̄ (electric),
    H = |P|² + |s² − σ²|²/t²   on the level H = 1.

Turned so that the wells lie on the real axis (ζ = x + iy = s σ̂*, p = P σ̂*, m = |σ|² = √(1 + t²)) the potential
is [r⁴ − 2m(x² − y²) + m²]/t², which is Garnier's integrable quartic: it separates in elliptic coordinates with
foci at ±2√m.  The second invariant, with L = x p_y − y p_x:

    I₂ = L² + 4m p_x² + (4m/t²) x² (r² − 2m).

So at every time the label sphere is filled by tori nested between two closed lines, the axial line and the cross
line, and on the torus K = −I₂/(4m) a line winds W(K) = T_ν/T_μ times, two periods of the curve
Y² = X(X − 1)Π(X), Π(X) = K + X − (4/ε²)X(X − ¾)², with ε = t/(2√(1 + t²)) the only parameter.  To leading order
W − 1 = (5/12)ε² = 5t²/(48(1 + t²)) on every torus.  The lines are never chaotic.

Recorded misses of MP22, kept: its forecast that the line through (0.6, −0.3, 0.5) at t = 0.8 runs more than 10
from the origin (it reaches 4.7); a count of crossings that included the starting point; and the older records'
"one turn in about 29 circuits" at t = 0.4, which is one turn in 63.9.
"""
from __future__ import annotations

import math

import mpmath as mp
import numpy as np

from . import fields

MP10_POINT = (0.6, -0.3, 0.5)                        # the point the lab's records MP10, MP11 and MP22 followed
TABLE_TIMES = (0.1, 0.2, 0.4, 0.8, 1.5, 3.0)
TABLE_FRACTIONS = (0.25, 0.5, 0.75, 0.999)


class Hopfion:
    """The reduced system at the time t > 0."""

    def __init__(self, t: float):
        if not t > 0:
            raise ValueError("the reduction is for t > 0; at t = 0 the lines of force are the hopfion's circles")
        self.t = float(t)
        self.sigma = complex(np.sqrt(1 - 1j * self.t))
        self.sigma_hat = self.sigma / abs(self.sigma)
        self.m = abs(self.sigma) ** 2
        self.focus = 2 * math.sqrt(self.m)
        self.eps = self.t / (2 * math.sqrt(1 + self.t ** 2))
        self.g = 4 / self.eps ** 2
        self.K_axial = 1 / self.t ** 2
        self.c0 = 0.5 - 0.25 * math.sqrt(1 + 4 * self.eps ** 2 / 3)
        self.K_cross = self.g * self.c0 * (self.c0 - 0.75) ** 2 - self.c0
        self.I2_axial, self.I2_cross = -4 * self.m * self.K_axial, -4 * self.m * self.K_cross
        self.I2_range = self.I2_axial - self.I2_cross

    # ---- labels and the plane
    def reduce(self, alpha, beta, kind: str = "B"):
        """(ζ, p) of a label point: ζ = s σ̂*, p = P σ̂*."""
        s = np.sqrt(1 - 1j * self.t * (1 - np.asarray(alpha, complex)))
        P = np.conj(beta) if kind == "B" else 1j * np.conj(beta)
        return s * np.conj(self.sigma_hat), P * np.conj(self.sigma_hat)

    def lift(self, zeta, p, kind: str = "B"):
        """The label point (α, β) of a phase point (ζ, p)."""
        s, P = zeta * self.sigma_hat, p * self.sigma_hat
        return (s * s - self.sigma ** 2) / (1j * self.t), (np.conj(P) if kind == "B" else np.conj(-1j * P))

    def potential(self, zeta):
        return np.abs(np.asarray(zeta) ** 2 - self.m) ** 2 / self.t ** 2

    def H(self, zeta, p):
        """|p|² + |ζ² − m|²/t²; it is |β|² + |α|², so 1 on the label sphere."""
        return np.abs(p) ** 2 + self.potential(zeta)

    def I2(self, zeta, p):
        """The second invariant, L² + 4m p_x² + (4m/t²) x² (r² − 2m)."""
        zeta, p = np.asarray(zeta, complex), np.asarray(p, complex)
        x, y, px, py = zeta.real, zeta.imag, p.real, p.imag
        return (x * py - y * px) ** 2 + 4 * self.m * px ** 2 + (4 * self.m / self.t ** 2) * x * x * (x * x + y * y - 2 * self.m)

    def I2_at(self, X, kind: str = "B"):
        """I₂ at points of space, through their labels."""
        return self.I2(*self.reduce(*fields.labels(self.t, X), kind))

    def place(self, I2):
        """Where a value of I₂ lies between the axial line (0) and the cross line (1): the fraction of ``torus``."""
        return (self.I2_axial - np.asarray(I2)) / self.I2_range

    def elliptic(self, zeta):
        """cosh μ and cos ν of ζ = 2√m cosh(μ + iν)."""
        a, b = np.abs(zeta + self.focus), np.abs(zeta - self.focus)
        return (a + b) / (2 * self.focus), (a - b) / (2 * self.focus)

    # ---- the tori
    def torus(self, fraction: float) -> float:
        """K on the torus at ``fraction`` of the way from the axial line (0) to the cross line (1)."""
        return self.K_axial + fraction * (self.K_cross - self.K_axial)

    def periods(self, K, dps: int = 25):
        """T_μ, T_ν and the three roots of Π on the torus K (mpmath numbers)."""
        return periods(K, self.eps, dps)

    def winding(self, K=None, fraction=None, dps: int = 25) -> float:
        """W = T_ν/T_μ on a torus given by K or by its fraction of the range; on the axial line with neither."""
        if K is None and fraction is None:
            return float(winding_axial(self.eps, dps))
        with mp.workdps(dps):
            if K is None:
                K = mp.mpf(self.K_axial) + mp.mpf(fraction) * (mp.mpf(self.K_cross) - mp.mpf(self.K_axial))
            Tm, Tn, _ = periods(K, self.eps, dps)
            return float(Tn / Tm)

    def leading_law(self) -> float:
        """(5/12)ε² = 5t²/(48(1 + t²)): W − 1 to leading order, the same on every torus."""
        return 5 * self.t ** 2 / (48 * (1 + self.t ** 2))

    # ---- the closed lines
    def axial_line(self, phi, kind: str = "B"):
        """The labels of the axial closed line: the great circle (−iσ̂² cos φ, σ̂* sin φ); iσ̂* sin φ for the electric type."""
        phi = np.asarray(phi, float)
        return -1j * self.sigma_hat ** 2 * np.cos(phi) + 0j, (np.conj(self.sigma_hat) if kind == "B" else 1j * np.conj(self.sigma_hat)) * np.sin(phi) + 0j

    def axial_tangent(self, phi, kind: str = "B"):
        """d/dφ of ``axial_line``."""
        phi = np.asarray(phi, float)
        return 1j * self.sigma_hat ** 2 * np.sin(phi) + 0j, (np.conj(self.sigma_hat) if kind == "B" else 1j * np.conj(self.sigma_hat)) * np.cos(phi) + 0j

    def cross_point(self, kind: str = "B"):
        """The labels of the point where the cross line meets the axis of the plane: ζ = 2√(m c₀), the momentum across it."""
        x0 = 2 * math.sqrt(self.m * self.c0)
        py = math.sqrt(1 - (x0 * x0 - self.m) ** 2 / self.t ** 2)
        return self.lift(x0 + 0j, 1j * py, kind)

    # ---- the motion in the plane
    def flow(self, y):
        """Hamilton's equations for (x, y, p_x, p_y), with the separated time τ′ as a fifth variable: dτ′ = 2dτ/|ζ² − 4m|."""
        x, yy, px, py = y[0], y[1], y[2], y[3]
        r2, c = x * x + yy * yy, 4 / self.t ** 2
        return np.array([2 * px, 2 * py, -c * (r2 - self.m) * x, -c * (r2 + self.m) * yy, 2 / abs(complex(x, yy) ** 2 - 4 * self.m)])

    def orbit(self, zeta, p, time: float, h: float | None = None, **kw):
        """The orbit of a phase point over ``time`` (see ``fields.integrate``); the state is (x, y, p_x, p_y, τ′)."""
        h = h or self.t / (40 * math.sqrt(self.m))
        return fields.integrate(self.flow, [zeta.real, zeta.imag, p.real, p.imag, 0.0], time, h, **kw)

    def measured_periods(self, zeta, p, count: int = 5) -> dict:
        """The two periods of an orbit in separated time, measured along it: the second route to the winding.

        Between successive upward crossings of y = 0 the separated time is T_μ/(2√m); between successive upward passages
        of cos ν through the middle of its range it is T_ν/(2√m).  Returns both, their spreads, and the K of the orbit."""
        zeta, p = complex(zeta), complex(p)
        K = float(-self.I2(zeta, p) / (4 * self.m))
        _, _, (X1, X2, _) = periods(K, self.eps, 20)
        mid = (math.sqrt(float(X1)) + math.sqrt(float(X2))) / 2
        span = (count + 2) * 2 * math.pi * self.t / (4 * math.sqrt(self.m)) * 1.3
        mu = self.orbit(zeta, p, span, event=lambda y: y[1], max_events=count + 1)["events"]
        nu_ = self.orbit(zeta, p, span, event=lambda y: self.elliptic(complex(y[0], y[1]))[1] - mid, max_events=count + 1)["events"]
        if len(mu) < 3 or len(nu_) < 3:
            raise ArithmeticError("too few crossings were found; lengthen the span")
        dm, dn = np.diff([e[1][4] for e in mu]), np.diff([e[1][4] for e in nu_])
        scale = 2 * math.sqrt(self.m)
        return {"K": K, "place": float(self.place(-4 * self.m * K)), "T_mu": float(dm.mean() * scale), "T_nu": float(dn.mean() * scale),
                "spread_mu": float(np.ptp(dm) * scale), "spread_nu": float(np.ptp(dn) * scale), "W": float(dn.mean() / dm.mean())}


# ----------------------------------------------------------------------------------------------------------------- the periods
def periods(K, eps, dps: int = 25):
    """T_μ, T_ν and the roots X₁ < X₂ < 1 < X₃ of Π(X) = K + X − (4/ε²)X(X − ¾)², at ``dps`` digits.

    T_ν = ∫ from X₁ to X₂ of dX/√R and T_μ = 2 ∫ from 1 to X₃ of dX/√R, R(X) = X(X − 1)Π(X)."""
    with mp.workdps(dps):
        g, K = 4 / mp.mpf(eps) ** 2, mp.mpf(K)
        X1, X2, X3 = sorted(mp.re(v) for v in mp.polyroots([-g, 3 * g / 2, 1 - 9 * g / 16, K], maxsteps=400, extraprec=400))

        def fn(th):
            X = (X1 + X2) / 2 + (X2 - X1) / 2 * mp.sin(th)
            return 1 / mp.sqrt(g * X * (1 - X) * (X3 - X))

        def fm(th):
            X = (1 + X3) / 2 + (X3 - 1) / 2 * mp.sin(th)
            return 1 / mp.sqrt(g * X * (X - X1) * (X - X2))
        Tn = mp.quad(fn, mp.linspace(-mp.pi / 2, mp.pi / 2, 5))
        Tm = 2 * mp.quad(fm, mp.linspace(-mp.pi / 2, mp.pi / 2, 5))
        return Tm, Tn, (X1, X2, X3)


def winding_axial(eps, dps: int = 25):
    """W on the axial line: T_μ = (4πε/3)(1 − 4ε²/9)^{−1/2} and T_ν an elementary integral."""
    with mp.workdps(dps):
        e = mp.mpf(eps)
        Tn = (e / 2) * mp.quad(lambda ph: 1 / ((mp.mpf(3) / 4 - e / 2 * mp.sin(ph)) * mp.sqrt(mp.mpf(1) / 4 + e / 2 * mp.sin(ph))), mp.linspace(-mp.pi / 2, mp.pi / 2, 5))
        return Tn / ((4 * mp.pi * e / 3) / mp.sqrt(1 - 4 * e * e / 9))


def table(times=TABLE_TIMES, fractions=TABLE_FRACTIONS, dps: int = 20) -> dict:
    """MP22's table: for each time ε, the leading law, and W − 1 on the axial line and on the tori at the given fractions."""
    out = {}
    for t in times:
        md = Hopfion(t)
        row = {"eps": md.eps, "leading": md.leading_law(), "axial": md.winding(dps=dps) - 1}
        for fr in fractions:
            row[fr] = md.winding(fraction=fr, dps=dps) - 1
        out[t] = row
    return out


# ----------------------------------------------------------------------------------------------------------------- lines in space
def line_field(t: float, kind: str = "B"):
    """dX/ds for a line of force of the squared hopfion at the time t, s being length on the label sphere.

    The same field as ``fields.line_field(t, kind, squared_field=True)``, written for one point in plain complex
    arithmetic, which is several times faster."""
    tc = t - 1j

    def f(X):
        x, y, z = float(X[0]), float(X[1]), float(X[2])
        a, b = complex(x, -y), tc - z
        d = x * x + y * y + z * z - tc * tc
        c = d ** 0.5 / d ** 3
        n = (4 * (b * b - a * a) * c, -4j * (a * a + b * b) * c, 8 * a * b * c)
        v = [w.imag for w in n] if kind == "B" else [w.real for w in n]
        nv = math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2])
        u = (v[0] / nv, v[1] / nv, v[2] / nv)
        xu = x * u[0] + y * u[1] + z * u[2]
        da = 2j * u[2] / d + 4j * b * xu / (d * d)
        db = 2 * complex(u[0], -u[1]) / d - 4 * a * xu / (d * d)
        sp = math.sqrt(abs(da) ** 2 + abs(db) ** 2)
        return np.array([u[0] / sp, u[1] / sp, u[2] / sp])
    return f


def line(t: float, X0, length: float, kind: str = "B", h: float = 0.1, **kw) -> dict:
    """A line of force of the squared hopfion through X0 at the time t, over ``length`` on the label sphere.
    One circuit is a length of about 2π."""
    return fields.integrate(line_field(t, kind), X0, length, h, **kw)


def invariant_drift(t: float, X0, circuits: float = 20, kind: str = "B", h: float = 0.1) -> dict:
    """I₂ and H along a line of force integrated in space: the drift of I₂ over its range, and of H."""
    md = Hopfion(t)
    run = line(t, X0, 2 * math.pi * circuits, kind, h)
    I = md.I2_at(run["y"], kind)
    Hs = md.H(*md.reduce(*fields.labels(t, run["y"]), kind))
    return {"place": float(md.place(I[0])), "drift": float(np.ptp(I) / md.I2_range), "H_off": float(np.abs(Hs - 1).max()),
            "farthest": float(np.linalg.norm(run["y"], axis=1).max()), "error": run["error"]}


def closed_line(t: float, X0, kind: str = "B", h: float = 0.02, max_length: float = 9.5) -> dict:
    """Follow the line through X0 until it first comes back through the plane across it at X0; for a closed line the
    return is X0.  Returns the label length of the circuit, the closure, and the points."""
    f = line_field(t, kind)
    X0 = np.asarray(X0, float)
    n0 = f(X0) / np.linalg.norm(f(X0))
    run = fields.integrate(f, X0, max_length, h, event=lambda X: float((X - X0) @ n0), max_events=1, min_s=1.0)
    if not run["events"]:
        raise ArithmeticError("the line did not return within the length allowed")
    s, Xr = run["events"][0]
    return {"length": s, "closure": float(np.linalg.norm(Xr - X0)), "points": run["y"][:-1], "error": run["error"]}


# ----------------------------------------------------------------------------------------------------------------- exact checks
def exact_checks() -> dict:
    """The derivation, exactly (sympy; a few seconds).  Every value should be True.

    The change of label; {H, I₂} = 0 for symbolic t and m; the second-order normal form at the 1 : 1 resonance, with
    the frequency difference (5/6)ε²; and W = 1 + (5/12)ε² + O(ε⁴) on the axial line."""
    import sympy as sp
    out = {}
    t, A = sp.symbols("t A")
    al = A + sp.I * t / 4 * (A - 1) ** 2
    out["alpha = A + (it/4)(A - 1)^2 gives (d alpha/dA)^2 = 1 - it(1 - alpha)"] = (sp.expand(sp.diff(al, A) ** 2 - (1 - sp.I * t * (1 - al))) == 0
                                                                               and sp.diff(al, A).subs(A, 1) == 1)
    x, y, px, py, m, tt = sp.symbols("x y p_x p_y m t", real=True)
    r2 = x * x + y * y
    V = (r2 ** 2 - 2 * m * (x * x - y * y) + m * m) / tt ** 2
    Hs = px ** 2 + py ** 2 + V
    I2s = (x * py - y * px) ** 2 + 4 * m * px ** 2 + (4 * m / tt ** 2) * x * x * (r2 - 2 * m)

    def pb(F, G):
        return sum(sp.diff(F, q) * sp.diff(G, p) - sp.diff(F, p) * sp.diff(G, q) for q, p in ((x, px), (y, py)))
    zc = x + sp.I * y
    out["{H, I_2} = 0"] = sp.expand(pb(Hs, I2s)) == 0 and sp.expand(V * tt ** 2 - sp.expand((zc ** 2 - m) * (sp.conjugate(zc) ** 2 - m))) == 0
    a, ab_, b, bb = sp.symbols("a abar b bbar")

    def br(F, G):
        return -2 * sp.I * (sp.diff(F, a) * sp.diff(G, ab_) - sp.diff(F, ab_) * sp.diff(G, a) + sp.diff(F, b) * sp.diff(G, bb) - sp.diff(F, bb) * sp.diff(G, b))
    u1, u2 = (a + ab_) / 2, (b + bb) / 2
    V1 = sp.expand(u1 * (u1 ** 2 + u2 ** 2))
    V2 = sp.expand((u1 ** 2 + u2 ** 2) ** 2 / 4)

    def weight(term):
        mm = sp.Poly(term, a, ab_, b, bb).monoms()[0]
        return mm[0] + mm[2] - mm[1] - mm[3]

    def harmonic(expr, k):
        return sum(term for term in sp.Add.make_args(sp.expand(expr)) if weight(term) == k)
    K2 = harmonic(V2, 0)
    for k in (1, 3):
        K2 += -sp.I * br(harmonic(V1, k), harmonic(V1, -k)) / (2 * k)
    na, nb = a * ab_, b * bb
    target = sp.expand((-12 * na ** 2 - sp.Rational(32, 3) * na * nb + sp.Rational(4, 3) * nb ** 2) / 32)
    Na, Nb = sp.symbols("n_a n_b", positive=True)
    K2n = (-12 * Na ** 2 - sp.Rational(32, 3) * Na * Nb + sp.Rational(4, 3) * Nb ** 2) / 32
    c1 = sp.expand(u1 ** 3)
    k1 = sp.expand(sum(-sp.I * br(harmonic(c1, k), harmonic(c1, -k)) / (2 * k) for k in (1, 3)))
    out["the normal form: -12, -32/3, 4/3 and no coupling of the phases"] = sp.expand(sp.expand(K2) - target) == 0
    out["the two frequencies differ by (5/6) eps^2 (n_a + n_b)"] = sp.simplify(2 * (sp.diff(K2n, Nb) - sp.diff(K2n, Na)) - sp.Rational(5, 6) * (Na + Nb)) == 0
    out["the method on the cubic oscillator: -15/32"] = sp.expand(k1 + sp.Rational(15, 32) * na ** 2) == 0
    es = sp.symbols("e", positive=True)
    ph = sp.symbols("phi", real=True)
    yv = es / 2 * sp.sin(ph)
    ser = sp.series(1 / ((sp.Rational(3, 4) - yv) * sp.sqrt(sp.Rational(1, 4) + yv)), es, 0, 3).removeO()
    Tn = es / 2 * sp.integrate(sp.expand(ser), (ph, -sp.pi / 2, sp.pi / 2))
    Wser = sp.series(Tn / (4 * sp.pi * es / 3) * sp.sqrt(1 - 4 * es ** 2 / 9), es, 0, 3).removeO()
    out["W = 1 + (5/12) eps^2 on the axial line"] = sp.simplify(Wser - 1 - sp.Rational(5, 12) * es ** 2) == 0
    return out
