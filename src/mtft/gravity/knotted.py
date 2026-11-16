r"""mtft.gravity.knotted — the (p, q) knotted fields: the tube round a core line, its rotation and its shear (MODEL; v1.1.0, 2026-11-16).

The fields F_pq = ∇α^p × ∇β^q, their invariant Ψ_B = Re(α^p β^q) and the tori it cuts out are Kedia, Bialynicki-Birula,
Peralta-Salas and Irvine's (Phys. Rev. Lett. 111, 150404, 2013).  On the label sphere write α = cos η e^{iξ₁},
β = sin η e^{iξ₂}, Θ = pξ₁ + qξ₂.  A magnetic line keeps H = cos^p η sin^q η cos Θ; the two core lines are
cos²η = p/(p + q), Θ = 0 or π, where H = ±h_max.  One circuit is ξ₁ → ξ₁ + 2πq, ξ₂ → ξ₂ − 2πp.

What this module adds to that, from the lab's records MP09 and MP21 (derived there, tested by three routes,
unreviewed; a search of the literature on 16 November 2026 did not find them, which does not make them new):

  * the lines next to a core line turn round it by 2πν per circuit in the torus frame, ν = √(2pq(p + q))
    (``nu``, ``generator``, ``tube_matrix``); the matrix is the identity exactly when 2pq(p + q) is a square;
  * on the torus H = h the winding is ι(h) = 2πq / ∮ dx/√(e^{λ(x) − 2 ln h} − 1) (``iota``), from ν at the core
    to 2pq at the edge;
  * the shear at the core line, with u = 1 − h/h_max (``shear``):
        dι/du = ν · [3pq(p + q) − 2(p² + pq + q²)] / (12pq(p + q)),
    zero for the hopfion only.

Three routes to the shear live here: the closed form; the profile integral at 30 digits; and the period of the
cross-section's Hamiltonian system integrated as a differential equation (``iota_ode``).  ``helicity`` and
``energy`` are the closed forms 2π² pq·p! q!/(p + q)! and (p + q) times it; their ratio 1/(p + q) is Kedia et al.'s.

Open, recorded and not tested (H-G7): for (1, 8) the circuit matrix is the identity while the shear is 35/36, so
the core line is degenerate and a perturbation should move the closed lines near it by its cube root.
"""
from __future__ import annotations

import math
from fractions import Fraction
from math import factorial, gcd, isqrt

import mpmath as mp
import numpy as np

from . import fields

PAIRS = ((1, 1), (1, 2), (2, 3), (3, 2), (2, 5), (1, 8))            # the lab's table (MP21)
OPEN = {"H-G7": "the core line of (1, 8) is degenerate (circuit matrix the identity, shear 35/36): closed lines near it should move "
                "by the cube root of a perturbation's strength; a hypothesis of MP21, not tested"}


# ----------------------------------------------------------------------------------------------------------------- the tube
def nu_squared(p: int, q: int) -> int:
    """2pq(p + q): the square of the rotation number at the core line."""
    return 2 * p * q * (p + q)


def nu(p: int, q: int) -> float:
    """ν = √(2pq(p + q)): turns of a neighbouring line round the core line per circuit, in the torus frame."""
    return math.sqrt(nu_squared(p, q))


def is_whole(p: int, q: int) -> bool:
    """True when ν is a whole number, which is when the matrix of one circuit is the identity."""
    r = isqrt(nu_squared(p, q))
    return r * r == nu_squared(p, q)


def whole_pairs(limit: int = 60, coprime: bool = True) -> list:
    """The pairs p ≤ q ≤ limit with ν a whole number; ``coprime`` keeps the pairs that give a knot or a single link."""
    return [(p, q) for p in range(1, limit + 1) for q in range(p, limit + 1) if is_whole(p, q) and (not coprime or gcd(p, q) == 1)]


def pell(count: int = 4) -> list:
    """For p = 1 the condition is Pell's equation (2q + 1)² − 2m² = 1: the first ``count`` values of q."""
    out, X, m = [], 3, 2
    for _ in range(count):
        out.append((X - 1) // 2)
        X, m = 3 * X + 4 * m, 2 * X + 3 * m
    return out


def generator(p: int, q: int) -> tuple:
    """The linear part of the line field at the core line, in the frame (δη, δΘ): L = [[0, −s], [2(p + q)s, 0]],
    s = √(pq)/(p + q).  det L = 2pq/(p + q); one circuit takes the parameter time 2π(p + q)."""
    s = math.sqrt(p * q) / (p + q)
    return ((0.0, -s), (2 * (p + q) * s, 0.0))


def tube_matrix(p: int, q: int) -> tuple:
    """The matrix of one circuit, exp(2π(p + q) L): a rotation by 2πν in the torus frame."""
    w, r = 2 * math.pi * nu(p, q), math.sqrt(2 * (p + q))
    if is_whole(p, q):
        return ((1.0, 0.0), (0.0, 1.0))
    return ((math.cos(w), -math.sin(w) / r), (r * math.sin(w), math.cos(w)))


def seifert_winding(p: int, q: int) -> float:
    """ν − pq: turns of the neighbours round the core line per circuit against the Seifert framing."""
    return nu(p, q) - p * q


def h_max(p: int, q: int) -> float:
    """The value of H on the core line: (p/(p + q))^{p/2} (q/(p + q))^{q/2} (Kedia et al.'s extreme of Ψ_B)."""
    return (p / (p + q)) ** (p / 2) * (q / (p + q)) ** (q / 2)


def H(alpha, beta, p: int, q: int):
    """Re(α^p β^q), constant along every magnetic line; Im(α^p β^q) is constant along every electric line."""
    return np.real(np.asarray(alpha) ** p * np.asarray(beta) ** q)


# ----------------------------------------------------------------------------------------------------------------- the shear
def shear_factor(p: int, q: int) -> Fraction:
    """[3pq(p + q) − 2(p² + pq + q²)] / (12pq(p + q)), exactly; for p = 1 the bracket is (q + 2)(q − 1)."""
    return Fraction(3 * p * q * (p + q) - 2 * (p * p + p * q + q * q), 12 * p * q * (p + q))


def shear(p: int, q: int) -> float:
    """dι/du at the core line, u = 1 − h/h_max: ν times ``shear_factor``."""
    return nu(p, q) * float(shear_factor(p, q))


def shear_dh(p: int, q: int) -> float:
    """dι/dh at the core line: −shear / h_max."""
    return -shear(p, q) / h_max(p, q)


def area_coefficients(p: int, q: int) -> tuple:
    """a, b, c of ε = a y² + b y³ + c y⁴ + Θ² + Θ⁴/6 near the core line (y = x − x₀, ε = 2 ln(h_max/h)), exactly."""
    n = p + q
    return (Fraction(p * n, 2 * q), Fraction(p * n * (2 * p + q), 6 * q * q), Fraction(p * n * (n * n + 4 * p * n + p * p), 24 * q ** 3))


def table() -> dict:
    """The lab's table (MP21): ν, dι/du, h_max and dι/dh at the core line for its six pairs."""
    return {(p, q): {"nu": nu(p, q), "shear": shear(p, q), "h_max": h_max(p, q), "shear_dh": shear_dh(p, q)} for p, q in PAIRS}


# ----------------------------------------------------------------------------------------------------------------- the profile
def _level_roots(pp, qq, lw, xc, l2h):
    """The two roots of λ(x) = 2 ln h, one each side of the core line x_c.  Since λ(x) < px and λ(x) < q ln(1 − eˣ), the
    points x = (2 ln h)/p − 1 and ln(1 − e^{(2 ln h − 1)/q}) lie outside them, for every level."""
    lo, hi = l2h / pp - 1, mp.log1p(-mp.exp((l2h - 1) / qq))
    x1 = mp.findroot(lambda x: lw(x) - l2h, (lo, xc), solver="anderson")
    x2 = mp.findroot(lambda x: lw(x) - l2h, (xc, hi), solver="anderson")
    if not (lo < x1 < xc < x2 < hi):
        raise ArithmeticError("the two roots of the level were not bracketed")
    return x1, x2


def iota(p: int, q: int, u, dps: int = 34):
    """ι on the torus H = h_max(1 − u), by MP09's integral at ``dps`` digits (an mpmath number).

    ι = 2πq / ∫ dx/√(e^{λ(x) − 2 ln h} − 1) between the two roots of λ(x) = 2 ln h, λ(x) = px + q ln(1 − eˣ)."""
    with mp.workdps(dps):
        pp, qq, u = mp.mpf(p), mp.mpf(q), mp.mpf(u)

        def lw(x):
            return pp * x + qq * mp.log1p(-mp.exp(x))
        xc = mp.log(pp / (pp + qq))
        l2h = lw(xc) + 2 * mp.log1p(-u)
        x1, x2 = _level_roots(pp, qq, lw, xc, l2h)
        mid, half, dmax = (x1 + x2) / 2, (x2 - x1) / 2, lw(xc) - l2h

        def slope(xe):
            return abs(pp - qq * mp.exp(xe) / (1 - mp.exp(xe)))

        def f(th):
            d = lw(mid + half * mp.sin(th)) - l2h
            if d < mp.mpf(10) ** (12 - dps) * dmax:                                      # the last sliver by its tangent
                d = slope(x2 if th > 0 else x1) * half * 2 * mp.sin((mp.pi / 2 - abs(th)) / 2) ** 2
            return half * mp.cos(th) / mp.sqrt(mp.expm1(d))
        return 2 * mp.pi * qq / mp.quad(f, mp.linspace(-mp.pi / 2, mp.pi / 2, 9))


def iota_edge(p: int, q: int) -> int:
    """The winding at the edge H = 0: 2pq."""
    return 2 * p * q


def richardson3(f1, f2, f4):
    """From f(u), f(2u), f(4u) with f = f₀ + c₁u + c₂u² + …: f₀ to O(u³)."""
    return (8 * f1 - 6 * f2 + f4) / 3


def shear_from_profile(p: int, q: int, u="1e-3", dps: int = 34) -> float:
    """The slope of ι at the core line from the profile integral at u, 2u and 4u, extrapolated: the second route."""
    with mp.workdps(dps):
        u, v = mp.mpf(u), mp.sqrt(mp.mpf(nu_squared(p, q)))
        s = [(iota(p, q, k * u, dps) - v) / (k * u) for k in (1, 2, 4)]
        return float(richardson3(*s))


def section_roots(p: int, q: int, u: float) -> tuple:
    """The two values of x = ln cos²η at which the torus H = h_max(1 − u) meets Θ = 0."""
    with mp.workdps(30):
        pp, qq = mp.mpf(p), mp.mpf(q)

        def lw(x):
            return pp * x + qq * mp.log1p(-mp.exp(x))
        xc = mp.log(pp / (pp + qq))
        x1, x2 = _level_roots(pp, qq, lw, xc, lw(xc) + 2 * mp.log1p(-mp.mpf(u)))
        return float(x1), float(x2)


def iota_ode(p: int, q: int, u: float, steps: int = 48) -> float:
    """ι on the torus H = h_max(1 − u) by integrating the cross-section's Hamiltonian system: the third route.

    With x = ln cos²η and ξ₁ as the time, dx/dξ₁ = 2 tan Θ and dΘ/dξ₁ = p − q eˣ/(1 − eˣ), which are Hamilton's
    equations for 2 ln H.  One turn round the core line takes ξ₁ the time P, and ι = 2πq/P."""
    _, x2 = section_roots(p, q, u)

    def f(y):
        ex = math.exp(y[0])
        return np.array([2 * math.tan(y[1]), p - q * ex / (1 - ex)])
    guess = 2 * math.pi * q / nu(p, q)
    run = fields.integrate(f, [x2, 0.0], 1.5 * guess, guess / steps, event=lambda y: -y[1], max_events=1)
    if not run["events"]:
        raise ArithmeticError("the section did not close within one and a half periods")
    return 2 * math.pi * q / run["events"][0][0]


# ----------------------------------------------------------------------------------------------------------------- helicity, energy
def helicity(p: int, q: int) -> float:
    """∫ A ∧ dA over the label sphere for A = ½ Re(α^p dβ^q − β^q dα^p), in closed form: 2π² pq·p! q!/(p + q)!."""
    return 2 * math.pi ** 2 * p * q * factorial(p) * factorial(q) / factorial(p + q)


def energy(p: int, q: int) -> float:
    """½ ∫ |F_pq|² d³x, in closed form: (p + q) times ``helicity`` (the ratio is Kedia et al.'s)."""
    return (p + q) * helicity(p, q)


def helicity_on_sphere(p: int, q: int, n_eta: int = 48, n_xi: int = 8) -> float:
    """The same helicity by quadrature of A ∧ dA on the label sphere (Gauss–Legendre in η, the trapezoid rule in ξ₁, ξ₂).

    The density does not depend on the two angles: its terms of frequency 2p and 2q are a product of three holomorphic
    1-forms in two variables, which is zero.  So a few points in each angle are enough."""
    g, w = np.polynomial.legendre.leggauss(n_eta)
    eta, w = (g + 1) * math.pi / 4, w * math.pi / 4
    xi = np.arange(n_xi) * 2 * math.pi / n_xi
    Et, X1, X2 = np.meshgrid(eta, xi, xi, indexing="ij")
    al, be = np.cos(Et) * np.exp(1j * X1), np.sin(Et) * np.exp(1j * X2)
    dal = np.array([1j * al, -np.sin(Et) * np.exp(1j * X1), 0 * al])                      # along ξ₁, η, ξ₂: the orientation space gives the sphere
    dbe = np.array([0 * be, np.cos(Et) * np.exp(1j * X2), 1j * be])
    f, g_ = al ** p, be ** q
    df, dg = p * al ** (p - 1) * dal, q * be ** (q - 1) * dbe
    A = 0.5 * np.real(f * dg - g_ * df)

    def dA(j, k):
        return np.real(df[j] * dg[k] - df[k] * dg[j])
    dens = A[0] * dA(1, 2) - A[1] * dA(0, 2) + A[2] * dA(0, 1)
    return float(np.sum(dens * w[:, None, None]) * (2 * math.pi / n_xi) ** 2)


def energy_in_space(p: int, q: int, t: float = 0.0, n: int = 96, L: float = 1.2) -> float:
    """½ ∫ |F_pq|² d³x on the slice of time t, by Gauss–Legendre quadrature on space mapped to a cube."""
    g, w = np.polynomial.legendre.leggauss(n)
    xs, ws = L * g / (1 - g * g), w * L * (1 + g * g) / (1 - g * g) ** 2
    X, Y, Z = np.meshgrid(xs, xs, xs, indexing="ij")
    W = (ws[:, None, None] * ws[None, :, None] * ws[None, None, :]).reshape(-1)
    Fv = fields.F(t, np.stack([X.reshape(-1), Y.reshape(-1), Z.reshape(-1)], -1), p, q)
    return float(0.5 * np.sum(np.sum(np.abs(Fv) ** 2, -1) * W))


# ----------------------------------------------------------------------------------------------------------------- exact checks
def exact_checks() -> dict:
    """The derivation, exactly, for symbolic p and q (sympy; about 2 s).  Every value should be True.

    From MP09's line field on the label sphere: H is kept; the cross-section is a Hamiltonian system for 2 ln H; its
    linear part at the core line is ``generator``, and one circuit is a rotation by 2π√(2pq(p + q)); the coefficients
    of ``area_coefficients``; the area of a level to second order; and the bracket of ``shear``."""
    import sympy as sp
    eta, Th, x = sp.symbols("eta Theta x", real=True)
    P, Q = sp.symbols("p q", positive=True)
    d_eta = -sp.sin(Th) * sp.sin(eta) * sp.cos(eta)
    d_xi1, d_xi2 = sp.cos(Th) * sp.sin(eta) ** 2, -sp.cos(Th) * sp.cos(eta) ** 2
    d_Th = P * d_xi1 + Q * d_xi2
    out = {}
    Hs = sp.cos(eta) ** P * sp.sin(eta) ** Q * sp.cos(Th)
    out["H is kept along the line field"] = sp.simplify((sp.diff(Hs, eta) * d_eta + sp.diff(Hs, Th) * d_Th) / Hs) == 0
    dx = sp.simplify(sp.diff(sp.log(sp.cos(eta) ** 2), eta) * d_eta / d_xi1)
    dT = sp.simplify(d_Th / d_xi1)
    lam = P * x + Q * sp.log(1 - sp.exp(x))
    Ham = lam + 2 * sp.log(sp.cos(Th))
    sub = {eta: sp.acos(sp.exp(x / 2))}
    rhs = P - Q * sp.exp(x) / (1 - sp.exp(x))
    out["the cross-section is Hamiltonian for 2 ln H"] = (sp.simplify((dx - 2 * sp.tan(Th)).subs(sub)) == 0 and sp.simplify(sp.simplify(dT.subs(sub)) - rhs) == 0
                                                          and sp.simplify(-sp.diff(Ham, Th) - 2 * sp.tan(Th)) == 0 and sp.simplify(sp.diff(Ham, x) - rhs) == 0)
    core = {eta: sp.acos(sp.sqrt(P / (P + Q))), Th: 0}
    J = sp.Matrix([[sp.diff(d_eta, eta), sp.diff(d_eta, Th)], [sp.diff(d_Th, eta), sp.diff(d_Th, Th)]]).subs(core)
    s = sp.sqrt(P * Q) / (P + Q)
    T = 2 * sp.pi * (P + Q)
    out["the linear part at the core line is [[0, -s], [2(p + q)s, 0]]"] = sp.simplify(J - sp.Matrix([[0, -s], [2 * (P + Q) * s, 0]])) == sp.zeros(2, 2)
    out["one circuit is a rotation by 2 pi sqrt(2pq(p + q))"] = (sp.simplify(J.det() - 2 * P * Q / (P + Q)) == 0 and sp.simplify(J.trace()) == 0
                                                                and sp.simplify(d_xi1.subs(core) * T - 2 * sp.pi * Q) == 0
                                                                and sp.simplify((T * sp.sqrt(2 * P * Q / (P + Q))) ** 2 - (2 * sp.pi) ** 2 * 2 * P * Q * (P + Q)) == 0)
    x0 = sp.log(P / (P + Q))
    l2, l3, l4 = [sp.simplify(sp.diff(lam, x, k).subs(x, x0)) for k in (2, 3, 4)]
    a, b, c = sp.simplify(-l2 / 2), sp.simplify(-l3 / 6), sp.simplify(-l4 / 24)
    n = P + Q
    out["the coefficients a, b, c"] = (sp.simplify(a - P * n / (2 * Q)) == 0 and sp.simplify(b - P * n * (2 * P + Q) / (6 * Q ** 2)) == 0
                                       and sp.simplify(c - P * n * (n ** 2 + 4 * P * n + P ** 2) / (24 * Q ** 3)) == 0
                                       and sp.simplify(sp.series(-2 * sp.log(sp.cos(Th)), Th, 0, 6).removeO() - (Th ** 2 + Th ** 4 / 6)) == 0)
    K = sp.Rational(15, 32) * b ** 2 / a ** 3 - sp.Rational(3, 8) * c / a ** 2 - sp.Rational(3, 8) * sp.Rational(1, 6)
    vv = sp.sqrt(2 * P * Q * n)
    slope = sp.simplify(2 * sp.pi * Q / (sp.pi / sp.sqrt(a)) * (-2 * K) * 2)
    closed = vv * (3 * P * Q * n - 2 * (P ** 2 + P * Q + Q ** 2)) / (12 * P * Q * n)
    out["the bracket of the shear"] = sp.simplify(slope - closed) == 0 and sp.simplify(2 * Q * sp.sqrt(a) - vv) == 0
    phi, eps = sp.symbols("phi epsilon", positive=True)
    aa, bb, cc, be = sp.symbols("a b c beta", positive=True)
    r1 = -bb / aa ** sp.Rational(3, 2) * sp.cos(phi) ** 3
    r2 = sp.Rational(3, 2) * (bb ** 2 / aa ** 3) * sp.cos(phi) ** 6 - (cc / aa ** 2) * sp.cos(phi) ** 4 - be * sp.sin(phi) ** 4
    rho2 = 1 + sp.sqrt(eps) * r1 + eps * r2
    lhs = rho2 + sp.sqrt(eps) * (bb / aa ** sp.Rational(3, 2)) * rho2 ** sp.Rational(3, 2) * sp.cos(phi) ** 3 + eps * rho2 ** 2 * ((cc / aa ** 2) * sp.cos(phi) ** 4 + be * sp.sin(phi) ** 4)
    res = sp.series(lhs - 1, eps, 0, sp.Rational(3, 2)).removeO()
    area2 = sp.integrate(r2, (phi, 0, 2 * sp.pi)) / 2
    out["the area of a level to second order"] = (sp.simplify(res) == 0 and sp.simplify(area2 / sp.pi - (sp.Rational(15, 32) * bb ** 2 / aa ** 3 - sp.Rational(3, 8) * cc / aa ** 2
                                                                                                    - sp.Rational(3, 8) * be)) == 0)
    return out
