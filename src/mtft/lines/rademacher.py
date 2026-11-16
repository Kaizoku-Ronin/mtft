r"""mtft.lines.rademacher — Dedekind sums, Rademacher's Ψ and the windings round the cusps (EXACT; LN-01, 2026-11-16).

Rademacher's function.  For γ = [[a, b], [c, d]] in SL(2, ℤ),

    Φ(γ) = b/d                                  (c = 0),
    Φ(γ) = (a + d)/c − 12 sign(c) s(d, |c|)     (c ≠ 0),        Ψ(γ) = Φ(γ) − 3 sign(c(a + d)),

s(h, k) the Dedekind sum.  Φ is the change of log η: log η(γτ) − log η(τ) = ½ log((cτ + d)/i) +
(πi/12) Φ(γ) for c > 0.  Ψ is an integer, a class function of SL(2, ℤ), and for a hyperbolic γ it is
the linking number of the closed orbit of the geodesic flow with the trefoil at level 1 (Ghys).  For
γ conjugate to a positive word in U = T and V = [[1, 0], [1, 1]], Ψ = #U − #V (``positive_word``).

In the 3-fold Γ₀(N)\PSL(2, ℝ), Ψ is the coordinate of a closed line along the fibre: the class of
a closed line is (Ψ; cls), cls being ``mtft.lines.darts.Graph.cls``, with Ψ ≡ χ(cls) modulo 6.

Windings.  A unit u = ∏ η(δτ)^{r_δ} (Σ r_δ = 0, with the usual congruences) is a function on
Γ₀(N)\ℍ without zeros or poles, and

    φ_u(γ) = (1/2πi) ∫_τ^{γτ} d log u = (1/24) Σ_δ r_δ Φ(γ_δ),      γ_δ = [[a, δ b], [c/δ, d]],

is a homomorphism from Γ₀(N) to ℤ: the number of times u winds round 0 along the loop.  On the
loop round a cusp it is the order of u at that cusp (``cusp_order``, Ligozat's formula: the two
are computed apart and compared in the tests).  At N = 143 three units give the rank-3 part of
H¹ carried by the cusps (``UNITS``):

    u₁ = η(τ)η(143τ)/(η(11τ)η(13τ)),    u_A = (η(τ)/η(11τ))¹²,    u_B = (η(τ)/η(13τ))².

u_A is a function on X₀(11) and u_B on X₀(13); only φ₁ belongs to level 143 alone.
"""
from __future__ import annotations

import math
from fractions import Fraction as Fr

UNITS = {"phi1": {1: 1, 11: -1, 13: -1, 143: 1}, "phiA": {1: 12, 11: -12}, "phiB": {1: 2, 13: -2}}
"""The three η-quotients of level 143, as exponents r_δ on η(δτ)."""

UNIT_ORDER = ("phi1", "phiA", "phiB")


def sign(v) -> int:
    return (v > 0) - (v < 0)


def dedekind(h: int, k: int) -> Fr:
    """The Dedekind sum s(h, k) for k > 0 and gcd(h, k) = 1, by the reciprocity law."""
    if k <= 0:
        raise ValueError("k must be positive")
    if math.gcd(h, k) != 1:
        raise ValueError("h and k must be coprime")
    h %= k
    s, sg = Fr(0), 1
    while h:
        s += sg * ((Fr(h, k) + Fr(k, h) + Fr(1, h * k)) / 12 - Fr(1, 4))
        sg = -sg
        h, k = k % h, h
    return s


def dedekind_by_definition(h: int, k: int) -> Fr:
    """s(h, k) = Σ_{r mod k} ((r/k)) ((hr/k)), straight from the definition (slow; the check of ``dedekind``)."""
    def saw(x):
        return Fr(0) if x.denominator == 1 else x - math.floor(x) - Fr(1, 2)
    return sum((saw(Fr(r, k)) * saw(Fr(h * r, k)) for r in range(1, k)), Fr(0))


def Phi(g) -> Fr:
    """Rademacher's Φ."""
    (a, b), (c, d) = g
    if c == 0:
        return Fr(b, d)
    return Fr(a + d, c) - 12 * sign(c) * dedekind(d, abs(c))


def Psi(g) -> int:
    """Rademacher's Ψ = Φ − 3 sign(c(a + d)), an integer and a class function of SL(2, ℤ)."""
    (a, _), (c, d) = g
    v = Phi(g) - 3 * sign(c * (a + d))
    if v.denominator != 1:
        raise ArithmeticError("Ψ is not an integer: the matrix is not in SL(2, ℤ)")
    return int(v)


def winding(g, r: dict) -> Fr:
    """φ_u(γ) for the unit u = ∏ η(δτ)^{r_δ}; γ must lie in Γ₀(N) for a level N that every δ divides."""
    (a, b), (c, d) = g
    if c < 0 or (c == 0 and d < 0):
        a, b, c, d = -a, -b, -c, -d
    tot = Fr(0)
    for dl, e in r.items():
        if c % dl:
            raise ValueError(f"the lower-left entry must be divisible by {dl}")
        tot += e * (Fr(dl * b, d) if c == 0 else Fr(dl * (a + d), c) - 12 * dedekind(d, c // dl))
    return tot / 24


def windings(g) -> tuple:
    """(φ₁, φ_A, φ_B) of an element of Γ₀(143), as integers."""
    out = []
    for k in UNIT_ORDER:
        v = winding(g, UNITS[k])
        if v.denominator != 1:
            raise ArithmeticError("a winding is not an integer: the matrix is not in Γ₀(143)")
        out.append(int(v))
    return tuple(out)


def cusp_order(r: dict, c: int, N: int) -> Fr:
    """The order of ∏ η(δτ)^{r_δ} at the cusp of denominator c | N of X₀(N) (Ligozat's formula).

    ord = (N/24) Σ_δ gcd(c, δ)² r_δ / (gcd(c, N/c) c δ), in the cusp's own uniformiser.
    """
    if N % c:
        raise ValueError("the denominator of the cusp must divide the level")
    return Fr(N, 24) * sum(Fr(math.gcd(c, dl) ** 2 * e, math.gcd(c, N // c) * c * dl) for dl, e in r.items())


def positive_word(A) -> tuple:
    """(#U, #V) of a positive word conjugate to the hyperbolic matrix A, U = T and V = [[1, 0], [1, 1]].

    The word is ``mtft.lines.lorenz.lorenz_word`` (x = U, y = V).  Ψ(A) = #U − #V: the Dedekind sums
    and the count of letters are two routes to one integer.
    """
    from .lorenz import lorenz_word
    w = lorenz_word(A)
    return w.count("x"), w.count("y")
