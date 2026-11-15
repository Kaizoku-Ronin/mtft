"""mtft.maass.ford — the Ford domain of Γ₀*(N) and integrals over it (CERTIFIED numerics).

For squarefree N the group Γ₀*(N) (Γ₀(N) with all its Atkin–Lehner matrices) has one cusp, and its
Ford domain F* is the part of the strip |x| ≤ ½ above every isometric circle.  The circle of
[[Q a, b], [N c, Q d]] (determinant Q, Q | N) has centre −Q d/(N c) and radius √Q/(N c).  The
domain is symmetric under x ↦ −x, so everything is computed on 0 ≤ x ≤ ½ and doubled.

Area.  F* has area (π/3)·ψ(N)/2^ω(N), ψ(N) = N ∏(1 + 1/p): 14π at N = 143.

Quadrature (the lab's record MP20).  The set of allowed x at height y changes only at the tops of
circles and at vertices, so the y-axis is cut into panels at those heights.  Inside a panel the
substitution y = y_b − (y_b − y_a)τ² removes the square-root behaviour below a top; above the
highest top the panels are logarithmic up to y = 4, and the strip above 4 is added in closed form
for the area (a cusp form is below e⁻²⁵ of its size there).  At N = 143 the area comes out as 14π
to 2e-11 and two forms of one sector are orthogonal to 2e-10.

``integrate`` takes forms as dictionaries with the keys ``R``, ``parity`` (0 even, 1 odd),
``twist`` and ``c`` (the coefficient vector, c[0] = 1), as ``mtft.maass.spectrum`` supplies them.

Symmetry.  An integral over F* means something only for an integrand that Γ₀*(N) and the mirror
x ↦ −x leave unchanged: then it does not depend on the choice of domain, and 2^ω(N) times it is the
integral over the curve X₀(N).  A product of forms qualifies when its parities, its Atkin–Lehner
signs and its twists multiply to the trivial character (``invariant``).  When they do not, the
integral over the curve is zero by symmetry, while the integral over F* is neither zero nor
meaningful; ``gram`` and ``triple`` therefore return 0 for such products without integrating.
"""
from __future__ import annotations

import math
from functools import lru_cache

import numpy as np

from .domain import is_squarefree, prime_factors
from .kbessel import k_bessel

Y_MAX = 4.0


def area(N: int = 143) -> float:
    """The hyperbolic area of Γ₀*(N)\\ℍ for squarefree N: (π/3) ψ(N) / 2^ω(N)."""
    if not is_squarefree(N):
        raise ValueError("the Ford domain of Gamma_0*(N) is implemented for squarefree N")
    ps = prime_factors(N)
    psi = 1
    for p in ps:
        psi *= p + 1
    return math.pi / 3 * psi / 2 ** len(ps)


def circles(N: int = 143, rmin: float = 0.004) -> list:
    """The isometric circles (centre, radius) of radius > rmin that meet the half-strip 0 ≤ x ≤ ½."""
    if not is_squarefree(N):
        raise ValueError("the Ford domain of Gamma_0*(N) is implemented for squarefree N")
    ps = prime_factors(N)
    Qs = sorted({math.prod(ps[i] for i in range(len(ps)) if k >> i & 1) for k in range(2 ** len(ps))})
    out = set()
    for Q in Qs:
        Nq, c = N // Q, 1
        while math.sqrt(Q) / (N * c) > rmin:
            r = math.sqrt(Q) / (N * c)
            for d in range(-N * c, N * c + 1):
                if math.gcd(Q * d, Nq * c) != 1:
                    continue
                xc = -Q * d / (N * c)
                if -r < xc < 0.5 + r:
                    out.add((xc, r))
            c += 1
    return sorted(out)


def intervals(y: float, circ: list) -> list:
    """The x in [0, ½] outside every circle at height y, as a list of intervals."""
    cut = sorted((xc - math.sqrt(r * r - y * y), xc + math.sqrt(r * r - y * y)) for xc, r in circ if r > y)
    out, a = [], 0.0
    for lo, hi in cut:
        if lo > a:
            out.append((a, min(lo, 0.5)))
        a = max(a, hi)
        if a >= 0.5:
            break
    if a < 0.5:
        out.append((a, 0.5))
    return [(p, q) for p, q in out if q - p > 1e-13]


def breakpoints(circ: list):
    """Heights where the set of allowed x changes: the tops of circles and the vertices.

    A top or a pairwise intersection counts when it lies in the closed half-strip and strictly
    inside no other circle.  Returns (tops, vertices), each sorted by height.
    """
    C = np.array(circ)
    xc_all, r2_all = C[:, 0], C[:, 1] ** 2 * (1 - 1e-10)

    def free(x, y, skip):
        if not (-1e-12 <= x <= 0.5 + 1e-12):
            return False
        inside = (x - xc_all) ** 2 + y * y < r2_all
        for k in skip:
            inside[k] = False
        return not inside.any()

    tops = [r for k, (xc, r) in enumerate(circ) if free(xc, r, (k,))]
    verts = []
    for i in range(len(circ)):
        x1, r1 = circ[i]
        for j in range(i + 1, len(circ)):
            x2, r2 = circ[j]
            d = x2 - x1
            if d < 1e-14 or d >= r1 + r2 or d <= abs(r1 - r2):
                continue
            x = x1 + (d * d + r1 * r1 - r2 * r2) / (2 * d)
            y2 = r1 * r1 - (x - x1) ** 2
            if y2 > 0 and free(x, math.sqrt(y2), (i, j)):
                verts.append((x, math.sqrt(y2)))
    return sorted(tops), sorted(verts, key=lambda v: v[1])


@lru_cache(maxsize=8)
def lowest_vertex(N: int = 143) -> tuple:
    """The lowest point (x, y) of the Ford domain of Γ₀*(N): its lowest vertex.

    Circles of radius below a height cannot reach it, so the search widens the family of circles
    until the lowest vertex found lies well above the smallest radius kept.
    """
    if N == 1:
        return (0.5, math.sqrt(3) / 2)
    rmin = 1.0 / N
    for _ in range(12):
        verts = breakpoints(circles(N, rmin))[1]
        if verts and verts[0][1] > 1.5 * rmin:
            return verts[0]
        rmin /= 2
    raise RuntimeError("the lowest vertex was not found")


class FordDomain:
    """The Ford domain of Γ₀*(N): its circles, tops, vertices, and a quadrature for dμ = dx dy/y²."""

    def __init__(self, N: int = 143, rmin: float | None = None):
        self.N = N
        self.area = area(N)
        low = lowest_vertex(N)
        self.rmin = 0.4 * low[1] if rmin is None else rmin
        if N == 1:
            self.circles = [(0.0, 1.0), (1.0, 1.0)]
            self.tops, self.vertices = [1.0], [low]
        else:
            self.circles = circles(N, self.rmin)
            self.tops, self.vertices = breakpoints(self.circles)
        self.lowest_vertex = self.vertices[0]

    def intervals(self, y: float) -> list:
        return intervals(y, self.circles)

    def heights(self) -> list:
        """The distinct breakpoint heights, from the lowest vertex to the highest top."""
        ys = sorted(set([round(t, 14) for t in self.tops] + [round(v[1], 14) for v in self.vertices]))
        bp = [ys[0]]
        for v in ys[1:]:
            if v - bp[-1] > 1e-11:
                bp.append(v)
        return bp

    def rows(self, n: int = 16, n_top: int = 40, panels_top: int = 12) -> list:
        """Quadrature rows (y, weight): the weight is dy/y² for both halves of the domain."""
        bp = self.heights()
        t0, w0 = np.polynomial.legendre.leggauss(n)
        tau, wt = 0.5 * (t0 + 1), 0.5 * w0
        rows = []
        for ya, yb in zip(bp[:-1], bp[1:]):
            for t, w in zip(tau, wt):
                y = yb - (yb - ya) * t * t
                rows.append((y, 2.0 * 2 * (yb - ya) * t * w / (y * y)))
        t1, w1 = np.polynomial.legendre.leggauss(n_top)
        edges = np.linspace(math.log(bp[-1]), math.log(Y_MAX), panels_top + 1)
        for a, b in zip(edges[:-1], edges[1:]):
            for t, w in zip(t1, w1):
                y = math.exp(0.5 * (b - a) * t + 0.5 * (a + b))
                rows.append((y, 2.0 * 0.5 * (b - a) * w / y))
        return rows

    def nodes(self, y: float, bandwidth: float, xfac: float = 1.0):
        """Gauss–Legendre nodes and weights in x on the allowed intervals at height y.

        ``bandwidth`` is the highest frequency of the integrand in x; the node count of an
        interval of length L is xfac·(2π·bandwidth·L + 30).
        """
        xs, ws = [], []
        for a, b in self.intervals(y):
            x0, w0 = _gl(int(xfac * (2 * math.pi * bandwidth * (b - a) + 30)))
            xs.append(0.5 * (b - a) * x0 + 0.5 * (a + b))
            ws.append(0.5 * (b - a) * w0)
        if not xs:
            return np.zeros(0), np.zeros(0)
        return np.concatenate(xs), np.concatenate(ws)

    def area_by_quadrature(self, n: int = 16) -> float:
        """The area from the quadrature, to compare with ``area``: a check of the domain itself."""
        tot = 0.0
        for y, wy in self.rows(n):
            tot += wy * sum(b - a for a, b in self.intervals(y))
        return tot + 1.0 / Y_MAX

    def summary(self) -> dict:
        return {"N": self.N, "area": self.area, "circles": len(self.circles), "tops": len(self.tops), "vertices": len(self.vertices),
                "distinct heights": len(self.heights()), "lowest vertex": self.lowest_vertex, "highest top": self.tops[-1]}


_GL = {}


def _gl(n: int):
    if n not in _GL:
        _GL[n] = np.polynomial.legendre.leggauss(n)
    return _GL[n]


@lru_cache(maxsize=4)
def domain(N: int = 143) -> FordDomain:
    return FordDomain(N)


def values_on_row(forms: list, x: np.ndarray, y: float) -> np.ndarray:
    """The forms' values at the points (x, y) of one height, inside the domain: an array forms × points."""
    U = np.zeros((len(forms), x.size))
    groups = {}
    for i, f in enumerate(forms):
        groups.setdefault((bool(f["twist"]), int(f["parity"])), []).append(i)
    for (twist, par), idx in groups.items():
        M = min(min(len(forms[i]["c"]) for i in idx), int(30.0 / (2 * math.pi * y)) + 2)
        fr = np.arange(1, M + 1) - (0.5 if twist else 0.0)
        Tm = (np.cos if par == 0 else np.sin)(2 * np.pi * np.outer(fr, x))
        for i in idx:
            amp = np.asarray(forms[i]["c"][:M]) * k_bessel(forms[i]["R"], 2 * np.pi * fr * y)
            U[i] = math.sqrt(y) * (amp @ Tm)
    return U


def character(forms: list) -> tuple:
    """The character of a product of forms: (parity, ε-signs…, twist), each multiplied up."""
    par = sum(int(f["parity"]) for f in forms) % 2
    tw = sum(bool(f["twist"]) for f in forms) % 2
    keys = sorted(k for k in forms[0] if k.startswith("eps") and k[3:].isdigit())
    return (par,) + tuple(int(np.prod([f[k] for f in forms])) for k in keys) + (tw,)


def invariant(forms: list) -> bool:
    """Whether the product of the forms is unchanged by Γ₀*(N) and the mirror: the trivial character."""
    ch = character(forms)
    return ch[0] == 0 and ch[-1] == 0 and all(e == 1 for e in ch[1:-1])


def integrate(forms: list, fun, N: int = 143, n: int = 16, xfac: float = 1.0, degree: int = 2) -> np.ndarray:
    """∫ over F* of fun(u₁, …, u_k) dμ, dμ = dx dy/y², for an integrand unchanged by Γ₀*(N) and the mirror.

    ``fun`` maps the array U (forms × points) of one row to an array (m × points), or to a vector
    of points; the result is the vector of the m integrals.  The half 0 ≤ x ≤ ½ is integrated and
    doubled, so an integrand that is odd under the mirror would come out wrong: see the module's
    note on symmetry.  ``degree`` is the degree of ``fun`` in the forms (2 for norms and overlaps,
    3 for triple products, 4 for u⁴): it sets the number of x-nodes.  The strip above y = 4 is left
    out: a cusp form is below e⁻²⁵ of its size there.
    """
    D = domain(N)
    total = None
    for y, wy in D.rows(n):
        band = degree * min(max(len(f["c"]) for f in forms), int(30.0 / (2 * math.pi * y)) + 2)
        x, w = D.nodes(y, band, xfac)
        if x.size == 0:
            continue
        vals = np.atleast_2d(fun(values_on_row(forms, x, y)))
        part = wy * (vals @ w)
        total = part if total is None else total + part
    return total


def gram(forms: list, N: int = 143, n: int = 16) -> np.ndarray:
    """The matrix of ∫ u_i u_j dμ over F*: norms on the diagonal, overlaps off it.

    Pairs from different sectors are orthogonal on the curve by symmetry and are entered as 0;
    pairs of one sector are integrated, and come out orthogonal only because their eigenvalues differ.
    """
    k = len(forms)
    pairs = [(i, j) for i in range(k) for j in range(i, k) if invariant([forms[i], forms[j]])]
    ii, jj = np.array([p[0] for p in pairs]), np.array([p[1] for p in pairs])
    v = integrate(forms, lambda U: U[ii] * U[jj], N, n, degree=2)
    G = np.zeros((k, k))
    G[ii, jj] = v
    G[jj, ii] = v
    return G


def participation_ratios(forms: list, N: int = 143, n: int = 16) -> np.ndarray:
    """I = area·∫u⁴/(∫u²)² for each form: 1 for a constant, larger the more the form is concentrated."""
    k = len(forms)
    v = integrate(forms, lambda U: np.vstack([U ** 2, U ** 4]), N, n, degree=4)
    return area(N) * v[k:] / v[:k] ** 2


def triple(f, g, h, N: int = 143, n: int = 16) -> float:
    """∫ f g h dμ over F*, divided by the three norms: the overlap of three forms.

    Zero, without integrating, when the product's character is not trivial (the selection rule).
    """
    if not invariant([f, g, h]):
        return 0.0
    v = integrate([f, g, h], lambda U: np.vstack([U[0] * U[1] * U[2], U[0] ** 2, U[1] ** 2, U[2] ** 2]), N, n, degree=3)
    return float(v[0] / math.sqrt(v[1] * v[2] * v[3]))
