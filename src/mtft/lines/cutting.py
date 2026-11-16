r"""mtft.lines.cutting — a closed geodesic of the modular surface cut along the fundamental domain, exactly (LN-01, 2026-11-16).

F is the standard domain of PSL(2, ℤ): |x| ≤ ½, |z| ≥ 1.  The axis of a hyperbolic matrix A, seen
in the frame z′ = m⁻¹z, is the axis of m⁻¹Am: a semicircle

    a(x² + y²) + bx + c = 0,     (a, b, c) = (lower left, lower right − upper left, −upper right),

with integer a, b, c.  It meets F in one arc whose ends have rational x and rational y², so the
whole cutting sequence is exact: ``walk`` follows the axis from frame to frame (across the right
edge by T, the left edge by T⁻¹, the bottom arc by S, and through a corner ρ or ρ + 1 into the
opposite sector) until it returns to the first frame, one period later.  Only the choice of the
first frame uses floating point (a generic point of the axis, reduced with 60 digits); every arc
is then exact, and the walk must close with the product of its steps equal to ±A.

``Cutting`` holds the arcs of one period with their frames.  Two semicircles cross at a point with
rational x and y² (``cross``), which is how ``mtft.lines.crossings`` counts every crossing of two
closed geodesics on the curve.

``mirror_relation`` finds whether the mirror image τ ↦ −τ̄ of the geodesic of A, which is the
geodesic of [[a, −b], [−c, d]], is the same closed geodesic of the modular surface, and in which
direction.

Provenance: the Modular Plane lab, records MP13 and MP14 (12 November 2026).
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from fractions import Fraction as Fr
from functools import lru_cache

import mpmath as mp

from .darts import I2, R2, TI, R, S, T, as_matrix, inv, is_hyperbolic, mm, mpow, neg, trace

HALF, Q34 = Fr(1, 2), Fr(3, 4)
_RHO, _RHO1 = complex(-0.5, math.sqrt(3) / 2), complex(0.5, math.sqrt(3) / 2)
_CORNER = {"rho": _RHO, "rho1": _RHO1}
_CAND = {"rho1": (T, mm(T, R), mm(T, R2), mm(mm(T, R), TI), mm(mm(T, R2), TI)),      # the other sectors round rho + 1
         "rho": (R, R2, TI, mm(R, TI), mm(R2, TI))}                                  # and round rho
_START_ANGLES = ("0.3141", "0.2718", "0.5772", "0.4142", "0.6931")


def pos(A):
    """The sign of ±A with positive trace (hyperbolic matrices only)."""
    return A if trace(A) > 0 else neg(A)


def same(A, B) -> bool:
    return A == B or A == neg(B)


def conj(g, A):
    """g⁻¹ A g with positive trace: the geodesic of A seen in the frame z′ = g⁻¹z."""
    return pos(mm(mm(inv(g), A), g))


def mirror_matrix(A):
    """[[a, −b], [−c, d]]: the matrix of the mirror image τ ↦ −τ̄ of the geodesic of A."""
    return ((A[0][0], -A[0][1]), (-A[1][0], A[1][1]))


def form(A) -> tuple:
    """(a, b, c) of the semicircle a(x² + y²) + bx + c = 0 that is the axis of A."""
    return (A[1][0], A[1][1] - A[0][0], -A[0][1])


def unoriented(f) -> tuple:
    """The form of a semicircle with a > 0, whatever its direction."""
    return f if f[0] > 0 else tuple(-v for v in f)


def kind(p) -> str:
    """Where a point (x, y²) of F lies: 'int', 'left', 'right', 'bot-', 'bot+', 'i', 'rho' or 'rho1'."""
    x, y2 = p
    on = x * x + y2 == 1
    if x == -HALF and on:
        return "rho"
    if x == HALF and on:
        return "rho1"
    if x == -HALF:
        return "left"
    if x == HALF:
        return "right"
    if on:
        return "i" if x == 0 else ("bot-" if x < 0 else "bot+")
    return "int"


@dataclass(frozen=True)
class Arc:
    """One arc of a closed geodesic inside F, in the frame ``m``."""
    N: tuple          # the matrix whose axis carries the arc: m⁻¹ A m with positive trace
    form: tuple       # (a, b, c)
    s: int            # +1 if the arc runs from left to right, −1 otherwise
    center: Fr        # −b/(2a)
    entry: tuple      # (x, y²), exact
    exit: tuple
    m: tuple          # the frame

    @property
    def length(self) -> float:
        """The hyperbolic length of the arc (60 digits inside, a float out)."""
        a, b, _ = self.form
        with mp.workdps(60):
            out = []
            for x, y2 in (self.entry, self.exit):
                out.append(mp.atan2(mp.sqrt(mp.mpf(y2.numerator) / y2.denominator),
                                    mp.mpf(x.numerator) / x.denominator - mp.mpf(-b) / (2 * a)))
            return float(abs(mp.log(mp.tan(out[1] / 2) / mp.tan(out[0] / 2))))


def _arc(Nm) -> dict | None:
    """The arc of the oriented axis of Nm inside F: exact ends (x, y²), or None if it has no length."""
    a, b, c = form(Nm)
    pts = set()
    yl = Fr(b, 2 * a) - Fr(c, a) - Fr(1, 4)                             # on x = -1/2
    if yl >= Q34:
        pts.add((-HALF, yl))
    yr = Fr(-b, 2 * a) - Fr(c, a) - Fr(1, 4)                            # on x = 1/2
    if yr >= Q34:
        pts.add((HALF, yr))
    if b != 0:                                                          # on |z| = 1
        xb = Fr(-(a + c), b)
        if abs(xb) <= HALF:
            pts.add((xb, 1 - xb * xb))
    pts = sorted(pts)
    if len(pts) != 2:
        return None
    s = 1 if a > 0 else -1                                              # the attracting end is to the right when a > 0
    en, ex = (pts[0], pts[1]) if s > 0 else (pts[1], pts[0])
    return {"N": Nm, "form": (a, b, c), "s": s, "center": Fr(-b, 2 * a), "entry": en, "exit": ex}


def _mob(g, z):
    return (g[0][0] * z + g[0][1]) / (g[1][0] * z + g[1][1])


def _step(ar: dict):
    """The step that carries the frame across the arc's exit, and the next arc."""
    k = kind(ar["exit"])
    x, y2 = ar["exit"]
    if k in ("right", "left", "i", "bot-", "bot+"):
        g = {"right": T, "left": TI}.get(k, S)
        nxt = _arc(conj(g, ar["N"]))
        if nxt is None or nxt["entry"] != (-x, y2):                     # T, T^-1 and S all send the exit to the next entry this way
            raise ArithmeticError("the cutting sequence broke at an edge")
        return g, nxt
    good = []                                                           # through a corner: into the opposite sector
    for g in _CAND[k]:
        nxt = _arc(conj(g, ar["N"]))
        if nxt is None or kind(nxt["entry"]) not in ("rho", "rho1"):
            continue
        if abs(_mob(g, _CORNER[kind(nxt["entry"])]) - _CORNER[k]) < 1e-9:
            good.append((g, nxt))
    if len(good) != 1:
        raise ArithmeticError("the cutting sequence broke at a corner")
    return good[0]


def _reduce_point(z):
    """m and z′ in F with z = m(z′), for a point z given with many digits."""
    m = I2
    for _ in range(100000):
        n = int(mp.nint(z.real))
        if n:
            z, m = z - n, mm(m, mpow(T, n))
        if abs(z) < 1:
            z, m = -1 / z, mm(m, S)
        else:
            return m, z
    raise ArithmeticError("the reduction of the starting point did not end")


def _first_arc(A):
    a, b, c = form(A)
    D = b * b - 4 * a * c
    with mp.workdps(60):
        xc, r = mp.mpf(-b) / (2 * a), mp.sqrt(D) / (2 * abs(a))
        for ang in _START_ANGLES:                                       # a generic point of the axis
            m0, _ = _reduce_point(xc + r * mp.expjpi(mp.mpf(ang)))
            ar = _arc(conj(m0, A))
            if ar is not None:
                return m0, ar
    raise ArithmeticError("no starting frame was found for the cutting sequence")


@dataclass(frozen=True)
class Cutting:
    """The cutting sequence of one period of the axis of a primitive hyperbolic matrix."""
    A: tuple
    arcs: tuple
    exits: str            # one letter per arc: R, L (edges), B (bottom arc), i, p (ρ), q (ρ + 1)

    @property
    def through_i(self) -> int:
        return self.exits.count("i")

    @property
    def through_rho(self) -> int:
        return self.exits.count("p") + self.exits.count("q")

    @property
    def length(self) -> float:
        return sum(ar.length for ar in self.arcs)


_LETTER = {"right": "R", "left": "L", "i": "i", "bot-": "B", "bot+": "B", "rho": "p", "rho1": "q"}


@lru_cache(maxsize=64)
def _walk(A) -> Cutting:
    m0, ar = _first_arc(A)
    arcs, m, first = [], m0, ar["N"]
    for _ in range(100000):
        arcs.append(Arc(ar["N"], ar["form"], ar["s"], ar["center"], ar["entry"], ar["exit"], m))
        g, ar = _step(ar)
        m = mm(m, g)
        if ar["N"] == first:
            break
    else:
        raise ArithmeticError("the cutting sequence did not close")
    if not same(mm(m, inv(m0)), A):
        raise ValueError("the walk closed before one period of the matrix: the matrix is a proper power "
                         "(only primitive hyperbolic matrices are cut here)")
    return Cutting(A, tuple(arcs), "".join(_LETTER[kind(a.exit)] for a in arcs))


def walk(A) -> Cutting:
    """The exact cutting sequence of the closed geodesic of a primitive hyperbolic matrix A of SL(2, ℤ)."""
    A = as_matrix(A)
    if not is_hyperbolic(A):
        raise ValueError("the matrix is not hyperbolic")
    return _walk(pos(A))


def cross(fa, fb):
    """The point (x, y²) of F where two semicircles cross, or None if they do not cross inside F (closed)."""
    a1, b1, c1 = fa
    a2, b2, c2 = fb
    den = a2 * b1 - a1 * b2
    if den == 0:
        return None
    x = Fr(a1 * c2 - a2 * c1, den)
    y2 = -(a1 * x * x + b1 * x + c1) / Fr(a1)
    if y2 <= 0 or abs(x) > HALF or x * x + y2 < 1:
        return None
    return (x, y2)


@lru_cache(maxsize=64)
def _mirror_relation(A):
    cut, cm = _walk(A), _walk(pos(mirror_matrix(A)))
    for ar_m in cm.arcs:
        for ar in cut.arcs:
            for e, Nm in ((1, ar.N), (-1, pos(inv(ar.N)))):
                if ar_m.N == Nm:
                    h = mm(ar_m.m, inv(ar.m))
                    if not same(mm(mm(h, A if e == 1 else inv(A)), inv(h)), pos(mirror_matrix(A))):
                        raise ArithmeticError("the mirror relation did not check")
                    return e, h
    return None


def mirror_relation(A):
    """(e, h) with h A^e h⁻¹ = ± the mirror matrix of A, e = ±1, or None.

    e = +1: the mirror image of the closed geodesic of A on the modular surface is itself, with its
    direction; e = −1: itself, reversed; None: another closed geodesic.
    """
    A = as_matrix(A)
    if not is_hyperbolic(A):
        raise ValueError("the matrix is not hyperbolic")
    return _mirror_relation(pos(A))
