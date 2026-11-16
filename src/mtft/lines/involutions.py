r"""mtft.lines.involutions — the Atkin–Lehner involutions and the mirror on the classes of Γ₀(143) (EXACT; LN-01, 2026-11-16).

The Atkin–Lehner matrices

    W₁₁ = [[66, 5], [143, 11]],    W₁₃ = [[78, 7], [143, 13]],    W₁₄₃ = [[0, −1], [143, 0]]

normalise Γ₀(143), so γ ↦ W γ W⁻¹ is an automorphism of the group and an integer matrix on its
classes (``al_matrix``: row e is the class of the image of generator e).  The three are commuting
involutions of ℤ²⁹ with W₁₁W₁₃ = W₁₄₃.  On the Hecke blocks they are signs: 143a1 (+, +, +),
f2 (−, +, −), f3 (+, −, −), and W₁₁ = −1 on the level-11 oldspace.  W₁₄₃ carries the Pell unit
[[12, 1], [143, 12]] to its inverse, so it negates that class.

The fibre number is not kept.  W_Q = g · diag(Q, 1) with g in SL(2, ℤ), so W_Q γ W_Q⁻¹ is conjugate
in SL(2, ℤ) to γ_Q = [[a, Q b], [c/Q, d]], and Ψ is a class function there.  Hence

    Ψ(W_Q γ W_Q⁻¹) = Ψ(γ) + Φ(γ_Q) − Φ(γ) = Ψ(γ) + the winding of (η(Qτ)/η(τ))²⁴ along γ:

    W₁₁: Ψ − 2φ_A,      W₁₃: Ψ − 12φ_B,      W₁₄₃: Ψ + 24φ₁ − 2φ_A − 12φ_B        (``psi_shift``).

On the thirty integers an involution is therefore (Ψ; cls) ↦ (Ψ + shift(cls); cls · W), and the
shift changes sign under the involution itself.  For the Pell unit, −21 + 0 + 30 + 12 = +21.

The mirror τ ↦ −τ̄ acts by γ ↦ [[a, −b], [−c, d]].  It is an involution of ℤ²⁹ that reverses Ψ and
the three windings, and on the block of 143a1 it is (n₁, n₂) ↦ (−n₁, n₁ + n₂): complex conjugation
of the period, because Re Ω₁ = Ω₂/2.

A recorded miss.  The lab conjectured, before its run and with low confidence (record MP12, P13),
that every closed line with no 143a1 component is negated in H₁ of the closed curve by one of the
three involutions.  It is false: ten of the 28 lines have no 143a1 component; W₁₄₃ negates eight of
them (the four systoles and the four Pell lines) and no involution negates the other two, which
are two of Sphereland 10's four lines of length 30.  ``conjecture_p13`` computes that; those two
zeros have no mechanism yet.
"""
from __future__ import annotations

from functools import lru_cache

from . import closed as _c
from . import darts as _d
from . import periods as _p
from . import rademacher as _r
from . import table as _t

WQ = {11: ((66, 5), (143, 11)), 13: ((78, 7), (143, 13)), 143: ((0, -1), (143, 0))}
"""Atkin–Lehner matrices of level 143 (determinants 11, 13, 143)."""


def al_conjugate(g, Q: int):
    """W_Q g W_Q⁻¹ for g in Γ₀(143), an element of Γ₀(143)."""
    W = WQ[Q]
    adj = ((W[1][1], -W[0][1]), (-W[1][0], W[0][0]))
    h = _d.mm(_d.mm(W, g), adj)
    if any(v % Q for row in h for v in row):
        raise ArithmeticError("the conjugate is not integral: the matrix is not in the group of level 143")
    return tuple(tuple(v // Q for v in row) for row in h)


PSI_SHIFT = {11: (0, -2, 0), 13: (0, 0, -12), 143: (24, -2, -12)}
"""Q → the coefficients of (φ₁, φ_A, φ_B) in Ψ(W_Q γ W_Q⁻¹) − Ψ(γ)."""


def psi_shift(g, Q: int) -> int:
    """Ψ(W_Q g W_Q⁻¹) − Ψ(g), from the windings of g: the winding of (η(Qτ)/η(τ))²⁴ along g."""
    w = _r.windings(g)
    return sum(k * v for k, v in zip(PSI_SHIFT[Q], w))


def mirror_conjugate(g):
    """[[a, −b], [−c, d]]: the element the mirror τ ↦ −τ̄ carries g to."""
    return ((g[0][0], -g[0][1]), (-g[1][0], g[1][1]))


@lru_cache(maxsize=4)
def al_matrix(Q: int) -> tuple:
    """29 × 29 integers: row e is the class of W_Q γ_e W_Q⁻¹.  The image of a class c is c · (this matrix)."""
    G = _d.graph(143)
    return tuple(G.cls(al_conjugate(g, Q)) for g in G.generators)


@lru_cache(maxsize=1)
def mirror_class_matrix() -> tuple:
    """29 × 29 integers: row e is the class of the mirror image of generator e."""
    G = _d.graph(143)
    return tuple(G.cls(mirror_conjugate(g)) for g in G.generators)


def apply(cls, M) -> tuple:
    """The image of a class under a matrix whose row e is the image of generator e."""
    n = len(M)
    return tuple(sum(cls[e] * M[e][j] for e in range(n)) for j in range(n))


def _matmul(A, B) -> tuple:
    n = len(A)
    return tuple(tuple(sum(A[i][k] * B[k][j] for k in range(n)) for j in range(n)) for i in range(n))


def group_law() -> dict:
    """The three matrices are involutions, commute, and W₁₁W₁₃ = W₁₄₃; the mirror is an involution."""
    A11, A13, A143, Mi = al_matrix(11), al_matrix(13), al_matrix(143), mirror_class_matrix()
    eye = tuple(tuple(int(i == j) for j in range(29)) for i in range(29))
    return {"involutions": all(_matmul(A, A) == eye for A in (A11, A13, A143)),
            "commute": _matmul(A11, A13) == _matmul(A13, A11),
            "product": _matmul(A11, A13) == A143,
            "mirror_involution": _matmul(Mi, Mi) == eye}


@lru_cache(maxsize=1)
def block_signs() -> dict:
    """Block name → (sign of W₁₁, W₁₃, W₁₄₃): +1, −1, or 0 where the involution is not a sign there."""
    C = _t.coord()
    out = {}
    for name in _t.BLOCKS:
        dm = _t.BLOCK_DIMS[name]
        sg = []
        for Q in (11, 13, 143):
            A = al_matrix(Q)
            img = tuple(tuple(sum(A[e][k] * C[name][k][i] for k in range(29)) for i in range(dm)) for e in range(29))
            sg.append(1 if img == C[name] else (-1 if img == tuple(tuple(-v for v in row) for row in C[name]) else 0))
        out[name] = tuple(sg)
    return out


def line_signs(cls) -> dict:
    """Q → +1, −1 or 0: whether W_Q keeps, negates or does neither to the class in H₁ of the closed curve (its 26 periods)."""
    me = _p.flat(_p.coordinates(cls))
    out = {}
    for Q in (11, 13, 143):
        im = _p.flat(_p.coordinates(apply(cls, al_matrix(Q))))
        out[Q] = 1 if im == me else (-1 if im == tuple(-v for v in me) else 0)
    return out


def mirror_sign(cls) -> int:
    """+1, −1 or 0: whether the mirror keeps, negates or does neither to a class of ℤ²⁹."""
    im = apply(cls, mirror_class_matrix())
    return 1 if im == tuple(cls) else (-1 if im == tuple(-v for v in cls) else 0)


def mirror_laws(samples) -> dict:
    """On a list of elements of Γ₀(143): the mirror reverses Ψ and the windings and acts on 143a1 by (n₁, n₂) ↦ (−n₁, n₁ + n₂)."""
    G = _d.graph(143)
    psi = wind = ell = lin = True
    Mi = mirror_class_matrix()
    for g in samples:
        h = mirror_conjugate(g)
        psi &= _r.Psi(h) == -_r.Psi(g)
        wind &= _r.windings(h) == tuple(-v for v in _r.windings(g))
        c, ch = G.cls(g), G.cls(h)
        lin &= ch == apply(c, Mi)
        n1, n2 = _p.coordinates(c)["ell"]
        ell &= _p.coordinates(ch)["ell"] == (-n1, n1 + n2)
    return {"psi_reversed": psi, "windings_reversed": wind, "linear": lin, "ell_conjugated": ell}


@lru_cache(maxsize=1)
def conjecture_p13() -> dict:
    """The lab's missed conjecture on the 28 lines: which lines have no 143a1 component, and which involution negates each."""
    L = _c.lines_143()
    zero = [n for n, l in enumerate(L, start=1) if not any(_p.coordinates(l.cls)["ell"])]
    signs = {n: line_signs(L[n - 1].cls) for n in zero}
    negated = [n for n in zero if -1 in signs[n].values()]
    return {"no_143a1_component": zero, "signs": signs, "negated_by_some": negated,
            "negated_by_none": [n for n in zero if n not in negated],
            "negated_by_W143": [n for n in zero if signs[n][143] == -1],
            "holds": len(negated) == len(zero),
            "converse_holds": all(not any(_p.coordinates(l.cls)["ell"]) for l in L if -1 in line_signs(l.cls).values())}


def al_image(line, Q: int) -> dict:
    """Where W_Q sends a closed line: the conjugate element, its class, and the lines of the 28 with ± that class in H₁ of the closed curve."""
    h = al_conjugate(line.gamma, Q)
    G = _d.graph(143)
    c = G.cls(h)
    me = _p.flat(_p.coordinates(c))
    hits = []
    for n, l in enumerate(_c.lines_143()):
        if l.motion != line.motion or l.ell != line.ell:
            continue
        other = _p.flat(_p.coordinates(l.cls))
        if other == me:
            hits.append((1, n))
        elif other == tuple(-v for v in me):
            hits.append((-1, n))
    return {"gamma": h, "cls": c, "psi": _r.Psi(h), "trace": abs(_d.trace(h)), "hits": hits}
