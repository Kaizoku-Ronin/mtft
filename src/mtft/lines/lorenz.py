r"""mtft.lines.lorenz — modular knots as Lorenz knots: words and linking numbers on the template (EXACT; LN-01, 2026-11-16).

Ghys: the closed orbits of the geodesic flow on PSL(2, ℤ)\PSL(2, ℝ), the complement of a trefoil
in the 3-sphere, are the periodic orbits of the Lorenz template.  A primitive hyperbolic matrix A
is conjugate to a positive word in U = [[1, 1], [0, 1]] and V = [[1, 0], [1, 1]], and its orbit is
the Lorenz orbit with that word read in x = U and y = V (``lorenz_word``; #x − #y is Rademacher's
Ψ).  On the template every crossing lies where the left strip passes over the right one, and two
strands cross there exactly when the one that starts on the left lands to the right of the one
that starts on the right; the linking number of two orbits is half the number of their mutual
crossings (``lorenz_linking``).  Nothing here uses geometry.

The transfer.  M_N = Γ₀(N)\PSL(2, ℝ) covers level 1 with ψ(N) sheets, so for two different motions
the total of crossings on X₀(N) between all the closed geodesics over one and all those over the
other is

    ψ(N) · Σ Lk(a, b),      a over the distinct words of A and A⁻¹, b over those of B and B⁻¹

(``transfer_sum``; a motion conjugate to its inverse has one word, counted once).  This is a second
route to the totals that ``mtft.lines.crossings`` finds by cutting geodesics, with no step shared.
For the lab's three motions the sums are 2, 2 and 6.

Provenance: the Modular Plane lab, record MP13, Route C (12 November 2026).
"""
from __future__ import annotations

from .darts import TI, S, T, as_matrix, inv, is_hyperbolic, mm, neg, trace

_U, _V = ((1, 1), (0, 1)), ((1, 0), (1, 1))


def positive_conjugate(A, depth: int = 40):
    """A conjugate of ±A with all entries ≥ 0, found breadth first over conjugation by S, T and T⁻¹."""
    A = as_matrix(A)
    if not is_hyperbolic(A):
        raise ValueError("the matrix is not hyperbolic")
    if trace(A) < 0:
        A = neg(A)
    frontier, seen = [A], {A}
    for _ in range(depth):
        for B in frontier:
            if min(B[0] + B[1]) >= 0:
                return B
        nxt = []
        for B in frontier:
            for h in (S, T, TI):
                C = mm(mm(h, B), inv(h))
                if C not in seen:
                    seen.add(C)
                    nxt.append(C)
        frontier = nxt
    raise ValueError("no conjugate with entries of one sign was found within the depth given")


def canon(w: str) -> str:
    """A cyclic word by its least rotation."""
    return min(w[i:] + w[:i] for i in range(len(w)))


def lorenz_word(A) -> str:
    """The cyclic word in x and y of the Lorenz orbit of a hyperbolic matrix, by its least rotation."""
    (a, b), (c, d) = B = positive_conjugate(A)
    w = ""
    while (a, b, c, d) != (1, 0, 0, 1):
        if a >= c and b >= d:
            a, b, w = a - c, b - d, w + "x"
        elif c >= a and d >= b:
            c, d, w = c - a, d - b, w + "y"
        else:
            raise ArithmeticError("the matrix does not peel into a positive word")
    P = ((1, 0), (0, 1))
    for ch in w:
        P = mm(P, _U if ch == "x" else _V)
    if P != B:
        raise ArithmeticError("the word does not spell the matrix")
    return canon(w)


def word_matrix(w: str):
    """The matrix a word in x and y spells."""
    P = ((1, 0), (0, 1))
    for ch in w:
        if ch not in "xy":
            raise ValueError("a Lorenz word uses the letters x and y only")
        P = mm(P, _U if ch == "x" else _V)
    return P


def _less(u: str, v: str) -> bool:
    """u^∞ < v^∞, lexicographically with x < y."""
    n = 2 * (len(u) + len(v))
    return (u * (n // len(u) + 1))[:n] < (v * (n // len(v) + 1))[:n]


def lorenz_linking(a: str, b: str) -> int:
    """The linking number of two different Lorenz orbits: half the crossings between their strands."""
    if canon(a) == canon(b):
        raise ValueError("the two words are the same orbit; its self-linking is not defined here")
    sa = [a[i:] + a[:i] for i in range(len(a))]
    sb = [b[i:] + b[:i] for i in range(len(b))]
    cr = 0
    for P, Q in ((sa, sb), (sb, sa)):
        for p in P:
            if p[0] != "x":
                continue
            for q in Q:
                if q[0] == "y" and _less(q[1:] + q[0], p[1:] + p[0]):
                    cr += 1
    if cr % 2:
        raise ArithmeticError("an odd number of crossings between two closed orbits")
    return cr // 2


def words(A) -> tuple:
    """The distinct Lorenz words of A and A⁻¹: one word if A is conjugate to its inverse, else two."""
    A = as_matrix(A)
    return tuple(sorted({lorenz_word(A), lorenz_word(inv(A))}))


def is_reciprocal(A) -> bool:
    """Whether A is conjugate to its inverse in PSL(2, ℤ) (its geodesic runs through a cone point and back)."""
    return len(words(A)) == 1


def transfer_sum(A, B) -> int:
    """Σ Lk(a, b) over the distinct words a of A^±1 and b of B^±1: the crossings of the two geodesics per sheet."""
    wa, wb = words(A), words(B)
    if set(wa) & set(wb):
        raise ValueError("the two motions give the same closed geodesic of the modular surface")
    return sum(lorenz_linking(a, b) for a in wa for b in wb)
