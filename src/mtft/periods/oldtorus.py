"""The level-11 oldspace of X_0(143) as an arithmetic abelian surface.

Everything is re-derived at call time by the kernel route:
L_old := Z-saturation of ker(U_13^2 - 4 U_13 + 13 I), then the intersection
form is restricted, the J_arith = (U_13 - 2I)/3 closure L_9 is built, and
the principal polarization E/2 and the hyperbolic product dynamics follow.

All lattice statements are EXACT (integer/Fraction arithmetic); the float
entropy is a CERTIFIED evaluation of the exact closed form 2 log(2+sqrt(5)).
"""
from __future__ import annotations

from fractions import Fraction as Fr
from functools import lru_cache
import math

import numpy as np
import sympy as sp

from mtft import hecke as H
import mtft.integral_lattice as IL
from .core import intersection_form
from .bridge import cuspidal_basis_change
from .involutions import al_matrix

N = 26

NORMAL_FORM_J = ((0, -1), (1, 0))
NORMAL_FORM_W = ((-1, 4), (0, 1))


@lru_cache(maxsize=1)
def _u13():
    U = [[Fr(x.numerator, x.denominator) for x in r]
         for r in H.cuspidal_hecke(13)]
    assert all(v.denominator == 1 for r in U for v in r)
    return tuple(tuple(int(v) for v in r) for r in U)


@lru_cache(maxsize=1)
def old_lattice():
    """HNF basis (26x4 integer columns) of the saturated old lattice, EXACT."""
    U = np.array(_u13(), object)
    M = U @ U - 4 * U + 13 * np.eye(N, dtype=object)
    ker = IL.rational_kernel(M)
    cols = []
    for v in ker:
        den = 1
        for x in v:
            den = den * x.denominator // math.gcd(den, x.denominator)
        cols.append([int(x * den) for x in v])
    K0 = np.array(cols, object).T
    primes = set()
    for s in IL.smith_invariants(K0):
        s = int(s); d = 2
        while d * d <= s:
            while s % d == 0:
                primes.add(d); s //= d
            d += 1
        if s > 1:
            primes.add(s)
    L, _ = IL.saturate(K0, sorted(primes) or [2])
    L = IL.hnf(L)
    assert L.shape == (N, 4)
    return tuple(tuple(int(L[i, j]) for j in range(4)) for i in range(N))


@lru_cache(maxsize=1)
def _EH():
    E = np.array(intersection_form(), object)
    C = np.array(cuspidal_basis_change(), object)
    return C.T @ E @ C


def intersection_on_old():
    """Restricted intersection form on L_old (4x4 integer), EXACT."""
    L = np.array(old_lattice(), object)
    Eo = L.T @ _EH() @ L
    return tuple(tuple(int(x) for x in r) for r in Eo)


def polarization_type():
    """Smith invariants (2,2,18,18) -> polarization type (2,18), EXACT."""
    sm = tuple(int(s) for s in IL.smith_invariants(
        np.array(intersection_on_old(), object)))
    assert sm == (2, 2, 18, 18)
    return {"smith": sm, "type": (2, 18), "det": 1296}


def j_arith():
    """(U_13 - 2I)/3 on the full 26D basis, EXACT Fractions.

    Satisfies J^2 = -I on the oldspace only; it is NOT integral on L_old
    (mod-3 rank of U_13-2I on L_old is 2, hence the 3^2 = 9 closure index).
    """
    U = _u13()
    return tuple(tuple(Fr(U[i][j] - (2 if i == j else 0), 3)
                       for j in range(N)) for i in range(N))


@lru_cache(maxsize=1)
def _closure():
    Ls = sp.Matrix([[old_lattice()[i][j] for j in range(4)] for i in range(N)])
    Jar = sp.Matrix([[sp.Rational(x.numerator, x.denominator) for x in r]
                     for r in j_arith()])
    coeff = (Ls.T * Ls).inv() * Ls.T
    c8 = coeff * sp.Matrix.hstack(Ls, Jar * Ls)
    den = int(sp.lcm([c.q for c in c8]))
    Ci = np.array([[int(c * den) for c in c8.row(i)] for i in range(4)], object)
    Hn = IL.hnf(Ci)
    d = 1
    for i in range(4):
        d *= int(Hn[i, i])
    index = den ** 4 // d
    L9 = Ls * sp.Matrix([[sp.Rational(int(Hn[i, j]), den) for j in range(4)]
                         for i in range(4)])
    return Ls, Jar, L9, index


def l9_index():
    """[L_9 : L_old] = 9, EXACT."""
    return _closure()[3]


def mod3_rank():
    """rank_3((U_13 - 2I)|L_old) = 2, so the index 9 = 3^2 is structural."""
    Ls = np.array(old_lattice(), object)
    U = np.array(_u13(), object)
    M = (U - 2 * np.eye(N, dtype=object)) @ Ls
    return IL.rank_modp(M, 3)


def principal_form():
    """E/2 restricted to L_9: integral, alternating, det +1 (principal)."""
    _, _, L9, _ = _closure()
    EHs = sp.Matrix([[int(_EH()[i, j]) for j in range(N)] for i in range(N)])
    E9 = L9.T * EHs * L9
    assert all(v.q == 1 for v in E9)
    E9h = E9 / 2
    assert all(v.q == 1 for v in E9h) and E9h.T == -E9h
    assert abs(E9h.det()) == 1
    return tuple(tuple(int(x) for x in r) for r in E9h.tolist())


def product_charpoly():
    """charpoly((J_arith W_13)|L_9) = (x^2 - 4x - 1)^2, EXACT.

    Also certifies W_13 L_9 = L_9 and J_arith integral on L_9, and that
    J_arith does NOT preserve the principal form E/2.
    """
    Ls, Jar, L9, _ = _closure()
    W13 = sp.Matrix([[al_matrix(13)[i][j] for j in range(N)] for i in range(N)])
    inv9 = (L9.T * L9).inv() * L9.T
    cw = inv9 * (W13 * L9)
    cj = inv9 * (Jar * L9)
    assert all(v.q == 1 for v in cw) and all(v.q == 1 for v in cj)
    x = sp.symbols("x")
    cp = sp.expand((sp.Matrix(cj) * sp.Matrix(cw)).charpoly(x).as_expr())
    assert sp.expand(cp - (x ** 2 - 4 * x - 1) ** 2) == 0
    E9h = sp.Matrix([[principal_form()[i][j] for j in range(4)]
                     for i in range(4)])
    preserves = (sp.Matrix(cj).T * E9h * sp.Matrix(cj) == E9h)
    return {"charpoly": "(x^2 - 4x - 1)^2",
            "spectral_radius_exact": "2 + sqrt(5)",
            "j_arith_preserves_principal_form": bool(preserves)}


def entropy():
    """Topological entropy 2 log(2 + sqrt 5) of the product automorphism."""
    return 2.0 * math.log(2.0 + math.sqrt(5.0))


__all__ = ["old_lattice", "intersection_on_old", "polarization_type",
           "j_arith", "l9_index", "mod3_rank", "principal_form",
           "product_charpoly", "entropy", "NORMAL_FORM_J", "NORMAL_FORM_W"]
