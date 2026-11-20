r"""mtft.lines.intersection — the intersection form of X₀(143) on its Hecke blocks, deg φ = 4 and the flux-free lattice.

EXACT; IF-01, 2026-11-20.

``periods`` gives every class 26 integer coordinates on the Hecke blocks (143a1, the level-11 oldspace, f2, f3:
2 + 4 + 8 + 12).  ``crossings`` counts, with signs, how two closed geodesics cross: their intersection number.  So the
intersection form J of the closed curve is a 26 × 26 rational matrix in block coordinates, with

    ι(γ, δ) = b(γ)ᵀ J b(δ).

``intersection_form`` solves J from 26 geodesics with independent coordinates and checks it on every pair of the
geodesics it computed.  What comes out (FF-19, 1,426 geodesics and 2,033,476 pairs; here from fewer):

  * J is block-diagonal: the Hecke blocks are orthogonal;
  * on the 143a1 block J = [[0, 1/4], [−1/4, 0]], so deg φ = |(J⁻¹)₁₂| = 4 (``modular_degree``), and the γ-flux,
    with γ² = −2 deg φ = −8, has ∫ tr F∧F ≠ 0 over X × E: no 5-manifold bounded by X × E carries it (FF-19);
  * on the lattice Λ that H₁(X₀(143), ℤ) makes in block coordinates (index 1,327,104 = 2¹⁴3⁴, LN-3) the form is
    unimodular (``homology_lattice``);
  * the flux-free classes, ker φ_* on Λ (φ_* = the 143a1 coordinates), are a rank-24 lattice whose form has
    determinant 16 = (deg φ)² and rank 22 mod 2 (``flux_free_lattice``).  Its isotropic sublattices have rank at most
    12, so no handlebody of the curve, which shrinks 13 meridians, can carry the flux.  Mod 2 the obstruction is
    invisible, because deg φ is even.

``cup_against_crossings`` ties this to ``twists``: the signed crossings mod 2 are the cup form's intersection of the classes'
26 bits, pair by pair.  ``pullbacks_on_flux_free`` checks that 143a1's twists vanish on every flux-free class.

Records: FF-19 in the MTFT Python project under claude/experiments/.  (Planned as ``mtft.lines.forms``; that name is taken twice,
by ``mtft.forms`` and ``mtft.periods.forms``.)
"""
from __future__ import annotations

from functools import lru_cache

import numpy as np
import sympy as sp

from . import periods as _per
from . import twists as _tw
from . import walls as _wl

DIM = 26
BLOCKS = ((0, 2), (2, 6), (6, 14), (14, 26))          # 143a1, the oldspace, f2, f3


def block_coordinates(cls) -> list:
    """The 26 block coordinates of a class of ℤ²⁹ (143a1 first)."""
    return list(_per.flat(_per.coordinates(tuple(cls))))


@lru_cache(maxsize=4)
def geodesics(tmax: int) -> dict:
    """Closed geodesics over the level-1 classes with trace ≤ tmax: block coordinates, lengths and signed crossings."""
    from .crossings import Crossings
    motions = {f"m{k}_t{A[0][0] + A[1][1]}": A for k, (blocks, A) in enumerate(_wl.level1_classes(tmax))}
    cr = Crossings(143, motions)
    lines = [cr.lines[g["rep"]] for g in cr.geodesics]
    return {"block coordinates": tuple(tuple(block_coordinates(l.cls)) for l in lines), "lengths": tuple(l.length for l in lines),
            "classes": tuple(tuple(l.cls) for l in lines), "iota": tuple(tuple(int(v) for v in row) for row in cr.iota)}


def intersection_form(tmax: int = 24) -> dict:
    """J on the block coordinates, solved from 26 geodesics with independent coordinates, checked on every pair."""
    g = geodesics(tmax)
    bc, iota = g["block coordinates"], g["iota"]
    pick, M = [], sp.zeros(0, DIM)
    for i, b in enumerate(bc):
        trial = M.col_join(sp.Matrix([list(b)]))
        if trial.rank() > M.rows:
            M, pick = trial, pick + [i]
            if len(pick) == DIM:
                break
    if len(pick) < DIM:
        raise ArithmeticError(f"the geodesics to trace {tmax} span only {len(pick)} of the 26 block coordinates")
    Minv = M.inv()
    J = Minv * sp.Matrix([[iota[a][b] for b in pick] for a in pick]) * Minv.T
    D = int(sp.ilcm(*[sp.fraction(x)[1] for x in J]))
    JD = np.array([[int(v) for v in row] for row in (J * D).tolist()], dtype=np.int64)
    B = np.array(bc, dtype=np.int64)
    if int(np.abs(B).max()) ** 2 * int(np.abs(JD).max()) * DIM * DIM >= 2 ** 62:
        raise OverflowError("block coordinates too large for the int64 check")
    pred = B @ JD @ B.T
    disagree = int(np.count_nonzero(pred != D * np.array(iota, dtype=np.int64)))
    off = max((abs(J[i, j]) for i in range(DIM) for j in range(DIM)
               if not any(a <= i < b and a <= j < b for a, b in BLOCKS)), default=0)
    return {"J": J, "geodesics": len(bc), "pairs checked": len(bc) ** 2, "pairs that disagree": disagree,
            "largest entry off the Hecke blocks": off, "common denominator": D}


def modular_degree(J) -> int:
    """|(J⁻¹)₁₂| on the 143a1 block: the degree of the modular parametrization X₀(143) → 143a1."""
    return int(abs(J.inv()[0, 1]))


@lru_cache(maxsize=1)
def homology_lattice() -> dict:
    """Λ = H₁(X₀(143), ℤ) in block coordinates (an HNF basis from the 29 generators) and its index in ℤ²⁶."""
    import mtft.integral_lattice as IL
    gens = [_per.flat(_per.coordinates(tuple(1 if k == e else 0 for k in range(29)))) for e in range(29)]
    H = sp.Matrix(IL.hnf(np.array(gens, dtype=object).T).tolist())
    return {"basis": H, "index": int(abs(H.det()))}


def flux_free_lattice(J) -> dict:
    """ker φ_* on Λ (φ_* = the 143a1 coordinates): its rank, the determinant of its form, its rank mod 2."""
    import mtft.integral_lattice as IL
    H = homology_lattice()["basis"]
    JL = H.T * J * H
    P = H[0:2, :]
    K = sp.Matrix(IL.int_kernel(np.array(P.tolist(), dtype=object)).tolist())
    JK = K.T * JL * K
    rows = [sum((int(JK[i, j]) % 2) << j for j in range(JK.cols)) for i in range(JK.rows)]
    return {"form on Λ integral": all(x.is_integer for x in JL), "det on Λ": int(JL.det()),
            "Smith invariants of φ_*": [int(x) for x in IL.smith_invariants(np.array(P.tolist(), dtype=object))],
            "rank": K.cols, "det": int(JK.det()), "rank over Q": JK.rank(), "rank mod 2": _tw.rank_f2(rows),
            "largest isotropic rank": JK.rank() // 2}                 # a nondegenerate skew form: its Lagrangians have half its rank


def cup_against_crossings(tmax: int = 24) -> dict:
    """Every pair of geodesics: the signed crossings mod 2 against the cup form's intersection of their class bits (C⁻¹)."""
    g = geodesics(tmax)
    Y = [_tw.class_bits(c) for c in g["classes"]]
    D = [_tw.apply_rows(_tw.class_form(), y) for y in Y]
    n = len(Y)
    odd = sum(1 for a in range(n) for b in range(n) if g["iota"][a][b] & 1)
    disagree = sum(1 for a in range(n) for b in range(n) if _tw.parity(D[a] & Y[b]) != (g["iota"][a][b] & 1))
    return {"geodesics": n, "pairs": n * n, "pairs with odd crossings": odd, "pairs that disagree mod 2": disagree}


def pullbacks_on_flux_free() -> dict:
    """143a1's twists on a basis of the flux-free classes of ℤ²⁹ (the kernel of the 143a1 coordinates): all zero."""
    import mtft.integral_lattice as IL
    K = np.array(IL.int_kernel(np.array(_tw.ell_matrix(), dtype=object)), dtype=object)
    cols = [tuple(int(K[i, j]) for i in range(K.shape[0])) for j in range(K.shape[1])]
    if any(sum(a * b for a, b in zip(row, c)) for row in _tw.ell_matrix() for c in cols):
        raise ArithmeticError("the kernel basis is not in the kernel")
    nonzero = sum(1 for c in cols for w in _tw.W()[1:] if _tw.parity(w & _tw.class_bits(c)))
    return {"kernel rank in Z^29": len(cols), "pullback values that are not 0": nonzero}
