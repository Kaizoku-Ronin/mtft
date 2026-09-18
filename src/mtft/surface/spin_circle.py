"""INT-02 (SC7-01, Astra 2026-09-16): the spin circle bundle P = S(S0) -> X0(143) and rank-one local systems on it.

Exact engine: the Seifert presentation of pi_1(P) for the unit circle bundle of Euler degree e over a genus-g surface,
the Fox-derivative cochain complex of a rank-one character, and Poincare duality with the inverse character for the
top degrees.  Theorem (fiber acyclicity): if the fiber holonomy lambda != 1 the fiber differential is multiplication
by the unit lambda - 1 and H^*(P; L) = 0.  For the degree-3 and degree-6 M1 sectors on e = 12, lambda = zeta_12^{-m}
!= 1, so the specified commuting twisted 7D complex has NO massless bulk modes for those sectors (a failed physical
gate of that construction, with the exact computation passing).  Pure functions; nothing is written to disk."""
from functools import lru_cache
import sympy as sp

def _clean(x): return sp.nsimplify(sp.simplify(x)) if isinstance(x, sp.Basic) else sp.nsimplify(x)

def commutator_word(a, b): return [a, b, -a, -b]

def seifert_relators(genus, euler):
    """pi_1(P) = < a_1, b_1, ..., a_g, b_g, t | [a_i, t] = [b_i, t] = 1, prod [a_i, b_i] = t^e >; generators one-based,
    negative index = inverse; the fiber generator is number 2g + 1."""
    fiber = 2 * genus + 1; central = [commutator_word(j, fiber) for j in range(1, fiber)]; surface = []
    for i in range(genus): surface += commutator_word(2 * i + 1, 2 * i + 2)
    surface += [-fiber] * euler
    return central + [surface]

def fox_row(word, character):
    """Left Fox derivative of a relator word evaluated in a rank-one character (tuple of unit scalars per generator)."""
    row = [sp.S.Zero] * len(character); prefix = sp.S.One
    for letter in word:
        j = abs(letter) - 1
        if letter > 0: row[j] += prefix; prefix = _clean(prefix * character[j])
        else: prefix = _clean(prefix / character[j]); row[j] -= prefix
    return [_clean(x) for x in row], _clean(prefix)

@lru_cache(maxsize=None)
def low_cohomology(genus, euler, character):
    n = 2 * genus + 1; d0 = sp.Matrix([_clean(x - 1) for x in character]); rows = []
    for rel in seifert_relators(genus, euler):
        row, value = fox_row(rel, character)
        if value != 1: raise ValueError("character does not define a representation of pi_1")
        rows.append(row)
    d1 = sp.Matrix(rows)
    if (d1 * d0).applyfunc(_clean) != sp.zeros(n, 1): raise ArithmeticError("d1 d0 != 0")
    r0, r1 = d0.rank(), d1.rank(); return {"rank_d0": r0, "rank_d1": r1, "h0": 1 - r0, "h1": n - r0 - r1}

def rank_one_cohomology(genus, euler, character):
    """H^0..H^3 dimensions of the rank-one local system on the closed oriented 3-manifold P; H^2, H^3 by Poincare duality
    with the INVERSE character (H-07/INT-02 duality rule)."""
    character = tuple(_clean(x) for x in character); low = low_cohomology(genus, euler, character)
    dual = low_cohomology(genus, euler, tuple(_clean(1 / x) for x in character))
    dims = [low["h0"], low["h1"], dual["h1"], dual["h0"]]
    if sum((-1) ** i * v for i, v in enumerate(dims)) != 0: raise ArithmeticError("Euler characteristic of a closed 3-manifold must vanish")
    return {**low, "dual_rank_d0": dual["rank_d0"], "dual_rank_d1": dual["rank_d1"], "cohomology_dimensions": dims}

def circle_bundle_topology(genus=13, euler=12):
    """H_1(P; Z) = Z^{2g} + Z/e, H^2(P; Z) = Z^{2g} + Z/e; a degree-m line bundle pulls back to torsion of order e/gcd(e, m)."""
    from math import gcd
    return {"H1": ("Z^%d" % (2 * genus), "Z/%d" % euler), "H2": ("Z^%d" % (2 * genus), "Z/%d" % euler),
            "pullback_torsion_order": lambda m: euler // gcd(euler, m % euler if m % euler else euler)}

def fiber_holonomy(degree, euler=12):
    """Holonomy of the flat pullback connection A_hat_m = pi^* A_m - (m/e) Theta around the fiber: exp(-2 pi i m/e)."""
    return _clean(sp.exp(-2 * sp.pi * sp.I * sp.Rational(degree, euler)))

def m1_sector_cohomology(degree, euler=12, genus=13, base=None):
    """Cohomology of the flat pullback of a degree-`degree` base bundle (with optional base character, default trivial)."""
    base = list(base) if base is not None else [sp.S.One] * (2 * genus)
    return rank_one_cohomology(genus, euler, [*base, fiber_holonomy(degree, euler)])

def spin_circle_metric(R, r, euler=12, genus=13):
    """Connection metric ds^2 = ds_X^2 + r^2 Theta^2 on P (base curvature -1/R^2, compact area 2 pi (2g-2) R^2):
    connection curvature F = f dvol with f = 2 pi e / area = e/((2g-2) R^2); Kaluza–Klein: R_P = R_X - (r^2/2) f^2 =
    -2/R^2 - r^2/(8 R^4) for e = 12, g = 13; volume 2 pi r * area = 96 pi^2 R^2 r.
    Ricci eigenvalues are unequal (not Einstein, not a hyperbolic 3-manifold; H-09)."""
    R, r = sp.nsimplify(R), sp.nsimplify(r); area = 2 * sp.pi * (2 * genus - 2) * R ** 2; c = sp.Rational(euler, 2 * genus - 2)
    return {"scalar_curvature": sp.simplify(-2 / R ** 2 - r ** 2 * c ** 2 / (2 * R ** 4)), "volume": sp.simplify(2 * sp.pi * r * area), "area_base": area, "einstein": False}

def radius_potential_7d(a, b, c):
    """Restricted Einstein-frame radius potential of the 7D reduction: V(rho) = a rho^-5 + b rho^-7 + c rho^-3 with a > 0,
    b, c >= 0 (negative curvature, flux, nonnegative bulk): every term is decreasing, so dV/d rho < 0 for all rho > 0 and
    there is no stationary point (SC7-01 §8).  Returns V and dV as functions of the symbol rho."""
    rho = sp.Symbol("rho", positive=True); V = a * rho ** -5 + b * rho ** -7 + c * rho ** -3
    return {"rho": rho, "V": sp.simplify(V), "dV": sp.simplify(sp.diff(V, rho)), "stationary_point_exists": False if (a > 0 and b >= 0 and c >= 0) else None}
