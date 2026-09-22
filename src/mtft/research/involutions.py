"""Oloid-to-Hodge handoff (Astra, 2026-09-19), exact algebra: one rational involution in three settings.

w(q) = -q/(1+q) is an involution (w(w(q)) = q).  It relates the oloid contact parameters cos u = -cos t/(1 + cos t)
(Dirnboeck–Stachel eq. 6, q = cos t), it is the Fricke parameter involution of the elliptic family in the T/H studies, and
the coordinate n(q) = -1/(q+2) conjugates it to the triangle reflection n -> -1-n, which preserves T(n) = n(n+1)/2
(and (2n+1)^2 = 8T + 1).  F_r(n) = 2 r T(n) is conjugate to z -> z^2 + c with z = r(n + 1/2), c = r(2 - r)/4 (T02).
A shared involution is an algebraic statement; no lift from the oloid's rolling dynamics to the Hodge transport exists
(O03), and two solid oloids intersected in R^3 cannot form S^4 (O05).  Real-parameter statements, not lattice claims."""
import sympy as sp

q, n, r, z = sp.symbols("q n r z")

def oloid_involution(x): return -x / (1 + x)

def triangle_number(x): return x * (x + 1) / 2

def conjugating_coordinate(x): return -1 / (x + 2)

def involution_checks():
    w = oloid_involution; nq = conjugating_coordinate(q)
    return {"w_is_involution": sp.simplify(w(w(q)) - q) == 0, "conjugates_to_triangle_reflection": sp.simplify(conjugating_coordinate(w(q)) - (-1 - nq)) == 0,
            "triangle_reflection_invariant": sp.expand(triangle_number(-1 - n) - triangle_number(n)) == 0, "square_identity": sp.expand((2 * n + 1) ** 2 - (8 * triangle_number(n) + 1)) == 0}

def logistic_conjugacy():
    """F_r(n) = 2 r T(n) = r n (n+1); with z = r (n + 1/2): z -> z^2 + c, c = r(2 - r)/4."""
    F = 2 * r * triangle_number(n); zn = r * (n + sp.Rational(1, 2)); c = r * (2 - r) / 4
    lhs = r * (F + sp.Rational(1, 2)); rhs = zn ** 2 + c
    return {"c": sp.simplify(c), "conjugacy_holds": sp.simplify(sp.expand(lhs - rhs)) == 0}

def oloid_intersection_dimension():
    return {"two_solid_oloids_in_R3_give_S4": False, "reason": "a convex intersection with interior is a 3-ball with S^2 boundary; a configuration space is a separate object (O05)"}
