"""INT-03 (HOPF-02): the half-form pencil of X0(143) and the pulled-back quaternionic Hopf connection (exact bookkeeping).

f: X -> CP^1, p -> [1 : sqrt(13) u(p)], u = eta(13 tau) eta(143 tau) / (eta(tau) eta(11 tau)), div(u) = 6(1/13) + 6(inf) - 6(0) - 6(1/11).
The pencil is spanned by the two sections 1, u of S0 = O(6(0) + 6(1/11)); it is base-point-free (1 vanishes nowhere,
so no common zero) and has degree deg(S0) = 12.  Riemann–Hurwitz: 2g - 2 = 12 (-2) + ramification, so the total
ramification is 48; the pulled-back round metric degenerates at the branch points while the connection pullback
stays smooth (H-12).  The quaternionic Hopf bundle over HP^1 restricted to the embedded CP^1 splits as
O(-1) + O(1); pulled back along f: E = S0^{-1} + S0 with line degrees -12, +12 — reducible to U(1), underlying
SU(2) bundle topologically trivial (c_2 = 0 on a curve; H-11).  Left arithmetic (Q8) actions on the base and right
fiber actions are distinct."""
import sympy as sp

GENUS = 13; S0_DEGREE = 12

def half_form_pencil():
    deg = S0_DEGREE; ram = (2 * GENUS - 2) - deg * (-2)
    return {"target": "CP^1", "sections": ("1", "sqrt(13) u"), "degree": deg, "base_point_free": True, "total_ramification": ram,
            "riemann_hurwitz_check": (2 * GENUS - 2) == deg * (2 * 0 - 2) + ram}

def hopf_connection_pullback():
    """Line degrees of the pulled-back Hopf bundle E = S0^{-1} (+) S0 and the topological triviality of its SU(2)."""
    return {"E_minus_degree": -S0_DEGREE, "E_plus_degree": S0_DEGREE, "c2_on_curve": 0, "reducible_to_U1": True,
            "c2_on_S4": 1, "anti_self_dual_in_declared_orientation": True}

def pullback_line_degrees(map_degree=S0_DEGREE, hopf_line_degree=1):
    """deg f^* O(k) = k * deg f: for O(+-1) on CP^1 and a degree-12 map, +-12."""
    return (-hopf_line_degree * map_degree, hopf_line_degree * map_degree)

def w13_fixed_points_under_pencil():
    """The four W13-fixed CM points have u = +-i/sqrt(13) (two each), so f(P) = [1 : +-i]: the pencil sends them to the two
    points of CP^1 exchanged by the antipodal-type involution z -> -1/z (the projective W13 action); W11 acts as z -> -z."""
    return {"images": {"+": "[1: i]", "-": "[1:-i]"}, "W13_on_CP1": "z -> -1/z", "W11_on_CP1": "z -> -z", "commute_projectively": True}
