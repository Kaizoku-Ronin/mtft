"""INT-11 (initial audit §7, SC7-01 §8, AXG-01/02 vacuum, AXG-04 §8): Einstein-frame radius potentials and the unwarped
product background of a 6D theory on X0(143)."""
import sympy as sp

def einstein_frame_potential(a, b, c, dim=6):
    """Restricted radius potential: 6D exponents (4, 6, 2), 7D (5, 7, 3): V = a rho^-p1 + b rho^-p2 + c rho^-p3 with a > 0
    (negative curvature), b >= 0 (flux), c >= 0 (nonnegative bulk): no stationary point.  Negative bulk terms can give
    conditional AdS minima."""
    p = {6: (4, 6, 2), 7: (5, 7, 3)}[dim]; rho = sp.Symbol("rho", positive=True); V = a * rho ** -p[0] + b * rho ** -p[1] + c * rho ** -p[2]
    return {"rho": rho, "V": V, "dV": sp.simplify(sp.diff(V, rho)), "stationary_point_exists": False if (a > 0 and b >= 0 and c >= 0) else None}

def product_background(M, U, rho_F):
    """Unwarped constant-field Einstein equations, M = M6^4, U = Lambda6 + V6(H0), rho_F > 0 flux energy: external
    R_mn = lambda4 g, internal R_mn = k g with lambda4 = (U - rho_F)/(2M), k = (U + 3 rho_F)/(2M)."""
    lam = (U - rho_F) / (2 * M); k = (U + 3 * rho_F) / (2 * M); return {"lambda4": sp.simplify(lam), "k": sp.simplify(k)}

def c3x_ads_control(M, R, rho_F):
    """Impose k = -1/R^2: U* = -2M/R^2 - 3 rho_F, lambda4 = -1/R^2 - 2 rho_F/M < 0 (AdS), ell4^2 = 3/(R^-2 + 2 rho_F/M) < 3 R^2;
    canonical radion mass^2 = 2/R^2 + 6 rho_F/M (control M = R = rho_F = 1: 8)."""
    U = -2 * M / R ** 2 - 3 * rho_F; bg = product_background(M, U, rho_F); ell2 = 3 / (1 / R ** 2 + 2 * rho_F / M)
    return {"U": sp.simplify(U), **bg, "ell4_squared": sp.simplify(ell2), "scale_separation": False, "radion_mass2": sp.simplify(2 / R ** 2 + 6 * rho_F / M), "classification": "classical control; not a vacuum, not a KK-stability result"}

def planck_reduction(M6_4, area): return {"M_P_squared": M6_4 * area, "note": "reduction relation, not a scale-selection mechanism (H-25)"}


def radion_brans_dicke(n_internal=2, cassini_bound=40000):
    """Kaluza–Klein reduction over n internal dimensions gives a 4D scalar–tensor theory with the volume modulus as the Brans–Dicke scalar,
    omega_BD = -(n - 1)/n (= -1/2 for a curve).  A massless radion therefore violates the Cassini bound omega > 4e4: radius stabilisation
    (a radion mass) is required by solar-system gravity, not merely desirable.  1/g_4^2 = A/g_6^2: the internal area is the vacuum
    'permittivity' of the gauge sector in the same reduction."""
    from fractions import Fraction
    w = Fraction(-(n_internal - 1), n_internal); return {"omega_BD": w, "passes_cassini_massless": w > cassini_bound, "requirement": "massive (stabilised) radion", "gauge_coupling_relation": "1/g4^2 = A/g6^2"}
