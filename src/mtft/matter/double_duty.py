"""DD-01 (2026-10-05): making the shape do double duty — a parity selection rule for the CKM matrix.

1. Parities (EXACT).  The torus involution w -> -w fixes the jet point P1 = (cusp, O) and both of the top's candidate homes (a point over a 2-torsion point of E, or
   another cusp over O), since 2-torsion points are their own negatives.  The constant-jet top is even, the curve-direction jet x is even, the torus-direction
   jet w is odd.
2. A decoupled family (EXACT).  In a background that respects the involution every Yukawa entry between opposite parities vanishes, in the up and the down sector
   alike: the w-family mixes with nothing.  The CKM matrix would then have a decoupled row and column, so either V_us = V_ub = 0 or V_cb = V_ub = 0 — both excluded.
   The data require the involution to be broken.
3. The hierarchy this gives (STRUCTURAL).  With breaking of size eps_P and the odd family as the first family: V_us ~ eps_P (x-w mixing), V_cb set by the top's tail
   (geometry, fixed by A_X), V_ub ~ c eps_P V_cb.  The data need eps_P ~ 0.225 and c = V_ub/(V_us V_cb) = 0.40 — a natural O(1).  Wolfenstein's lambda, lambda^2,
   lambda^3 pattern reads as parity breaking, geometry, and their product.
4. The source (STRUCTURAL).  The gamma-flux is odd under the involution, so gamma-charged matter blocks (TH-3's c_L = c_R, which also supply the vector-like pairs)
   break it: the Cabibbo angle measures the matter blocks' gamma-charge at the jet point.
5. A refinement of SP-2 (ESTIMATE).  The top couples to the degree-1 x-jet through the gradient of its tail: V_cb ~ 3 (r1/d)(r2/d)^3.  The required separation is
   0.953 (top area A_E/3) to 1.094 (0.55 A_E), and over a 2-torsion point the curve's area becomes 18–28 A_E (was 18–36).
OPEN: the O(1) coefficient c from the overlaps; eps_P from the gamma-flux at the jet point."""
import math

V_US, V_CB, V_UB = 0.2250, 0.0410, 0.0037


def parities():
    return {"top": +1, "x_jet": +1, "w_jet": -1}

def allowed(a, b):
    p = parities(); return p[a] * p[b] == +1

def decoupled_family():
    fams = ("top", "x_jet", "w_jet")
    return [f for f in fams if all(not allowed(f, g) for g in fams if g != f)]

def required_coefficient():
    return V_UB / (V_US * V_CB)

def separation(eps2, V=V_CB):
    r2 = math.sqrt(eps2 / math.pi); r1 = math.sqrt((1 - eps2) / 2 / math.pi)
    return (3 * r1 * r2 ** 3 / V) ** 0.25

def curve_area(eps2, path=0.1781, offset=0.5629):
    d = separation(eps2); return (d * d - offset * offset) / path ** 2

def report():
    return {"parities": parities(), "decoupled": decoupled_family(), "c_required": required_coefficient(),
            "d": {e: separation(e) for e in (1 / 3, 0.55)}, "A_X": {e: curve_area(e) for e in (1 / 3, 0.55)}}
