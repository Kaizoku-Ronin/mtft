"""TA-01 (2026-10-05): auditing the two-point design — the top's blow-up is too big to fit beside the jets on one torus.

1. How far apart two points of the fibre can be (NUMERICAL).  On the unit-area torus of 143a1 the deep hole sits 0.626 from O; the farthest torsion points are the
   3-torsion points at 0.598 (2-torsion: 0.494, 0.563).  No point of the cusp fibre is farther from (cusp, O) than 0.626.
2. How big the blow-ups must be (EXACT constraints, ESTIMATE target).  D-flatness of the Q block: 2 eps1 + eps2 = a_Q A_E (a_Q >= 1).  Both points on one fibre:
   eps1 + eps2 < A_E.  The top (TP-2) needs (2 r2)^2 ~ y_t = 0.70, so eps2 ~ 0.55 A_E (radius 0.418); the low end eps2 = A_E/3 (radius 0.326) needs an O(1) boost.
   With a_Q = 1: eps1 = 0.225 (radius 0.268) or 0.333 (0.326).
3. The mixing (ESTIMATE).  The harmonic tail of a degree-k jet falls like (r/d)^(3+k) (the Bochner–Martinelli kernel of C^2 and its derivatives).  At the best
   torsion separation 0.598 the top's tail at the jet point is 0.34 (eps2 = 0.55) or 0.16 (eps2 = 1/3); the jets' tail at the top is 0.04 or 0.09.  An
   unsuppressed off-diagonal entry of that size gives V_cb ~ 0.16–0.34 and m_c/m_t ~ its square — against 0.041 and 0.0036.  The same-fibre design fails by 4–8x.
4. Where the top must sit (STRUCTURAL).  V_cb = 0.041 needs d ~ r2 / 0.041^(1/3) = 2.9 r2 ~ 0.95–1.21, beyond the fibre's 0.626: the top's point must lie on a
   different curve fibre.  Another cusp is too far (CC-37: mixing <= 1/143), so it must be a non-cusp point at about three blow-up radii from the cusp.
   Natural candidates: special points of X_0(143) such as the preimages of torsion points of E (NU-02's points).  The jets' three units are odd, so with
   torus parity (LB-1) at least one unit must be curve-direction — the global condensate's curve Taylor data can supply it.
OPEN: a natural non-cusp point at the right distance; PA-01 stability with two points; the condensate's global fraction."""
import math
import numpy as np


def fibre_geometry(grid=401):
    from mtft.research.theta_torus import tau_143a1
    tau = complex(tau_143a1()); s = 1 / math.sqrt(tau.imag)
    L = [s * (m + n * tau) for m in range(-3, 4) for n in range(-3, 4)]
    dist = lambda z: min(abs(z - l) for l in L)
    us = np.linspace(0, 1, grid)
    hole = max(dist(s * (u + v * tau)) for u in us for v in us)
    tors = {n: sorted({round(dist(s * (i + j * tau) / n), 4) for i in range(n) for j in range(n) if (i, j) != (0, 0)}) for n in (2, 3)}
    return {"deep_hole": hole, "torsion": tors}

def radius(eps):
    return math.sqrt(eps / math.pi)

def blowups(eps2, a_Q=1, A_E=1.0):
    eps1 = (a_Q * A_E - eps2) / 2
    return {"eps1": eps1, "eps2": eps2, "r1": radius(eps1), "r2": radius(eps2), "ample_on_fibre": eps1 + eps2 < A_E, "d_flat": abs(2 * eps1 + eps2 - a_Q * A_E) < 1e-12}

def tail(r, d, k):
    return (r / d) ** (3 + k)

def mixing(eps2, d=0.598):
    b = blowups(eps2)
    return {"top_tail_at_jets": tail(b["r2"], d, 0), "jet_tail_at_top": tail(b["r1"], d, 1)}

def required_separation(eps2, V=0.041):
    return radius(eps2) / V ** (1 / 3)

def report():
    g = fibre_geometry()
    return {"fibre": g, "eps2=0.55": (blowups(0.55), mixing(0.55)), "eps2=1/3": (blowups(1 / 3), mixing(1 / 3)),
            "required_d": {0.55: required_separation(0.55), 1 / 3: required_separation(1 / 3)}}
