"""GC-01 (2026-10-05): what the gauge couplings fix — not the curve's shape, but the Pati–Salam scale, and they disagree with b–tau.

1. Shape versus volume (EXACT).  The volume of the blown-up surface is J^2/2 = A_X A_E - sum eps_k^2 / 2.  A stack's gauge kinetic function is
   Re f_i = J^2/2 - kappa F_i^2/2 with F_i^2 = 2 a_i b_i - 8 c_i^2 - sum_k e_ik^2 (kappa = 0 for pure 8D SYM, > 0 with the usual flux correction).  With no torus
   flux, A_X enters every coupling only through the common volume: coupling ratios fix V/kappa, never the shape A_X/A_E.  The couplings cannot turn V_cb
   into a prediction; V_cb remains the measure of the shape.
2. The matching condition (NUMERICAL, SM one loop).  With alpha_4 = alpha_3, alpha_L = alpha_2 and alpha_Y^-1 = alpha_R^-1 + (2/3) alpha_4^-1, the flux-only splitting
   predicts R = (alpha_R^-1 - alpha_4^-1)/(alpha_L^-1 - alpha_4^-1) = (sum e_R^2 - sum e_4^2)/(sum e_L^2 - sum e_4^2), independent of V/kappa.  The data give R = 3.99
   at m_Z, 3.678 at the b–tau crossing (1.12e6 GeV), 1 at 5.08e13 GeV (g_L = g_R); alpha_2 = alpha_3 at 9.71e16 GeV.
3. Quantised options (EXACT given f).  Two-point designs (family charges (2, 1) in both Q and u^c) give R = (5 + 2m)/(5 - 2m): 1/9, 3/7, 1, 7/3, 9; the single-point
   |e| = 3 design gives 1/5, 1, 5.  None reaches 3.678.  The matching scales are 8.9e11 (7/3), 5.1e13 (1), 1.4e14 (3/7), 2.1e14 GeV (1/9); 5 and 9 have none
   above m_Z.
4. The tension (ESTIMATE).  Every allowed matching scale lies 6–8 orders above the b–tau crossing, where y_b/y_tau = 0.71 (8.9e11) to 0.66: b–tau would need a
   30–35 per cent correction there.  The scale used since PS-01 (1.1e6 GeV) is not where the gauge couplings put Pati–Salam; PS1-4's reading as M_PS is a caveat
   (CC-48), and everything evaluated at 1.1e6 GeV (LEP-01's ratios, LB-01's y_t = 0.70) inherits it.
5. Which home (EXACT by symmetry).  P1 = (cusp, O) is fixed by the torus involution [-1], so the potential for the top's point is [-1]-invariant and its torus
   coordinate is critical at the four 2-torsion points (O included).  No automorphism of X_0(143) fixes a cusp, so the curve coordinate has no symmetry selection.
OPEN: threshold corrections (vector-like pairs, adjoints, KK towers) that could reconcile the two scales; the shape enters the couplings only there."""
import math
from scipy.optimize import brentq


def _inv(k, mu):
    from mtft.matter.ps_breaking import INPUTS_MZ as I
    mZ = I["m_Z"]; aem = 1 / I["alpha_inv"]; s2 = I["sin2w"]
    inv0 = {"Y": (1 - s2) / aem, "2": s2 / aem, "3": 1 / I["alpha_s"]}; b = {"Y": 41 / 6, "2": -19 / 6, "3": -7.0}
    return inv0[k] - b[k] / (2 * math.pi) * math.log(mu / mZ)

def ratio(mu):
    a4 = _inv("3", mu); aL = _inv("2", mu); aR = _inv("Y", mu) - 2 / 3 * a4
    return (aR - a4) / (aL - a4)

def matching_scale(R, lo=91.1876, hi=1e15):
    f = lambda t: ratio(math.exp(t)) - R
    ts = [math.log(lo) + k * (math.log(hi) - math.log(lo)) / 400 for k in range(401)]
    for t1, t2 in zip(ts, ts[1:]):
        if f(t1) * f(t2) < 0: return math.exp(brentq(f, t1, t2))
    return None

def alpha23_scale():
    return math.exp(brentq(lambda t: _inv("2", math.exp(t)) - _inv("3", math.exp(t)), math.log(1e10), math.log(1e19)))

def design_ratios(charges=(2, 1), mmax=2):
    n = sum(c * c for c in charges)
    return sorted({(n + 2 * m) / (n - 2 * m) for m in range(-mmax, mmax + 1) if n - 2 * m > 0 and n + 2 * m > 0})

def single_point_ratios():
    return sorted({(9 + 6 * m) / (9 - 6 * m) for m in (-1, 0, 1)})

def b_tau_at(mu):
    import numpy as np
    from mtft.matter.ps_breaking import b_tau_running
    c = b_tau_running()["curve_log10mu_R"]
    return float(np.interp(math.log10(mu), [e for e, _ in c], [r for _, r in c]))

def volume(A_X, A_E, eps):
    return A_X * A_E - sum(e * e for e in eps) / 2

def report():
    bt = 1.1236e6
    return {"R_mZ": ratio(91.1876), "R_btau": ratio(bt), "LR_scale": matching_scale(1.0), "alpha23": alpha23_scale(),
            "two_point": {R: matching_scale(R) for R in design_ratios()}, "single_point": {R: matching_scale(R) for R in single_point_ratios()},
            "b_tau_at_matching": {R: b_tau_at(s) for R, s in ((7 / 3, matching_scale(7 / 3)), (1.0, matching_scale(1.0)))}}
