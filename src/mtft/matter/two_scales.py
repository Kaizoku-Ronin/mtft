"""GC-02 (2026-10-17): two scales, not one — the flux-split couplings matched at the curve's compactification scale M_C, above the
Pati–Salam scale M_PS, with Pati–Salam running between them.

GC-01 matched at M_PS, i.e. took M_C = M_PS.  VE-5 puts the curve's Kaluza–Klein tower at M_C = 1/(2 R_C) = 0.063-0.070 M_KK, and the
breaking is a 4D condensate below the KK scales (SC-01's reading R1), so M_PS <= M_C.  Chain: SM one loop to M_PS; alpha_4 = alpha_3,
alpha_L = alpha_2, alpha_Y^-1 = alpha_R^-1 + (2/3) alpha_4^-1 there; Pati–Salam one loop to M_C; the tree-level flux splitting as the boundary
condition at M_C.  With PS-01's content b_4 = -32/3 and b_L = b_R = 2/3 + n_H/3 for any bidoublet Higgs content, so both differences
alpha_{L,R}^-1 - alpha_4^-1 shrink at the same rate and R = (alpha_R^-1 - alpha_4^-1)/(alpha_L^-1 - alpha_4^-1) RISES above M_PS: the designs
below the SM value (1/9, 3/7, 1, 7/3, 1/5) need M_C = M_PS (GC-01's scales), while R = 5 and 9 — unreachable in GC-01 — are met at a
computable M_C(M_PS).  With M_PS at the b–tau crossing (1.12e6 GeV): M_C = 1.7e7 (R = 5) and 2.8e8 GeV (R = 9), M_KK = 2.5e8-4.4e9 GeV,
M_SS = 1.2e8-2.2e9 GeV, and GC-4's 30-35 per cent correction is not needed.  Convention d(1/alpha)/d ln mu = -b/(2 pi).
Status: GC-6 EXACT (one loop), GC-7 NUMERICAL; not included: thresholds at M_C, which background realises R = 5 or 9, the K_L -> mu e bound."""
import math

from scipy.optimize import brentq

from mtft.matter.gauge_couplings import _inv, ratio, design_ratios, single_point_ratios, b_tau_at

PS_CONTENT = {"(4,2,1)": 3, "(4bar,1,2)": 3, "(1,2,2)": 6}        # Weyl fermions (PS-01's net indices)
CURVE_RUNG = (0.063, 0.070)                                        # M_C / M_KK (VE-5)


def ps_b(n_bidoublet_fermions=6, n_higgs=1, n_delta=0):
    """One-loop Pati–Salam coefficients (b_4, b_L, b_R) for 3 (4,2,1) + 3 (4bar,1,2) + n_bf (1,2,2) Weyl fermions, n_higgs complex (1,2,2)
    scalars and n_delta complex (4bar,1,2) scalars.  T(fund) = 1/2; Weyl 2/3, complex scalar 1/3, vector -11/3 C_A."""
    f4 = 6 + 6                                                      # Weyl SU(4) (anti)fundamentals
    fL = 12 + 2 * n_bidoublet_fermions; fR = 12 + 2 * n_bidoublet_fermions
    b4 = -11 / 3 * 4 + (2 / 3) * 0.5 * f4 + (1 / 3) * 0.5 * (2 * n_delta)
    bL = -11 / 3 * 2 + (2 / 3) * 0.5 * fL + (1 / 3) * 0.5 * (2 * n_higgs)
    bR = -11 / 3 * 2 + (2 / 3) * 0.5 * fR + (1 / 3) * 0.5 * (2 * n_higgs + 4 * n_delta)
    return {"4": b4, "L": bL, "R": bR}


def ps_inverse_couplings(mu, M_PS, **content):
    """alpha_4^-1, alpha_L^-1, alpha_R^-1 at mu >= M_PS: SM one loop to M_PS, the PS relation there, then PS one loop."""
    a4 = _inv("3", M_PS); aL = _inv("2", M_PS); aR = _inv("Y", M_PS) - 2 / 3 * a4
    b = ps_b(**content); l = math.log(mu / M_PS)
    return {"4": a4 - b["4"] / (2 * math.pi) * l, "L": aL - b["L"] / (2 * math.pi) * l, "R": aR - b["R"] / (2 * math.pi) * l}


def two_scale_ratio(mu, M_PS, **content):
    a = ps_inverse_couplings(mu, M_PS, **content)
    return (a["R"] - a["4"]) / (a["L"] - a["4"])


def compactification_scale(R, M_PS, hi=1e19, **content):
    """The M_C >= M_PS at which the flux-split ratio R is met after Pati–Salam running from M_PS; None if never (R below the SM value at M_PS)."""
    f = lambda t: two_scale_ratio(math.exp(t), M_PS, **content) - R
    ts = [math.log(M_PS) + k * (math.log(hi) - math.log(M_PS)) / 600 for k in range(601)]
    for t1, t2 in zip(ts, ts[1:]):
        a = ps_inverse_couplings(math.exp(t2), M_PS, **content)
        if a["L"] - a["4"] <= 0: break                                  # the denominator's zero: R has run through infinity
        if f(t1) * f(t2) < 0: return math.exp(brentq(f, t1, t2))
    return None


def two_scale_report(M_PS=1.1236e6, **content):
    """GC-01's quantised ratios met at M_C after PS running from M_PS, with the scales VE-5 and SC-01 tie to M_C."""
    out = {"M_PS": M_PS, "R_at_M_PS": ratio(M_PS), "ps_b": ps_b(**content), "designs": {}}
    for R in sorted(set(design_ratios()) | set(single_point_ratios())):
        Mc = compactification_scale(R, M_PS, **content)
        out["designs"][R] = None if Mc is None else {"M_C": Mc, "M_KK": tuple(Mc / c for c in CURVE_RUNG),
                                                     "M_SS": tuple(0.4943 * Mc / c for c in CURVE_RUNG), "y_b_over_y_tau_at_M_PS": b_tau_at(M_PS)}
    return out
