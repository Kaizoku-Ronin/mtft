"""DS-01 (2026-10-07): the down sector is parity-symmetric — and V_ub then has a floor the data sit below.

1. What the condensate touches (EXACT).  The nu^c condensate sits at the (A1, A2) intersection, so its first-order D-term is felt by the fields sharing those stacks
   — u^c (+1), e^c (+1), nu^c (+2), L (-1), H_u (-1) — and by none of Q, d^c, H_d (0): the same pattern as CD-1's dressing.  EPS-01's eps_d/eps_u split is not a
   split: the down sector has no first-order parity breaking at all, and the Cabibbo angle is an up-sector quantity, V_us = sin theta_12^u.
2. The floor (EXACT).  With the down sector parity-symmetric and its (2,2) entry negligible (TO-01), m_s/m_b = tan theta_L^d tan theta_R exactly, and TL-3 gives
   V_ub = V_us sin theta_L^d (the up sector's small local arm shifts this by <= 5%).  Hence
        c = V_ub/(V_us V_cb) >= (m_s/m_b) / (V_cb tan theta_R),     |V_ub| >= V_us (m_s/m_b) / tan theta_R.
   For a right-handed down 2–3 angle at or below 45 degrees: c >= 0.47 (0.42–0.52 with +-10% on m_s/m_b) and |V_ub| >= 0.0043 (0.0039–0.0048), against the
   measured c = 0.40 and |V_ub| = 0.0037 (exclusive) to 0.0042 (inclusive).  The pre-registered window [0.35, 0.45] is missed at the central inputs and grazed
   within their uncertainties.
3. What the data then say (NUMERICAL).  Solving the structure with c = 0.40, V_cb, m_s/m_b, m_c/m_t and gamma = 65.7: tan theta_R = 1.17 (49.5 degrees, the
   right-handed down 2–3 mixing beyond maximal), sin theta_L^u = 0.050 (three times the down angle), 2–3 phase 48 degrees, u^c dressing z = 0.063.
4. The tails cannot do it (ESTIMATE).  With FL-01's assignment (theta_R from the Q block's tail, theta_L from the right-handed block's), tan theta_R = 1.17 needs
   T(Q)/T(d^c) = 71.  Curve flux reaches at most 13.8 (FL-01), which caps c at >= 0.91 and |V_ub| at >= 0.0084: 2.3 times the measured value.  A Q-block tail of
   order one is not a tail — the left-handed top would not be localised.
5. Verdict (STRUCTURAL).  A two-point design whose down sector is parity-symmetric predicts V_ub too large by a factor >= 2.3.  First-order parity breaking in
   the down sector — which the nu^c condensate cannot supply — or a different home for the left-handed top is required.  Recorded as the fifth open failure.
CC-59 supersedes TL-5."""
import math
import numpy as np

STACKS = ("c", "A1", "L", "B", "A2")
CHARGES = {"Q": (-1, 0, 1, 0, 0), "L": (0, -1, 1, 0, 0), "u^c": (1, 0, 0, 0, -1), "d^c": (1, 0, 0, -1, 0), "nu^c": (0, 1, 0, 0, -1),
           "e^c": (0, 1, 0, -1, 0), "H_u": (0, 0, -1, 0, 1), "H_d": (0, 0, -1, 1, 0)}
V_US, V_CB, V_UB = 0.2250, 0.0410, 0.0037
MS_MB, MC_MT, GAMMA = 0.01924, 0.0037, 65.7
FLUX_RATIO_MAX = 13.8                       # FL-01


def touch(field):
    return int(np.dot(CHARGES[field], CHARGES["nu^c"]))

def touched():
    return {f: touch(f) for f in CHARGES}

def angles_2x2(b, bp):
    """Exact: m2/m3 and tan(thetaL) tan(thetaR) for [[0, b], [bp, 1]]."""
    Y = np.array([[0, b], [bp, 1.0]]); U, s, Vt = np.linalg.svd(Y)
    thL = math.atan2(abs(U[0, 0]), abs(U[1, 0])); thR = math.atan2(abs(Vt[0, 0]), abs(Vt[0, 1]))
    return s[1] / s[0], math.tan(thL) * math.tan(thR)

def c_floor(tan_thetaR=1.0, ms_mb=MS_MB):
    return ms_mb / (V_CB * tan_thetaR)

def vub_floor(tan_thetaR=1.0, ms_mb=MS_MB):
    return V_US * ms_mb / tan_thetaR

def solve(c=V_UB / (V_US * V_CB), gamma=GAMMA):
    thLd = math.asin(c * V_CB); tanR = MS_MB / math.tan(thLd)
    vcb = V_CB * np.exp(-1j * math.radians(180 - gamma)); e = math.sin(thLd) - vcb
    s23u = abs(e); z = MC_MT / (math.tan(math.asin(s23u)) * tanR)
    return {"sin_thL_d": math.sin(thLd), "tan_thR": tanR, "thR_deg": math.degrees(math.atan(tanR)), "sin_thL_u": s23u, "phase_deg": math.degrees(np.angle(e)), "z": z}

def tail_ratio_needed(c=V_UB / (V_US * V_CB)):
    thLd = math.asin(c * V_CB); return (MS_MB / math.tan(thLd)) / math.tan(thLd)

def c_with_flux_cap(ratio=FLUX_RATIO_MAX):
    return math.sqrt(MS_MB / ratio) / V_CB
