"""FL-01 (2026-10-05): which tail sets which angle, what the data demand, and why curve flux cannot supply it.

1. Which tail sets which angle (EXACT).  In Y_ij = int Q_i d^c_j H with the second family at P1 (degree-1 jets) and the third at P2 (degree-0), each off-diagonal
   entry is dominated by the degree-0 field's tail reaching P1 (it falls like (r/d)^3, the jets' like (r/d)^4): Y_23 ~ T(d^c), Y_32 ~ T(Q).  So the left-handed
   angle that enters the CKM matrix, theta_L = Y_23/Y_33, is set by the right-handed block's tail, and theta_R by the Q block's, which up and down share.  TO-01's
   item 4 had the direction backwards (CC-51).
2. What the data demand (EXACT, induced texture plus measured inputs, no geometry).  m_c/m_t = theta_L^u theta_R, m_s/m_b = theta_L^d theta_R, V_cb = |theta_L^d -
   e^{i phi} theta_L^u|.  Hence theta_L^d/theta_L^u = 5.31 and the tails must stand T(Q) : T(d^c) : T(u^c) = 40 : 5.3 : 1 (aligned phases) to 86 : 5.3 : 1 (opposite):
   the Q block's tail longest, u^c's shortest.
3. Curve flux cannot supply it (ESTIMATE).  A block of curve degree a has magnetic length sqrt(A_X / 2 pi a); beyond it tails fall like exp(-d sqrt(2 pi a / A_X)).
   For a separation along the curve d = sqrt(A_X) l (l in area-1 units) the ratio exp(sqrt(2 pi) l (sqrt(a_2) - sqrt(a_1))) does not depend on A_X.  With the
   unique-Higgs condition a_Q + a_u = 24 the most it reaches (a_Q = 0, a_u = 24, l = 0.214, cusp to cusp) is 13.8, short of 40 — and VS-3's stability ordering
   a_Q > a_u points it the wrong way.
4. The remaining source (STRUCTURAL).  The condensate lives in the right-handed block: Pati–Salam-scale effects act on u^c and d^c and not on Q, and its SU(2)_R
   D-term splits u^c from d^c.  That is the one ingredient with the right shape — shorter right-handed tails, u^c shorter than d^c if the sign is right.
OPEN: right-handed tail ranges with the condensate's D-terms; the sign of the u^c–d^c splitting."""
import math

V_CB = 0.041


def required_angles(phase_sign=+1):
    from mtft.matter.third_observable import sector_ratios
    mc_mt, ms_mb = sector_ratios(); k = ms_mb / mc_mt
    thd = V_CB / (1 - phase_sign / k); thu = thd / k; thR = ms_mb / thd
    return {"theta_L_d": thd, "theta_L_u": thu, "theta_R": thR, "ratio_d_over_u": k}

def required_tails(phase_sign=+1):
    a = required_angles(phase_sign)
    return a["theta_R"] / a["theta_L_u"], a["theta_L_d"] / a["theta_L_u"], 1.0

def flux_ratio(a1, a2, l):
    """T(block with a1)/T(block with a2) for a curve-direction separation of l (area-1 units); independent of A_X."""
    return math.exp(math.sqrt(2 * math.pi) * l * (math.sqrt(a2) - math.sqrt(a1)))

def best_flux_ratio(total=24, l=0.214):
    return max(flux_ratio(aq, total - aq, l) for aq in range(0, total + 1))

def report():
    return {"aligned": required_tails(+1), "opposite": required_tails(-1), "best_curve_flux": best_flux_ratio(),
            "stable_pair_13_11": flux_ratio(13, 11, 0.214)}
