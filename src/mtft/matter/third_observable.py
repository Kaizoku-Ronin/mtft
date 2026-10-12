"""TO-01 (2026-10-05): a pre-registered third observable — and why V_cb never measured the curve's size.

Pre-registration (fixed before computing).  Observable: V_cb predicted from the two sectors' own mass ratios, using only the design's structure (no fitted
parameter, independent of where the top sits).  Pass window: [0.030, 0.055] (the measured 0.041 with estimate-level slack).
1. The second family is induced (EXACT within the design).  Its direct Yukawa entry costs three units of curve-direction breaking (TF-1), payable only by the
   global condensate's curve Taylor data: at most (r1/sqrt(r))^3 = 2.4e-4 (r = 18.6) to 1.8e-5 (r = 105), far below m_c/m_t = 0.0036 and m_s/m_b = 0.019.  So each
   sector's second-family mass comes from the third family's tail: m2/m3 = theta_L theta_R, and with Q and u^c/d^c carrying the same exceptional charges at both
   points the tails are left–right symmetric: theta_L = theta_R = sqrt(m2/m3).
2. The prediction (NUMERICAL).  The CKM element is a difference of rotations: V_cb = |sqrt(m_s/m_b) - e^{i phi} sqrt(m_c/m_t)|, between 0.078 and 0.198 for every phase.
   The minimum misses the window by a factor 1.9.  Fail — wherever the top sits.
3. What this corrects (STRUCTURAL; CC-50).  A geometric angle common to both sectors cancels in the CKM matrix; with an induced second family V_cb is fixed by mass
   ratios, not by the top's distance.  SP-2, DD-4 and CA-2's r-from-V_cb become caveats, and CA-3's Fricke candidate loses its premise (superseded).
4. What the design needs (STRUCTURAL/ESTIMATE).  Left–right asymmetric tails: with V_cb = 0.041 mostly from the down sector, theta_L^d ~ 0.05 and theta_R^d =
   (m_s/m_b)/theta_L^d ~ 0.38, a ratio near 7–8.  [Corrected by CC-51 / FL-01: theta_L is set by the right-handed block's tail, so the Q block's tail must be the
   longer one, and VS-3's ordering a_Q > a_u points the wrong way.]
OPEN: tail ranges as a function of the blocks' curve degrees; the up–down asymmetry from the SU(2)_R-breaking fluxes."""
import math

WINDOW = (0.030, 0.055)


def sector_ratios():
    from mtft.matter.leptons import UP_MZ, DOWN_MZ
    return UP_MZ[1] / UP_MZ[2], DOWN_MZ[1] / DOWN_MZ[2]

def vcb_range():
    mc_mt, ms_mb = sector_ratios(); a, b = math.sqrt(ms_mb), math.sqrt(mc_mt)
    return abs(a - b), a + b

def direct_entry_bound(r, r1=0.268):
    return (r1 / math.sqrt(r)) ** 3

def passes():
    lo, hi = vcb_range(); return not (hi < WINDOW[0] or lo > WINDOW[1])

def lr_asymmetry(V=0.041, theta_u=0.009):
    mc_mt, ms_mb = sector_ratios(); thL = V + theta_u; thR = ms_mb / thL
    return thL, thR, thR / thL

def report():
    return {"ratios": sector_ratios(), "vcb_range": vcb_range(), "passes": passes(), "direct_bound": {r: direct_entry_bound(r) for r in (18.6, 105.0)},
            "lr": lr_asymmetry()}
