"""CA-01 (2026-10-05): the Cabibbo angle from the gamma-flux, and two measurements of the curve's size.

1. The Cabibbo angle measures r (ESTIMATE).  EPS-01's gamma-tilt, eps ~ (2/3 ... 2 pi/3) |c| phi', phi' = 3.216/sqrt(r), r = A_X/A_E, is the parity breaking DD-3 needs.
   V_us = 0.225 then requires r = 90.8-896 for gamma-charge c = 1 and 363-3584 for c = 2.
2. Two measurements disagree (STRUCTURAL/ESTIMATE).  V_cb with the top over a 2-torsion point needs r = 18.6-27.8 (DD-4).  The two ranges miss by a factor 3.3 at
   their closest ends; EPS-01's down-sector tilt (230-2300) and NU-02's neutrino suppression (96-280, 380-1080) both side with the larger window.  The 2-torsion
   home is the odd one out.
3. A post-hoc candidate (MODEL, spotted after the tension).  The Fricke fixed point i/sqrt(143) — 0.107 from the cusp, fixed by W_143, and with torus point O fixed
   by the torus involution — would give V_cb = 0.041 at r = 79-104, overlapping the Cabibbo window's lower edge and NU-02's second-family window.  Chosen after
   seeing the numbers, so it is a hypothesis for a third observable to test, not evidence.
4. V_ub's coefficient (STRUCTURAL).  The data's c = 0.40 lies between the parity texture's down-sector value (about 1) and the up-sector Gatto–Sartori–Tonin value
   sqrt(m_u/m_c)/V_us = 0.20: consistent, not sharp.
5. The vector-like pairs, again (ESTIMATE).  c = 1 gives 4c^2 + 1 = 5 pairs per block; TH-01's n l = 29.2 then puts the right-handed pairs at M_PS/e^5.84, about
   3 TeV — within collider reach.  c = 2's 17 pairs would need r far from every other measurement.
OPEN: the overlap coefficients behind EPS-01's (2/3 ... 2 pi/3) band; a third observable at the Fricke point."""
import math

V_US, V_CB, V_UB = 0.225, 0.041, 0.0037
PHI_COEFF = 3.216
BAND = (2 / 3, 2 * math.pi / 3)


def r_for_cabibbo(c=1, V=V_US):
    return tuple((b * PHI_COEFF * c / V) ** 2 for b in BAND)

def r_from_vcb(path, offset=0.0):
    from mtft.matter.double_duty import separation
    return tuple((separation(e) ** 2 - offset ** 2) / path ** 2 for e in (1 / 3, 0.55))

def gap(lo_range, hi_range):
    return min(hi_range) / max(lo_range)

def c_bracket(m_u=0.00127, m_c=0.619):
    return (math.sqrt(m_u / m_c) / V_US, 1.0, V_UB / (V_US * V_CB))

def pairs_and_scale(c=1, M_PS=1.1236e6):
    from mtft.matter.thresholds import required_nl
    n = 4 * c * c + 1; l = required_nl(7 / 3) / n
    return n, M_PS / math.exp(l)

def report():
    two_torsion = r_from_vcb(0.1781, 0.5629); fricke = r_from_vcb(0.107, 0.0)
    return {"cabibbo_c1": r_for_cabibbo(1), "cabibbo_c2": r_for_cabibbo(2), "vcb_two_torsion": two_torsion, "vcb_fricke": fricke,
            "gap": gap(two_torsion, r_for_cabibbo(1)), "c_bracket": c_bracket(), "pairs_c1": pairs_and_scale(1)}
