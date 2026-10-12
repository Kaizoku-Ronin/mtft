"""DP-01 (2026-10-05): how strong the dressing must be, and why one dressing cannot serve both quarks and leptons.

1. The light state stays massless (EXACT).  Gauge invariance pairs u^c with the conjugate partner of the leptoquark gaugino, which has the Dirac Kaluza–Klein mass M:
   one Dirac combination takes mass sqrt(M^2 + m^2), the orthogonal one stays massless with weight Z = M / sqrt(M^2 + m^2) on the original family state.  No seesaw
   mass appears.
2. Character selection (EXACT).  The gauginos are constant along the curve at lowest order, so the charm's u^c — an x-jet at P1 — mixes only through the condensate's
   own x-jet component at P1, and the top's u^c — degree 0 at P2 — through the condensate's value at P2.  z_u = Z_2/Z_3 compares those two local amplitudes.
3. How strong (ESTIMATE).  z_u = 0.25 (CD-3's symmetric-tail point) needs m_2/M = 3.9 with an undressed top, 4.9 if Z_3 = 0.8 (m_3/M = 0.75): the condensate's x-jet
   at P1 several times its value at P2, and the gaugino tower near M_PS/2 — the scale VS-4 already wanted.  Demanding, not unnatural.
4. Quarks and leptons conflict (EXACT given the structure).  With Pati–Salam-symmetric tails (e^c with d^c, L with Q), (m_mu/m_tau)/(m_s/m_b) = z_e = 3.09.  The same
   condensate components and the same gaugino towers dress e^c exactly as they dress u^c, so z_e = z_u: the up sector's 0.25 and the leptons' 3.09 differ by a
   factor 12.  Taking z = 3.09 instead brings back a u^c/d^c tail ratio of 16, beyond what flux can give (FL-3).  One dressing cannot serve both sectors unless the
   leptoquark and W_R towers differ; SU(2)_R-breaking flux on the (B, A2) block can split them, but it would have to reverse the spatial pattern, not rescale it.
OPEN: whether the (B, A2) flux can make the W_R tower dress the tau's e^c more than the muon's while the leptoquark tower does the opposite for u^c."""
import math


def weight(m_over_M):
    return 1 / math.sqrt(1 + m_over_M ** 2)

def masses_for_z(z, Z3=1.0):
    Z2 = z * Z3
    return math.sqrt(1 / Z2 ** 2 - 1), math.sqrt(max(1 / Z3 ** 2 - 1, 0.0))

def z_leptons():
    from mtft.matter.leptons import LEPTON_POLE, DOWN_MZ
    return (LEPTON_POLE[1] / LEPTON_POLE[2]) / (DOWN_MZ[1] / DOWN_MZ[2])

def conflict(z_up=0.25):
    return z_leptons() / z_up

def tail_ratio_if_shared(z):
    from mtft.matter.third_observable import sector_ratios
    mc_mt, ms_mb = sector_ratios()
    return (mc_mt / ms_mb) / z

def report():
    return {"masses_z025_Z3_1": masses_for_z(0.25, 1.0), "masses_z025_Z3_08": masses_for_z(0.25, 0.8), "z_e": z_leptons(), "conflict": conflict(),
            "u_over_d_tail_at_ze": tail_ratio_if_shared(z_leptons())}
