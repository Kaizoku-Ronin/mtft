"""FH-01 (2026-10-07): the top at the Fricke point — parity holds, four r measurements meet, and V_ub tests the parity breaking, not the home.

1. Where the Fricke point lands on 143a1 (EXACT/NUMERICAL).  143a1 = [0,-1,1,-1,-2] has root number -1, rank 1 and trivial torsion, so its newform has
   W_143-eigenvalue +1 and phi(W tau) = phi(0) + phi(tau) with phi(0) = O.  The involution therefore pins nothing on the torus: the Fricke point maps to
   phi(i/sqrt 143) = 0.22215 omega_1, on the real component and not a 2-torsion point.  A home on the graph of phi sits there; CA-3 placed it over O.
2. Parity holds anyway (NUMERICAL).  A tail's torus-nonuniform part decays along the curve as exp(-2 pi |k| s), |k| the shortest dual vector of the unit-area
   lattice (0.989).  At the Fricke separation s = 0.107 sqrt(r) = 1.62–1.79 (r = 230–280) that is 1.5e-5 to 4.2e-5: the top's torus position is invisible to the
   tails, the w-family's direct coupling to the top is suppressed by the same factor, and DD-01's selection rule holds whether the home sits over O or over phi(i/sqrt 143).
3. Four measurements of r meet (NUMERICAL).  V_cb at the Fricke home 167–306 (AG-4), V_us 91–896 (CA-1), EPS-01's down-sector tilt 230–2300 and NU-02's first
   window 96–280 share r = 230–280.  With the tau/2 home V_cb's r (51–94) was disjoint from EPS-01's.  The windows are broad: supportive, not decisive.
4. V_ub tests the breaking, not the home (STRUCTURAL).  With a torus-uniform tail the w-family reaches the top only through the x-component parity breaking gives it.
   If the gamma-flux acts as a rotation of the jet basis, the down-sector part of V_ub cancels and |V_ub/V_cb| = theta_12^u, the texture-zero relation: symmetric
   up-sector 1–2 entries give sqrt(m_u/m_c) = 0.045 against the data's 0.090, a factor 2.0; it works if the up-sector left angle is twice its symmetric value
   (entry ratio 4).  A non-rotation breaking gives c ~ 1 (DD-01).  The data's c = 0.40 sits between, as CA-4 said.
OPEN: the form of the gamma-flux's action at the jet point (rotation or independent), which V_ub now measures."""
import math
import numpy as np

PHI_FRICKE = 0.222153                      # lattice coordinate of phi(i/sqrt 143) along omega_1 (PARI, 1200 terms)
CURVE = {"model": (0, -1, 1, -1, -2), "root_number": -1, "rank": 1, "torsion": 1, "w143": +1}
ELL_FRICKE = 0.107
R_WINDOWS = {"V_cb at the Fricke home": (166.8, 305.8), "V_us (CA-1)": (90.8, 896.0), "EPS-01 tilt": (230.0, 2300.0), "NU-02 first window": (96.0, 280.0)}
V_UB, V_CB, V_US = 0.0037, 0.041, 0.225
MU_OVER_MC = 0.0020                        # m_u/m_c at a common scale (about 2e-3)


def tau_143a1():
    from mtft.matter.bl_wilson_line import tau_143a1 as t
    return t()

def shortest_dual():
    t = tau_143a1(); s = 1 / math.sqrt(t.imag)
    B = np.array([[1.0, t.real], [0.0, t.imag]]) * s
    D = np.linalg.inv(B).T
    return min(np.linalg.norm(D @ np.array([m, n])) for m in range(-4, 5) for n in range(-4, 5) if (m, n) != (0, 0))

def separation(r, l=ELL_FRICKE):
    return l * math.sqrt(r)

def torus_suppression(r):
    return math.exp(-2 * math.pi * shortest_dual() * separation(r))

def concordance():
    lo = max(w[0] for w in R_WINDOWS.values()); hi = min(w[1] for w in R_WINDOWS.values()); return (lo, hi)

def texture_zero_ratio():
    """|V_ub/V_cb| from the texture-zero relation with symmetric up-sector 1-2 entries, against the data."""
    return math.sqrt(MU_OVER_MC), V_UB / V_CB

def is_two_torsion(x=PHI_FRICKE, tol=1e-3):
    return min(abs(x - c) for c in (0.0, 0.5, 1.0)) < tol

def report():
    return {"phi_fricke": PHI_FRICKE, "two_torsion": is_two_torsion(), "dual": shortest_dual(), "suppression": [torus_suppression(r) for r in (230, 280)],
            "concordance": concordance(), "texture_zero": texture_zero_ratio()}
