"""PAR-01 (2026-10-08): the parity audit (filed as PA-01, which is Path A's label; relabelled by CC-62; claims PA-1 to PA-4) — every source odd under the torus involution, what it reaches at first order, and the verdict on the two-point design.

The selection rule of DD-01 is the torus involution w -> -w on 143a1.  For the down sector to mix its first family it needs an odd source that reaches Q, d^c or
H_d at first order with more than a basis rotation.  The design's ingredients, sorted:

  odd sources                      first-order reach                                               effect on 1–2 mixing
  gamma mixed flux (Poincare type) every gamma-charged block (Q, u^c, d^c, L, e^c, nu^c), condensate basis rotation of each block — cancels in U_u^dagger U_d (TL-2);
                                                                                                    acts physically only through the condensate's profile
  condensate's odd components      u^c (+1), e^c (+1), nu^c (+2), L (-1), H_u (-1); not Q, d^c, H_d  real — but up sector and leptons only (DS-1)
  P2's torus displacement, and the 13 odd complex-structure moduli H^1(O_X) (x) H^0(T_E)
                                   the third families of every block                                through tails only: exp(-2 pi |k| s) ~ 1e-5 at the Fricke
                                                                                                    separation (FH-2)
  torus Wilson lines               gamma-charged blocks                                            continuous ones have positive masses and sit at even values
                                                                                                    (VM-01); the Z/2 ones are saddles (BL-3)
  even (no entry)                  curve fluxes a_Q, a_u, a_d; exceptional charges; torus flux densities; the lattice's six-fold structure and its two-unit
                                   breaking (torus charge 2 = 0 mod 2); the three even spin structures (ι acts trivially on H_1(E, Z/2)); the Phi-type Higgs
                                   modes (torus-constant); the curve involutions.

1. No odd source reaches Q, d^c or H_d at first order with a physical effect (EXACT, given the ingredient list).  The down sector is exactly parity-symmetric
   at first order; DS-2's floor applies: |V_ub| >= V_us (m_s/m_b) / tan theta_R.
2. The one loophole (NUMERICAL).  The floor used the up sector's local arm as m_u fixes it (v_u = 0.0044).  Evading the floor through that arm needs
   |v_u| >= 3.1 z V_us — 0.18 at z = 0.25 — which makes m_u/m_c = |v_u (2 u_u - v_u)| between 0.048 and 0.11, 24 to 56 times the measured 0.002, unless tuned.
3. The 13 odd deformations (EXACT).  The product's complex-structure deformations split into even (H^1(T_X), 36; H^1(T_E), 1) and odd (H^1(O_X) (x) H^0(T_E), 13).
   The odd ones vary the torus origin along the curve — the graph-of-phi structure FH-1 met — and reach the jet point only through tails.
4. Verdict (STRUCTURAL).  With tails capped at 13.8 (FL-01) the two-point design predicts |V_ub| >= 0.0084, 2.3 times the measured value: it is ruled out as a
   flavour model.  Two redesigns remain: a down-sector odd source, which needs a field on the B stack with an odd profile that the spectrum does not contain; or a
   different home for the left-handed top, so that the right-handed 2–3 angle is of order one — whose first question is whether y_t = 0.70 survives H_u bridging
   two points."""
import math

SOURCES = {
    "gamma mixed flux": {"odd": True, "reach": ("Q", "u^c", "d^c", "L", "e^c", "nu^c", "condensate"), "effect": "basis rotation; physical only via the condensate"},
    "condensate odd components": {"odd": True, "reach": ("u^c", "e^c", "nu^c", "L", "H_u"), "effect": "real, up sector and leptons"},
    "P2 torus displacement / odd moduli": {"odd": True, "reach": ("Q3", "u^c3", "d^c3"), "effect": "tails only, 1e-5"},
    "torus Wilson lines": {"odd": True, "reach": (), "effect": "even values (VM-01) or saddles (BL-3)"},
    "curve fluxes": {"odd": False}, "exceptional charges": {"odd": False}, "torus flux densities": {"odd": False},
    "lattice six-fold and two-unit breaking": {"odd": False}, "even spin structures": {"odd": False}, "Phi-type Higgs modes": {"odd": False},
}
DOWN = ("Q", "d^c", "H_d")
V_US, V_CB, V_UB = 0.2250, 0.0410, 0.0037
MS_MB, MC_MT, MU_MC = 0.01924, 0.0037, 0.0020
THETA_R_NATURAL = 0.43                       # FL-01-type tails
GENUS = 13


def down_sector_odd_sources():
    """Odd sources whose first-order reach includes a down-sector field with a physical (non-rotation, non-tail) effect."""
    return [k for k, v in SOURCES.items() if v["odd"] and any(f in v["reach"] for f in DOWN) and "rotation" not in v["effect"] and "tails" not in v["effect"]]

def odd_moduli():
    return {"even": (3 * GENUS - 3) + 1, "odd": GENUS}

def floor_vub(tan_theta_R=1.0):
    return V_US * MS_MB / tan_theta_R

def loophole(z=0.25, theta_R=THETA_R_NATURAL, c=V_UB / (V_US * V_CB)):
    """Up-sector local arm needed to bring c down to the data with a parity-symmetric down sector, and the m_u/m_c it implies (both signs)."""
    s23d = MS_MB / math.tan(theta_R); s23u = MC_MT / (z * math.tan(theta_R))
    kappa_minus_1 = (s23d - c * V_CB) / s23u; v = kappa_minus_1 * V_US
    out = {}
    for sign in (-1, +1):
        vv = sign * v; u = V_US + vv; out[sign] = {"v_u": vv, "mu_mc": abs(vv * (2 * u - vv)), "ratio_to_data": abs(vv * (2 * u - vv)) / MU_MC}
    return {"v_needed": v, "cases": out}

def verdict(flux_ratio_max=13.8):
    c_floor = math.sqrt(MS_MB / flux_ratio_max) / V_CB
    return {"c_floor": c_floor, "vub_floor": c_floor * V_US * V_CB, "factor_over_data": c_floor * V_US * V_CB / V_UB}
