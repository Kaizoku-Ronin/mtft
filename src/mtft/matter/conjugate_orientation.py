"""CO-01 (2026-10-05): is there a conjugate condensate?  No — and the one thing that tells right from left.

1. The D-term vacuum in closed form (EXACT).  For a rank-one Sigma in the (4bar,1,2) block with FI xi and equal couplings, V = (1/2)[(v^2 - xi)^2 + (5/8) v^4]
   (SU(2)_R contributes v^4/4, SU(4) 3 v^4/8), so v^2 = (8/13) xi: the residual U(1) D-term is -(5/13) xi and the nonabelian ones are set by v^2 — PB-01's x_*
   with equal couplings.  An aligned conjugate Sigma-bar enters only through v^2 - v-bar^2: the D-term potential is exactly degenerate along that family, and any
   positive mass (the one-loop masses of Matter-05) lifts it to v-bar = 0.  A conjugate condensate would need a negative mass of its own, which no FI can give
   (conjugate fields carry opposite U(1) charges) and no other Pati–Salam-scale condensate exists to drive it.  One orientation is a theorem of the field content,
   not an artefact of PB-01's scan; VM-4's candidate is closed.
2. Gauge loops do not tell right from left (NUMERICAL).  Radiative Dirac masses scale with C_4 alpha_4 + C_2 alpha_{L,R} (Casimirs 15/8 and 3/4): at the b–tau scale
   m_R / m_L = 0.90 — ten per cent, not the factor five TH-01 needs.
3. What does (EXACT at tree level).  The right-handed block hosts the condensate, so its other scalars sit in the residual D-terms: the e^c-type partner gets
   m^2 = (5/13 + 2/13 + 3/13) g^2 xi = (10/13) g^2 xi, the coloured partners (5/13 - 2/13 - 1/13) = (2/13) g^2 xi — Pati–Salam-scale masses — while the left-handed
   block feels no SU(2)_R D-term.  Fermion pairs stay massless at tree level on both sides (VM-2).
4. The mechanism this leaves (STRUCTURAL).  Radiative Dirac masses need the scalar partners in the loop; heavy partners suppress them.  The condensate therefore
   lowers the right-handed fermion pairs' one-loop masses relative to the left-handed ones — the direction TH-01 needs — by a factor set by M_PS over the
   Scherk–Schwarz scale.  Its size is the next computation.
OPEN: the one-loop Dirac masses with Scherk–Schwarz gauginos and D-term-split partners; the sign of the (4,2,1) scalars' residual U(1)_4 D-term."""
import math


def d_term_potential(v2, xi=1.0, vbar2=0.0):
    u = v2 - vbar2
    return 0.5 * ((u - xi) ** 2 + (5 / 8) * u * u)

def d_term_minimum(xi=1.0):
    x = 8 / 13
    return {"v2_over_xi": x, "residual_U1": x - 1, "residual_SU2R": x / 2, "SU4_sum_D2_over_xi2": (3 / 8) * x * x}

def lifted_minimum(xi=1.0, m2=0.05, n=400):
    """Minimise V_D + m2 (v^2 + vbar^2) on a grid: the conjugate component vanishes."""
    best = None
    for i in range(n + 1):
        for j in range(n + 1):
            v2, vb2 = 1.5 * xi * i / n, 1.5 * xi * j / n
            val = d_term_potential(v2, xi, vb2) + m2 * (v2 + vb2)
            if best is None or val < best[0]: best = (val, v2, vb2)
    return {"v2": best[1], "vbar2": best[2]}

def scalar_masses():
    """m^2 / (g^2 xi) of the right-handed block's non-condensate scalars at the D-term minimum."""
    return {"e_c_partner": 5 / 13 + 2 / 13 + 3 / 13, "coloured_partner": 5 / 13 - 2 / 13 - 1 / 13}

def gauge_loop_asymmetry(mu=1.1236e6):
    from mtft.matter.gauge_couplings import _inv
    a4 = 1 / _inv("3", mu); aL = 1 / _inv("2", mu); aR = 1 / (_inv("Y", mu) - 2 / 3 * _inv("3", mu))
    C4, C2 = 15 / 8, 3 / 4
    return (C4 * a4 + C2 * aR) / (C4 * a4 + C2 * aL)

def report():
    return {"minimum": d_term_minimum(), "lifted": lifted_minimum(), "flat_check": d_term_potential(1.0, 1.0, 0.4) - d_term_potential(0.6, 1.0, 0.0),
            "scalars": scalar_masses(), "gauge_loop_asymmetry": gauge_loop_asymmetry()}
