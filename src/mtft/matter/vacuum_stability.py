"""VS-01 (2026-10-05): the condensate's feedback on the (4,2,1) scalars — a D-flat Q block breaks colour.

1. The feedback (EXACT, tree level, PB-01's single relative U(1)).  Sigma in the 4bar carries U(1)_4 charge -1, so D_4 = xi_4 - v^2, and its SU(4) D-term is
   D^A = -v^2 (T^A)_44.  A (4,2,1) scalar in SU(4) direction b has m^2 = g^2 [xi_4L - v^2] + Sum_A D^A (T^A)_bb with Sum_A (T^A)_44 (T^A)_bb = (delta_b4 - 1/4)/2:
   squark doublets  m^2 = g^2 [xi_4L - (7/8) v^2],   slepton doublet  m^2 = g^2 [xi_4L - (11/8) v^2],   with v^2 = (8/13) xi_R4 (CO-1).
2. Consequence (EXACT).  A D-flat Q block (xi_4L = 0 — TP-4's relation, FS-1's design criterion) makes every (4,2,1) scalar tachyonic once the condensate forms:
   colour and the electroweak group break at the Pati–Salam scale.  Stability needs xi_4L >= (11/13) xi_R4; CC-49 caveats TP-4 and FS-1.
3. What it selects (STRUCTURAL).  With exceptional terms E entering both FIs alike, xi_4L / xi_R4 = (a_Q A_E - E)/(a_u A_E - E): the unique-Higgs pair (a_Q, a_u) =
   (11, 13) sits below 11/13 for any E > 0 and is unstable; (13, 11) is stable.  Stability orders the curve degrees: a_Q > a_u.
4. The asymmetry this fixes (STRUCTURAL/ESTIMATE).  Left-handed scalar partners sit near the stability edge (light); right-handed ones at (2/13)-(10/13) g^2 xi
   (CO-3, Pati–Salam scale).  Radiative Dirac masses with heavy partners fall like (M_SS / m_s)^2, so m_R / m_L ~ (M_SS / M_PS)^2 x O(1): the factor five of
   TH-01 needs the Scherk–Schwarz scale within about a factor two of the Pati–Salam scale, and then the left-handed pairs, at (alpha C/4 pi) M_SS, are the lightest
   new states — near the TeV scale.
OPEN: the loop function itself; the Green–Schwarz relaxation of the anomalous U(1)s; the two-U(1) coefficients (v^2 = (8/21) xi with both D-terms active)."""


def su4_pair(b_is_lepton):
    """Sum_A (T^A)_44 (T^A)_bb for the fundamental, b the lepton direction or a colour."""
    return (1.0 - 0.25) / 2 if b_is_lepton else (0.0 - 0.25) / 2

def scalar_masses_4L(xi_4L, xi_R4):
    v2 = (8 / 13) * xi_R4
    return {"squark": xi_4L - v2 - v2 * su4_pair(False), "slepton": xi_4L - v2 - v2 * su4_pair(True)}

def stable(xi_4L, xi_R4):
    m = scalar_masses_4L(xi_4L, xi_R4); return min(m.values()) >= -1e-12

def stability_bound():
    return 11 / 13

def ordering(a_Q, a_u, A_E=1.0, E=0.3):
    return stable(a_Q * A_E - E, a_u * A_E - E)

def radiative_asymmetry(M_SS_over_M_PS):
    return M_SS_over_M_PS ** 2

def report():
    return {"d_flat_Q_block": scalar_masses_4L(0.0, 1.0), "bound": stability_bound(), "pair_11_13": ordering(11, 13), "pair_13_11": ordering(13, 11),
            "edge_E0": scalar_masses_4L(11.0, 13.0), "M_SS_for_factor_5": (1 / 5) ** 0.5}
