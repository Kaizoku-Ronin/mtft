"""QL-01 (2026-10-06): the quark–lepton split — the b–tau bound closes the dressing route, and names what the vertex must supply.

1. The b–tau bound (EXACT).  Kinetic dressing only lowers a Yukawa (Z <= 1).  With tree-level b–tau equality (the SU(4) index is contracted alike for the quark and the
   lepton component) and d^c undressed (CD-1), explaining (m_mu/m_tau)/(m_s/m_b) = 3.09 by dressing forces y_tau/y_b = Z_3^e <= 0.32 at the matching scale.  GC-01's
   running gives y_tau/y_b = 1.40 (8.9e11 GeV) to 1.53 (2.1e14): the dressing route misses by a factor >= 4.4 (>= 3.1 even with b–tau unification).  CD-4 is retracted
   (CC-52); GT-4's reversal and DP-3/DP-4's conflict are moot.  The Wilson-line split of GT-4 hits the ratio in the massless-gaugino limit (z_e = 2.7-3.3 at r ~ 100,
   from sqrt(0.72/0.0024)/(5.2-6.5)) but there every lepton weight vanishes, so it fails the same bound.
2. The targets (NUMERICAL).  The lepton/quark Yukawa ratio per family: f_3 = 1.40-1.53 (SM running to the matching scales), f_2 = 4.33-4.72, f_1 = 0.42-0.45.  Relative
   to the third family, independent of running: f_2/f_3 = 3.09, f_1/f_3 = 0.30 — the Georgi–Jarlskog pattern (3, 1/3, 1) within 3% and 10%.
3. What the vertex must supply (STRUCTURAL).  In the down sector the only SU(4)-breaking at the Yukawa vertex is the condensate in the nu^c direction, whose insertion acts
   on the lepton component alone: f_j = 1 + kappa_j with kappa = (-0.58, +3.3, +0.40) at 8.9e11 GeV (-0.55, +3.7, +0.53 at 2.1e14).  The two families at P1 need opposite signs — the insertion must be odd
   under the jet parity (x against w), the role the gamma-flux already plays (DD-01) — and the third family a mild positive one.  A kinetic weight cannot do this; a
   condensate factor inside the Yukawa amplitude (XQ-01, YA-01) can, and it is the only candidate left.
OPEN: the condensate insertion's sign and size at the x-jet, w-jet and constant jet; whether f_1 f_2 = 1.9 (not 1) is reproduced."""
import math

MASS = {"e_over_mu": 0.004836, "mu_over_tau": 0.0595, "d_over_s": 0.0503, "s_over_b": 0.0595 / 3.092}


def z_e():
    return MASS["mu_over_tau"] / MASS["s_over_b"]

def first_over_third():
    return (MASS["e_over_mu"] * MASS["mu_over_tau"]) / (MASS["d_over_s"] * MASS["s_over_b"])

def dressing_bound():
    """Largest y_tau/y_b that dressing allows once z_e is explained by dressing."""
    return 1 / z_e()

def tau_over_b(mu):
    from mtft.matter.gauge_couplings import b_tau_at
    return 1 / b_tau_at(mu)

def shortfall(mu):
    return tau_over_b(mu) * z_e()

def factors(mu):
    f3 = tau_over_b(mu); return (first_over_third() * f3, z_e() * f3, f3)

def insertions(mu):
    return tuple(f - 1 for f in factors(mu))

def wilson_line_limit(share_home=0.7239, share_jets=0.002367, c_ratio=(5.2, 6.5)):
    """z_e in the massless-W_R-gaugino limit: the mixing ratio m_3^e/m_2^e = (c_2/c_1) sqrt(share_home/share_jets)."""
    s = math.sqrt(share_home / share_jets); return (s / c_ratio[1], s / c_ratio[0])

def report():
    return {"z_e": z_e(), "f1_over_f3": first_over_third(), "bound": dressing_bound(), "tau_over_b": [tau_over_b(m) for m in (8.9e11, 2.1e14)],
            "shortfall": [shortfall(m) for m in (8.9e11, 2.1e14)], "factors": [factors(m) for m in (8.9e11, 2.1e14)], "insertions": insertions(8.9e11),
            "wilson_line_limit": wilson_line_limit()}
