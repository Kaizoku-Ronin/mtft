"""LB-01 (2026-10-04): how the lattice of 143a1 breaks the blown-up point's rotation, and whether jet families can carry an O(1) Yukawa.

1. Parity (EXACT).  j(143a1) = -262144/1859 (tau = 1/2 + 1.02327 i) is neither 0 nor 1728, so the lattice's only rotations are +-1: torus characters are conserved
   mod 2.  Of YA-3's surviving channels (n, a, k - l), the torus mismatch 1 + (k - l) + n is 6 for (3, 0, 2) and (4, 0, 1) and 7 for (4, 1, 2), which dies.  So the
   Higgs's slope drops out: Y3 ~ h0 c2 c0*, Y4 ~ h0 (c1 c0* + c2 c1*), both needing exactly six units of torus-rotation breaking.  The texture shape of YA-3 stands.
2. The lattice multipoles (NUMERICAL).  For the unit-area lattice, G4 = 1.4016 and G6 = 3.8595; the six shortest vectors sit at 0, +-64, +-116, 180 degrees, and their
   six-fold sum is 5.66 against 6 for a hexagonal lattice — 143a1's lattice is nearly as six-fold as a lattice can be.  The torus Green's function's regular part
   carries Re(G_2k w^2k)/(2k) (log sigma(w) = log w - sum G_2k w^2k/(2k)), so the six-unit breaking felt across a disc of radius r is |G6| r^6/6 at the edge and at
   most |G6| (2r)^6/6 between opposite edges.
3. The size at the blow-up scale (ESTIMATE).  The minimal family's exceptional area is A_E/3, a disc of radius 0.326: the six-unit breaking is 7.7e-4 (edge) to 0.049
   (opposite edges).  Massive messengers only lower it.  The heaviest jet entry also carries the condensate's tilt (Y4 ~ sqrt 2 t), so jet Yukawas are at most about
   0.05 x sqrt 2 t <~ 0.01 for t <~ 0.14.
4. Verdict (ESTIMATE).  The top Yukawa at the b–tau crossing scale (1.1e6 GeV) is 0.70 (SM one loop): a jet top falls short by a factor of about 73
   even at the bound; y_b there is 0.0104, right at it.  The first two
   families (and b, tau, at the edge) are within reach of jet couplings; the top is not.  In backgrounds of this class the third family cannot be a jet.
OPEN: a background whose third family comes from an unsuppressed (non-jet) channel."""
import itertools
import math
import mpmath as mp


def lattice_data(dps=30):
    from mtft.research.theta_torus import tau_143a1
    mp.mp.dps = dps
    tau = mp.mpc(tau_143a1(dps)); q = mp.exp(2j * mp.pi * tau)
    def E(k, terms=60):
        c = {4: 240, 6: -504}[k]
        return 1 + c * mp.fsum(mp.mpf(sum(d ** (k - 1) for d in range(1, n + 1) if n % d == 0)) * q ** n for n in range(1, terms))
    E4, E6 = E(4), E(6); s = 1 / mp.sqrt(tau.imag)
    vecs = sorted((s * (m + n * tau) for m, n in itertools.product(range(-3, 4), repeat=2) if (m, n) != (0, 0)), key=abs)[:6]
    return {"tau": complex(tau), "j": float(mp.re(1728 * E4 ** 3 / (E4 ** 3 - E6 ** 2))), "G4": float(mp.re(2 * mp.zeta(4) * E4 / s ** 4)),
            "G6": float(mp.re(2 * mp.zeta(6) * E6 / s ** 6)), "shell_angles": sorted(round(float(mp.arg(v)) * 180 / math.pi, 1) for v in vecs),
            "sixfold_sum": float(mp.re(mp.fsum(mp.exp(6j * mp.arg(v)) for v in vecs)))}

def blowup_radius(eps_over_AE=1 / 3):
    return math.sqrt(eps_over_AE / math.pi)

def breaking(k, eps_over_AE=1 / 3, across=False, data=None):
    d = data or lattice_data(); G = d["G4"] if k == 4 else d["G6"]
    r = blowup_radius(eps_over_AE) * (2 if across else 1)
    return abs(G) * r ** k / k

def parity_channels():
    from mtft.matter.yukawa_audit import curve_allowed
    return [(n, a, kl, 1 + kl + n) for (n, a, kl) in curve_allowed()]

def jet_yukawa_bound(t=0.14, eps_over_AE=1 / 3):
    from mtft.matter.yukawa_audit import coherent_amplitudes
    c = coherent_amplitudes(t)
    return breaking(6, eps_over_AE, across=True) * abs(c[1] * c[0] + c[2] * c[1])

def top_yukawa_at_crossing():
    """y_t at the b–tau equality scale, from the package's SM one-loop running (PS-01)."""
    from scipy.integrate import solve_ivp
    from mtft.matter.ps_breaking import INPUTS_MZ, _rge, b_tau_running
    I = dict(INPUTS_MZ); mu = b_tau_running()["mu_equal"]
    g1 = math.sqrt(4 * math.pi / (I["alpha_inv"] * (1 - I["sin2w"]))); g2 = math.sqrt(4 * math.pi / (I["alpha_inv"] * I["sin2w"])); g3 = math.sqrt(4 * math.pi * I["alpha_s"])
    r2 = math.sqrt(2) / I["v"]
    sol = solve_ivp(_rge, [math.log(I["m_Z"]), math.log(mu)], [g1, g2, g3, r2 * I["m_t"], r2 * I["m_b"], r2 * I["m_tau"]], rtol=1e-10, atol=1e-13)
    return {"mu": mu, "y_t": float(sol.y[3, -1]), "y_b": float(sol.y[4, -1])}

def report():
    d = lattice_data(); top = top_yukawa_at_crossing()
    return {"lattice": d, "radius": blowup_radius(), "delta6_edge": breaking(6, data=d), "delta6_across": breaking(6, across=True, data=d),
            "parity_channels": parity_channels(), "jet_bound": jet_yukawa_bound(), "top": top, "shortfall": top["y_t"] / jet_yukawa_bound()}
