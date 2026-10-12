"""TH-01 (2026-10-05): can thresholds reconcile the gauge couplings with the b–tau scale?

1. Left–right pairs cancel (EXACT).  Write the matching ratio as R = N/D with N = alpha_Y^-1 - (5/3) alpha_3^-1 and D = alpha_2^-1 - alpha_3^-1.  A vector-like (Dirac)
   (4,2,1) pair adds (b_3, b_2, b_Y) = (4/3, 8/3, 8/9), a (4bar,1,2) pair (4/3, 0, 32/9).  Per unit log below the matching scale they shift (N, D) by
   (+4/3, -4/3)/(2 pi) and (-4/3, +4/3)/(2 pi) respectively: equal numbers of left and right pairs at equal masses leave R untouched.  Only an asymmetry moves
   the matching, and light right-handed pairs move it the right way (down from 3.68).
2. What it takes (NUMERICAL).  At the b–tau crossing N = 56.6, D = 15.4.  Reaching the two-point design's R = 7/3 there needs n l = 29.2 (right-handed pairs times log);
   R = 1 needs 96.9 (n pairs a factor e^l below the scale; scalar partners would add half again).  With n = 19 that is a factor 4.7 (l = 1.54) for 7/3 and e^5.1
   (about 7 TeV) for R = 1.  Colour adjoints alone (Delta b_3 only) need Delta b_3 l = 195 for 7/3: implausible.
3. The design's own pairs (EXACT counts).  The gamma-neutral two-point design has one vector-like pair per block (jets: 3 vs 1 at P1, 1 vs 0 at P2) — n l = 29.3
   would need l = 29: impossible.  Matter blocks with gamma-charge c_L = c_R = c (bidoublet still gamma-neutral, Higgs still Phi-type) add 4c^2 non-chiral modes
   to each block: c = 2 gives 16 + 1 = 17 pairs per block, enough at a factor 5.6 if the right-handed ones sit below the left-handed ones.
4. The shape (STRUCTURAL).  Torus Kaluza–Klein towers enter through the DKL factor ln(A_E Im tau |eta(tau)|^4); for 143a1 Im tau |eta|^4 = 0.3527 (log -1.042).  The curve's
   towers bring the analogous spectral invariant of X_0(143) with ln A_X; the shape A_X/A_E enters the matching only through the difference of the two towers'
   beta coefficients — still to be computed per block.
OPEN: the masses of the vector-like pairs (which condensate gives the right-handed ones their mass); the towers' beta coefficients."""
from fractions import Fraction as F
import math


def pair_b(kind):
    """(b_3, b_2, b_Y) of one vector-like Dirac pair (fermions only)."""
    if kind == "421":
        return (F(4, 3) * F(1, 2) * 2, F(4, 3) * F(1, 2) * 4, F(4, 3) * (6 * F(1, 36) + 2 * F(1, 4)))
    if kind == "412":
        return (F(4, 3) * F(1, 2) * 2, F(0), F(4, 3) * (3 * F(4, 9) + 3 * F(1, 9) + 1))
    raise ValueError(kind)

def slopes(kind):
    b3, b2, bY = pair_b(kind)
    return (-(bY - F(5, 3) * b3), -(b2 - b3))                     # times 1/(2 pi) per unit log

def N_D(mu=1.1236e6):
    from mtft.matter.gauge_couplings import _inv
    a3, a2, aY = _inv("3", mu), _inv("2", mu), _inv("Y", mu)
    return aY - 5 / 3 * a3, a2 - a3

def required_nl(R_target, mu=1.1236e6):
    N, D = N_D(mu); x = (N - R_target * D) / (1 + R_target)
    return x / (float(F(12, 9)) / (2 * math.pi))

def adjoint_requirement(R_target, mu=1.1236e6):
    N, D = N_D(mu)                                                 # (N + 5/3 y)/(D + y) = R  with y = Delta b_3 l / (2 pi)
    y = (N - R_target * D) / (R_target - F(5, 3).__float__())
    return y * 2 * math.pi

def dkl_factor():
    import mpmath as mp
    from mtft.research.theta_torus import tau_143a1
    tau = mp.mpc(tau_143a1(30)); q = mp.exp(2j * mp.pi * tau)
    eta = mp.exp(2j * mp.pi * tau / 24) * mp.nprod(lambda n: 1 - q ** n, [1, mp.inf])
    return float(tau.imag * abs(eta) ** 4)

def vector_like_pairs(c=0):
    from mtft.matter.background import _jets
    return (min(_jets(-2), _jets(2)) + min(_jets(-1), _jets(1))) + 4 * c * c

def report():
    return {"slopes_421": slopes("421"), "slopes_412": slopes("412"), "N_D": N_D(), "nl_7/3": required_nl(7 / 3), "nl_1": required_nl(1.0),
            "factor_n19_7/3": math.exp(required_nl(7 / 3) / 19), "adjoint_7/3": adjoint_requirement(7 / 3), "dkl": dkl_factor(),
            "pairs": {c: vector_like_pairs(c) for c in (0, 1, 2)}}
