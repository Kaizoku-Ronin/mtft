"""GT-01 (2026-10-05): the two gaugino towers — the semi-Dirac-like split exists, with the right shape and the wrong strength.

1. Same mass, different homes (EXACT).  Each tower's lowest modes are zero modes along the curve (Dolbeault forms) times the lowest twisted torus mode, so towers
   twisted alike share their lowest mass M; they differ only in where those modes live.  The leptoquark tower carries no flux: its curve modes are the constant and
   the thirteen differentials, whose density in the canonical (Arakelov) metric is the Bergman kernel divided by itself — uniform.  It dresses every point equally.
2. Where the W_R tower lives (NUMERICAL).  The W_R gauginos sit in the (B, A2) block, which carries the SU(2)_R-breaking flux.  If that flux sits holomorphically at the
   blown-up cusp (L = O(k cusp)), the tower's curve modes are the differentials vanishing there to order k, and their share of the canonical density falls like
   |q|^2k at the jets — k = 2: 0.0024 at |q| = 0.17, 0.055 at |q| = 0.39 — while staying 0.72 at the points over tau/2 and 0.20 at the Fricke point.  The split has the
   shape e^c needs: the W_R tower avoids the jets and reaches the top's homes.
3. The strength bound (EXACT).  A subspace's Bergman density never exceeds the full one, so with equal masses and couplings the W_R tower dresses no point more strongly
   than the leptoquark tower: z_e <= 1/Z_3^u.  Reaching z_e = 3.09 forces Z_3^u <= 0.32, and z_u = 0.25 then forces m_2^u/M = 12.3 — every u^c family mostly gaugino.
4. What the reversal needs (STRUCTURAL).  A lighter W_R tower, M_LQ/M_WR >= 3.9 (4.6 with the k = 2 share at the homes) — for instance different Z/2 Wilson lines on B and
   A2, leaving the W_R gauginos untwisted.  That breaks SU(2)_R at the compactification scale and changes GC-01's matching.  Shape is available; strength is the obstacle.
OPEN: the Wilson-line split and its effect on the gauge matching; the holomorphic structure of the SU(2)_R flux."""
import math
import numpy as np


def _setup():
    from mtft.matter.special_points import density_setup
    A, Yi = density_setup(); return A, Yi, np.linalg.inv(Yi)

def share(tau, k, data=None):
    """Bergman density of the differentials vanishing to order k at the cusp, as a fraction of the full Bergman density at tau."""
    A, Yi, G = data or _setup()
    w = (np.exp(2j * np.pi * tau) ** np.arange(A.shape[0])) @ A
    if k == 0: return 1.0
    _, _, vh = np.linalg.svd(A[1:k + 1, :]); V = vh[k:].conj().T
    phi = V.T @ w; GV = V.conj().T @ G @ V
    return float(np.real(phi.conj() @ np.linalg.inv(GV.conj()) @ phi) / np.real(np.conj(w) @ Yi @ w))

def tau_at(qabs):
    return complex(0, -math.log(qabs) / (2 * math.pi))

def lepton_bound(Z3u):
    return 1 / Z3u

def reversal_requirements(z_e=3.092, z_u=0.25):
    Z3u = 1 / z_e; Z2u = z_u * Z3u
    return Z3u, math.sqrt(1 / Z2u ** 2 - 1)

def mass_ratio_needed(share_home=1.0, m3u_over_M=0.75, Z3e=1 / 3.092):
    m3e = math.sqrt(1 / Z3e ** 2 - 1)
    return m3e / (m3u_over_M * math.sqrt(share_home))

def report():
    d = _setup()
    return {"k2_jets": (share(tau_at(0.17), 2, d), share(tau_at(0.39), 2, d)), "k2_homes": (share(complex(0.18363, 0.03346), 2, d), share(complex(0, 1 / math.sqrt(143)), 2, d)),
            "reversal": reversal_requirements(), "mass_ratio": (mass_ratio_needed(1.0), mass_ratio_needed(share(complex(0.18363, 0.03346), 2, d)))}
