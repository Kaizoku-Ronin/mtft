"""RT-01 (2026-10-02): the ratio-3 working background of FS-01 — does eps_d/eps_u = 3 survive its order-one factors, and what do its neutrinos do?

Background: SU(4) kept, Higgs gamma-charges (Delta_u, Delta_d) = (1, 3) (cheapest: a = (-1, -1, 1, 0, 2), c = (0, 0, -1, 2, 0), exceptional curve area
eps_D = 2 A_E / 3, 672 fermion zero modes), nu^c uncharged.

1. Leading order (STRUCTURAL, exact in the local model).  The up and down Higgs entries (L, A2) and (L, B) have the same exceptional charge (-6) at the same point
   over the same cusp; their gamma-localised zero modes differ only through Delta (the curve degrees do not enter localised modes).  Every factor in eps except
   the plane-wave momentum k ~ Delta phi' is common, so eps_d/eps_u = |Delta_d/Delta_u| = 3 exactly at leading order.
2. First correction (ESTIMATE).  A charge-Delta Higgs is a Gaussian of width^2 ~ A_E/(pi |Delta| phi') around the cusp; the jets fill the exceptional region
   (area eps_D).  The momentum is sampled over the overlap of the two Gaussians, so eps_Delta ~ K |Delta| phi' / sqrt(1 + 2 u_Delta) with
   u_Delta = |Delta| (eps_D/A_E) phi' and phi' = 3.2166/sqrt(r) (EPS-01), K in [4/3, 4 pi/3] for eps_D = 2 A_E/3.  The narrower down-type Higgs lowers the ratio.
   Fixing r by eps_u = 0.060 gives r ~ 4,800 - 49,000 and eps_d/eps_u = 2.84 - 2.95; fixing it by eps_d gives r ~ 7,500 - 81,000 and 2.87 - 2.96.  Against 2.30: a 25 per cent
   excess that the first correction does not remove and the estimate's precision cannot resolve.
3. Jet-sector neutrinos (STRUCTURAL).  nu^c is uncharged, so nu^c and its partners N are jets at the blown-up point.  In the one-orientation condensate vacuum
   (PB-01): nu^c has the condensate's charges, so its Majorana mass needs the conjugate orientation (zero); N's Majorana mass mu needs Sigma^2 and arises at tree
   level through heavy exchange; the B-L-conserving nu^c–N pairing M needs Sigma^dagger Sigma and only comes through the suppressed D-term route, like the Yukawas.
   That is the reverse of the inverse-seesaw ordering: integrating out N leaves nu^c with M^2/mu, and with M and m_D sharing the D-term suppression the light mass is
   m_D^2 mu / M^2 ~ v^2 / v_R — a type-I seesaw at the Pati–Salam scale, 27 MeV for the third family whatever the suppression, about 10^9 too heavy (NU-01).  The jet sector alone
   does not rescue the neutrinos in this background; FS-5 already showed a charged nu^c and ratio 3 cannot be had cheaply together.
OPEN: sector-dependent factors beyond the overlap model; a neutrino mechanism compatible with the ratio-3 background."""
import math
import numpy as np

PHI_COEFF = 3.2166                      # phi' = PHI_COEFF / sqrt(r)   (EPS-01, Arakelov cusp density from the frozen periods)
EPS_D = 2 / 3                            # exceptional area / A_E in the cheapest ratio-3 representative
K_BAND = (4 / 3, 4 * math.pi / 3)         # EPS-01's (2/3 ... 2 pi/3) scaled by eps_D / (A_E/3) = 2
TARGET = {"eps_u": math.sqrt(0.619 / 171.0), "eps_d": math.sqrt(0.055 / 2.86)}


def leading_ratio(delta_u=1, delta_d=3):
    return abs(delta_d) / abs(delta_u)

def eps_sector(r, delta, K):
    phip = PHI_COEFF / math.sqrt(r); u = abs(delta) * EPS_D * phip
    return K * abs(delta) * phip / math.sqrt(1 + 2 * u)

def _solve_r(target, delta, K):
    lo, hi = 1.0, 1e9
    for _ in range(200):
        mid = math.sqrt(lo * hi)
        if eps_sector(mid, delta, K) > target: lo = mid
        else: hi = mid
    return math.sqrt(lo * hi)

def corrected_ratio_window():
    out = {}
    for K in K_BAND:
        ru = _solve_r(TARGET["eps_u"], 1, K); rd = _solve_r(TARGET["eps_d"], 3, K)
        out[round(K, 4)] = {"r_from_eps_u": ru, "ratio_at_r_u": eps_sector(ru, 3, K) / eps_sector(ru, 1, K),
                            "r_from_eps_d": rd, "ratio_at_r_d": eps_sector(rd, 3, K) / eps_sector(rd, 1, K)}
    return out


# ------------------------------------------------------------------ the jet-sector neutrino structure
def neutral_matrix(mD, M, mu):
    """(nu_L, nu^c, N): Dirac m_D, nu^c–N pairing M, N Majorana mu, no nu^c Majorana mass (needs the conjugate orientation)."""
    return np.array([[0.0, mD, 0.0], [mD, 0.0, M], [0.0, M, mu]])

def lightest(m):
    return float(np.min(np.abs(np.linalg.eigvalsh(m))))

def jet_neutrino_estimate(theta, v=174.0, v_R=1.1e6, y=0.7, mu_over_vR=1.0):
    """m_D = theta y v and M = theta y v_R share the D-term suppression theta; mu ~ v_R is holomorphic and unsuppressed.  Returns the light mass in GeV."""
    return lightest(neutral_matrix(theta * y * v, theta * y * v_R, mu_over_vR * v_R))

def report():
    w = corrected_ratio_window()
    return {"leading_ratio": leading_ratio(), "targets": TARGET, "window": w,
            "nu_light_GeV": {t: jet_neutrino_estimate(t) for t in (0.3, 1.0, 3.0)}, "nu_seesaw_v2_over_vR_GeV": 0.7 ** 2 * 174.0 ** 2 / 1.1e6}
