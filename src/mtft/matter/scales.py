"""SC-01 (2026-10-17): the scale ladder of the twisted compactification — what is fixed, by what, and what the measured background implies.

Two numbers set every scale of the product X0(143) x 143a1 with the even twist: the Kaluza–Klein scale M_KK = 2 pi / sqrt(A_E) and the
dimensionless 8D coupling at that scale, lambda = g_8^2 / A_E^2 = g_4^2 r with r = A_X / A_E.

1. The ladder (EXACT in PA-01's product approximation; one convention, C1, stated).  In units of M_KK:
     messenger (Scherk–Schwarz) scale   M_SS = 1 / (2 sqrt(Im tau)) = 0.4943          (PA-01's M sqrt(A_E) = 3.106 = 2 pi x this)
     messenger-pair threshold           2 M_SS = 0.9886 = the next boson level
     the saddles' gap                   |tau| / (2 sqrt(Im tau)) = 0.5629              (PG-01, 2026-10-16: both towers, both saddles)
     Wilson-line masses                 omega_hat = sqrt(96 T_eta / pi^6) sqrt(lambda) = 0.2175 sqrt(lambda), 0.4039 sqrt(lambda)   (Matter-05 closed form)
     Phi-adjoint masses                 m_hat_Phi = sqrt(8 S_eta / pi^6) sqrt(lambda) = 0.1324 sqrt(lambda)
   C1: Matter-05's closed forms with n = 8 real bosonic degrees of freedom per charged direction and sum_alpha |alpha(e)|^2 = 8 (Tr e^2 = 1/2);
   the repo leaves the group factor open, and it rescales lambda by an O(1) number.  Everything else is a ratio of lattice sums on 143a1.
2. The measured background (NUMERICAL).  FH-3 puts r = 230–280; GC-01 puts the Pati–Salam matching at 8.9e11–2.1e14 GeV, a lower bound on
   M_KK (the breaking is a 4D effect below the KK scale; reading R1).  With g_4^2 = 4 pi alpha_3(M_KK) = 0.31–0.37 there, lambda = 72–104.  Then
   omega_hat = 1.85–2.2 (light mode) and 3.4–4.1 (heavy), m_hat_Phi = 1.12–1.35: every one-loop 'modulus' is HEAVIER than the KK scale.  The
   light-modulus picture behind PA-01's path and PG-01's cosmological run does not hold on the measured background.
3. Dissipation (EXACT given 2 and PG-01's thresholds; width an ESTIMATE).  Both Wilson-line modes lie above the messenger-pair threshold
   (lambda > 20.7 and 6.0): the jet vacuum damps its Wilson line by linear decay into messenger pairs; broad parametric resonance (lambda < 1.3)
   and moduli trapping are excluded.  Width: Gamma / omega = (N g_4^2 / 8 pi) beta^3 with N the number of Dirac messengers at the threshold
   (25 pairs per matter entry times their multiplicity, 25 <= N <= 300) gives 0.2–3: at best a broad resonance, at worst no particle.
4. The expansion parameters (STRUCTURAL).  Matter-05's closed forms are leading order in A_E / R_C^2 = 56 pi / r = 0.63–0.77 — not small —
   and the 8D loop parameter at the KK scale, g_8^2 M_KK^4 / (4 pi)^4 = lambda / 16, is 4.5–6.5.  The one-loop numbers above are O(1) uncertain;
   the full heat kernel on X0(143) (the geodesic sum) is what would fix them.
5. What fixes M_KK itself (STRUCTURAL).  A lower bound from GC-01; an upper bound would need the volume stabilised (I.6), which the programme
   does not have (PV-01, MS-6, PG-4).  The absolute scale stays unassigned; every ratio is fixed.

Status: 1 EXACT (given C1), 2 NUMERICAL, 3 EXACT/ESTIMATE, 4–5 STRUCTURAL."""
import math

import numpy as np

R_WINDOW = (230.0, 280.0)                     # FH-3
M_PS_WINDOW = (8.9e11, 2.1e14)                # GC-01: allowed matching scales, GeV
N_DOF, GROUP_FACTOR = 8, 8.0                  # convention C1
LAMBDA_THRESHOLDS = {"messenger_pair": (20.66, 5.99), "jet_messenger": (5.17, 1.50), "broad_resonance": (1.32, 0.38)}   # PG-01 (2026-10-16), light/heavy mode


def tau():
    from mtft.research.theta_torus import tau_143a1
    return complex(tau_143a1())


def lattice_sums(eps=(0, 1), N=300):
    from mtft.matter.one_loop_adjoints import lattice_sums as ls
    S, T = ls(eps, N)
    return S, tuple(sorted(float(x) for x in np.linalg.eigvalsh(T)))


def ladder(eps=(0, 1)):
    """Every scale in units of M_KK; the lambda-dependent ones as coefficients of sqrt(lambda)."""
    t = tau(); y = t.imag; S, (T1, T2) = lattice_sums(eps)
    return {"M_SS": 1 / (2 * math.sqrt(y)), "pair_threshold": 1 / math.sqrt(y), "saddle_gap": abs(t) / (2 * math.sqrt(y)),
            "omega_hat_over_sqrt_lambda": (math.sqrt(96 * T1 / math.pi ** 6), math.sqrt(96 * T2 / math.pi ** 6)),
            "m_phi_hat_over_sqrt_lambda": math.sqrt(8 * S / math.pi ** 6), "S_eta": S, "T_eta": (T1, T2),
            "pa01_messenger_M_sqrtA": 2 * math.pi / (2 * math.sqrt(y))}


def g4_squared(mu):
    from mtft.matter.gauge_couplings import _inv
    return 4 * math.pi / _inv("3", mu)


def background(M_KK=None, r=None):
    """The measured background: lambda = g_4^2(M_KK) r and every scale, for M_KK in GeV and r (defaults: the GC-01 and FH-3 windows)."""
    L = ladder()
    Ms = (M_KK,) if M_KK else M_PS_WINDOW; rs = (r,) if r else R_WINDOW
    rows = []
    for M in Ms:
        for rr in rs:
            g2 = g4_squared(M); lam = g2 * rr
            rows.append({"M_KK_GeV": M, "r": rr, "g4_squared": g2, "lambda": lam, "sqrt_A_E_GeV^-1": 2 * math.pi / M,
                         "omega_hat": tuple(c * math.sqrt(lam) for c in L["omega_hat_over_sqrt_lambda"]),
                         "m_phi_hat": L["m_phi_hat_over_sqrt_lambda"] * math.sqrt(lam),
                         "expansion_parameter_56pi_over_r": 56 * math.pi / rr, "loop_parameter_lambda_over_16": lam / 16})
    return rows


def dissipation(lam, N_range=(25, 300)):
    """Which of PG-01's regimes the Wilson line is in at this lambda, and the estimated width of the light mode into messenger pairs."""
    L = ladder(); w = L["omega_hat_over_sqrt_lambda"][0] * math.sqrt(lam)
    above = lam > LAMBDA_THRESHOLDS["messenger_pair"][0]
    regime = ("linear decay into messenger pairs" if above else "below the messenger-pair threshold" if lam > LAMBDA_THRESHOLDS["broad_resonance"][0]
              else "broad parametric resonance")
    beta = math.sqrt(max(0.0, 1 - L["pair_threshold"] ** 2 / w ** 2)) if above else 0.0
    g2 = g4_squared(M_PS_WINDOW[0])
    return {"lambda": lam, "regime": regime, "both_modes_above_pair_threshold": lam > LAMBDA_THRESHOLDS["messenger_pair"][0],
            "broad_resonance_possible": lam < LAMBDA_THRESHOLDS["broad_resonance"][0], "velocity_beta": beta,
            "width_over_mass_light_mode": tuple(N * g2 / (8 * math.pi) * beta ** 3 for N in N_range)}


def report():
    L = ladder(); B = background(); lam = [b["lambda"] for b in B]
    return {"ladder": L, "background": B, "lambda_range": (min(lam), max(lam)),
            "omega_hat_range": ((min(b["omega_hat"][0] for b in B), max(b["omega_hat"][0] for b in B)), (min(b["omega_hat"][1] for b in B), max(b["omega_hat"][1] for b in B))),
            "m_phi_hat_range": (min(b["m_phi_hat"] for b in B), max(b["m_phi_hat"] for b in B)),
            "all_moduli_above_M_KK": all(b["omega_hat"][0] > 1 and b["m_phi_hat"] > 1 for b in B),
            "dissipation": {k: dissipation(v) for k, v in (("lambda_min", min(lam)), ("lambda_max", max(lam)))},
            "convention": "C1: n = 8, sum |alpha(e)|^2 = 8", "reading_R1": "M_KK >= the GC-01 matching scale (Pati–Salam breaks in 4D below the KK scale)"}
