"""VE-01 (2026-10-17): the vacuum energy of the twisted compactification — the programme's own terms, their signs, and the one point
where the one-loop part vanishes.  Addresses puzzle I.6 (a vacuum energy from the parent action); corrects SC-4 (CC-70).

1. The perfect square (EXACT, one loop, product approximation, through the a_1 term of the curve heat kernel).  For the U(8) multiplet
   with torus Wilson lines U (a flat connection on 143a1, holonomy U^w around the winding w), the Scherk–Schwarz supertrace factorises
   block by block as (1 - eta(w)) e^{i (theta_i - theta_j) . w}, so
       V(U) = -c sum_{eta(w) = -1} K(w) |Tr U^w|^2,     c > 0,  K(w) = |w|^-8 + (2 pi/3r) |w|^-6 + ... > 0  (summed: S0 (1 + 2.78/r)).
   V <= 0 everywhere; V = 0 exactly where Tr U^w = 0 on every odd winding, and that locus is the global MAXIMUM.  MS-5 makes the
   block coefficient flux-independent through a_1, so the square survives to relative order 1/r; flux enters at a_2 ((B t)^2).
2. The Z/2 table (EXACT).  For the Z/2 Wilson lines of the three Pati–Salam groups, N(bosons untwisted) - N(fermions untwisted) =
   (Tr W)^2 for all eight assignments: 64 (jet vacuum, W = 1), 16 (one SU(2) flipped), 0 (W = diag(-1_4, 1_2, 1_2) or its complement).
   PA-01's matter-untwisted background is the 0 point — the Bose–Fermi-balanced ('super no-scale', Kounnas–Partouche) point of the
   compactification, where the one-loop vacuum energy vanishes.  V(jet) = -64 c S0, V(saddles) = -32 c S0, V(untwisted) = 0.
3. Its stability (EXACT).  At the balanced point the 16 x 16 Hessian in the Cartan Wilson lines has signature (2 negative, 14 zero): the two
   negative directions are the relative U(1) T = diag(1_4, -1_2, -1_2) on the two torus cycles, with curvatures exactly minus the jet
   vacuum's Wilson-line masses (PG-4); of the 14 flat directions, 2 (the overall U(1)) are exactly flat and 12 (the SU(4), SU(2)_L, SU(2)_R
   Cartans and theta_L - theta_R) are maxima at quartic order, V = -c K u^4.  Downhill from 0 lie the 16-points and then the jet vacuum.
   So the balanced point is a vacuum only in a construction with no continuous torus Wilson lines at all.  T is anomalous (SU(2)_L^2 T = 12,
   SU(2)_R^2 T = -12, SU(4)^2 T = grav^2 T = 0): a Green–Schwarz mass M_St for it (I.2) lifts the two tachyons iff M_St exceeds the jet
   vacuum's Wilson-line masses, 1.85-2.2 and 3.4-4.1 M_KK on the measured background (SC-2); the twelve non-abelian directions have no such
   lift.  Atkin–Lehner left-right parity (the two-cusps-dparity background) freezes theta_4 at a half-period and leaves theta_L = theta_R
   continuous, whose projection on T is -4 b: still tachyonic.
4. The curve corrections to Matter-05's closed forms (EXACT at order a_1; the rest ESTIMATE) — CC-70.  The torus images weight the proper
   time at <t> = 0.1 A_E, so the curve's expansion parameter is <t>/R_C^2 = 0.1 x 56 pi/r, not 56 pi/r.  Term by term, relative to the leading
   closed form: a_1 (topological, MS-5's -chi/12 per block) gives +(pi/3)(a_1/N)(S1/S0)/r = +2.78/r = +1.0 to +1.2% to V and +1.5 to +1.8% to the
   masses at r = 230-280; a_2 (curvature^2, scalar-Laplacian value as the size) about 2 x 10^-4; the shortest closed geodesics of Gamma_0(143)
   (length 2 arccosh 2) 3 x 10^-5 per primitive class.  SC-4's 'O(1) uncertain' was wrong for the curve: the closed forms are good to about
   a per cent.  The O(1) uncertainty that remains is lambda/16 = 4.5-6.5, and no geodesic sum fixes that.
5. The curve's own KK rung (EXACT in the cusped metric; DIAGNOSTIC for the first discrete mode).  The hyperbolic Laplacian of Gamma_0(143)
   has no eigenvalue in (0, 1/4): Selberg's conjecture is a theorem at this level (Booker–Lee–Strömbergsson 2018, Gamma_1(N) for N <= 880), and
   1/4 is also the Breitenlohner–Freedman bound of the AdS_2 the curve is a quotient of.  So the curve's KK tower starts at 1/(2 R_C) =
   sqrt(56 pi/r)/(4 pi) M_KK = 0.063-0.070 M_KK, seven to eight times below M_SS (the ladder's lowest rung, missing from SC-01), with the first
   Maass form at sqrt(0.39)/R_C = 0.08 M_KK (SPEC-01, DIAGNOSTIC).  In the smooth blown-up metric the lowest eigenvalue is not computed.
6. The ledger (STRUCTURAL; numbers ESTIMATES under convention C1).  On the measured background the classical curvature term of the product
   compactification, +M_P^2/R_C^2 (negative internal curvature is positive 4D energy; Douglas–Kallosh's no-go is PV-01's runaway), is
   10^105-10^110 rho_Lambda; the one-loop energy of the jet vacuum, -(6 n A_C/pi^4 A_E^3) S0 x 32 complex directions = -3.7 to -4.5 M_KK^4, is
   10^95-10^105 rho_Lambda and 10^-5 to 10^-10 of the curvature term.  Nothing here is a vacuum (I.6 stays open); what is new is that the
   AdS the programme has — the twisted Casimir sink — is exactly cancelled at one point of the Wilson-line torus, and that point is its top.

Status: 1-3 EXACT (product approximation), 4 EXACT/ESTIMATE, 5 EXACT/DIAGNOSTIC, 6 STRUCTURAL."""
import math
from fractions import Fraction

import numpy as np

GROUPS = ("4", "L", "R")
RANKS = {"4": 4, "L": 2, "R": 2}
CHI = -24                                      # chi(X0(143))
A1_OVER_N = 2.0                                # MS-5: -chi/12 per block, any flux, twist or Wilson line
R_WINDOW = (230.0, 280.0)                      # FH-3


def tau():
    from mtft.research.theta_torus import tau_143a1
    return complex(tau_143a1())


def odd_windings(eps=(0, 1), N=40):
    """Odd windings (m, n) of the unit-area 143a1 lattice with their vectors w = (m + n tau)/sqrt(Im tau)."""
    t = tau(); s = 1 / math.sqrt(t.imag)
    mn, w = [], []
    for m in range(-N, N + 1):
        for n in range(-N, N + 1):
            if (eps[0] * m + eps[1] * n) % 2 == 1:
                z = s * (m + n * t); mn.append((m, n)); w.append((z.real, z.imag))
    return np.array(mn), np.array(w)


def lattice_sums(eps=(0, 1), N=40):
    """S_k = sum over odd windings of |w|^-(8 - 2k): S0 (leading), S1 (a_1 term), S2 (a_2 term)."""
    _, w = odd_windings(eps, N); r2 = (w ** 2).sum(1)
    return {"S0": float((r2 ** -4).sum()), "S1": float((r2 ** -3).sum()), "S2": float((r2 ** -2).sum())}


# ------------------------------------------------------------------ 1. the perfect square
def trace_holonomy(theta, mn):
    """Tr U^w = sum_i exp(2 pi i theta_i . (m, n)) for every winding; theta is (8, 2) in torus coordinates."""
    return np.exp(2j * math.pi * (np.asarray(theta) @ mn.T)).sum(0)


def potential(theta, eps=(0, 1), N=40):
    """V(U) = -sum_odd |w|^-8 |Tr U^w|^2 in PG-6's units (16 of PG-6's units per unit here along the relative U(1))."""
    mn, w = odd_windings(eps, N); K = ((w ** 2).sum(1)) ** -4
    return float(-(K * np.abs(trace_holonomy(theta, mn)) ** 2).sum())


def block_sum(theta, eps=(0, 1), N=40):
    """The same potential summed block by block, -sum_odd K(w) sum_{i,j} cos(2 pi (theta_i - theta_j) . w): the identity behind the square."""
    mn, w = odd_windings(eps, N); K = ((w ** 2).sum(1)) ** -4; th = np.asarray(theta)
    tot = 0.0
    for i in range(th.shape[0]):
        for j in range(th.shape[0]):
            tot += (K * np.cos(2 * math.pi * ((th[i] - th[j]) @ mn.T))).sum()
    return float(-tot)


def z2_theta(wilson, eps=(0, 1)):
    """Cartan Wilson lines (8, 2) for Z/2 assignments wilson = {'4': 0/1, 'L': 0/1, 'R': 0/1}: 1 means the half-period eps/2."""
    half = np.array(eps, float) / 2
    rows = []
    for g in GROUPS:
        rows += [half * wilson[g]] * RANKS[g]
    return np.array(rows)


# ------------------------------------------------------------------ 2. the Z/2 table
def z2_table():
    """All eight Z/2 Wilson lines: untwisted-boson and untwisted-fermion counts from Background.twist, against (Tr W)^2."""
    import itertools
    from mtft.matter.background import Background
    rows = []
    for bits in itertools.product((0, 1), repeat=3):
        wl = dict(zip(GROUPS, bits))
        bg = Background("ve", ("inf",), {"4": (0,), "L": (-3,), "R": (3,)}, wilson=wl)
        nB = sum(RANKS[g] ** 2 for g in GROUPS); nF = 0
        for i, j in itertools.combinations(GROUPS, 2):
            tw = bg.twist(i, j); n = 2 * RANKS[i] * RANKS[j]
            if tw["fermions"] == "untwisted": nF += n
            elif tw["bosons"] == "untwisted": nB += n
        trW = sum((-1) ** wl[g] * RANKS[g] for g in GROUPS)
        rows.append({"wilson": wl, "bosons_untwisted": nB, "fermions_untwisted": nF, "TrW_squared": trW ** 2,
                     "identity": nB - nF == trW ** 2, "V_over_cS0": -trW ** 2, "stability": bg.wilson_line_stability()})
    return rows


# ------------------------------------------------------------------ 3. stability of the balanced point
def hessian(theta, eps=(0, 1), N=40):
    """Analytic 16 x 16 Hessian of potential() in the Cartan coordinates theta_{i a}."""
    mn, w = odd_windings(eps, N); K = ((w ** 2).sum(1)) ** -4; th = np.asarray(theta); n = th.shape[0]
    phase = np.exp(2j * math.pi * (th @ mn.T))                       # (n, W)
    S = phase.sum(0)                                                   # (W,)
    H = np.zeros((2 * n, 2 * n))
    d = 2j * math.pi * mn.T                                            # (2, W): d phi_i / d theta_{i a}
    for i in range(n):
        for a in range(2):
            dS_ia = d[a] * phase[i]
            for j in range(n):
                for b in range(2):
                    dS_jb = d[b] * phase[j]
                    term = np.conj(dS_jb) * dS_ia
                    if i == j:
                        term = term + np.conj(S) * d[a] * d[b] * phase[i]
                    H[2 * i + a, 2 * j + b] = -(K * 2 * term.real).sum()
    return H


def signature(H, tol=1e-6):
    ev = np.linalg.eigvalsh(H); scale = max(1.0, np.abs(ev).max())
    return {"negative": int((ev < -tol * scale).sum()), "zero": int((np.abs(ev) <= tol * scale).sum()), "positive": int((ev > tol * scale).sum()),
            "eigenvalues": ev}


def balanced_point(eps=(0, 1)):
    return z2_theta({"4": 1, "L": 0, "R": 0}, eps)


def quartic_profile(direction, amplitudes=(0.02, 0.04, 0.08), eps=(0, 1), N=40):
    """V along theta = balanced + u * direction, with the fitted power of u (4 for the intra-block Cartan directions)."""
    base = balanced_point(eps); dirn = np.asarray(direction, float)
    vals = [potential(base + u * dirn, eps, N) - potential(base, eps, N) for u in amplitudes]
    p = math.log(abs(vals[-1] / vals[0])) / math.log(amplitudes[-1] / amplitudes[0]) if vals[0] else float("nan")
    return {"values": vals, "fitted_power": p, "sign": "maximum" if vals[-1] < 0 else "minimum"}


def su4_cartan_direction():
    d = np.zeros((8, 2)); d[0, 1] = 1; d[1, 1] = -1; return d


def t_direction():
    d = np.zeros((8, 2)); d[:4, 1] = 1; d[4:, 1] = -1; return d


def t_anomalies():
    """Anomaly coefficients of T = Q4 - QL - QR from PS-01's exact table."""
    from mtft.matter.ps_breaking import anomalies
    A = anomalies(); lin = lambda key: A[f"{key} Q4"] - A[f"{key} QL"] - A[f"{key} QR"]
    out = {k: lin(k) for k in ("SU4^2", "SU2L^2", "SU2R^2", "grav^2")}
    out["anomalous"] = any(v != 0 for v in out.values())
    return out


def dparity_projection(b=1.0):
    """Atkin–Lehner left–right parity: theta_4 -> -theta_4, theta_L <-> theta_R.  Its fixed locus keeps theta_L = theta_R = b free; the
    projection of that direction on T, sum_i eps_i delta theta_i = -4 b, is nonzero: the direction stays tachyonic."""
    return {"fixed_locus": "theta_4 in {0, half-periods}, theta_L = theta_R", "projection_on_T": -4 * b, "tachyon_removed": False}


def stueckelberg_criterion():
    """The two T tachyons are lifted by a Green–Schwarz mass M_St iff M_St > the jet vacuum's Wilson-line masses (PG-4: the two Hessians are
    negatives of each other); on the measured background those are SC-2's omega_hat."""
    from mtft.matter.scales import report
    R = report(); (lo1, hi1), (lo2, hi2) = R["omega_hat_range"]
    return {"M_St_over_M_KK_light_direction": (lo1, hi1), "M_St_over_M_KK_heavy_direction": (lo2, hi2),
            "lifts_nonabelian_directions": False, "status": "I.2 (the Green–Schwarz construction) is open"}


# ------------------------------------------------------------------ 4. curve corrections (CC-70)
def t_moments(eps=(0, 1), N=40):
    """<t^k> of the leading integrand t^-5 e^{-|w|^2/4t} summed over odd windings (A_E = 1): the torus images' proper-time weight."""
    _, w = odd_windings(eps, N); r2 = (w ** 2).sum(1)
    mom = lambda k: (math.gamma(4 - k) * (4 / r2) ** (4 - k)).sum() / (math.gamma(4) * (4 / r2) ** 4).sum()
    return {"<t>": float(mom(1)), "<t^2>": float(mom(2))}


def curve_corrections(r, eps=(0, 1), N=40, geodesic_classes=1):
    """Relative corrections to Matter-05's leading closed forms at shape r = A_C/A_E: a_1 (EXACT, MS-5), a_2 (ESTIMATE: scalar value per
    degree of freedom), shortest geodesics (per primitive class, trace 4 in Gamma_0(143)); for V and for the Wilson-line mass^2."""
    _, w = odd_windings(eps, N); r2 = (w ** 2).sum(1); R2 = r / (56 * math.pi)
    S0, S1, S2 = (r2 ** -4).sum(), (r2 ** -3).sum(), (r2 ** -2).sum()
    # K(w) = (A_C/4pi) Gamma(4) (4/|w|^2)^4 + a_1 Gamma(3) (4/|w|^2)^3 + a_2 Gamma(2) (4/|w|^2)^2 + ...; ratios of lattice sums follow
    a1 = (math.pi / 3) * A1_OVER_N * S1 / S0 / r
    a2 = (S2 / S0) / (1440 * R2 ** 2)                 # a_2 = (1/60 pi) int K^2 dA per degree of freedom (scalar Laplacian, curvature -1/R_C^2)
    # masses: the Hessian weights are w w^T K(w); the two eigen-directions of T_k = sum w w^T |w|^-(8-2k) coincide (Re tau = 1/2)
    T0 = (w[:, :, None] * w[:, None, :] * (r2 ** -4)[:, None, None]).sum(0); T1 = (w[:, :, None] * w[:, None, :] * (r2 ** -3)[:, None, None]).sum(0)
    e0, e1 = np.linalg.eigvalsh(T0), np.linalg.eigvalsh(T1)
    a1_mass = tuple(float((math.pi / 3) * A1_OVER_N * (b / a) / r) for a, b in zip(e0, e1))
    # geodesics: sum_gamma l/(2 sinh(l/2)) e^{-s/4} e^{-l^2/4s}/sqrt(4 pi s), s = t/R_C^2, against A_C/(4 pi t), weighted by t^-4 e^{-|w|^2/4t}
    ell = 2 * math.acosh(2.0)
    ts = np.logspace(-4, 2, 2000); lt = np.log(ts)
    P = (np.exp(-r2[None, :] / (4 * ts[:, None]))).sum(1) / ts
    lead = np.trapezoid(ts ** -3 * (r / (4 * math.pi * ts)) * P * ts, lt)
    s = ts / R2
    geo = np.trapezoid(ts ** -3 * (ell / (2 * math.sinh(ell / 2)) * np.exp(-s / 4) * np.exp(-ell ** 2 / (4 * s)) / np.sqrt(4 * math.pi * s)) * P * ts, lt)
    return {"r": r, "expansion_parameter_<t>/R_C^2": t_moments(eps, N)["<t>"] / R2, "56pi_over_r": 56 * math.pi / r,
            "a1_relative_V": float(a1), "a1_relative_mass2": a1_mass, "a2_relative_V_estimate": float(a2),
            "geodesic_relative_V_per_class": float(geodesic_classes * geo / lead), "systole": ell}


# ------------------------------------------------------------------ 5. the curve's KK rung
def curve_kk_rung(r, lambda_1=0.39):
    """1/(2 R_C) and sqrt(lambda_1)/R_C in units of M_KK = 2 pi/sqrt(A_E), with R_C^2 = r A_E/(56 pi)."""
    inv_RC = math.sqrt(56 * math.pi / r) / (2 * math.pi)
    from mtft.matter.scales import ladder
    M_SS = ladder()["M_SS"]
    return {"continuum_threshold_1/(2R_C)": inv_RC / 2, "first_maass_sqrt(lambda_1)/R_C": math.sqrt(lambda_1) * inv_RC,
            "M_SS": M_SS, "M_SS_over_threshold": M_SS / (inv_RC / 2), "selberg_bound": "theorem at level 143 (Booker–Lee–Strömbergsson 2018)",
            "BF_bound_AdS2": 0.25}


# ------------------------------------------------------------------ 6. the ledger
M_P_REDUCED = 2.435e18                                  # GeV
RHO_LAMBDA = 2.52e-47                                   # GeV^4 (Omega_Lambda = 0.685, H0 = 67.4)


def energy_ledger(n=8, eps=(0, 1)):
    from mtft.matter.scales import M_PS_WINDOW
    S0 = lattice_sums(eps)["S0"]; rows = []
    for M in M_PS_WINDOW:
        for r in R_WINDOW:
            AE = (2 * math.pi / M) ** 2; RC2 = r * AE / (56 * math.pi)
            V_curv = M_P_REDUCED ** 2 / RC2
            V_loop = -(6 * n * r / math.pi ** 4) * S0 / AE ** 2 * 32
            rows.append({"M_KK_GeV": M, "r": r, "V_curvature_GeV4": V_curv, "V_loop_jet_GeV4": V_loop, "V_loop_jet_in_M_KK4": V_loop / M ** 4,
                         "curvature_over_loop": V_curv / abs(V_loop), "log10_curvature_over_rho_Lambda": math.log10(V_curv / RHO_LAMBDA),
                         "log10_loop_over_rho_Lambda": math.log10(abs(V_loop) / RHO_LAMBDA)})
    return rows


def report(N=40):
    base = balanced_point(); jet = z2_theta({"4": 0, "L": 0, "R": 0})
    sig = signature(hessian(base, N=N))
    return {"lattice_sums": lattice_sums(N=N), "z2_table": z2_table(),
            "V_jet_over_cS0": potential(jet, N=N) / lattice_sums(N=N)["S0"], "V_balanced": potential(base, N=N),
            "square_vs_block_sum_at_random": None,
            "hessian_signature_at_balanced_point": {k: sig[k] for k in ("negative", "zero", "positive")},
            "hessian_signature_at_jet": {k: v for k, v in signature(hessian(jet, N=N)).items() if k != "eigenvalues"},
            "quartic_su4_cartan": quartic_profile(su4_cartan_direction(), N=N), "t_anomalies": t_anomalies(), "dparity": dparity_projection(),
            "stueckelberg": stueckelberg_criterion(), "curve_corrections": [curve_corrections(r) for r in R_WINDOW],
            "kk_rung": [curve_kk_rung(r) for r in R_WINDOW], "ledger": energy_ledger(), "t_moments": t_moments()}
