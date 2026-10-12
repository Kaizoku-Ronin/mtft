"""PS-01 (2026-10-02): the Pati–Salam branch of the minimal Scherk–Schwarz flux — how it breaks, what is anomalous, and the first comparison
with data (b–tau).  Builds on `mtft.matter.pati_salam` (first look) and `mtft.surface.blowup` (lattice U + <-8> + <-1>, classes (a, b, c, e)).

1. Hypercharge (EXACT).  With T_3R(B) = +1/2, T_3R(A2) = -1/2, (B-L)(colour) = 1/3, (B-L)(A1) = -1, the stack hypercharges Y = T_3R + (B-L)/2 are
   (c 1/6, A1 -1/2, L 0, B 1/2, A2 -1/2) — exactly the PV-02 assignment (y_L = 0) — and the blocks branch as Pati–Salam requires:
   (4,2,1) -> (3,2)_{1/6} + (1,2)_{-1/2};  (4bar,1,2) -> (3bar,1)_{1/3} + (3bar,1)_{-2/3} + (1,1)_1 + (1,1)_0;  (1,2,2) -> (1,2)_{+-1/2}.
2. Anomalies of the net chiral content (EXACT): 3 (4, 2bar_L) with stack charges (Q4, QL) = (1, -1); 3 (4bar, 2_R) with (Q4, QR) = (-1, 1);
   6 (2bar_L, 2_R) with (QL, QR) = (-1, 1).  SU(4)^3, SU(4)^2 Q4 and the Witten counts (24 doublets each) cancel; the stack U(1)s carry mixed anomalies
   (SU(2)_L^2 QL = -12, SU(2)_R^2 QR = 12, SU(2)_L^2 Q4 = 6, ...) and must become massive through Green–Schwarz axions; Y is a nonabelian generator and
   anomaly-free.  The axion completion is the Matter-02/03 photon gate.
3. Wilson lines on 143a1 (EXACT): flat, so they keep every index and every slope, and keep the rank.  N_A1 != N_c breaks SU(4) -> SU(3) x U(1)_{B-L};
   N_A2 != N_B breaks SU(2)_R -> U(1)_{3R}; the Matter-03 Higgs pair (both orientations of the (L, B) block) survives iff N_L = N_B.  They cannot give the
   Standard Model: one more rank must go.
4. The Kahler FI (EXACT D-term algebra).  Slopes of the classes c (0,0,0,0), L (1,0,0,-3), B = A2 (-1,0,0,3), A1 (0,0,0,0) in
   omega = A_X f1 + A_E f2 + kappa gamma - eps E are (0, delta, -delta, 0, -delta) with delta = A_E - 3 eps (kappa unused).  Stack FI terms
   xi_k = n_k mu_k = (0, 2 delta, -delta, 0, -delta) are cancelled by the (1,2,2) block alone (its B and A2 components): delta drives ELECTROWEAK breaking,
   not Pati–Salam breaking.
5. The gamma direction (EXACT on the product).  gamma (gamma^2 = -8, the primitive part of the graph of phi: X0(143) -> 143a1, degree 4) is invisible to
   the canonical class, so gamma-flux leaves every chirality unchanged, and it enters slopes as -8 c kappa.  The Pati–Salam-symmetric flux
   (gamma-charge 1 on the quark–lepton group c + A1, 0 on L, -2 on B + A2; traceless) at delta = 0 gives xi = (4-group: -32 kappa, R-group: +32 kappa):
   the (4bar,1,2) block alone cancels it — a Pati–Salam-breaking condensate driven by kappa, separate from the electroweak delta.  O(gamma) restricted to
   each fibre {x} x E is the flat bundle of the point phi(x); so a block with gamma-charge c and exceptional charge e has chi = -4c^2 - e(e-1)/2: its
   extra 4c^2 modes localise where c phi(x) is the fibre's twist — for bosons at phi^-1(E[c]); for c = 1 that is exactly the four cusps of X0(143)
   (all cusps are rational, Manin–Drinfeld makes their images torsion, 143a1 has no rational torsion, and deg phi = 4).  Under the Scherk–Schwarz twist
   the fermionic ones sit at phi^-1(P_eta) instead.  Cost: 8c^2 vector-like modes per rank pair.  The cusp infinity carries the blow-up point: that
   interaction, the nonabelian D-terms of the condensate and the F-terms are OPEN.
6. b–tau (NUMERICAL, SM one loop).  Quark and lepton blocks are related by the unbroken SU(4) of the flux, so Y_d = Y_e (and Y_u = Y_nu) at the matching
   scale.  Running the SM from m_Z, y_b = y_tau at mu = 1.1e6 GeV (0.6e6 - 2.3e6 for m_b(m_Z) = 2.79 - 2.93 GeV) — at the PeV scale where K_L -> mu e
   already bounds the Pati–Salam vector leptoquark (2100 TeV for maximal mixing; Smirnov, arXiv:1210.3697).  At that scale m_s/m_mu ~ 0.33 and
   m_d/m_e ~ 3.4: the Georgi–Jarlskog pattern, which the minimal relation does not give.  And Y_nu = Y_u with M_R below 1e6 GeV gives a type-I seesaw mass
   m_t^2/M_R ~ 0.03 GeV, about 10^8 times too heavy: an inverse seesaw (the spectrum has extra neutral vector-like states) or suppressed Dirac couplings
   is required.  Two-loop running, thresholds and the extra light states are not included."""
from fractions import Fraction
import math
import itertools

# ------------------------------------------------------------------ 1. hypercharge and branching
T3R = {"c": Fraction(0), "A1": Fraction(0), "L": Fraction(0), "B": Fraction(1, 2), "A2": Fraction(-1, 2)}
BML = {"c": Fraction(1, 3), "A1": Fraction(-1), "L": Fraction(0), "B": Fraction(0), "A2": Fraction(0)}
RANK = {"c": 3, "A1": 1, "L": 2, "B": 1, "A2": 1}

def stack_hypercharge():
    return {k: T3R[k] + BML[k] / 2 for k in RANK}

def block_hypercharge(i, j):
    """Hom(j -> i) = V_i (x) V_j^*: Y_i - Y_j."""
    y = stack_hypercharge(); return y[i] - y[j]

def branching():
    """The SM content of the three Pati–Salam blocks, from stack hypercharges (EXACT)."""
    y = lambda i, j: block_hypercharge(i, j)
    return {"(4,2,1)": {"Q (c<-L)": y("c", "L"), "l (A1<-L)": y("A1", "L")},
            "(4bar,1,2)": {"d^c (B<-c)": y("B", "c"), "u^c (A2<-c)": y("A2", "c"), "e^c (B<-A1)": y("B", "A1"), "nu^c (A2<-A1)": y("A2", "A1")},
            "(1,2,2)": {"H (B<-L)": y("B", "L"), "H' (A2<-L)": y("A2", "L")},
            "pv02_consistent": stack_hypercharge() == {"c": Fraction(1, 6), "A1": Fraction(-1, 2), "L": Fraction(0), "B": Fraction(1, 2), "A2": Fraction(-1, 2)}}


# ------------------------------------------------------------------ 2. anomalies
CONTENT = [  # (mult, dims {SU4, SU2L, SU2R}, conj flags, U(1) charges)
    {"mult": 3, "rep": {"SU4": 4, "SU2L": 2, "SU2R": 1}, "su4_sign": +1, "Q": {"Q4": 1, "QL": -1, "QR": 0}, "name": "(4, 2bar_L)"},
    {"mult": 3, "rep": {"SU4": 4, "SU2L": 1, "SU2R": 2}, "su4_sign": -1, "Q": {"Q4": -1, "QL": 0, "QR": 1}, "name": "(4bar, 2_R)"},
    {"mult": 6, "rep": {"SU4": 1, "SU2L": 2, "SU2R": 2}, "su4_sign": 0, "Q": {"Q4": 0, "QL": -1, "QR": 1}, "name": "(2bar_L, 2_R)"}]

def anomalies():
    """Exact anomaly coefficients of the net chiral content (Dynkin index T(fund) = 1/2, A(4) = 1)."""
    out = {}
    dim = lambda f: f["rep"]["SU4"] * f["rep"]["SU2L"] * f["rep"]["SU2R"]
    out["SU4^3"] = sum(f["mult"] * f["su4_sign"] * f["rep"]["SU2L"] * f["rep"]["SU2R"] for f in CONTENT if f["rep"]["SU4"] == 4)
    for g, other in (("SU2L", ("SU4", "SU2R")), ("SU2R", ("SU4", "SU2L"))):
        out[f"Witten_{g}_doublets"] = sum(f["mult"] * f["rep"][other[0]] * f["rep"][other[1]] for f in CONTENT if f["rep"][g] == 2)
    for U in ("Q4", "QL", "QR"):
        out[f"SU4^2 {U}"] = sum(Fraction(1, 2) * f["mult"] * f["rep"]["SU2L"] * f["rep"]["SU2R"] * f["Q"][U] for f in CONTENT if f["rep"]["SU4"] == 4)
        out[f"SU2L^2 {U}"] = sum(Fraction(1, 2) * f["mult"] * f["rep"]["SU4"] * f["rep"]["SU2R"] * f["Q"][U] for f in CONTENT if f["rep"]["SU2L"] == 2)
        out[f"SU2R^2 {U}"] = sum(Fraction(1, 2) * f["mult"] * f["rep"]["SU4"] * f["rep"]["SU2L"] * f["Q"][U] for f in CONTENT if f["rep"]["SU2R"] == 2)
        out[f"grav^2 {U}"] = sum(f["mult"] * dim(f) * f["Q"][U] for f in CONTENT)
    for U, V, W in itertools.combinations_with_replacement(("Q4", "QL", "QR"), 3):
        out[f"{U} {V} {W}"] = sum(f["mult"] * dim(f) * f["Q"][U] * f["Q"][V] * f["Q"][W] for f in CONTENT)
    return out


# ------------------------------------------------------------------ 3. Wilson lines
def wilson_lines(N):
    """N: stack -> label of its flat line bundle on 143a1.  Returns the unbroken nonabelian group, the twisted blocks and the Higgs-pair verdict."""
    groups = []
    if N["c"] == N["A1"]: groups.append("SU(4)")
    else: groups += ["SU(3)", "U(1)_{B-L}"]
    groups.append("SU(2)_L")
    if N["B"] == N["A2"]: groups.append("SU(2)_R")
    else: groups.append("U(1)_{3R}")
    twisted = sorted({tuple(sorted((i, j))) for i in N for j in N if i < j and N[i] != N[j]})
    return {"unbroken": groups, "rank": 5, "twisted_blocks": twisted, "higgs_pair_survives": N["L"] == N["B"],
            "reaches_SM": False, "reason": "flat bundles commute with the Cartan torus: the rank stays 5 (Standard Model rank 4)"}


# ------------------------------------------------------------------ 4./5. slopes, FI terms, single-block D-flat vacua
CLASSES = {"c": (0, 0, 0, 0), "A1": (0, 0, 0, 0), "L": (1, 0, 0, -3), "B": (-1, 0, 0, 3), "A2": (-1, 0, 0, 3)}

def slopes(delta=0.0, kappa=0.0, gamma_charge=None, A_E=3.0):
    """mu = a A_E + b A_X - 8 c kappa + e eps with eps = (A_E - delta)/3 (A_X does not enter: b = 0).  gamma_charge adds c to a stack."""
    eps = (A_E - delta) / 3; g = gamma_charge or {}
    return {k: a * A_E - 8 * (c + g.get(k, 0)) * kappa + e * eps for k, (a, b, c, e) in CLASSES.items()}

def fi_terms(**kw):
    mu = slopes(**kw); return {k: RANK[k] * mu[k] for k in mu}

def single_block_vacua(xi, tol=1e-12):
    """Blocks Hom(j -> i) whose scalar alone cancels every stack FI term: xi supported on {i, j} with xi_i = -xi_j > 0 (group-level for coinciding stacks)."""
    groups = {"4": ("c", "A1"), "L": ("L",), "R": ("B", "A2")}
    gx = {g: sum(xi[s] for s in st) for g, st in groups.items()}
    support = [g for g, v in gx.items() if abs(v) > tol]
    same_within = all(abs(xi[s] / RANK[s] - xi[st[0]] / RANK[st[0]]) < tol for st in groups.values() for s in st)
    out = []
    if len(support) == 2 and abs(gx[support[0]] + gx[support[1]]) < tol and same_within:
        i, j = (support[0], support[1]) if gx[support[0]] > 0 else (support[1], support[0])
        out.append({"Hom": f"{j} -> {i}", "rep": {("4", "R"): "(4,1,2)/(4bar,1,2)", ("R", "4"): "(4,1,2)/(4bar,1,2)", ("4", "L"): "(4,2,1)", ("L", "4"): "(4,2,1)",
                                                   ("L", "R"): "(1,2,2)", ("R", "L"): "(1,2,2)"}[(i, j)], "vev_squared": gx[i]})
    return {"group_FI": gx, "vacua": out}

PS_GAMMA_FLUX = {"c": 1, "A1": 1, "L": 0, "B": -2, "A2": -2}      # Pati–Salam symmetric, traceless: 3 + 1 - 2 - 2 = 0

def gamma_facts():
    """gamma^2 = -8; chirality neutral; chi of a block with gamma-charge c, exceptional charge e (b = 0): -4c^2 - e(e-1)/2; boson localisation set for
    c = 1 = phi^-1(O) = the four cusps (rational cusps, Manin–Drinfeld, E(Q)_tors = 0 checked by point counts, deg phi = 4)."""
    from mtft.surface.blowup import dot, K_TILDE, chi
    from math import gcd
    a1, a2, a3, a4, a6 = (0, -1, 1, -1, -2)
    def count(p): return 1 + sum(1 for x in range(p) for y in range(p) if (y * y + a1 * x * y + a3 * y - (x ** 3 + a2 * x * x + a4 * x + a6)) % p == 0)
    tors_bound = 0
    for p in (3, 5, 7, 17, 19, 23):
        tors_bound = gcd(tors_bound, count(p))
    chis = {(c, e): chi((0, 0, c, e)) for c in (0, 1, 2) for e in (-3, 0, 3)}
    return {"gamma_sq": dot((0, 0, 1, 0), (0, 0, 1, 0)), "K_dot_gamma": dot((0, 0, 1, 0), K_TILDE), "chi_table": chis,
            "chi_formula_holds": all(chis[(c, e)] == -4 * c * c - Fraction(e * (e - 1), 2) for (c, e) in chis),
            "E_Q_torsion_order_divides": tors_bound, "cusps": 4, "modular_degree": 4, "cusps_fill_fibre_over_O": tors_bound == 1}


# ------------------------------------------------------------------ 6. b–tau against data (SM one loop)
def _rge(t, y):
    g1, g2, g3, yt, yb, ytau = y; k = 1 / (16 * math.pi ** 2); Y2 = 3 * yt ** 2 + 3 * yb ** 2 + ytau ** 2
    return [k * 41 / 6 * g1 ** 3, -k * 19 / 6 * g2 ** 3, -k * 7 * g3 ** 3,
            k * yt * (1.5 * (yt ** 2 - yb ** 2) + Y2 - 17 / 12 * g1 ** 2 - 9 / 4 * g2 ** 2 - 8 * g3 ** 2),
            k * yb * (1.5 * (yb ** 2 - yt ** 2) + Y2 - 5 / 12 * g1 ** 2 - 9 / 4 * g2 ** 2 - 8 * g3 ** 2),
            k * ytau * (1.5 * ytau ** 2 + Y2 - 15 / 4 * g1 ** 2 - 9 / 4 * g2 ** 2)]

INPUTS_MZ = {"alpha_inv": 127.95, "sin2w": 0.23122, "alpha_s": 0.1180, "m_t": 171.0, "m_b": 2.86, "m_tau": 1.7463, "m_s": 0.055, "m_mu": 0.1027,
             "m_d": 0.0027, "m_e": 0.000487, "v": 246.22, "m_Z": 91.1876}

def b_tau_running(m_b=None, alpha_s=None, mu_max=1e17):
    """One-loop SM running of y_b / y_tau from m_Z; returns the equality scale, the ratio curve and the light-generation ratios there."""
    from scipy.integrate import solve_ivp
    from scipy.optimize import brentq
    I = dict(INPUTS_MZ); I["m_b"] = m_b or I["m_b"]; I["alpha_s"] = alpha_s or I["alpha_s"]
    g1 = math.sqrt(4 * math.pi / (I["alpha_inv"] * (1 - I["sin2w"]))); g2 = math.sqrt(4 * math.pi / (I["alpha_inv"] * I["sin2w"])); g3 = math.sqrt(4 * math.pi * I["alpha_s"])
    r2 = math.sqrt(2) / I["v"]
    sol = solve_ivp(_rge, [math.log(I["m_Z"]), math.log(mu_max)], [g1, g2, g3, r2 * I["m_t"], r2 * I["m_b"], r2 * I["m_tau"]], dense_output=True, rtol=1e-10, atol=1e-13)
    R = lambda t: sol.sol(t)[4] / sol.sol(t)[5]
    t_star = brentq(lambda t: R(t) - 1, math.log(I["m_Z"]), math.log(mu_max / 10)); mu_star = math.exp(t_star)
    shrink = R(t_star) / R(math.log(I["m_Z"]))                     # the common QCD/EW factor for down-type over lepton ratios
    curve = [(e, R(e * math.log(10))) for e in range(2, 17)]
    return {"mu_equal": mu_star, "R_mZ": R(math.log(I["m_Z"])), "curve_log10mu_R": curve,
            "ms_over_mmu_at_mu": I["m_s"] / I["m_mu"] * shrink, "md_over_me_at_mu": I["m_d"] / I["m_e"] * shrink}

def seesaw_estimate(M_R=1e6, m_D=171.0):
    """Type-I seesaw with a Dirac mass equal to the top mass (Y_nu = Y_u): m_nu ~ m_D^2 / M_R in GeV."""
    return m_D ** 2 / M_R

def report():
    bt = b_tau_running()
    return {"branching": branching(), "anomalies": anomalies(),
            "wilson_lines": wilson_lines({"c": "O", "A1": "N1", "L": "N2", "B": "N2", "A2": "N3"}),
            "delta_vacua": single_block_vacua(fi_terms(delta=0.1)), "kappa_vacua": single_block_vacua(fi_terms(kappa=-0.01, gamma_charge=PS_GAMMA_FLUX)),
            "gamma": gamma_facts(), "b_tau": {k: v for k, v in bt.items() if k != "curve_log10mu_R"}, "seesaw_GeV": seesaw_estimate()}
