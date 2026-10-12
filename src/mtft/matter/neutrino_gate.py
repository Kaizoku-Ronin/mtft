"""NU-01 (2026-10-02): the neutrino gate of the Pati–Salam branch — what neutral fermions the construction provides, where they live, which masses
the gauge symmetry allows, and what an inverse seesaw would need.

Background: PS-01 (minimal Scherk–Schwarz flux, Pati–Salam groups 4 = c + A1, 2_L = L, 2_R = B + A2).  Pati–Salam breaking uses gamma-flux; the
cheapest Pati–Salam-symmetric choice that keeps the bidoublet block free of gamma-flux is gamma-charges (1, 1, -1, -1, -1) on (c, A1, L, B, A2) with
delta = -8 kappa (FI only on the 4 and R groups: the (4,1,2) block alone is D-flat).  Its block gamma-differences are (4,L) 2, (4,R) 2, (L,R) 0
(localised extras 256 per orientation, against 384 for PS-01's (1, 1, 0, -2, -2)).

1. Inventory (EXACT, index level; fermion zero modes under the twist sit in one cohomological degree, so each orientation carries -chi).
   Entries Hom(j -> i) per rank-one component, classes from `mtft.surface.blowup` (a, b, c, e); "3 + 16" = 3 exceptional jets + 16 gamma-localised:
     lepton doublet entry (A1 <- L):  3 + 16 modes against 6 + 16 in the conjugate orientation — net 3, plus 3 + 16 = 19 vector-like pairs;
     nu^c entry (A2 <- A1):           3 + 16 against 6 + 16 — net 3 (B-L = +1), plus 19 vector-like (nu^c, N) pairs;
     bidoublet entries (L <- B), (L <- A2): no gamma-flux, Delta e = 6: 21 against 15 per entry, B-L = 0, Y = -+1/2.
   (Which orientation is left-handed is a convention; the counts are not.)
   There are no B-L-neutral SINGLET fermion zero modes: adjoint fermions are removed by the twist (Matter-04).  The only B-L-neutral neutral fermions
   are the bidoublet ("higgsino-like") components.
2. Localisation (EXACT on the product).  A block with gamma-charge c restricted to a fibre {x} x 143a1 is the flat bundle of c phi(x); bosons see it
   untwisted, fermions twisted by the spin structure eta (2-torsion point P_eta).  So the 4c^2 extra bosonic modes sit on phi^-1(E[c]) and the 4c^2
   fermionic ones on phi^-1({Q : cQ = P_eta}); the two sets are DISJOINT (c phi(x) cannot be both O and P_eta).  For c = 1 the bosonic set is exactly
   the four cusps.  Consequence: a holomorphic coupling whose fermion classes are supported on the fermionic set and whose bosonic vev is supported on the
   bosonic set vanishes identically (classes with disjoint supports have zero cup product).  This is a selection rule, not yet a mass: it removes the
   leading holomorphic Majorana coupling of the gamma-localised N through the gamma-localised Pati–Salam condensate; jet-sector couplings at the blown-up
   point are not covered.
3. Atkin–Lehner at the cusps (EXACT).  143a is the (+,+) newform (rank 1, root number -1 = -epsilon_143), so phi o W_Q = phi for Q in {11, 13, 143}
   (all cusps map to O), and W_Q x id_E preserves the graph, gamma, K and the spin structure.  The Atkin–Lehner group (Z/2)^2 permutes the four cusps
   simply transitively: the four cusp-localised modes of a c = +-1 block carry its regular representation — one mode in each sector (+,+), (+,-), (-,+),
   (-,-), including the (-,-) sector that contains no weight-2 newform (Paper 37's residue sector).  The single-point blow-up and the cusp line bundle
   O(infinity) break this symmetry; an Atkin–Lehner-symmetric background would blow up the whole orbit of (infinity, O).
4. Selection rules (EXACT).  B-L is gauged inside SU(4); the condensate carries |Delta(B-L)| = 1, so Majorana masses (Delta(B-L) = 2) for nu^c or N need
   two insertions: dimension five, M ~ v_R^2 / M_*.  Dirac pairings of the vector-like (nu^c, N) conserve B-L.
5. Requirements (NUMERICAL).  Y_nu = Y_u at M_PS = 1.1e6 GeV (PS-01).  With one-loop SM running, m_D(3) = y_t(M_PS) v / sqrt 2.  Normal ordering with
   m_3 = 50 meV (the scale Paper 37 targets), m_2 = 8.7 meV:
     type I needs M_R(3) ~ m_D^2 / m_3 ~ 3e14 GeV (m_D(3) = 121 GeV) — eight orders above M_PS: excluded in this branch;
     an inverse seesaw with M = M_PS needs mu(3) = m_3 M^2 / m_D(3)^2 = 4.3 MeV and mu(2) = 65 GeV, i.e. lepton-number suppressions mu / M of 3.8e-9 and 5.7e-5
     (third family) — far below the dimension-five estimate mu ~ M.  Item 2 is the first structural source of such a suppression found in the construction.
OPEN: the B-L-conserving Dirac pairing M (no singlet scalar vev at tree level), jet-sector Majorana couplings, the PMNS structure (Y_nu = Y_u has small
mixing), Majorana versus Dirac, and the Atkin–Lehner-symmetric background."""
from fractions import Fraction
import math

from mtft.surface.blowup import chi

STACK_CLASS = {"c": (0, 0, 0, 0), "A1": (0, 0, 0, 0), "L": (1, 0, 0, -3), "B": (-1, 0, 0, 3), "A2": (-1, 0, 0, 3)}
GAMMA_CHEAP = {"c": 1, "A1": 1, "L": -1, "B": -1, "A2": -1}
GAMMA_PS01 = {"c": 1, "A1": 1, "L": 0, "B": -2, "A2": -2}

def stack_class(s, gamma):
    a, b, c, e = STACK_CLASS[s]; return (a, b, c + gamma.get(s, 0), e)

def entry(i, j, gamma=None):
    """Fermion zero modes of Hom(j -> i) and of its conjugate Hom(i -> j) (one degree each under the twist): (-chi(L), -chi(L^-1)), net, pairs."""
    gamma = gamma if gamma is not None else GAMMA_CHEAP
    ci, cj = stack_class(i, gamma), stack_class(j, gamma); L = tuple(x - y for x, y in zip(ci, cj)); Linv = tuple(-x for x in L)
    f, fbar = -chi(L), -chi(Linv)
    return {"class": L, "modes": f, "conjugate_modes": fbar, "net": f - fbar, "vector_like_pairs": min(f, fbar), "gamma_charge": L[2], "localised_per_orientation": 4 * L[2] ** 2}

def fi_check(gamma=None, kappa=-0.01):
    """D-flatness of the Pati–Salam condensate for a gamma assignment with delta chosen to zero the L-group FI."""
    from mtft.matter.ps_breaking import fi_terms, single_block_vacua
    gamma = gamma if gamma is not None else GAMMA_CHEAP
    delta = 8 * kappa * gamma["L"]
    return single_block_vacua(fi_terms(delta=delta, kappa=kappa, gamma_charge=gamma))

def localised_extras(gamma):
    groups = {"(4,L)": ("c", "L"), "(4,R)": ("c", "B"), "(L,R)": ("L", "B")}
    dims = {"(4,L)": 4 * 2, "(4,R)": 4 * 2, "(L,R)": 2 * 2}
    return {g: 4 * (gamma[i] - gamma[j]) ** 2 * dims[g] for g, (i, j) in groups.items()}

def inventory(gamma=None):
    gamma = gamma if gamma is not None else GAMMA_CHEAP
    from mtft.matter.ps_breaking import block_hypercharge, BML
    rows = {"nu_L (A1 <- L)": entry("A1", "L", gamma), "nu^c (A2 <- A1)": entry("A2", "A1", gamma),
            "bidoublet (L <- B)": entry("L", "B", gamma), "bidoublet (L <- A2)": entry("L", "A2", gamma)}
    charges = {"nu_L (A1 <- L)": ("A1", "L"), "nu^c (A2 <- A1)": ("A2", "A1"), "bidoublet (L <- B)": ("L", "B"), "bidoublet (L <- A2)": ("L", "A2")}
    for k, (i, j) in charges.items():
        rows[k]["B_minus_L"] = BML[i] - BML[j]; rows[k]["Y"] = block_hypercharge(i, j)
    return {"rows": rows, "BL_neutral_singlets": 0, "localised_extras_total": sum(localised_extras(gamma).values())}


# ------------------------------------------------------------------ localisation and Atkin–Lehner
def localisation_sets(c):
    """Counts (with multiplicity) of the bosonic and fermionic localisation sets of a block with gamma-charge c on the product, and disjointness."""
    if c == 0: return {"bosons": 0, "fermions": 0, "disjoint": True}
    return {"bosons": 4 * c * c, "fermions": 4 * c * c, "disjoint": True, "boson_condition": f"{c} phi(x) = O", "fermion_condition": f"{c} phi(x) = P_eta",
            "bosons_are_the_cusps": abs(c) == 1}

def atkin_lehner_on_cusps():
    """Cusps of X0(143) <-> divisors d of 143 (cusp 1/d); W_Q sends d to d Q / gcd(d, Q)^2.  Returns the permutation table and the multiplicity of each
    character of (Z/2)^2 in the permutation representation (regular: each once)."""
    N = 143; divs = [1, 11, 13, 143]; Qs = [11, 13, 143]
    act = {Q: {d: d * Q // math.gcd(d, Q) ** 2 for d in divs} for Q in Qs}
    elems = [(0, 0), (1, 0), (0, 1), (1, 1)]                                  # exponents of (W_11, W_13)
    def perm(e):
        def f(d):
            if e[0]: d = act[11][d]
            if e[1]: d = act[13][d]
            return d
        return f
    fixed = {e: sum(1 for d in divs if perm(e)(d) == d) for e in elems}
    chars = {(s11, s13): sum((s11 ** e[0]) * (s13 ** e[1]) * fixed[e] for e in elems) // 4 for s11 in (1, -1) for s13 in (1, -1)}
    transitive = len({perm(e)(143) for e in elems}) == 4
    return {"action": act, "simply_transitive": transitive and all(fixed[e] == 0 for e in elems[1:]), "character_multiplicities": chars,
            "phi_commutes": "143a is (+,+), root number -1 = -epsilon_143: phi o W_Q = phi + (phi(W_Q infinity) - O) = phi"}


# ------------------------------------------------------------------ requirements
def up_yukawas_at(mu_target=1.1236e6):
    """One-loop SM running of y_t, y_b, y_tau (as in ps_breaking) with y_c, y_u passive (own squares dropped)."""
    from scipy.integrate import solve_ivp
    from mtft.matter.ps_breaking import INPUTS_MZ, _rge
    I = INPUTS_MZ
    def rhs(t, y):
        base = _rge(t, y[:6]); g1, g2, g3, yt, yb, ytau = y[:6]; k = 1 / (16 * math.pi ** 2); Y2 = 3 * yt ** 2 + 3 * yb ** 2 + ytau ** 2
        light = k * (Y2 - 17 / 12 * g1 ** 2 - 9 / 4 * g2 ** 2 - 8 * g3 ** 2 - 1.5 * yb ** 2)
        return base + [y[6] * light, y[7] * light]
    g1 = math.sqrt(4 * math.pi / (I["alpha_inv"] * (1 - I["sin2w"]))); g2 = math.sqrt(4 * math.pi / (I["alpha_inv"] * I["sin2w"])); g3 = math.sqrt(4 * math.pi * I["alpha_s"])
    r2 = math.sqrt(2) / I["v"]
    y0 = [g1, g2, g3, r2 * I["m_t"], r2 * I["m_b"], r2 * I["m_tau"], r2 * 0.619, r2 * 0.00127]       # m_c(m_Z), m_u(m_Z) MS-bar (approximate)
    sol = solve_ivp(rhs, [math.log(I["m_Z"]), math.log(mu_target)], y0, rtol=1e-10, atol=1e-14)
    yt, yc, yu = sol.y[3][-1], sol.y[6][-1], sol.y[7][-1]
    return {"y_t": yt, "y_c": yc, "y_u": yu, "m_D_GeV": {"3": yt * I["v"] / math.sqrt(2), "2": yc * I["v"] / math.sqrt(2), "1": yu * I["v"] / math.sqrt(2)}}

NU_MASSES_EV = {"3": 0.050, "2": 0.0087, "1": 0.0}

def seesaw_requirements(M=1.1236e6):
    mD = up_yukawas_at(M)["m_D_GeV"]; out = {}
    for f in ("3", "2"):
        m = NU_MASSES_EV[f] * 1e-9                                               # GeV
        out[f] = {"m_D_GeV": mD[f], "typeI_M_R_GeV": mD[f] ** 2 / m, "inverse_mu_GeV": m * M ** 2 / mD[f] ** 2, "mu_over_M": m * M / mD[f] ** 2}
    return out

def report():
    return {"inventory": inventory(), "fi_cheap": fi_check(), "localised_cheap": localised_extras(GAMMA_CHEAP), "localised_ps01": localised_extras(GAMMA_PS01),
            "loc_c1": localisation_sets(1), "loc_c2": localisation_sets(2), "atkin_lehner": atkin_lehner_on_cusps(), "requirements": seesaw_requirements()}
