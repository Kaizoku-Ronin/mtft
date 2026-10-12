"""CK-01 (2026-10-02): quark mixing and the cusps of X0(143).

1. Inter-cusp distances (EXACT).  The horoball {Im tau > 1} at infinity and its Atkin–Lehner images are embedded; the image under W_Q (bottom row (143, d),
   determinant Q) has Euclidean diameter Q / 143^2, so the distance between the horoballs of infinity and W_Q(infinity) along the vertical geodesic is
   log(143^2 / Q):  d(inf, 0) = log 143 = 4.963,  d(inf, 1/11) = log(143^2/13) = 7.361,  d(inf, 1/13) = log(143^2/11) = 7.528; the other pairs follow by the
   Atkin–Lehner isometries, so any three cusps form the same triangle.  Any smooth metric that agrees with the hyperbolic one below height 1 can only lengthen
   these distances.
2. The bound (EXACT within the stated assumptions).  Families localised at different cusps mix only through the tails of their harmonic jet representatives,
   which decay at least like the massless hyperbolic Green function, e^{-d} (the twist mass only steepens it), so
       |V| <= e^{-d}  :  1/143 = 0.0070 for the closest pair, 13/143^2 = 6.4e-4 and 11/143^2 = 5.4e-4 for the others.
   |V_us| = 0.225 and |V_cb| = 0.041 both exceed 1/143: no family may sit at a cusp of its own.  Every split assignment (1+1+1 of AL-01, 2+1, 1+2) is excluded by
   quark mixing; all three families share one blown-up point — the Matter-04 background — and the family index is the jet degree 0, 1, 2 (CC-37).
3. The family triplet and the jet-charge texture (EXACT structure, STRUCTURAL magnitudes).  At one blown-up point the (4,2,1) and (4bar,1,2) families are
   the exceptional jets of a bundle with charge -3: H^0(F|_{3E}) = H^0(O_E) + H^0(O_E(1)) + H^0(O_E(2)) (1 + 2 + 3 modes) against 2 + 1 + 0 for the conjugate
   orientation, so the three net families are the degree-2 jets H^0(P^1, O(2)) = Sym^2 of the tangent plane at the point: the monomials x^2, x w, w^2 in the
   curve direction x and the torus direction w.  The torus rotation at the point acts on them with charges q = 0, 1, 2 (powers of w); it is only an approximate
   symmetry, broken by the lattice of 143a1.  Heavy bidoublet modes that are constant along the torus (the only ones the torus-constant D-term mixes with the Higgs)
   couple families (q, q') only for q + q' = 0: at leading order the Yukawa matrix is rank one — the x^2 family alone is heavy, as in local F-theory Yukawas.
   Torus-dependent components (E-momentum of the heavy modes and of the condensate) feed the lighter families with one power of a small parameter eps per unit
   charge: Y_{qq'} = X_{qq'} eps^{q+q'} with order-one X — a Froggatt–Nielsen texture whose charge is the torus degree of the jet.  Masses ~ (eps^4, eps^2, 1),
   mixing ~ (eps, eps, eps^2).  Fitting eps_u^2 = m_c/m_t and eps_d^2 = m_s/m_b (0.060, 0.139), the texture reproduces m_u/m_t and m_d/m_b within a factor of two
   to three and |V_us| ~ 0.15, but overpredicts |V_cb| (~0.15 against 0.041) and |V_ub| (~0.02 against 0.0037): the known weakness of symmetric textures, since
   the left- and right-handed jets carry the same charges.  `monte_carlo` quantifies the spread with independent order-one coefficients.
4. Consequences.  Mixing is local, so the size ratio r is no longer constrained by the CKM matrix: the large r the neutrino sector wants (NU-02) is compatible.
   The three-cusp background of AL-01 is retired for phenomenology; its index-level statements stand.
QUALIFICATION (EPS-01): 'rank one at leading order' holds when the D-term the Higgs mixes through is torus-constant at the point — true for product-type and
single-band gamma-localised condensates; epsilon is the condensate's torus-charge misalignment there.  See `mtft.matter.epsilon`.
OPEN: the Taylor coefficients x_n of the actual heavy modes (which also fix the relative sizes between up, down and lepton sectors), the phases, and the
Georgi–Jarlskog problem (Y_d = Y_e still holds entry by entry)."""
import math
import numpy as np

N = 143
DATA = {"mu_over_mt": 1.27e-3 / 171.0, "mc_over_mt": 0.619 / 171.0, "md_over_mb": 2.7e-3 / 2.86, "ms_over_mb": 0.055 / 2.86,
        "me_over_mtau": 0.487e-3 / 1.746, "mmu_over_mtau": 0.1027 / 1.746, "Vus": 0.2250, "Vcb": 0.0410, "Vub": 0.00367}


def cusp_distances():
    d = {("inf", "0"): math.log(N), ("inf", "1/11"): math.log(N * N / 13), ("inf", "1/13"): math.log(N * N / 11)}
    d[("1/11", "1/13")] = d[("inf", "0")]; d[("0", "1/13")] = d[("inf", "1/11")]; d[("0", "1/11")] = d[("inf", "1/13")]   # Atkin–Lehner isometries
    return d

def mixing_bounds():
    return {pair: math.exp(-dist) for pair, dist in cusp_distances().items()}

def split_assignments():
    """Which arrangements of three families over the cusps survive |V_us| and |V_cb|?"""
    bmax = 1 / N
    out = {}
    for name, needs in (("1+1+1 (one per cusp)", ("Vus", "Vcb", "Vub")), ("2+1 (third family alone)", ("Vcb", "Vub")), ("1+2 (first family alone)", ("Vus", "Vub")), ("3 (one point)", ())):
        violated = [v for v in needs if DATA[v] > bmax]
        out[name] = {"inter-cusp elements": needs, "exceed 1/143": violated, "allowed": not violated}
    return out


# ------------------------------------------------------------------ the jet-degree (Hankel) texture
def texture(eps, X):
    """Y_{qq'} = X_{qq'} eps^{q+q'}, q the torus charge of the jet (0, 1, 2)."""
    return np.array([[X[a, b] * eps ** (a + b) for b in range(3)] for a in range(3)], dtype=complex)

def angles(V):
    """Standard-parametrisation magnitudes from a unitary V: |V_us|, |V_cb|, |V_ub|."""
    return abs(V[0, 1]), abs(V[1, 2]), abs(V[0, 2])

def singular(Y):
    U, s, Vh = np.linalg.svd(Y)                       # Y = U diag(s) Vh, s descending
    return U, s, Vh

def monte_carlo(eps_u, eps_d, n=4000, seed=1, spread=2.0):
    """Independent order-one X entries (log-uniform magnitude in [1/spread, spread], random phase) for up and down; medians and 16-84 bands."""
    rng = np.random.default_rng(seed); out = []
    for _ in range(n):
        xu = np.exp(rng.uniform(-math.log(spread), math.log(spread), (3, 3))) * np.exp(2j * math.pi * rng.uniform(0, 1, (3, 3)))
        xd = np.exp(rng.uniform(-math.log(spread), math.log(spread), (3, 3))) * np.exp(2j * math.pi * rng.uniform(0, 1, (3, 3)))
        Uu, su, _ = singular(texture(eps_u, xu)); Ud, sd, _ = singular(texture(eps_d, xd))
        V = Uu.conj().T @ Ud                           # left-handed rotations (Y acts as Q Y u^c)
        # order singular values ascending (u, c, t) and build the mixing in that basis
        Pu = Uu[:, ::-1]; Pd = Ud[:, ::-1]; V = Pu.conj().T @ Pd
        out.append([su[2] / su[0], su[1] / su[0], sd[2] / sd[0], sd[1] / sd[0], *angles(V)])
    a = np.array(out); med = np.median(a, axis=0); lo = np.percentile(a, 16, axis=0); hi = np.percentile(a, 84, axis=0)
    keys = ["mu_over_mt", "mc_over_mt", "md_over_mb", "ms_over_mb", "Vus", "Vcb", "Vub"]
    return {k: {"median": float(m), "p16": float(l), "p84": float(h), "data": DATA[k]} for k, m, l, h in zip(keys, med, lo, hi)}

def fit_eps(ratio_second_to_third):
    """Leading order: m_2 / m_3 ~ eps^2."""
    return math.sqrt(ratio_second_to_third)

def report():
    eu, ed = fit_eps(DATA["mc_over_mt"]), fit_eps(DATA["ms_over_mb"])
    return {"distances": cusp_distances(), "bounds": mixing_bounds(), "splits": split_assignments(), "eps_u": eu, "eps_d": ed,
            "eps_e_from_mmu": fit_eps(DATA["mmu_over_mtau"]), "eps_e_from_me": DATA["me_over_mtau"] ** 0.25, "monte_carlo": monte_carlo(eu, ed)}
