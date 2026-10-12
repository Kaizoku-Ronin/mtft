"""EPS-01 (2026-10-02): what sets epsilon, the family-hierarchy parameter of CK-01.

Local model at the blown-up point (infinity, O): coordinates x along X0(143), w along 143a1.  The product metric, the blow-up of the point and the exceptional
fluxes are invariant under the separate rotations of x and of w; the three families are the degree-2 jets x^2, x w, w^2 with torus charges q = 0, 1, 2.
The Yukawa entry (q, q') reads the coefficient of torus charge q + q' in the effective coupling, i.e. in the D-term profile the Higgs mixes through (PB-01).

1. Rank one (EXACT in the local model).  If the D-term is torus-constant at the point, only q = q' = 0 couples: Y = rank one, one massive family.
2. The gamma-localised condensate keeps it rank one (EXACT).  A mode built from one band of the fibre family O(c gamma)|_{x x E} = N_{c phi(x)} is a plane wave
   along the torus with x-dependent momentum; its density |Sigma|^2 has no torus dependence, so its D-term is torus-constant (1. applies).
3. Rotating the families does not help (EXACT).  Any effect that only mixes the jets x^2, x w, w^2 among themselves keeps Y an outer product (rank one).
4. What gives rank three (EXACT structure).  Y_{qq'} is proportional to the torus-charge-(q + q') component D_{q+q'} of the D-term at the point.  A condensate
   spread over jets of degree <= 2 has D-charges in [-2, 2], so the texture is
        [[1, eps, eps^2], [eps, eps^2, 0], [eps^2, 0, 0]],   eps = |D_1| / |D_0|,
   determinant -eps^6, masses ~ (1, eps^2, eps^4), mixing ~ (eps, eps, eps^2): epsilon is the condensate's torus-charge MISALIGNMENT at the blown-up point —
   a vacuum-alignment quantity, not a pure geometric constant.
5. Up versus down (EXACT structure).  A rank-one aligned condensate gives every D^a the same spatial profile, hence the same eps for up and down.  The data need
   eps_d / eps_u = 2.3 (from m_s/m_b and m_c/m_t): the condensate's SU(2)_R components must have different jet profiles, or H_u and H_d must see different
   heavy spectra.
6. Parity (EXACT).  The odd D-charges are odd under w -> -w; the eps-linear entries (V_us, V_cb) need that symmetry broken — the gamma-flux breaks it (gamma -> -gamma).
7. Size (ESTIMATE).  The gamma-flux's torus momentum across the exceptional region (area A_E/3) at the blown-up cusp bounds the misalignment it can induce.  In the
   Arakelov metric scaled to area r A_E, the cusp density from the frozen periods is rho_Ar(infinity) = K(inf, inf) / (13 * 4 pi^2) = 0.02434 (K the cusp Bergman
   kernel of PB-01), so phi' = 1/sqrt(A_f rho_Ar r) = 3.216/sqrt(r) and eps_gamma ~ (2/3 ... 2 pi/3) |c| phi' = (2.1 ... 6.7) |c| / sqrt(r).
   eps_d = 0.14 needs r ~ 230 - 2300 (c = 1): the large-r window the neutrino sector also prefers (NU-02).
CK-01's 'rank one at leading order' holds under 1-2; this module supersedes its informal wording."""
import math
import numpy as np

JETS = {"x^2": (2, 0), "x w": (1, 1), "w^2": (0, 2)}           # (x-degree, w-degree)

def allowed_entries(d_charges):
    """Yukawa entries (q, q') allowed when the D-term carries torus charges d_charges (q + q' must be one of them)."""
    return {(a, b) for a in range(3) for b in range(3) if (a + b) in set(d_charges)}

def d_charges_from_condensate(jet_degrees):
    """Torus charges present in |Sigma|^2 when the condensate spreads over jets with the given w-degrees: all differences."""
    return sorted({abs(a - b) for a in jet_degrees for b in jet_degrees})

def texture(eps, X=None):
    X = np.ones((3, 3)) if X is None else X
    return np.array([[X[a, b] * eps ** (a + b) if a + b <= 2 else 0.0 for b in range(3)] for a in range(3)], dtype=complex)

GENERIC_X = np.array([[1.0, 0.8, 1.3], [1.1, 2.0, 0.0], [0.7, 0.0, 0.0]])  # order-one, non-degenerate (all-ones is a special cancelling point)

def texture_properties(eps=0.1, X=GENERIC_X):
    Y = texture(eps, X); s = np.linalg.svd(Y, compute_uv=False)
    return {"det": complex(np.linalg.det(Y)), "singular_values": s, "ratios": (s[1] / s[0], s[2] / s[0])}

def plane_wave_density_is_constant(k=(0.7, -0.3), n=64):
    """A single fibre band is a plane wave along the torus: its density has no torus dependence (checked on a grid)."""
    x = np.linspace(0, 1, n, endpoint=False); X, Wy = np.meshgrid(x, x, indexing="ij")
    psi = np.exp(2j * math.pi * (k[0] * X + k[1] * Wy)) * 0.37
    return float(np.ptp(np.abs(psi) ** 2))

def outer_product_rank(rotation_seed=0):
    """Rotating the families by any unitary keeps a coupling to the charge-0 component an outer product (rank one)."""
    rng = np.random.default_rng(rotation_seed)
    U = np.linalg.qr(rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3)))[0]
    V = np.linalg.qr(rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3)))[0]
    e0 = np.array([1, 0, 0]); Y = np.outer(U.T @ e0, V.T @ e0)
    return int(np.linalg.matrix_rank(Y, tol=1e-12))

def monte_carlo(eps_u, eps_d, n=3000, seed=5, spread=2.0):
    """The EPS-01 texture (zeros where q + q' > 2) with independent order-one coefficients: median mass ratios and CKM magnitudes."""
    from mtft.matter.ckm_cusps import DATA, angles
    rng = np.random.default_rng(seed); out = []
    for _ in range(n):
        Xu = np.exp(rng.uniform(-math.log(spread), math.log(spread), (3, 3))) * np.exp(2j * math.pi * rng.uniform(0, 1, (3, 3)))
        Xd = np.exp(rng.uniform(-math.log(spread), math.log(spread), (3, 3))) * np.exp(2j * math.pi * rng.uniform(0, 1, (3, 3)))
        Uu, su, _ = np.linalg.svd(texture(eps_u, Xu)); Ud, sd, _ = np.linalg.svd(texture(eps_d, Xd))
        V = Uu[:, ::-1].conj().T @ Ud[:, ::-1]
        out.append([su[2] / su[0], su[1] / su[0], sd[2] / sd[0], sd[1] / sd[0], *angles(V)])
    a = np.array(out); keys = ["mu_over_mt", "mc_over_mt", "md_over_mb", "ms_over_mb", "Vus", "Vcb", "Vub"]
    return {k: {"median": float(np.median(a[:, i])), "p16": float(np.percentile(a[:, i], 16)), "p84": float(np.percentile(a[:, i], 84)), "data": DATA[k]} for i, k in enumerate(keys)}

def eps_ratio_needed():
    from mtft.matter.ckm_cusps import DATA
    return math.sqrt(DATA["ms_over_mb"]) / math.sqrt(DATA["mc_over_mt"])

def arakelov_cusp_density():
    from mtft.matter.messenger_yukawa import cusp_kernel
    K = cusp_kernel()["kernel"]; return K[0, 0] / (13 * 4 * math.pi ** 2)

A_F = 1.9699231645720703 ** 2 * 1.0232745926964613           # area of the period lattice of 143a1 (NU-02)

def eps_gamma(r, c=1, rho=None):
    rho = arakelov_cusp_density() if rho is None else rho
    phip = 1 / math.sqrt(A_F * rho * r)
    return (2 / 3 * abs(c) * phip, 2 * math.pi / 3 * abs(c) * phip)

def r_for_eps(target, c=1, rho=None):
    rho = arakelov_cusp_density() if rho is None else rho
    k = 1 / math.sqrt(A_F * rho)
    return ((2 / 3 * abs(c) * k / target) ** 2, (2 * math.pi / 3 * abs(c) * k / target) ** 2)

def report():
    rho = arakelov_cusp_density()
    return {"rank_one_entries": sorted(allowed_entries([0])), "texture": texture_properties(), "plane_wave_ptp": plane_wave_density_is_constant(),
            "rotated_rank": outer_product_rank(), "eps_ratio_needed": eps_ratio_needed(), "rho_Ar_cusp": rho, "phi_prime_coeff": 1 / math.sqrt(A_F * rho),
            "eps_gamma_r176": eps_gamma(176, rho=rho), "r_for_eps_d": r_for_eps(0.14, rho=rho), "r_for_eps_u": r_for_eps(0.06, rho=rho)}
