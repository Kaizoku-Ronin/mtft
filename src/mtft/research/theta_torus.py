"""R2C-09 (2026-09-19): the torus factor of the surface Yukawa — theta functions on 143a1.

143a1: y^2 + y = x^3 - x^2 - x - 2, c4 = 64, Delta = -1859 = -11 * 13^2, j = -262144/1859 (exact).  tau lies on Re tau = 1/2 (negative real j):
tau = 1/2 + 1.02327459...i, found from mpmath's kleinj (normalised so that 1728 kleinj(i) = 1728).  H^0(E, N_k) has the standard basis
theta[j/k, 0](k z, k tau), orthogonal with equal norms for the metric h_k = exp(-2 pi k (Im z)^2 / Im tau) (classical; reproduced here to 1e-26).
Products theta_1 theta_2 lie in the degree-3 span (theta multiplication; residual 1e-25).  Torus factors of the mixed-origin solutions:
  S1 (Q curve (3,1), u^c torus (1,-3), Higgs (-4,2)):  B_{j alpha} = <theta3_j, theta1 theta2_alpha>, 3 x 2, singular values (1.103, 0.859), rank 2.
  S2 (Q curve (-3,1), u^c torus (1,3), Higgs (2,-4)):  B'_{alpha j} = <theta4_alpha, theta1 theta3_j>, 4 x 3, singular values (1.124, 0.979, 0.843), rank 3.
Structural consequence (EXACT): in both solutions the u^c curve factor is H^0(S0(P)) with h^0 = 2, so the up mass matrix M_{i,(a,j)} =
sum v A_{i a ...} B_{... j} has rank <= 2 at leading order: the lightest up-type quark is massless at leading order (a Riemann–Roch number).
The curve factor A of S2 is the Serre pairing of the family space with its 2-dimensional subspace of sections vanishing at two CM points;
m_c/m_t then depends only on the four torus Higgs directions — the next computation.
v0.33.0 (2026-09-23, theorem compendium VI.4–VI.5): the rank bound is qualified — rank M <= min(3, 2 * kunneth_rank(v)); unconditional
(every VEV) in S2, where the curve Higgs factor is one-dimensional; in S1 a Kuenneth-rank-two VEV reaches rank 3 (`up_mass_rank_examples`).
`torus_factor_closed_form` evaluates B through the theta multiplication formula (no quadrature) and agrees with `torus_factor` to 1e-25;
`tau_143a1_qseries` is an independent route to tau (q-series of E4^3/Delta)."""
import mpmath as mp
import numpy as np

J_143A1 = mp.mpf(64 ** 3) / mp.mpf(-1859)

def tau_143a1(dps=25):
    mp.mp.dps = dps; g = lambda y: mp.re(1728 * mp.kleinj(mp.mpc(0.5, mp.re(y))) - J_143A1); y = mp.re(mp.findroot(g, mp.mpf(1.08)))
    return mp.mpc(0.5, y)

def theta_basis(k, tau, terms=14):
    def make(jc):
        def f(z):
            s = mp.mpc(0)
            for n in range(-terms, terms + 1): m = n + mp.mpf(jc) / k; s += mp.exp(mp.pi * 1j * k * tau * m * m + 2 * mp.pi * 1j * k * m * z)
            return s
        return f
    return [make(jc) for jc in range(k)]

def inner(F, G, k, tau, N=24):
    tot = mp.mpc(0)
    for i in range(N):
        for l in range(N):
            z = (i + mp.mpf(0.5)) / N + (l + mp.mpf(0.5)) / N * tau; tot += F(z) * mp.conj(G(z)) * mp.exp(-2 * mp.pi * k * (mp.im(z)) ** 2 / mp.im(tau))
    return tot * mp.im(tau) / N ** 2

def torus_factor(k1, k2, tau=None, N=24):
    """Projection matrix of the products theta_{k1} theta_{k2} onto the degree-(k1+k2) basis, in orthonormal bases: rows = target basis,
    columns = product index (only k1 = 1 used: one theta of degree 1).  Returns matrix, singular values, closure residual, Gram off-diagonal."""
    tau = tau if tau is not None else tau_143a1(); k3 = k1 + k2; T1, T2, T3 = theta_basis(k1, tau), theta_basis(k2, tau), theta_basis(k3, tau)
    n1 = [inner(f, f, k1, tau, N) for f in T1]; n2 = [inner(f, f, k2, tau, N) for f in T2]; G3 = mp.matrix(k3, k3)
    for i in range(k3):
        for j in range(k3): G3[i, j] = inner(T3[i], T3[j], k3, tau, N)
    off = max(abs(G3[i, j]) for i in range(k3) for j in range(k3) if i != j); B = mp.matrix(k3, k1 * k2); resid = []
    for a in range(k1):
        for b in range(k2):
            prod = lambda z, a=a, b=b: T1[a](z) * T2[b](z); col = a * k2 + b; coeffs = [inner(prod, T3[j], k3, tau, N) / G3[j, j] for j in range(k3)]
            rec = lambda z, coeffs=coeffs: sum(coeffs[j] * T3[j](z) for j in range(k3)); d = lambda z, prod=prod, rec=rec: prod(z) - rec(z)
            resid.append(mp.sqrt(abs(inner(d, d, k3, tau, N) / inner(prod, prod, k3, tau, N))))
            for j in range(k3): B[j, col] = coeffs[j] * mp.sqrt(G3[j, j]) / mp.sqrt(n1[a] * n2[b])
    U, S, V = mp.svd_c(B); return {"tau": tau, "B": B, "singular_values": [float(abs(s)) for s in S], "rank": sum(1 for s in S if abs(s) > 1e-10), "closure_residual": float(max(resid)), "gram_offdiag": float(off)}

def up_mass_rank_bound():
    """Rank of the leading-order up mass matrix M_{i,(a,j)} = sum_{r,s} v_{rs} A_{iar} B_{js} (compendium VI.4, v0.33.0 qualification).
    The u^c curve factor H^0(S0(P)) is 2-dimensional in both solutions, so each Kuenneth component of the Higgs VEV v contributes a
    3 x 2 curve matrix: rank M <= min(3, 2 * kunneth_rank(v)).  S2: the curve Higgs factor H^0(O(sum P - P)) is one-dimensional, every
    VEV has Kuenneth rank one, rank M <= 2 for ALL VEVs (m_u = 0 at leading order, unconditional).  S1: the bound holds for VEVs of
    Kuenneth rank one (a single Higgs mode, or a curve mode times a torus mode); a VEV of Kuenneth rank two can give rank 3
    (`up_mass_rank_examples`), so R2C-09's "m_u = 0 in both solutions" is restricted to Kuenneth-rank-one VEVs in S1."""
    return {"u_c_curve_factor": "H^0(S0(P)), h^0 = 2", "rank_bound": 2, "rank_bound_hypothesis": "Higgs VEV of Kuenneth rank one",
            "rank_bound_general": "min(3, 2 * kunneth_rank(v))",
            "S2": "curve Higgs factor 1-dimensional: rank <= 2 for every VEV; m_u = 0 at leading order (unconditional)",
            "S1": "curve Higgs factor 16-dimensional, torus factor 2-dimensional: rank 3 is reached by a Kuenneth-rank-two VEV",
            "consequence": "m_u = 0 at leading order: unconditional in S2; in S1 for Kuenneth-rank-one VEVs"}


# ------------------------------------------------ v0.33.0 (2026-09-23): second route to the torus factor (theta multiplication formula)
def tau_143a1_qseries(dps=30, terms=40):
    """tau of 143a1 on Re tau = 1/2 from j = E4^3/Delta summed as q-series (independent of mpmath's kleinj): y = 1.02327459269646120559956631..."""
    import sympy as _sp
    mp.mp.dps = dps
    def j_of(tau):
        q = mp.exp(2j * mp.pi * tau); E4 = 1 + 240 * sum(int(_sp.divisor_sigma(n, 3)) * q ** n for n in range(1, terms))
        D = q * mp.fprod((1 - q ** n) ** 24 for n in range(1, terms)); return E4 ** 3 / D
    y = mp.findroot(lambda t: mp.re(j_of(mp.mpc(0.5, t))) - J_143A1, mp.mpf("1.02")); return mp.mpc(0.5, y)

def theta_norm_closed_form(k, tau):
    """||theta_{k,j}||_k = (Im tau / 2k)^(1/4) for every j (orthogonal basis with equal norms; compendium VI.5(c))."""
    return (mp.im(tau) / (2 * k)) ** mp.mpf(0.25)

def theta_product_coefficient(j, b, k2, tau, terms=12):
    """c_j(b) in theta_{1,0} theta_{k2,b} = sum_j c_j(b) theta_{k2+1,j}: a theta constant, sum_l exp(pi i tau k2 (k2+1) (l + (j - b - b/k2)/(k2+1))^2)
    (completing the square in the double sum; compendium VI.5(d))."""
    k3 = k2 + 1; sh = mp.mpf(j - b) - mp.mpf(b) / k2
    return sum(mp.exp(mp.pi * 1j * tau * k2 * k3 * (l + sh / k3) ** 2) for l in range(-terms, terms + 1))

def torus_factor_closed_form(k1, k2, tau=None, dps=30):
    """The torus factor B_{j,b} = <theta_1 theta_{k2,b}, theta_{k2+1,j}> / (||theta_1|| ||theta_{k2,b}|| ||theta_{k2+1,j}||) in closed form,
    B_{j,b} = c_j(b) (2 k2 / ((k2+1) Im tau))^(1/4), with no quadrature (only k1 = 1 is implemented, as in `torus_factor`).
    Agrees with `torus_factor` entrywise to ~1e-25 (compendium App. Q); singular values S1 (1.10302625437, 0.859234662985),
    S2 (1.12362454225, 0.978985971006, 0.843100702144)."""
    if k1 != 1: raise ValueError("closed form implemented for k1 = 1")
    mp.mp.dps = dps; tau = tau if tau is not None else tau_143a1(dps); k3 = k2 + 1
    pref = theta_norm_closed_form(k3, tau) / (theta_norm_closed_form(1, tau) * theta_norm_closed_form(k2, tau))
    B = mp.matrix([[theta_product_coefficient(j, b, k2, tau) * pref for b in range(k2)] for j in range(k3)])
    S = mp.svd_c(B, compute_uv=False)
    return {"tau": tau, "B": B, "singular_values": [float(abs(x)) for x in S], "rank": sum(1 for x in S if abs(x) > 1e-10), "route": "theta multiplication formula (no quadrature)"}

def up_mass_rank_examples(seed=143, n_random=50):
    """Numerical witnesses for `up_mass_rank_bound` with the closed-form torus factors: S2 is rank 2 for random VEVs; in S1 the
    Kuenneth-rank-two VEV with curve matrices A(w1) = E11, A(w2) = E22 + E31 (both in the image {A : A_12 = A_21} of w -> A(w)) gives rank 3,
    while Kuenneth-product VEVs give rank <= 2."""
    B1 = torus_factor_closed_form(1, 2)["B"]; B2 = torus_factor_closed_form(1, 3)["B"]
    nB1 = np.array([[complex(B1[i, j]) for j in range(B1.cols)] for i in range(B1.rows)]); nB2 = np.array([[complex(B2[i, j]) for j in range(B2.cols)] for i in range(B2.rows)])
    rng = np.random.default_rng(seed); A2 = rng.standard_normal((3, 2)) + 1j * rng.standard_normal((3, 2))
    s2 = {int(np.linalg.matrix_rank(np.kron(A2, (nB2.T @ (rng.standard_normal(4) + 1j * rng.standard_normal(4)))[None, :]), tol=1e-9)) for _ in range(n_random)}
    E11 = np.array([[1, 0], [0, 0], [0, 0]], complex); E22_31 = np.array([[0, 0], [0, 1], [1, 0]], complex)
    r3 = int(np.linalg.matrix_rank(np.kron(E11, nB1[:, 0][None, :]) + np.kron(E22_31, nB1[:, 1][None, :]), tol=1e-9))
    prod = set()
    for _ in range(n_random):
        Aw = rng.standard_normal((3, 2)) + 1j * rng.standard_normal((3, 2)); Aw[0, 1] = Aw[1, 0]
        prod.add(int(np.linalg.matrix_rank(np.kron(Aw, (nB1 @ rng.standard_normal(2))[None, :]), tol=1e-9)))
    return {"S2_ranks_random_VEVs": sorted(s2), "S1_rank_kunneth_rank_two_VEV": r3, "S1_ranks_product_VEVs": sorted(prod)}
