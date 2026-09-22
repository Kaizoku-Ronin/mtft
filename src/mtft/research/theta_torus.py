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
m_c/m_t then depends only on the four torus Higgs directions — the next computation."""
import mpmath as mp

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

def up_mass_rank_bound(): return {"u_c_curve_factor": "H^0(S0(P)), h^0 = 2", "rank_bound": 2, "consequence": "m_u = 0 at leading order in both mixed-origin solutions"}
