"""The mass gap, three ways: classically there is no scale; quantum-mechanically the scale is Lambda from the running coupling; on the
lattice the gap is a theorem at strong coupling.  2026-10-19.  EXACT unless marked.

Checks: (1) Derrick scaling of a static Yang-Mills energy, E[A_lambda] = lambda^{4-d} E[A] (numerically, on an abelian test field);
(2) the two-loop Lambda is renormalisation-group invariant to the order it is defined; (3) the SU(2) strong-coupling observables of the
Wilson action and of CC-69's repaired action by Weyl integration, reproducing the R&D notes section 10 (u_f shift 1.8 per cent at kappa = 1,
beta = 2.3; Delta beta = 0.051; every character coefficient non-negative); (4) the Bhanot-Creutz bulk-transition threshold against the
Monte Carlo bracket of section 12; (5) the lattice glueball numbers quoted in the lab, in units."""
import math
import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import iv

# ---------------------------------------------------------------- 1. Derrick scaling: no static lumps in d = 3
def energy(F):                                              # F = field strength samples on a grid of spacing h in d dims; E = sum F^2 h^d
    return float(np.sum(F ** 2))
for d in (1, 2, 3, 4):
    n = 64; L = 8.0; h = L / n
    ax = [np.linspace(-L / 2, L / 2, n, endpoint=False)] * d; X = np.meshgrid(*ax, indexing='ij'); r2 = sum(x ** 2 for x in X)
    A = np.exp(-r2); F = np.gradient(A, h)[0] if d > 1 else np.gradient(A, h)          # one component of dA, enough for the scaling
    E0 = energy(F) * h ** d
    lam = 2.0; Al = lam * np.exp(-lam ** 2 * r2); Fl = np.gradient(Al, h)[0] if d > 1 else np.gradient(Al, h)   # A_lambda(x) = lambda A(lambda x)
    E1 = energy(Fl) * h ** d
    assert abs(E1 / E0 - lam ** (4 - d)) < 0.05 * lam ** (4 - d), (d, E1 / E0, lam ** (4 - d))
print("Derrick: E[A_lambda] / E[A] = lambda^(4-d) for d = 1..4 (numerically, lambda = 2): in d = 3 the energy of a static lump scales as lambda, so it can be spread to zero; no classical threshold")

# ---------------------------------------------------------------- 2. dimensional transmutation: Lambda is RG-invariant
def beta_fn(N):
    b0 = 11 * N / 3; b1 = 34 * N ** 2 / 3
    return b0, b1
def Lambda_over_mu(g, N):                                   # two-loop: Lambda = mu (b0 g^2 / 16 pi^2)^{-b1/(2 b0^2)} exp(-8 pi^2 / (b0 g^2))
    b0, b1 = beta_fn(N); x = b0 * g * g / (16 * math.pi ** 2)
    return x ** (-b1 / (2 * b0 ** 2)) * math.exp(-8 * math.pi ** 2 / (b0 * g * g))
for N in (2, 3):
    b0, b1 = beta_fn(N)
    rhs = lambda t, g: [-(b0 * g[0] ** 3 / (16 * math.pi ** 2) + b1 * g[0] ** 5 / (16 * math.pi ** 2) ** 2)]   # dg/d ln mu
    mism = []
    for g0 in (1.0, 0.7, 0.5):
        sol = solve_ivp(rhs, [0, 3.0], [g0], rtol=1e-12, atol=1e-14)
        L0 = Lambda_over_mu(g0, N); L1 = Lambda_over_mu(sol.y[0, -1], N) * math.exp(3.0)   # Lambda = (Lambda/mu) mu, mu' = e^3 mu
        mism.append(abs(L1 / L0 - 1))
    assert mism[0] < 1e-2 and mism[1] < mism[0] * (0.7 / 1.0) ** 3 and mism[2] < mism[1] * (0.5 / 0.7) ** 3, (N, mism)   # the mismatch is O(g^4): it falls faster than g^3
print("two-loop Lambda is RG-invariant up to O(g^4): running g from mu to e^3 mu by the two-loop beta function changes Lambda by %.2e, %.2e, %.2e at g = 1, 0.7, 0.5 (SU(3)); b0 = 11N/3, b1 = 34N^2/3; Lambda/mu at g = 1: SU(2) %.3e, SU(3) %.3e" % (*mism, Lambda_over_mu(1, 2), Lambda_over_mu(1, 3)))

# ---------------------------------------------------------------- 3. SU(2) strong coupling: Wilson and the repaired action by Weyl integration
Y = 0.18174; NMAX = 24
W = np.zeros(NMAX + 1)
for n in range(1, NMAX + 1): W[n] = sum(math.log(d) / d for d in range(1, n + 1) if n % d == 0)
A_N = np.array([W[n] * math.exp(-2 * math.pi * Y * n) for n in range(NMAX + 1)]); A_N[1] = 0
th = (np.arange(4000) + 0.5) * math.pi / 4000; mu = (2 / math.pi) * np.sin(th) ** 2 * (math.pi / 4000)
chi = lambda j2: np.sin((j2 + 1) * th) / np.sin(th)        # spin j = j2/2: character sin((2j+1)theta)/sin theta
def action(beta, kappa):
    s = beta * (1 - np.cos(th))
    for n in range(2, NMAX + 1): s += (kappa / 2) * A_N[n] * ((n + 1) - chi(n))
    return s
def coefficients(beta, kappa, jmax2=16):                    # c_j = (1/d_j) int f chi_j dmu, u_j = c_j / c_0
    f = np.exp(-action(beta, kappa)); c0 = np.sum(f * mu)
    return [np.sum(f * chi(j2) * mu) / (j2 + 1) / c0 for j2 in range(0, jmax2 + 1)]
# the Wilson action's u_f is I_2(beta)/I_1(beta)
for b in (0.5, 1.0, 2.3):
    assert abs(coefficients(b, 0)[1] - iv(2, b) / iv(1, b)) < 1e-6
# every coefficient non-negative on a grid (positive-type class function): the Osterwalder-Seiler hypothesis
mn = min(min(coefficients(b, k)) for b in (0, 1, 2.5, 6) for k in (0, 1, 10, 100, 1000))
assert mn > -1e-12, mn
uW = coefficients(2.3, 0)[1]; uR = coefficients(2.3, 1)[1]
shift = uR / uW - 1
from scipy.optimize import brentq
beta_eff = brentq(lambda b: coefficients(b, 0)[1] - uR, 2.0, 3.0); dbeta = beta_eff - 2.3
assert abs(shift - 0.018) < 0.004 and abs(dbeta - 0.051) < 0.006, (shift, dbeta)
act = lambda beta, kappa: sum((j2 + 1) ** 2 * u for j2, u in enumerate(coefficients(beta, kappa)) if j2 >= 1)
a1 = act(1.0, 1); a23 = act(2.3, 1)                       # the repaired action at kappa = 1, as in section 10
assert abs(a1 - 1.7) < 0.1, a1
print("SU(2) Weyl integration: Wilson u_f = I2/I1 (1e-6); all character coefficients >= 0 up to kappa = 1000 (min %.1e); at beta = 2.3 the tower (kappa = 1) shifts u_f by %.1f%%, worth Delta beta = %.3f -- section 10 reproduced" % (mn, 100 * shift, dbeta))
print("   strong-coupling leading order: sigma a^2 = -ln u_f, m(0++) a = -4 ln u_f; character activity sum d^2 u = %.2f at beta = 1 (the theorem's regime), %.1f at beta = 2.3 (repaired action, kappa = 1)" % (a1, a23))

# ---------------------------------------------------------------- 4. the bulk transition: Bhanot-Creutz threshold against the Monte Carlo bracket
beta_A = lambda kappa: 0.053 * kappa                        # the Sym^2 (adjoint) coupling of the repaired SU(2) action in Bhanot-Creutz's convention (section 12)
kappa_BC = 1.25 / 0.053
assert 10 < kappa_BC < 30
print("bulk transition: beta_A = 0.053 kappa reaches Bhanot-Creutz's first-order region beta_A >= 1.25 at kappa = %.1f, inside the Monte Carlo bracket (10, 30) of section 12; the package default kappa = 1 is far below it" % kappa_BC)

# ---------------------------------------------------------------- 5. the numbers quoted in the lab
MP = {'0++': (1730, 50, 80), '2++': (2400, 25, 120), '0-+': (2590, 40, 130)}      # Morningstar-Peardon 1999, SU(3), MeV (stat)(sys)
sqrt_sigma = 440.0
for k, (m, s1, s2) in MP.items(): print("   SU(3) glueball %s: %d(%d)(%d) MeV = %.2f sqrt(sigma) at sqrt(sigma) = 440 MeV" % (k, m, s1, s2, m / sqrt_sigma))
print("ALL CHECKS PASSED")
