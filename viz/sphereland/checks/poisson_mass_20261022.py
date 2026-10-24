#!/usr/bin/env python3
"""Poisson's ratio as a mass mechanism — checks for Sphereland Instrument 20.

Detour of 22 Oct 2026 (Roger: "exploring Poisson ratio as a mass mechanism with the
insights we have from modular time and knots"; attached: Wikipedia's Poisson's ratio,
Kronecker delta, Degree of a continuous mapping, Planck units).

Thesis checked here: Poisson's ratio is dimensionless, so it never makes a mass scale; it
sets ratios, signs and which modes exist.  Every isotropic quadratic energy splits exactly as
    W = (M/2)(tr e)^2 - 2 mu sigma_2(e),      M = lambda + 2 mu (P-wave modulus),
with sigma_2 the generalized Kronecker delta contraction (sum of 2x2 principal minors), and
nu fixes the share 2 mu / M = (1 - 2 nu)/(1 - nu).  The sigma_2 part is a null Lagrangian on
gradients, Gauss's curvature on surfaces (topological by Gauss-Bonnet), the ADM kinetic term of
gravity and the Fierz-Pauli mass; the labile aether of Cauchy (1839) and Kelvin (1888), M = 0,
is pure sigma_2.

Sections:  1 Kronecker deltas and the ratio   2 Lamb waves: masses from thickness
           3 the elastic torus: modular time is a mass ratio   4 the labile aether: graviton
           mass, DeWitt supermetric, Horava   5 rods and knots: twist/bend = 1/(1+nu), Michell,
           Lk = Tw + Wr as a degree   6 topology's rest energy: shells, Willmore, dislocations
           7 the programme's bridge under stress (model)

Every numbered block prints what it checks; ends with ALL CHECKS PASSED.  Runtime ~1-2 min.
Needs numpy, mpmath, sympy.  Layers: [structure] theorem/exact identity; [model] assumption;
[audit] a number the programme has to answer for.
"""
from __future__ import annotations

import math
import sys
import time
from fractions import Fraction

import mpmath as mp
import numpy as np
import sympy as sp

T0 = time.time()
FAILS: list[str] = []


def check(name: str, ok: bool, detail: str = "") -> None:
    print(f"  [{'ok  ' if ok else 'FAIL'}] {name}" + (f" — {detail}" if detail else ""))
    if not ok:
        FAILS.append(name)


rng = np.random.default_rng(20261022)

# =====================================================================================
print("§1  Two Kronecker deltas, one ratio  [structure]")
# =====================================================================================


def rand_orth(d, reflect=False):
    Q, R = np.linalg.qr(rng.standard_normal((d, d)))
    Q = Q @ np.diag(np.sign(np.diag(R)))
    if np.linalg.det(Q) < 0:
        Q[:, 0] *= -1
    if reflect:
        Q = Q @ np.diag([-1.0] + [1.0] * (d - 1))
    return Q


def perm_matrix(d, perm):
    """Matrix of T_{i0 i1 i2 i3} -> T_{i_perm} on row-major vec(T)."""
    n = d ** 4
    P = np.zeros((n, n))
    for idx in range(n):
        ijkl = np.unravel_index(idx, (d,) * 4)
        src = tuple(ijkl[perm[a]] for a in range(4))
        P[idx, np.ravel_multi_index(src, (d,) * 4)] = 1.0
    return P


def invariant_dim(d, group="O", symmetric=False):
    rows = []
    mats = [rand_orth(d) for _ in range(3)]
    if group == "O":
        mats.append(rand_orth(d, reflect=True))
    for R in mats:
        R4 = np.kron(np.kron(R, R), np.kron(R, R))
        rows.append(R4 - np.eye(d ** 4))
    if symmetric:
        for perm in [(1, 0, 2, 3), (0, 1, 3, 2), (2, 3, 0, 1)]:
            rows.append(perm_matrix(d, perm) - np.eye(d ** 4))
    s = np.linalg.svd(np.vstack(rows), compute_uv=False)
    return int(np.sum(s < 1e-9))


dims = {d: invariant_dim(d, "O") for d in (2, 3, 4)}
check("O(d)-invariant 4-tensors are spanned by the three products δ_ijδ_kl, δ_ikδ_jl, δ_ilδ_jk (d = 2, 3, 4)", all(v == 3 for v in dims.values()), str(dims))
check("SO(2) alone allows 6 (chiral ones), but none survive the symmetries of an elastic energy", invariant_dim(2, "SO") == 6 and invariant_dim(2, "SO", True) == 2)
symd = {d: invariant_dim(d, "O", True) for d in (2, 3, 4)}
check("with minor and major symmetry two survive: C_ijkl = λ δ_ijδ_kl + μ(δ_ikδ_jl + δ_ilδ_jk)", all(v == 2 for v in symd.values()), str(symd))


def stiffness(d, lam, mu):
    I = np.eye(d)
    return lam * np.einsum("ij,kl->ijkl", I, I) + mu * (np.einsum("ik,jl->ijkl", I, I) + np.einsum("il,jk->ijkl", I, I))


def sym_basis(d):
    B = []
    for i in range(d):
        for j in range(i, d):
            E = np.zeros((d, d))
            if i == j:
                E[i, i] = 1.0
            else:
                E[i, j] = E[j, i] = 1 / math.sqrt(2)
            B.append(E)
    return B


ok = True
for d in (2, 3, 4):
    lam, mu = rng.uniform(-0.5, 2), rng.uniform(0.3, 2)
    C = stiffness(d, lam, mu)
    B = sym_basis(d)
    Mat = np.array([[np.einsum("ij,ijkl,kl->", Ba, C, Bb) for Bb in B] for Ba in B])
    ev = np.sort(np.linalg.eigvalsh(Mat))
    expect = np.sort([2 * mu] * (len(B) - 1) + [d * lam + 2 * mu])
    ok &= np.allclose(ev, expect, atol=1e-12)
check("on symmetric strains C has eigenvalues 2μ (shear, d(d+1)/2 − 1 times) and dλ + 2μ = d·K (bulk, once)", ok)

lam_s, mu_s, d_s, nu_s = sp.symbols("lambda mu d nu", real=True)
nu_d = lam_s / ((d_s - 1) * lam_s + 2 * mu_s)
# invert: lambda/mu as a function of nu
x_of_nu = sp.solve(sp.Eq(nu_s, nu_d.subs(lam_s, sp.Symbol("x") * mu_s)), sp.Symbol("x"))[0]
cL2_cT2 = sp.simplify((x_of_nu + 2))
check("ν_d = λ/((d−1)λ + 2μ); c_L²/c_T² = (λ + 2μ)/μ = 2(1 − (d−2)ν)/(1 − (d−1)ν)", sp.simplify(cL2_cT2 - 2 * (1 - (d_s - 2) * nu_s) / (1 - (d_s - 1) * nu_s)) == 0, f"3D: {sp.simplify(cL2_cT2.subs(d_s, 3))}")
E_young = 2 * mu_s * (d_s * lam_s + 2 * mu_s) / ((d_s - 1) * lam_s + 2 * mu_s)
check("E = 2μ(1 + ν) and E = d·K_d(1 − (d−1)ν) in every dimension", sp.simplify(E_young - 2 * mu_s * (1 + nu_d)) == 0 and sp.simplify(E_young - (d_s * lam_s + 2 * mu_s) * (1 - (d_s - 1) * nu_d)) == 0)

ok = True
for d in (2, 3, 4):
    for nu in np.linspace(-1.4, 1.4, 281):
        if abs(1 - (d - 1) * nu) < 1e-9 or abs(nu + 1) < 1e-9:
            continue
        x = 2 * nu / (1 - (d - 1) * nu)             # lambda/mu at mu = 1
        posdef = (1 > 0) and (x + 2 / d > 0)          # mu > 0 and K_d = lambda + 2mu/d > 0
        ok &= posdef == (-1 < nu < 1 / (d - 1))
check("positive-definite (free body: μ > 0, K > 0) ⇔ −1 < ν < 1/(d−1)", ok, "2D: (−1, 1); 3D: (−1, ½); 4D: (−1, ⅓)")
ok = True
for d in (2, 3, 4):
    n = rng.standard_normal(d)
    n /= np.linalg.norm(n)
    lam, mu = rng.uniform(-3, 2), rng.uniform(0.3, 2)
    Q = np.einsum("ijkl,j,l->ik", stiffness(d, lam, mu), n, n)
    ev = np.sort(np.linalg.eigvalsh(Q))
    ok &= np.allclose(ev, np.sort([mu] * (d - 1) + [lam + 2 * mu]))
check("acoustic tensor C_ijkl n_j n_l has eigenvalues μ (×(d−1)) and M = λ + 2μ: strong ellipticity ⇔ μ > 0, M > 0", ok)
check("M = 0 (zero longitudinal speed) ⇔ ν_d = 1/(d−2): ν = 1 in 3D, ½ in 4D", all(sp.simplify(nu_d.subs({lam_s: -2 * mu_s, d_s: d}) - sp.Rational(1, d - 2)) == 0 for d in (3, 4, 5)))

# the split W = (M/2)(tr e)^2 - 2 mu sigma_2(e)
e = sp.Matrix(3, 3, lambda i, j: sp.Symbol(f"e{min(i, j)}{max(i, j)}"))
tr = e.trace()
sigma2 = sum(e[i, i] * e[j, j] - e[i, j] ** 2 for i in range(3) for j in range(i + 1, 3))
W = mu_s * sum(e[i, j] ** 2 for i in range(3) for j in range(3)) + lam_s / 2 * tr ** 2
check("every isotropic energy is W = (M/2)(tr ε)² − 2μ σ₂(ε), σ₂ = sum of 2×2 principal minors", sp.expand(W - ((lam_s + 2 * mu_s) / 2 * tr ** 2 - 2 * mu_s * sigma2)) == 0)
gkd = sum((1 if (i == k and j == l) else 0) * 0 for i in range(1) for j in range(1) for k in range(1) for l in range(1))  # placeholder for clarity
delta2 = lambda i, j, k, l: (1 if (i == k and j == l) else 0) - (1 if (i == l and j == k) else 0)
gk_contraction = sum(delta2(i, j, k, l) * e[k, i] * e[l, j] for i in range(3) for j in range(3) for k in range(3) for l in range(3))
check("σ₂(ε) is the generalized Kronecker delta: δ^{ij}_{kl} ε^k_i ε^l_j = (tr ε)² − tr ε² = 2σ₂(ε)", sp.expand(gk_contraction - 2 * sigma2) == 0)
share = sp.simplify((2 * mu_s / (lam_s + 2 * mu_s)).subs(lam_s, x_of_nu.subs(d_s, 3) * mu_s))
check("Poisson's ratio is the Kronecker share: 2μ/M = (1 − 2ν)/(1 − ν) in 3D", sp.simplify(share - (1 - 2 * nu_s) / (1 - nu_s)) == 0, "ν = ½: 0 (all uniaxial); ν = 0: 1; ν = −1: 3/2; ν = 1: ∞ (all Kronecker)")
check("σ₂ vanishes on uniaxial strain (rank one) — so at M = 0 a uniaxial strain costs nothing", sp.simplify(sigma2.subs({sp.Symbol('e11'): 0, sp.Symbol('e22'): 0, sp.Symbol('e01'): 0, sp.Symbol('e02'): 0, sp.Symbol('e12'): 0})) == 0)

# null Lagrangian and the MacCullagh identity
X = sp.symbols("x0:3")
u = [sp.Function(f"u{i}")(*X) for i in range(3)]
grad = sp.Matrix(3, 3, lambda i, j: sp.diff(u[j], X[i]))       # grad[i,j] = d_i u_j
div = sum(grad[i, i] for i in range(3))
A_ = sum(grad[i, j] * grad[j, i] for i in range(3) for j in range(3))
two_sigma2 = div ** 2 - A_
flux = [u[i] * div - sum(u[j] * grad[j, i] for j in range(3)) for i in range(3)]
check("null Lagrangian: 2σ₂(∇u) = ∂_i(u_i div u − u_j ∂_j u_i) for any field u", sp.simplify(sp.expand(two_sigma2 - sum(sp.diff(flux[i], X[i]) for i in range(3)))) == 0)
curl = [sp.diff(u[2], X[1]) - sp.diff(u[1], X[2]), sp.diff(u[0], X[2]) - sp.diff(u[2], X[0]), sp.diff(u[1], X[0]) - sp.diff(u[0], X[1])]
eps = sp.Matrix(3, 3, lambda i, j: (grad[i, j] + grad[j, i]) / 2)
labile = mu_s * (sum(eps[i, j] ** 2 for i in range(3) for j in range(3)) - div ** 2)
check("the labile aether (λ = −2μ) is MacCullagh's rotational aether plus a divergence: W = (μ/2)|curl u|² − 2μσ₂(∇u)", sp.expand(labile - (mu_s / 2 * sum(c ** 2 for c in curl) - mu_s * two_sigma2)) == 0, "so with a fixed boundary its equations are Maxwell's in temporal gauge (u ↔ A)")
print("    Roger's page: (1 + ΔL/L)^(−ν) = 1 + ΔL'/L is Hencky's law λ_t = λ_a^(−ν): ν^Hencky = −ln λ_t / ln λ_a, used exactly in §3.")

# =====================================================================================
print("§2  Lamb waves: thickness makes masses, ν sets them  [structure]")
# =====================================================================================
mp.mp.dps = 30


def kappa_of(nu):
    return mp.sqrt(2 * (1 - nu) / (1 - 2 * nu))


def _cs(x2):
    return mp.cos(mp.sqrt(x2)) if x2 >= 0 else mp.cosh(mp.sqrt(-x2))


def _sx(x2):
    if x2 == 0:
        return mp.mpf(1)
    return mp.sin(mp.sqrt(x2)) / mp.sqrt(x2) if x2 > 0 else mp.sinh(mp.sqrt(-x2)) / mp.sqrt(-x2)


def _xs(x2):
    return mp.sqrt(x2) * mp.sin(mp.sqrt(x2)) if x2 >= 0 else -mp.sqrt(-x2) * mp.sinh(mp.sqrt(-x2))


def F_lamb(Om, K, kap, sym):
    """Rayleigh–Lamb, entire form, half-thickness h = 1, c_T = 1: Om = ωh/c_T, K = kh."""
    q2, p2 = Om ** 2 - K ** 2, Om ** 2 / kap ** 2 - K ** 2
    if sym:
        return (q2 - K ** 2) ** 2 * _cs(p2) * _sx(q2) + 4 * K ** 2 * _xs(p2) * _cs(q2)
    return 4 * K ** 2 * _xs(q2) * _cs(p2) + (q2 - K ** 2) ** 2 * _sx(p2) * _cs(q2)


def root_near(Om0, K, kap, sym, span=mp.mpf("0.02")):
    f = lambda O: F_lamb(O, K, kap, sym)
    a, b = Om0 - span, Om0 + span
    xs = [a + (b - a) * i / 400 for i in range(401)]
    vals = [f(x) for x in xs]
    best = None
    for i in range(400):
        if vals[i] == 0 or vals[i] * vals[i + 1] < 0:
            r = mp.findroot(f, (xs[i], xs[i + 1]), solver="anderson")
            if best is None or abs(r - Om0) < abs(best - Om0):
                best = r
    return best


nu = mp.mpf("0.25")
kap = kappa_of(nu)
cut_s = sorted([mp.pi * n for n in (1, 2)] + [(2 * n + 1) * mp.pi * kap / 2 for n in (0, 1)])
cut_a = sorted([(2 * n + 1) * mp.pi / 2 for n in (0, 1, 2)] + [n * mp.pi * kap for n in (1,)])
ok = True
for c in cut_s[:3]:
    r = root_near(c, mp.mpf("1e-4"), kap, True, mp.mpf("0.01"))
    ok &= r is not None and abs(r - c) < 1e-6
for c in cut_a[:3]:
    r = root_near(c, mp.mpf("1e-4"), kap, False, mp.mpf("0.01"))
    ok &= r is not None and abs(r - c) < 1e-6
check("cutoffs (masses) at k → 0: symmetric nπc_T/h and (2n+1)πc_L/2h, antisymmetric (2n+1)πc_T/2h and nπc_L/h", ok, f"ν = 0.25, c_L/c_T = {mp.nstr(kap, 8)}")
check("the two lowest symmetric cutoffs coincide at ν = 1/3 (c_L = 2c_T)", abs(kappa_of(mp.mpf(1) / 3) - 2) < 1e-25)


def s1_cutoff(nu):
    k_ = kappa_of(nu)
    return min(mp.pi, mp.pi * k_ / 2)


def s1_curvature(nu, K=mp.mpf("0.01")):
    k_ = kappa_of(nu)
    Oc = s1_cutoff(nu)
    vals = []
    for KK in (K, 2 * K):
        r = root_near(Oc, KK, k_, True, mp.mpf("0.004"))
        vals.append(2 * (r - Oc) / KK ** 2)
    curv = (4 * vals[0] - vals[1]) / 3
    return curv * Oc        # = 1 for a relativistic particle with speed c_T


curvs = {float(v): s1_curvature(mp.mpf(v)) for v in ("-0.5", "0", "0.25", "0.29", "0.4", "0.47")}
print("    S1 at its cutoff: ω''·ω_c/c_T² (1 = relativistic with speed c_T; negative = negative effective mass):")
print("      " + ", ".join(f"ν={k}: {mp.nstr(v, 5)}" for k, v in curvs.items()))
nu_hi = mp.findroot(lambda v: s1_curvature(v), (mp.mpf("0.44"), mp.mpf("0.46")), solver="anderson")
nu_lo = mp.findroot(lambda v: s1_curvature(v), (mp.mpf("-0.56"), mp.mpf("-0.52")), solver="anderson")
check("S1's effective mass ħ/ω'' is negative (a backward wave, with a zero-group-velocity point) exactly for ν between two roots", curvs[0.25] < 0 and curvs[0.47] > 0 and curvs[-0.5] < 0, f"roots ν = {mp.nstr(nu_lo, 6)} and {mp.nstr(nu_hi, 6)}")
u_hi = mp.findroot(lambda u: u * mp.tan(u) - mp.pi ** 2 / 8, mp.mpf("0.93"))           # u = π c_T/c_L, shear-cutoff branch
v_lo = mp.findroot(lambda v: v ** 3 * mp.tan(v) + 2 * mp.pi ** 2, mp.mpf("1.91"))        # v = π c_L/2c_T, stretch-cutoff branch
nu_of_kappa = lambda k_: (k_ ** 2 - 2) / (2 * k_ ** 2 - 2)
nu_hi_exact, nu_lo_exact = nu_of_kappa(mp.pi / u_hi), nu_of_kappa(2 * v_lo / mp.pi)
check("…the edges in closed form, from expanding Rayleigh–Lamb at the cutoffs: ω''ω_c = 1 − (8/πκ)tan(π/κ) (shear cutoff) and κ² + (16/πκ)cot(πκ/2) (stretch cutoff), κ = c_L/c_T",
      abs(nu_hi - nu_hi_exact) < 1e-6 and abs(nu_lo - nu_lo_exact) < 1e-6
      and abs(curvs[0.4] - (1 - 8 / (mp.pi * kappa_of(mp.mpf('0.4'))) * mp.tan(mp.pi / kappa_of(mp.mpf('0.4'))))) < 1e-4
      and abs(curvs[0.0] - (kappa_of(0) ** 2 + 16 / (mp.pi * kappa_of(0)) * mp.cot(mp.pi * kappa_of(0) / 2))) < 1e-4,
      f"u·tan u = π²/8 → ν = {mp.nstr(nu_hi_exact, 6)};  v³·tan v = −2π² → ν = {mp.nstr(nu_lo_exact, 6)}; Prada et al. (2008) found (−0.540, 0.451) by tracking the minimum")
# the Dirac point at nu = 1/3
k3 = kappa_of(mp.mpf(1) / 3)
slopes = []
for K in (mp.mpf("0.004"), mp.mpf("0.008")):
    r_lo = root_near(mp.pi - K / 3, K, k3, True, mp.mpf("0.02"))
    xs_ = [mp.pi - mp.mpf("0.03") + mp.mpf("0.06") * i / 1200 for i in range(1201)]
    vals = [F_lamb(x, K, k3, True) for x in xs_]
    roots = [mp.findroot(lambda O: F_lamb(O, K, k3, True), (xs_[i], xs_[i + 1]), solver="anderson") for i in range(1200) if vals[i] * vals[i + 1] < 0]
    roots = sorted(roots, key=lambda r: abs(r - mp.pi))[:2]
    slopes.append(sorted([(r - mp.pi) / K for r in roots]))
vD = [(sl[1] - sl[0]) / 2 for sl in slopes]
check("at ν = 1/3 the S1 and S2 branches leave the common cutoff linearly, ω − πc_T/h = ±(2c_T/π)k + O(k²): a Dirac cone, effective mass zero", slopes[0][0] < 0 < slopes[0][1] and abs(vD[0] - 2 / mp.pi) < 1e-5 and abs(vD[1] - 2 / mp.pi) < 1e-5,
      f"slopes {mp.nstr(slopes[1][0], 5)}, {mp.nstr(slopes[1][1], 5)}; half-difference {mp.nstr(vD[1], 8)} vs 2/π = {mp.nstr(2 / mp.pi, 8)} (degenerate perturbation: π³δ²/2 = 2πk²)")
# A1, the lightest massive mode (thickness shear, cutoff π c_T/2h for every ν)
def a1_curvature(nu, K=mp.mpf("0.01")):
    k_ = kappa_of(nu)
    Oc = mp.pi / 2
    vals = [2 * (root_near(Oc, KK, k_, False, mp.mpf("0.01")) - Oc) / KK ** 2 for KK in (K, 2 * K)]
    return (4 * vals[0] - vals[1]) / 3 * Oc


ok = True
for v in ("-0.5", "0", "0.3", "0.45"):
    k_ = kappa_of(mp.mpf(v))
    ok &= abs(a1_curvature(mp.mpf(v)) - (1 + 16 / (mp.pi * k_) * mp.cot(mp.pi / (2 * k_)))) < 1e-4
check("A1 (thickness shear, the lightest massive mode, ω_c = πc_T/2h for every ν): ω''ω_c = 1 + (16/πκ)cot(π/2κ) > 1", ok, f"ν = 0.3: {mp.nstr(a1_curvature(mp.mpf('0.3')), 6)} — inertia below the relativistic value m_c")
# A0 and S0 at small k, Rayleigh at large k
nu = mp.mpf("0.3")
kap = kappa_of(nu)
K = mp.mpf("0.005")
a0_pred = K ** 2 * mp.sqrt(mp.mpf(2) / (3 * (1 - nu)))
a0 = mp.findroot(lambda O: F_lamb(O, K, kap, False), a0_pred * mp.mpf("1.001"))
check("A0 flexural: ω h/c_T → (kh)²√(2/(3(1−ν))), i.e. a Schrödinger particle with m = (ħ/2c_T h)√(3(1−ν)/2)", abs(a0 / a0_pred - 1) < 1e-4, f"ν = 0.3: ratio {mp.nstr(a0 / a0_pred, 8)}")
K = mp.mpf("0.002")
s0 = mp.findroot(lambda O: F_lamb(O, K, kap, True), K * mp.sqrt(2 / (1 - nu)))
check("S0 extensional: ω/k → c_P = c_T√(2/(1−ν)) = √(E/ρ(1−ν²))", abs(s0 / K / mp.sqrt(2 / (1 - nu)) - 1) < 1e-5)


def rayleigh_xi(nu):
    k_ = kappa_of(nu)
    f = lambda x: (2 - x ** 2) ** 2 - 4 * mp.sqrt(1 - x ** 2) * mp.sqrt(1 - x ** 2 / k_ ** 2)
    return mp.findroot(f, (mp.mpf("0.6"), mp.mpf("0.999")), solver="bisect")


xiR = rayleigh_xi(mp.mpf("0.25"))
check("Rayleigh speed c_R/c_T = 0.919402 at ν = ¼ (the massless surface wave both lowest modes approach)", abs(xiR - mp.mpf("0.9194")) < 1e-4, mp.nstr(xiR, 7))
K = mp.mpf(12)
xi3 = rayleigh_xi(nu)
a0_hi = root_near(K * xi3, K, kap, False, mp.mpf("0.05"))
s0_hi = root_near(K * xi3, K, kap, True, mp.mpf("0.05"))
check("A0 and S0 phase speeds → c_R at large kh (ν = 0.3; the split is exponentially small, ~e^(−2kh√(1−c_R²/c_T²)))", abs(a0_hi / K - xi3) < 2e-4 and abs(s0_hi / K - xi3) < 2e-4, f"c_R/c_T = {mp.nstr(xi3, 7)}; at kh = 12: A0 {mp.nstr(a0_hi / K, 7)}, S0 {mp.nstr(s0_hi / K, 7)}")

# =====================================================================================
print("§3  The elastic torus: modular time is a mass ratio  [structure]")
# =====================================================================================
mp.mp.dps = 30


def stretch(theta, la, nu):
    """Uniaxial Hencky stretch: λ_a along angle θ, λ_t = λ_a^(−ν) across."""
    c, s = math.cos(theta), math.sin(theta)
    R = np.array([[c, -s], [s, c]])
    return R @ np.diag([la, la ** (-nu)]) @ R.T


def apply(F, z):
    v = F @ np.array([z.real, z.imag])
    return complex(v[0], v[1])


def hdist(z, w):
    return math.acosh(1 + abs(z - w) ** 2 / (2 * z.imag * w.imag))


ok_d = ok_A = True
worst = 0.0
for _ in range(200):
    tau = complex(rng.uniform(-0.5, 0.5), rng.uniform(0.3, 3))
    th, la, nu = rng.uniform(0, math.pi), math.exp(rng.uniform(-0.4, 0.4)), rng.uniform(-0.95, 0.95)
    F = stretch(th, la, nu)
    tau2 = apply(F, tau) / apply(F, 1)
    dH = hdist(tau, tau2)
    worst = max(worst, abs(dH - (1 + nu) * abs(math.log(la))))
    ok_A &= abs(math.log(abs(np.linalg.det(F))) - (1 - nu) * math.log(la)) < 1e-12
check("Hencky stretch λ_a along any direction: the shape point moves a hyperbolic distance (1 + ν)|ln λ_a| exactly, for any τ", worst < 1e-9, f"200 random (τ, θ, λ_a, ν): worst {worst:.1e} (the affine map is Teichmüller-extremal)")
check("…and the area by ln(A'/A) = (1 − ν) ln λ_a exactly: shape : size = (1 + ν) : (1 − ν)", ok_A, "μ/K₂ = (1 − ν)/(1 + ν) in 2D")
ok = True
for _ in range(20):
    tau = complex(rng.uniform(-0.5, 0.5), rng.uniform(0.3, 3))
    th, nu, h = rng.uniform(0, math.pi), rng.uniform(-0.9, 0.9), 1e-7
    F = stretch(th, math.exp(h), nu)
    dtau = (apply(F, tau) / apply(F, 1) - tau) / h
    ok &= abs(dtau - (-1j * tau.imag * (1 + nu) * complex(math.cos(2 * th), math.sin(2 * th)))) < 1e-5
check("first order: δτ = −i y (1 + ν) e^{2iθ} δ ln λ_a", ok)
ok = True
for _ in range(20):
    tau = complex(rng.uniform(-0.5, 0.5), rng.uniform(0.4, 2))
    b = complex(*rng.standard_normal(2)) * 1e-5
    Fb = np.eye(2) + np.array([[b.real, b.imag], [b.imag, -b.real]])     # u = b z̄, pure shear
    tau2 = apply(Fb, tau) / apply(Fb, 1)
    mu_, A_ = 1.7, 2.3
    E_exact = A_ * 2 * mu_ * abs(b) ** 2
    E_hyp = mu_ * A_ / 2 * abs(tau2 - tau) ** 2 / tau.imag ** 2
    ok &= abs(E_hyp / E_exact - 1) < 1e-4
check("the elastic energy of a shape change is (μA/2)|δτ|²/y²: the Poincaré metric is the torus's shear modulus", ok)


def kk_masses(tau, A, nmax=6):
    """√eigenvalues of −Δ on C/Λ, Λ = L(Z + τZ), area A, from the dual lattice."""
    y = tau.imag
    L = math.sqrt(A / y)
    w1, w2 = np.array([L, 0.0]), np.array([L * tau.real, L * y])
    Bm = np.array([w1, w2])
    D = np.linalg.inv(Bm).T                # dual basis rows
    out = {}
    for n in range(-nmax, nmax + 1):
        for m in range(-nmax, nmax + 1):
            if n == 0 and m == 0:
                continue
            k = n * D[0] + m * D[1]
            out[(n, m)] = 2 * math.pi * np.linalg.norm(k)
    return out


tau, A = complex(0.31, 1.27), 2.4
ms = kk_masses(tau, A)
ok = all(abs(ms[(n, m)] - 2 * math.pi * abs(m - n * tau) / math.sqrt(A * tau.imag)) < 1e-12 for (n, m) in ms)
L = math.sqrt(A / tau.imag)
ok2 = all(abs(ms[(n, m)] - 2 * math.pi * L * abs(m - n * tau) / A) < 1e-12 for (n, m) in ms)
check("KK masses m_(n,m) = 2π|m − nτ|/√(A y) = (2π/A)·ℓ, ℓ the length of a closed geodesic (a torus knot on the Clifford torus)", ok and ok2)
y0 = 2.3
ms = kk_masses(complex(0, y0), 1.0)
check("rectangular torus: modular time is a mass ratio, y = m_(1,0)/m_(0,1)", abs(ms[(1, 0)] / ms[(0, 1)] - y0) < 1e-12)
ok = True
for nu in (-0.5, 0.0, 0.3, 0.49):
    for la in (1.01, 1.3, 0.8):
        F = stretch(0.0, la, nu)
        tau2 = apply(F, complex(0, y0)) / apply(F, 1)
        A2 = abs(np.linalg.det(F))
        m2 = kk_masses(tau2, A2)
        r = -math.log(m2[(0, 1)] / ms[(0, 1)]) / math.log(m2[(1, 0)] / ms[(1, 0)])
        ok &= abs(r - nu) < 1e-10
check("Hencky's ν is minus the ratio of the log-shifts of the two KK masses: ν = −δln m_(0,1) / δln m_(1,0), exactly", ok)


def sigma_k(n, k):
    return mp.fsum(mp.mpf(d) ** k for d in range(1, n + 1) if n % d == 0)


def xi(s):
    return mp.pi ** (-s / 2) * mp.gamma(s / 2) * mp.zeta(s)


def f_epstein(tau, s, nmax=None):
    """2ζ(2s)E(τ,s) = Σ'(y/|mτ+n|²)^s, analytic near s = 0 in this form."""
    x, y = mp.re(tau), mp.im(tau)
    nmax = nmax or int(mp.ceil(80 / (2 * mp.pi * y))) + 5
    acc = mp.fsum(mp.mpf(n) ** (s - mp.mpf(1) / 2) * sigma_k(n, 1 - 2 * s) * mp.besselk(s - mp.mpf(1) / 2, 2 * mp.pi * n * y) * 2 * mp.cos(2 * mp.pi * n * x) for n in range(1, nmax + 1))
    return 2 * mp.zeta(2 * s) * y ** s + 2 * mp.pi ** s * mp.rgamma(s) * (xi(2 * s - 1) * y ** (1 - s) + 2 * mp.sqrt(y) * acc)


def log_eta(tau):
    q = mp.exp(2j * mp.pi * tau)
    return 1j * mp.pi * tau / 12 + mp.fsum(mp.log(1 - q ** n) for n in range(1, 400))


ok = True
for tau, A in ((mp.mpc(0, 1), mp.mpf(1)), (mp.mpc("0.3", "1.1"), mp.mpf("2.5"))):
    zp = -mp.log(A / (4 * mp.pi ** 2)) + mp.diff(lambda s: f_epstein(tau, s), 0)
    target = mp.log(A * mp.im(tau)) + 4 * mp.re(log_eta(tau))
    ok &= abs(-zp - target) < 1e-15 and abs(f_epstein(tau, mp.mpf("1e-25")) + 1) < 1e-20
check("the KK spectrum's zeta is Instrument 19's Eisenstein series, ζ_Δ(s) = (A/4π²)^s·2ζ(2s)E(τ,s), and det'Δ = A·y·|η(τ)|⁴", ok, "−ζ'_Δ(0) by numerical differentiation at 30 digits; ζ_Δ(0) = −1")
# 143a1's own torus
c4, c6, Delta = 64, 1864, -1859
check("143a1: c4 = 64, c6 = 1864, Δ = −1859 = −11·13², 1728Δ = c4³ − c6²", 1728 * Delta == c4 ** 3 - c6 ** 2 and Delta == -11 * 13 ** 2)
jE = mp.mpf(c4) ** 3 / Delta
t143 = mp.findroot(lambda t: mp.re(mp.kleinj(mp.mpc(mp.mpf(1) / 2, t))) - jE / 1728, mp.mpf("1.08"))
tau143 = mp.mpc(mp.mpf(1) / 2, t143)
check("its period lattice has τ = ½ + i·t on the edge of the fundamental domain (Δ < 0: a rhombic torus)", abs(mp.kleinj(tau143) - jE / 1728) < 1e-20, f"j = −2¹⁸/1859 = {mp.nstr(jE, 10)}, τ = ½ + {mp.nstr(t143, 12)} i; mass ratio |τ| = {mp.nstr(abs(tau143), 8)}, the pair (0,1),(1,1) degenerate")

# =====================================================================================
print("§4  The labile aether: the graviton's only mass, gravity's own kinetic term  [structure]")
# =====================================================================================
G = 0.5 * (np.einsum("ik,jl->ijkl", np.eye(3), np.eye(3)) + np.einsum("il,jk->ijkl", np.eye(3), np.eye(3))) - np.einsum("ij,kl->ijkl", np.eye(3), np.eye(3))
B3 = sym_basis(3)
Gm = np.array([[np.einsum("ij,ijkl,kl->", Ba, G, Bb) for Bb in B3] for Ba in B3])
evG = np.sort(np.linalg.eigvalsh(Gm))
uni = np.zeros((3, 3))
uni[0, 0] = 1
check("DeWitt's supermetric G = ½(δδ + δδ) − δδ: shear +1 (×5), bulk −2 (×1), uniaxial 0 — μ = ½, λ = −1, M = 0, ν = 1", np.allclose(evG, [-2, 1, 1, 1, 1, 1]) and abs(np.einsum("ij,ijkl,kl->", uni, G, uni)) < 1e-15, "the conformal mode's negative stiffness is the conformal-factor problem")
ok = True
for Kd, Rexp in ((np.eye(3) / 2.0, 6 / 4.0), (np.diag([1 / 2.0, 1 / 2.0, 0]), 2 / 4.0), (np.diag([1 / 2.0, 0, 0]), 0.0)):
    ok &= abs(np.trace(Kd) ** 2 - np.sum(Kd * Kd) - Rexp) < 1e-15
check("Gauss's equation R = K² − K_ijK^ij = 2σ₂(K): S³ (6/r²), S²×ℝ (2/r²), S¹×ℝ² (0) — uniaxial bending leaves space flat", ok)
# Gauss for an ellipsoid point: K_G from the shape operator vs the closed formula
aa, bb, cc = 1.3, 0.9, 0.6
uu, vv = 0.7, 1.1
Ssym = sp.symbols("U V")
Pt = sp.Matrix([aa * sp.cos(Ssym[0]) * sp.sin(Ssym[1]), bb * sp.sin(Ssym[0]) * sp.sin(Ssym[1]), cc * sp.cos(Ssym[1])])
Pu, Pv = Pt.diff(Ssym[0]), Pt.diff(Ssym[1])
Nn = Pu.cross(Pv)
Nn = Nn / sp.sqrt(Nn.dot(Nn))
I1 = sp.Matrix([[Pu.dot(Pu), Pu.dot(Pv)], [Pu.dot(Pv), Pv.dot(Pv)]])
II = sp.Matrix([[Pt.diff(Ssym[0], 2).dot(Nn), Pt.diff(Ssym[0], Ssym[1]).dot(Nn)], [Pt.diff(Ssym[0], Ssym[1]).dot(Nn), Pt.diff(Ssym[1], 2).dot(Nn)]])
Sh = (I1.inv() * II).subs({Ssym[0]: uu, Ssym[1]: vv}).evalf(30)
KG = Sh.det()
xyz = [aa * math.cos(uu) * math.sin(vv), bb * math.sin(uu) * math.sin(vv), cc * math.cos(vv)]
KG_formula = 1 / (aa * bb * cc) ** 2 / (xyz[0] ** 2 / aa ** 4 + xyz[1] ** 2 / bb ** 4 + xyz[2] ** 2 / cc ** 4) ** 2
check("…and on a surface σ₂ of the shape operator is Gauss's curvature (ellipsoid point vs the closed formula)", abs(float(KG) - KG_formula) < 1e-12, f"K = {KG_formula:.10f}")

# Fierz–Pauli with a general trace coefficient
w, p, m = sp.symbols("omega p m", positive=True)
a = sp.Symbol("a")
eta = sp.diag(-1, 1, 1, 1)
hs = {}
for i in range(4):
    for j in range(i, 4):
        hs[(i, j)] = sp.Symbol(f"h{i}{j}")
h = sp.Matrix(4, 4, lambda i, j: hs[(min(i, j), max(i, j))])
hup = eta * h * eta


def fp_lagrangian(kup, a_val, m_val):
    kdn = [sum(eta[i, j] * kup[j] for j in range(4)) for i in range(4)]
    htr = sum(eta[i, i] * h[i, i] for i in range(4))
    T1 = -sp.Rational(1, 2) * sum(kdn[l] * kup[l] for l in range(4)) * sum(h[i, j] * hup[i, j] for i in range(4) for j in range(4))
    T2 = sum(sum(kdn[i] * hup[i, l] for i in range(4)) * sum(kup[j] * h[j, l] for j in range(4)) for l in range(4))
    T3 = -sum(sum(kdn[i] * hup[i, j] for i in range(4)) * kdn[j] * htr for j in range(4))
    T4 = sp.Rational(1, 2) * sum(kdn[l] * kup[l] for l in range(4)) * htr ** 2
    M_ = -sp.Rational(1, 2) * sp.sympify(m_val) ** 2 * (sum(h[i, j] * hup[i, j] for i in range(4) for j in range(4)) - a_val * htr ** 2)
    return sp.expand(T1 + T2 + T3 + T4 + M_)


syms = list(hs.values())
L0 = fp_lagrangian([w, 0, 0, 0], a, m)
# zero momentum: spin 2, h0i, and the scalar sector
s1_, t_, n_ = sp.symbols("s t n")
sub = {hs[(1, 2)]: s1_, hs[(0, 0)]: n_, hs[(1, 1)]: t_ / 3, hs[(2, 2)]: t_ / 3, hs[(3, 3)]: t_ / 3}
sub.update({hs[k]: 0 for k in hs if k not in ((1, 2), (0, 0), (1, 1), (2, 2), (3, 3))})
Lr = sp.expand(L0.subs(sub))
L_spin2 = sp.expand(Lr.coeff(s1_, 2))
check("zero momentum: the EH kinetic term is ½(ḣ_ijḣ_ij − ḣ_ii²) — DeWitt's labile form; spin 2 has L = (ω² − m²)s² (healthy, mass m)", sp.simplify(L_spin2 - (w ** 2 - m ** 2)) == 0)
Lsc = sp.expand(Lr - L_spin2 * s1_ ** 2)
n_star = sp.solve(sp.diff(Lsc, n_), n_)[0]
Lt = sp.simplify(Lsc.subs(n_, n_star))
mg2 = m ** 2 * (4 * a - 1) / (2 * (1 - a))
check("a ≠ 1: h₀₀ is algebraic; eliminating it leaves L_t = −⅓(ω² − m_g²)t², a ghost with m_g² = m²(4a − 1)/(2(1 − a))", sp.simplify(Lt + sp.Rational(1, 3) * (w ** 2 - mg2) * t_ ** 2) == 0, "Hinterbichler's (3 − 4a')/(2a') with a' = 1 − a")
L_fp = sp.expand(Lsc.subs(a, 1))
check("a = 1 (Fierz–Pauli): h₀₀ enters linearly — a Lagrange multiplier imposing t = 0, the ghost removed", sp.degree(L_fp, n_) == 1 and sp.simplify(sp.diff(L_fp, n_) + m ** 2 * t_) == 0)
# full momentum, Lorentz-invariant dispersion
a_val = sp.Rational(3, 10)
Lk = fp_lagrangian([w, 0, 0, p], a_val, 1)
Mk = sp.hessian(Lk, syms)
detk = sp.factor(Mk.det())
mg2_val = (4 * a_val - 1) / (2 * (1 - a_val))
fac = sp.factor_list(detk)[1]
facs = {sp.expand(f): k for f, k in fac}
want2 = sp.expand(w ** 2 - p ** 2 - 1)
want0 = sp.expand(7 * w ** 2 - 7 * p ** 2 - 1)        # m_g^2 = 1/7 at a = 3/10
got2 = any(sp.simplify(f - want2) == 0 or sp.simplify(f + want2) == 0 for f in facs if facs[f] == 5)
got0 = any(sp.simplify(f - want0) == 0 or sp.simplify(f + want0) == 0 for f in facs)
check("at any momentum (a = 3/10): det = c·(ω² − p² − m²)⁵(ω² − p² − m_g²), m_g² = 1/7 — five spin-2 modes and one ghost, both Lorentz-invariant", got2 and got0, str([f"{sp.factor(f)}^{k}" for f, k in fac]))
Mmod = sp.Symbol("M")
check("in elastic terms (μ = 1, λ = −2a, M = 2(1 − a)): m_g² = m²(3μ/M − 2), and in the 3D Poisson ratio m_g² = −m²(1 + 2ν)/(2(1 − ν))",
      sp.simplify(mg2 - m ** 2 * (3 / (2 * (1 - a)) - 2)) == 0 and sp.simplify(mg2.subs(a, nu_s / (2 * nu_s - 1)) + m ** 2 * (1 + 2 * nu_s) / (2 * (1 - nu_s))) == 0,
      "M → 0 (the labile aether) sends the ghost to infinite mass")
lH = sp.Symbol("lambda_H")
lamL, muL = -lH, sp.Rational(1, 2)
K3 = lamL + 2 * muL / 3
M3 = lamL + 2 * muL
nu3 = sp.simplify(lamL / (2 * (lamL + muL)))
check("Hořava's K_ijK^ij − λK²: the scalar graviton's kinetic factor (3λ−1)/(λ−1) = 3K/M = (1 + ν)/(1 − ν)", sp.simplify((3 * lH - 1) / (lH - 1) - 3 * K3 / M3) == 0 and sp.simplify(3 * K3 / M3 - ((1 + nu3) / (1 - nu3))) == 0,
      "ghost-free ⇔ −1 < ν < 1 (λ < ⅓ or λ > 1); GR λ = 1 is ν = 1; the Weyl point λ = ⅓ is ν = −1")
check("…and the IR sound speed c_s² = −(λ−1)/(3λ−1) = −(1 − ν)/(1 + ν): kinetic factor × c_s² = −1 for every ν", sp.simplify(-(lH - 1) / (3 * lH - 1) + (1 - nu3) / (1 + nu3)) == 0)
hbar, c_, Gn, kB = 1.054571817e-34, 299792458.0, 6.67430e-11, 1.380649e-23
FP = c_ ** 4 / Gn
mP = math.sqrt(hbar * c_ / Gn)
check("Planck force c⁴/G = 1.2103e44 N is spacetime's stiffness (ADM: (c⁴/16πG)(R + K_ijK^ij − K²)); m_P = √(ħ F_P/c³)", abs(FP / 1.2103e44 - 1) < 2e-4 and abs(math.sqrt(hbar * FP / c_ ** 3) / mP - 1) < 1e-12, f"F_P = {FP:.5e} N, m_P = {mP:.5e} kg")

# =====================================================================================
print("§5  Rods and knots: twist/bend = 1/(1 + ν), Michell's threshold, Lk = Tw + Wr as a degree  [structure]")
# =====================================================================================
nu_r = sp.Symbol("nu")
Ey = sp.Symbol("E", positive=True)
r_ = sp.Symbol("r", positive=True)
I_ = sp.pi * r_ ** 4 / 4
J_ = sp.pi * r_ ** 4 / 2
check("circular cross-section: C/B = GJ/(EI) = 2G/E = 1/(1 + ν) — Kirchhoff's top has I₃/I₁ = 1/(1 + ν)", sp.simplify((Ey / (2 * (1 + nu_r)) * J_) / (Ey * I_) - 1 / (1 + nu_r)) == 0)


def ring_curve(params, n, N):
    al, be, ga, de = params
    th = np.linspace(0, 2 * np.pi, N, endpoint=False)
    r = 1 + al * np.cos(n * th) + be * np.sin(n * th)
    r1 = n * (-al * np.sin(n * th) + be * np.cos(n * th))
    r2 = -n * n * (al * np.cos(n * th) + be * np.sin(n * th))
    z = ga * np.cos(n * th) + de * np.sin(n * th)
    z1 = n * (-ga * np.sin(n * th) + de * np.cos(n * th))
    z2 = -n * n * (ga * np.cos(n * th) + de * np.sin(n * th))
    c, s = np.cos(th), np.sin(th)
    P = np.stack([r * c, r * s, z], 1)
    D1 = np.stack([r1 * c - r * s, r1 * s + r * c, z1], 1)
    D2 = np.stack([r2 * c - 2 * r1 * s - r * c, r2 * s + 2 * r1 * c - r * s, z2], 1)
    return P, D1, D2, 2 * np.pi / N


def writhe(P, D1, dth):
    N = len(P)
    tot = 0.0
    for i0 in range(0, N, 256):
        Pi, Ti = P[i0:i0 + 256, None, :], D1[i0:i0 + 256, None, :]
        R = Pi - P[None, :, :]
        cr = np.cross(Ti, D1[None, :, :])
        num = np.einsum("ijk,ijk->ij", cr, R)
        den = np.linalg.norm(R, axis=2) ** 3
        idx = np.arange(i0, min(i0 + 256, N))
        den[np.arange(len(idx)), idx] = np.inf
        tot += np.sum(num / den)
    return tot * dth * dth / (4 * np.pi)


def bend_energy(P, D1, D2, dth):
    sp_ = np.linalg.norm(D1, axis=1)
    kap2 = np.sum(np.cross(D1, D2) ** 2, axis=1) / sp_ ** 6
    L = np.sum(sp_) * dth
    return 0.5 * np.sum(kap2 * sp_) * dth * L / (2 * np.pi)     # rescaled to length 2π, B = 1


def hessians(n, N=480, hstep=2e-3):
    def ev(pv):
        P, D1, D2, dth = ring_curve(pv, n, N)
        return bend_energy(P, D1, D2, dth), writhe(P, D1, dth)
    e0 = np.zeros(4)
    Hb, Hw = np.zeros((4, 4)), np.zeros((4, 4))
    f0 = ev(e0)
    for i in range(4):
        ei = np.eye(4)[i] * hstep
        fp_, fm_ = ev(ei), ev(-ei)
        Hb[i, i] = (fp_[0] - 2 * f0[0] + fm_[0]) / hstep ** 2
        Hw[i, i] = (fp_[1] - 2 * f0[1] + fm_[1]) / hstep ** 2
        for j in range(i + 1, 4):
            ej = np.eye(4)[j] * hstep
            fpp, fpm, fmp, fmm = ev(ei + ej), ev(ei - ej), ev(-ei + ej), ev(-ei - ej)
            Hb[i, j] = Hb[j, i] = (fpp[0] - fpm[0] - fmp[0] + fmm[0]) / (4 * hstep ** 2)
            Hw[i, j] = Hw[j, i] = (fpp[1] - fpm[1] - fmp[1] + fmm[1]) / (4 * hstep ** 2)
    return Hb, Hw


def michell_threshold(n):
    Hb, Hw = hessians(n)
    Lc = np.linalg.cholesky(Hb)
    Li = np.linalg.inv(Lc)
    mu_max = np.max(np.abs(np.linalg.eigvalsh(Li @ Hw @ Li.T)))
    t_c = 1 / mu_max                       # Hb − t·Hw singular; t = 4π²C·Lk/L₀, L₀ = 2π, B = 1
    return t_c / (2 * np.pi)               # = (C/B)·Lk_c


thr2, thr3 = michell_threshold(2), michell_threshold(3)
check("twisted ring, reduced energy (B/2)∫κ² + (C/2L)(2π(Lk − Wr))²: the circle loses stability to the n = 2 writhing mode at (C/B)·Lk = √3 (Michell 1889)", abs(thr2 - math.sqrt(3)) < 2e-3, f"computed {thr2:.5f} vs √3 = {math.sqrt(3):.5f}; with C/B = 1/(1+ν): Lk_c = √3(1 + ν)")
check("…the n = 3 mode at (C/B)·Lk = √8 = √(n² − 1)", abs(thr3 - math.sqrt(8)) < 4e-3, f"computed {thr3:.5f} vs √8 = {math.sqrt(8):.5f}")


def torus_knot(pp, qq, Rt, rt, N, eps=0.0):
    t = np.linspace(0, 2 * np.pi, N, endpoint=False)
    nvec = np.stack([np.cos(qq * t) * np.cos(pp * t), np.cos(qq * t) * np.sin(pp * t), np.sin(qq * t)], 1)
    P = np.stack([(Rt + rt * np.cos(qq * t)) * np.cos(pp * t), (Rt + rt * np.cos(qq * t)) * np.sin(pp * t), rt * np.sin(qq * t)], 1)
    D = np.stack([-pp * (Rt + rt * np.cos(qq * t)) * np.sin(pp * t) - rt * qq * np.sin(qq * t) * np.cos(pp * t),
                  pp * (Rt + rt * np.cos(qq * t)) * np.cos(pp * t) - rt * qq * np.sin(qq * t) * np.sin(pp * t),
                  rt * qq * np.cos(qq * t)], 1)
    dn = np.stack([-qq * np.sin(qq * t) * np.cos(pp * t) - pp * np.cos(qq * t) * np.sin(pp * t),
                   -qq * np.sin(qq * t) * np.sin(pp * t) + pp * np.cos(qq * t) * np.cos(pp * t), qq * np.cos(qq * t)], 1)
    return P + eps * nvec, D + eps * dn, nvec, dn, 2 * np.pi / N


def gauss_link(P1, D1, P2, D2, dt1, dt2):
    tot = 0.0
    for i0 in range(0, len(P1), 400):
        R = P1[i0:i0 + 400, None, :] - P2[None, :, :]
        cr = np.cross(D1[i0:i0 + 400, None, :], D2[None, :, :])
        tot += np.sum(np.einsum("ijk,ijk->ij", cr, R) / np.linalg.norm(R, axis=2) ** 3)
    return tot * dt1 * dt2 / (4 * np.pi)


def crossing_degree(P1, P2):
    """lk = deg of Ψ = (γ₂ − γ₁)/|γ₂ − γ₁|: signed preimages of ẑ, i.e. crossings with γ₂ over γ₁, sign of (γ₂' × γ₁')·ẑ."""
    A1, B1 = P1, np.roll(P1, -1, 0)
    A2, B2 = P2, np.roll(P2, -1, 0)
    total, count = 0, 0
    for i0 in range(0, len(A1), 300):
        a1, b1 = A1[i0:i0 + 300, None, :], B1[i0:i0 + 300, None, :]
        d1 = b1 - a1
        d2 = (B2 - A2)[None, :, :]
        rr = A2[None, :, :] - a1
        den = d1[..., 0] * d2[..., 1] - d1[..., 1] * d2[..., 0]
        with np.errstate(divide="ignore", invalid="ignore"):
            s = (rr[..., 0] * d2[..., 1] - rr[..., 1] * d2[..., 0]) / den
            u_ = (rr[..., 0] * d1[..., 1] - rr[..., 1] * d1[..., 0]) / den
        hit = (s >= 0) & (s < 1) & (u_ >= 0) & (u_ < 1)
        ii, jj = np.nonzero(hit)
        for i, j in zip(ii, jj):
            z1 = a1[i, 0, 2] + s[i, j] * d1[i, 0, 2]
            z2 = A2[j, 2] + u_[i, j] * d2[0, j, 2]
            if z2 > z1:
                count += 1
                total += int(np.sign(d2[0, j, 0] * d1[i, 0, 1] - d2[0, j, 1] * d1[i, 0, 0]))
    return total, count


Pk, Dk, nk, dnk, dtk = torus_knot(2, 3, 2.0, 0.8, 3000)
Pe, De, _, _, _ = torus_knot(2, 3, 2.0, 0.8, 3000, eps=0.06)
Lk = gauss_link(Pk, Dk, Pe, De, dtk, dtk)
Wr = writhe(Pk, Dk, dtk)
Tw = np.sum(np.einsum("ij,ij->i", np.cross(Dk / np.linalg.norm(Dk, axis=1)[:, None], nk), dnk)) * dtk / (2 * np.pi)
deg, ncross = crossing_degree(Pk, Pe)
check("trefoil T(2,3) with its torus framing: |Lk| (Gauss integral) = pq = 6, the sign is the handedness", abs(abs(Lk) - 6) < 2e-3, f"Lk = {Lk:.5f} (this parametrisation is the left-handed trefoil)")
check("Călugăreanu–White–Fuller: Lk = Tw + Wr", abs(Lk - (Tw + Wr)) < 2e-3, f"Tw = {Tw:.5f}, Wr = {Wr:.5f}, sum {Tw + Wr:.5f}")
check("Lk is the degree of the Gauss map (γ₂ − γ₁)/|γ₂ − γ₁|: signed preimages of a regular value give the integer directly", deg == round(Lk), f"{ncross} preimages of ẑ, signed sum {deg}")

# =====================================================================================
print("§6  Topology has a rest energy: shells, Willmore tori, dislocations  [structure]")
# =====================================================================================
k1, k2, Dp = sp.symbols("kappa1 kappa2 D", real=True)
Wplate = Dp / 2 * ((k1 + k2) ** 2 - 2 * (1 - nu_r) * k1 * k2)
Hq = sp.hessian(Wplate, (k1, k2))
check("Kirchhoff plate: (D/2)[(κ₁+κ₂)² − 2(1−ν)κ₁κ₂] = (D/2)(tr κ)² − D(1−ν)σ₂(κ): Gaussian modulus κ̄ = −D(1 − ν)", sp.expand(Wplate - (Dp / 2 * (k1 + k2) ** 2 - Dp * (1 - nu_r) * k1 * k2)) == 0)
check("…positive-definite ⇔ |ν| < 1 — the 2D Poisson bounds are Helfrich's window −2κ < κ̄ < 0", sp.simplify(Hq.det() - Dp ** 2 * (1 - nu_r ** 2)) == 0)


def surface_integrals(Xfun, N=400):
    u = np.linspace(0, 2 * np.pi, N, endpoint=False)
    v = np.linspace(0, 2 * np.pi, N, endpoint=False) + np.pi / N      # midpoints: no grid point on a pole
    U, V = np.meshgrid(u, v, indexing="ij")
    Xs, Xu, Xv, Xuu, Xuv, Xvv = Xfun(U, V)
    E_, F_, G_ = (Xu * Xu).sum(0), (Xu * Xv).sum(0), (Xv * Xv).sum(0)
    nrm = np.cross(Xu, Xv, axis=0)
    W_ = np.linalg.norm(nrm, axis=0)
    nrm = nrm / W_
    e_, f_, g_ = (Xuu * nrm).sum(0), (Xuv * nrm).sum(0), (Xvv * nrm).sum(0)
    H = (e_ * G_ - 2 * f_ * F_ + g_ * E_) / (2 * (E_ * G_ - F_ ** 2))
    K = (e_ * g_ - f_ ** 2) / (E_ * G_ - F_ ** 2)
    dA = W_ * (2 * np.pi / N) ** 2
    return np.sum(H ** 2 * dA), np.sum(K * dA), np.sum(dA)


def homogeneous_torus(yv):
    r1, r2 = 1 / math.sqrt(1 + yv * yv), yv / math.sqrt(1 + yv * yv)
    U_, V_ = sp.symbols("U V")
    Xs4 = sp.Matrix([r1 * sp.cos(U_), r1 * sp.sin(U_), r2 * sp.cos(V_), r2 * sp.sin(V_)])
    P3 = sp.Matrix([Xs4[0], Xs4[1], Xs4[2]]) / (1 - Xs4[3])
    exprs = [P3, P3.diff(U_), P3.diff(V_), P3.diff(U_, 2), P3.diff(U_, V_), P3.diff(V_, 2)]
    fns = [sp.lambdify((U_, V_), list(ex), "numpy") for ex in exprs]

    def Xfun(U, V):
        return [np.array([np.broadcast_to(c, U.shape) for c in f(U, V)], dtype=float) for f in fns]
    return Xfun


ok = True
rows = []
for yv in (1.0, 1.5, 2.0, 0.5, 3.0):
    W_, Kint, _ = surface_integrals(homogeneous_torus(yv))
    rows.append(f"y={yv}: {W_:.6f} vs {math.pi ** 2 * (yv + 1 / yv):.6f}")
    ok &= abs(W_ - math.pi ** 2 * (yv + 1 / yv)) < 1e-6 and abs(Kint) < 1e-6
check("stereographic image of the homogeneous torus of modular class τ = iy: ∫H²dA = π²(y + 1/y), ∫K dA = 0", ok, "; ".join(rows[:3]))
print("    so a toroidal shell's bending energy along this family is 2D·π²(y + 1/y): S-invariant, least at the square torus (Clifford, 4π²D; Marques–Neves: W ≥ 2π² for every genus ≥ 1).")


def ellipsoid(U, V, a=1.4, b=1.0, c=0.7):
    # polar angle v/2 in (0, π): the ellipsoid once
    X = np.array([a * np.cos(U) * np.sin(V / 2), b * np.sin(U) * np.sin(V / 2), c * np.cos(V / 2)])
    Xu = np.array([-a * np.sin(U) * np.sin(V / 2), b * np.cos(U) * np.sin(V / 2), 0 * U])
    Xv = np.array([a * np.cos(U) * np.cos(V / 2) / 2, b * np.sin(U) * np.cos(V / 2) / 2, -c * np.sin(V / 2) / 2])
    Xuu = np.array([-a * np.cos(U) * np.sin(V / 2), -b * np.sin(U) * np.sin(V / 2), 0 * U])
    Xuv = np.array([-a * np.sin(U) * np.cos(V / 2) / 2, b * np.cos(U) * np.cos(V / 2) / 2, 0 * U])
    Xvv = np.array([-a * np.cos(U) * np.sin(V / 2) / 4, -b * np.sin(U) * np.sin(V / 2) / 4, -c * np.cos(V / 2) / 4])
    return X, Xu, Xv, Xuu, Xuv, Xvv


Wel, Kel, _ = surface_integrals(ellipsoid, N=600)
check("Gauss–Bonnet / Hopf: ∫K dA = 2πχ = 4π·deg(Gauss map) — ellipsoid 4π (χ = 2), torus 0 (χ = 0)", abs(Kel - 4 * math.pi) < 1e-3, f"ellipsoid ∫K = {Kel:.6f} = 4π·{Kel / (4 * math.pi):.6f}")
gg = sp.Symbol("g", integer=True)
U_shell = 2 * Dp * sp.Symbol("W") - 2 * sp.pi * Dp * (1 - nu_r) * (2 - 2 * gg)
check("closed shell: U = 2D·W + 4πD(1−ν)(g − 1); sphere (W = 4π, g = 0): 4πD(1 + ν); each handle costs 4πD(1 − ν)", sp.simplify(U_shell.subs({sp.Symbol("W"): 4 * sp.pi, gg: 0}) - 4 * sp.pi * Dp * (1 + nu_r)) == 0 and sp.simplify(U_shell.subs(gg, gg + 1) - U_shell - 4 * sp.pi * Dp * (1 - nu_r)) == 0,
      "genus 13 (X₀(143)): the topological part is 48πD(1 − ν) — the same 48π as its hyperbolic area, both −2πχ")
# dislocations
xs_, ys_, ts_, b_, v_, cT = sp.symbols("x y t b v c_T", positive=True)
gam = 1 / sp.sqrt(1 - v_ ** 2 / cT ** 2)
uz = b_ / (2 * sp.pi) * sp.atan2(ys_, gam * (xs_ - v_ * ts_))
wave = sp.diff(uz, xs_, 2) + sp.diff(uz, ys_, 2) - sp.diff(uz, ts_, 2) / cT ** 2
check("moving screw dislocation u_z = (b/2π)·arctan(y/γ(x − vt)), γ = (1 − v²/c_T²)^(−½), solves the wave equation exactly (Frank 1949): Lorentz-contracted along the motion at c_T", sp.simplify(wave) == 0)
thg = np.linspace(0, 2 * np.pi, 20001)
uz0 = lambda th: 0.37 / (2 * np.pi) * np.unwrap(np.arctan2(np.sin(th), np.cos(th)))
check("its Burgers vector ∮du = b is a winding number — the degree of the angle map S¹ → S¹, times b", abs((uz0(thg)[-1] - uz0(thg)[0]) - 0.37) < 1e-12)


def edge_energy_per_lnr(nu, mu=1.0, b=1.0, N=200000):
    th = np.linspace(0, 2 * np.pi, N, endpoint=False)
    x, y = np.cos(th), np.sin(th)
    Dd = mu * b / (2 * np.pi * (1 - nu))
    r4 = (x * x + y * y) ** 2
    sxx = -Dd * y * (3 * x * x + y * y) / r4
    syy = Dd * y * (x * x - y * y) / r4
    sxy = Dd * x * (x * x - y * y) / r4
    szz = nu * (sxx + syy)
    tr_ = sxx + syy + szz
    exx = (sxx - nu / (1 + nu) * tr_) / (2 * mu)
    eyy = (syy - nu / (1 + nu) * tr_) / (2 * mu)
    ezz = (szz - nu / (1 + nu) * tr_) / (2 * mu)
    exy = sxy / (2 * mu)
    Wd = 0.5 * (sxx * exx + syy * eyy + szz * ezz + 2 * sxy * exy)
    return np.sum(Wd) * 2 * np.pi / N


ok = True
for nu in (0.0, 0.25, 1 / 3, 0.45):
    ok &= abs(edge_energy_per_lnr(nu) - 1 / (4 * np.pi * (1 - nu))) < 1e-9
check("edge dislocation: dE/d ln r = μb²/(4π(1 − ν)) vs screw μb²/4π — the edge particle is 1/(1 − ν) times heavier", ok, "Poisson's ratio is the mass ratio of the two topological particles of a crystal")

# =====================================================================================
print("§7  The programme's bridge under stress  [model / audit]")
# =====================================================================================
mp.mp.dps = 25
yc = mp.mpf("0.181256")
alpha_inv = mp.mpf("137.035999084")
kexp = 2 * alpha_inv / mp.pi
E2 = 1 - 24 * mp.fsum(sigma_k(n, 1) * mp.exp(-2 * mp.pi * n * yc) for n in range(1, 500))
E2_mod = (6 * yc / mp.pi - (1 - 24 * mp.fsum(sigma_k(n, 1) * mp.exp(-2 * mp.pi * n / yc) for n in range(1, 20)))) / yc ** 2
check("E₂(i·y_c) two ways (q-series; quasi-modular transform to i/y_c)", abs(E2 - E2_mod) < 1e-15, f"E₂ = {mp.nstr(E2, 10)}")
coef = kexp * mp.pi * yc * E2 / 12
for nu in (0.0, 0.3):
    hh = 1e-7
    F = stretch(0.0, math.exp(hh), nu)
    tau2 = apply(F, complex(0, float(yc))) / apply(F, 1)
    fd = kexp * (mp.re(log_eta(mp.mpc(tau2.real, tau2.imag))) - mp.re(log_eta(mp.mpc(0, yc)))) / hh
    check(f"d ln(mₑ/M_P)/d ln λ_a = k(πy/12)E₂(iy)(1 + ν) at ν = {nu} (finite difference with the Hencky stretch)", abs(fd / (coef * (1 + nu)) - 1) < 1e-4, f"{mp.nstr(fd, 6)} = {mp.nstr(coef, 6)}·(1 + ν)")
print(f"    a uniaxial Hencky strain of 10⁻³ along the space circle moves Paper 22's mₑ/M_P by {mp.nstr(100 * coef * 1e-3, 3)}(1 + ν) per cent; the area channel (1 − ν) does not enter — the bridge is pure modular time.")
inv_coef = -(kexp / 2) * (mp.mpf(1) / 2 - mp.pi * yc * E2 / 6)
print(f"    the covariant version (√y|η|²)^(k/2) = y^(k/4)|η|^k responds as {mp.nstr(inv_coef, 6)}·(1 + ν): the η part plus the zero mode's −k/4.")

print()
print(f"runtime {time.time() - T0:.0f} s")
if FAILS:
    print("FAILED:", *FAILS, sep="\n  ")
    sys.exit(1)
print("ALL CHECKS PASSED")
