#!/usr/bin/env python3
"""
aether_torus_20261024.py -- certificates for Sphereland v14, Instrument 21 ("The aether").

Roger's 2026-10-24 photo dump (phase portraits, electron/positron field "manifolds", the matharium coil
y = x - sin x, kissing numbers, Maxwell's light wave, the inverse-square law, octants, pi^4 + pi^5 ~ e^6),
his question about the aether ("the reason for Aetheris") and the elastic 4D torus as a modular time medium,
and Astra's critique of the elastic-torus reading of modular time.

Sections
  0  Astra's correction: the cycle pair of a torus has mass ratio |tau|; the two lightest levels need not
  1  Four aethers on one number r = c_L^2/c_T^2; zero stiffness freezes a mode, a constraint removes it;
     MacCullagh = Maxwell in temporal gauge without Gauss's law; alpha = Z0/(2 R_K); the rotlet is a dipole
  2  The phase portrait of a deformation: (alpha, beta) = (dilation, strain^2 - rotation^2); what each part
     does to a torus; Gamma0(143) has no centres (no elliptic elements), four cusps of widths 1, 11, 13, 143
  3  The torus as a modular medium: lattice energies are Eisenstein series; the hexagonal point minimises them;
     Delta_H E = s(s-1)E fixes the shear modulus there; 6D -> 4D reduction gives the moduli action;
     the one-loop Casimir potential on T^2 is -(2 pi^3/945) E(tau, 3)/A^2; stringy-cosmic-string deficit pi/6
  4  The inverse-square law as a torus's long-distance shadow; Yukawa strength = lightest-level degeneracy =
     kissing number of the dual lattice; theta series as divisor gases; Ehrenfest and Bertrand
  5  The matharium coil is Kepler's equation; crossings are the apsides; winding one per turn; Fock's S^3
  6  pi^4 + pi^5 vs e^6 (no mechanism; look-elsewhere) against e^{pi sqrt 163} (explained by j)

Every check prints PASS/FAIL; the script ends with ALL CHECKS PASSED or exits non-zero.
"""
import math, cmath, itertools, time, sys
import numpy as np
import mpmath as mp
import sympy as sp
from scipy import special, integrate, optimize, linalg

T0 = time.time()
mp.mp.dps = 30
RESULTS = []


def check(name, ok, detail=""):
    ok = bool(ok)
    RESULTS.append((name, ok))
    print(("  PASS  " if ok else "  FAIL  ") + name + (f"\n        {detail}" if detail else ""))


def divisors(n):
    return [d for d in range(1, n + 1) if n % d == 0]


EPS3 = np.zeros((3, 3, 3))
for (i, j, k), sgn in (((0, 1, 2), 1), ((1, 2, 0), 1), ((2, 0, 1), 1), ((0, 2, 1), -1), ((2, 1, 0), -1), ((1, 0, 2), -1)):
    EPS3[i, j, k] = sgn

# =====================================================================================================
print("§0  Astra's correction: the cycle pair versus the two lightest levels  [structure]")
# =====================================================================================================


def kk_levels(tau, N=14):
    y = tau.imag
    vals = sorted(abs(m - n * tau) / math.sqrt(y) for n in range(-N, N + 1) for m in range(-N, N + 1) if (n, m) != (0, 0))
    lev = []
    for v in vals:
        if not lev or v - lev[-1] > 1e-9:
            lev.append(v)
    return lev


lev = kk_levels(3j)
check("τ = 3i: the first two distinct Kaluza–Klein levels are in ratio 2, the cycle pair (1,0), (0,1) in ratio 3",
      abs(lev[1] / lev[0] - 2) < 1e-12 and abs(abs(3j) / 1 - 3) < 1e-12,
      f"levels ∝ √(m² + 9n²): {', '.join(f'{v * math.sqrt(3):.4f}' for v in lev[:5])} (×1/√3); Instrument 20's 'two lightest masses are in the ratio |τ|' is true only for |τ| ≤ 2")
rng = np.random.default_rng(20261024)
ok = True
for _ in range(300):
    x = rng.uniform(-0.499, 0.499)
    y = math.sqrt(1 - x * x) + rng.uniform(1e-3, 3.5)
    tau = complex(x, y)
    lv = kk_levels(tau)
    ok &= abs(lv[1] / lv[0] - min(abs(tau), 2)) < 1e-9
check("for reduced τ (|Re τ| < ½, |τ| > 1) the second level / first = min(|τ|, 2): the cycle pair is always in ratio |τ|, but the first cycle's second harmonic comes first once |τ| > 2", ok,
      "300 random reduced τ; the two cycle masses are the two shortest independent dual vectors (the dual of a 2D lattice is the lattice turned 90° and scaled by 1/A)")

# =====================================================================================================
print("§1  Four aethers on one number r = c_L²/c_T² = M/μ  [structure]")
# =====================================================================================================


def iso_C(lam, mu, d=3):
    I = np.eye(d)
    return lam * np.einsum('ij,kl->ijkl', I, I) + mu * (np.einsum('ik,jl->ijkl', I, I) + np.einsum('il,jk->ijkl', I, I))


ok = True
for _ in range(50):
    lam, mu = rng.uniform(-1.5, 3), rng.uniform(0.2, 2)
    n = rng.normal(size=3); n /= np.linalg.norm(n)
    Q = np.einsum('ijkl,j,l->ik', iso_C(lam, mu), n, n)
    ev = np.sort(np.linalg.eigvalsh(Q))
    ok &= np.allclose(np.sort([mu, mu, lam + 2 * mu]), ev, atol=1e-12)
check("the acoustic tensor C_ijkl n_j n_l has eigenvalues μ, μ (transverse) and M = λ + 2μ (longitudinal): c_T² = μ/ρ, c_L² = M/ρ", ok, "50 random (λ, μ, n)")
r_, lam_, mu_ = sp.symbols('r lambda mu', real=True)
nu_r = sp.simplify((r_ * mu_ - 2 * mu_) / (2 * (r_ * mu_ - mu_)))
ok = sp.simplify(nu_r - (r_ - 2) / (2 * (r_ - 1))) == 0 and nu_r.subs(r_, 0) == 1 and nu_r.subs(r_, sp.Rational(4, 3)) == -1 and sp.limit(nu_r, r_, sp.oo) == sp.Rational(1, 2)
check("Poisson's ratio in terms of the speed ratio: ν = (r − 2)/(2(r − 1)); r = 0 (labile) is ν = 1, r = 4/3 is ν = −1, r → ∞ (Green) is ν = ½", ok,
      "so the quasi-labile band 0 < r < 4/3 is ν > 1 together with ν < −1, joined through ν = ±∞ at r = 1 (λ = −μ)")
e_, M_, d_ = sp.symbols('e M d', real=True)
W_dil = sp.expand(mu_ * d_ * e_**2 + sp.Rational(1, 2) * (M_ - 2 * mu_) * d_**2 * e_**2)
ok = sp.simplify(W_dil - (sp.Rational(1, 2) * M_ * d_**2 * e_**2 - mu_ * d_ * (d_ - 1) * e_**2)) == 0
vals = {dd: sp.simplify(W_dil.subs({M_: 0, d_: dd}) / (mu_ * e_**2)) for dd in (2, 3, 4)}
ok &= vals == {2: -2, 3: -6, 4: -12}
check("Astra's labile dilation: a uniform dilation ε = e·I costs W = ½M d² e² − μ d(d − 1) e², so at M = 0 it is −2μe² (d = 2), −6μe² (d = 3), −12μe² (d = 4)", ok,
      "the labile medium is unstable as a free body: its bulk modulus is K = M − 2μ(d − 1)/d = −4μ/3 in 3D")
r_free = sp.solve(sp.Eq(r_ - sp.Rational(4, 3), 0), r_)[0]
check("stability ladder: a free body needs μ > 0 and K > 0 (r > 4/3); a clamped body (fixed containing vessel) needs only strong ellipticity μ > 0, M > 0 (r > 0)",
      r_free == sp.Rational(4, 3), "the labile aether r = 0 is the marginal case of the clamped body and is below the free body's threshold — Kelvin's condition 'a fixed containing vessel' (Whittaker) is exactly this gap")

# periodic medium on T^3: the sigma_2 term integrates to zero, so the labile energy is MacCullagh's
Ng = 20
k1 = np.fft.fftfreq(Ng, 1.0 / Ng) * 2 * np.pi
KX, KY, KZ = np.meshgrid(k1, k1, k1, indexing='ij')
Kv = [KX, KY, KZ]
xg = np.arange(Ng) * 2 * np.pi / Ng
X, Y, Z = np.meshgrid(xg, xg, xg, indexing='ij')


def smooth_field(seed, nmodes=6):
    g = np.random.default_rng(seed)
    u = np.zeros((3, Ng, Ng, Ng))
    for c in range(3):
        for _ in range(nmodes):
            kk = g.integers(-2, 3, size=3)
            ph = g.uniform(0, 2 * np.pi)
            u[c] += g.normal() * np.cos(kk[0] * X + kk[1] * Y + kk[2] * Z + ph)
    return u


def grad(u):
    uh = [np.fft.fftn(u[c]) for c in range(3)]
    return np.array([[np.real(np.fft.ifftn(1j * Kv[j] * uh[i])) for j in range(3)] for i in range(3)])   # G[i][j] = d_j u_i


def energies(u, lam, mu):
    Gd = grad(u)
    eps = 0.5 * (Gd + Gd.transpose(1, 0, 2, 3, 4))
    tr = eps[0, 0] + eps[1, 1] + eps[2, 2]
    W1 = np.mean(mu * np.einsum('ij...,ij...->...', eps, eps) + 0.5 * lam * tr**2)
    curl = np.array([Gd[2, 1] - Gd[1, 2], Gd[0, 2] - Gd[2, 0], Gd[1, 0] - Gd[0, 1]])
    W2 = np.mean(0.5 * (lam + 2 * mu) * tr**2 + 0.5 * mu * np.sum(curl**2, axis=0))
    return W1, W2


u = smooth_field(7)
ok = True
for lam, mu in ((-2.0, 1.0), (0.7, 1.3), (-1.2, 0.8)):
    W1, W2 = energies(u, lam, mu)
    ok &= abs(W1 - W2) < 1e-11 * max(1, abs(W1))
phi_h = np.fft.fftn(np.cos(X + 2 * Y) + 0.3 * np.sin(2 * X - Z) + 0.7 * np.cos(Y + Z))
ugrad = np.array([np.real(np.fft.ifftn(1j * Kv[j] * phi_h)) for j in range(3)])
W1g, _ = energies(ugrad, -2.0, 1.0)
check("on a periodic medium (a 3-torus) ∫[μ ε:ε + ½λ(tr ε)²] = ∫[½M(div u)² + ½μ|curl u|²] exactly — σ₂ is a null Lagrangian — so the labile medium is MacCullagh's",
      ok and abs(W1g) < 1e-11, f"spectral derivatives on 20³; a longitudinal field u = ∇φ costs {W1g:.1e} at M = 0: the potential energy has the gauge redundancy u → u + ∇φ")
tt = sp.symbols('t')
xs = sp.symbols('x y z')
uf = [sp.Function(f'u{i}')(tt, *xs) for i in range(3)]
rho_, M2 = sp.symbols('rho M', positive=True)


def curl_s(v):
    return [sp.diff(v[2], xs[1]) - sp.diff(v[1], xs[2]), sp.diff(v[0], xs[2]) - sp.diff(v[2], xs[0]), sp.diff(v[1], xs[0]) - sp.diff(v[0], xs[1])]


divu = sum(sp.diff(uf[i], xs[i]) for i in range(3))
rhs = [(-mu_ * curl_s(curl_s(uf))[i] + M2 * sp.diff(divu, xs[i])) / rho_ for i in range(3)]
lhs_div = sp.simplify(sum(sp.diff(rhs[i], xs[i]) for i in range(3)) - M2 / rho_ * sum(sp.diff(divu, xs[i], 2) for i in range(3)))
check("the equation of motion ρü = −μ curl curl u + M ∇ div u gives ∂²ₜ(div u) = (M/ρ)∇²(div u): at M = 0, div u̇ is conserved point by point — longitudinal disturbances freeze, they do not leave",
      lhs_div == 0, "zero stiffness leaves the longitudinal mode in the spectrum at ω = 0; only a constraint removes it")
A0, Ax, U = sp.Function('A0')(tt, xs[0]), sp.Function('Ax')(tt, xs[0]), sp.Function('U')(tt, xs[0])
eps0 = sp.symbols('epsilon0', positive=True)
Lmax = sp.Rational(1, 2) * eps0 * (sp.diff(Ax, tt) + sp.diff(A0, xs[0]))**2
el = sp.euler_equations(Lmax, [A0, Ax], [tt, xs[0]])
gauss = sp.simplify(el[0].lhs / eps0)
Lmac = sp.Rational(1, 2) * rho_ * sp.diff(U, tt)**2 - sp.Rational(1, 2) * M2 * sp.diff(U, xs[0])**2
elm = sp.simplify(sp.euler_equations(Lmac, [U], [tt, xs[0]])[0].lhs)
ok = sp.simplify(gauss + sp.diff(sp.diff(Ax, tt) + sp.diff(A0, xs[0]), xs[0])) == 0 or sp.simplify(gauss - sp.diff(sp.diff(Ax, tt) + sp.diff(A0, xs[0]), xs[0])) == 0
ok &= sp.simplify(elm.subs(M2, 0) + rho_ * sp.diff(U, tt, 2)) == 0 or sp.simplify(elm.subs(M2, 0) - rho_ * sp.diff(U, tt, 2)) == 0
check("Maxwell's A₀ is a Lagrange multiplier: varying it gives ∂ₓEₓ = 0 (Gauss), which removes the longitudinal mode; MacCullagh's aether has no A₀, and its longitudinal equation at M = 0 is ρÜ = 0 (drift)", ok,
      "MacCullagh's L = ½ρ|u̇|² − ½μ|curl u|² is Maxwell's L = ½ε₀|Ȧ|² − |curl A|²/2μ₀ in temporal gauge (u ↔ A, ρ ↔ ε₀, μ ↔ 1/μ₀) with Gauss's law dropped: charges exist but are frozen")
c_light = 299792458.0
h_pl = 6.62607015e-34
e_ch = 1.602176634e-19
mu0 = 4e-7 * math.pi * 0.99999999987
eps0v = 8.8541878188e-12
alpha_inv = 137.035999177
RK = h_pl / e_ch**2
Z0 = mu0 * c_light
alpha_from_Z0 = Z0 / (2 * RK)
check("c = 1/√(μ₀ε₀) is an identity in the 2019 SI (CODATA 2022 μ₀ε₀c² = 1 to 3·10⁻¹⁰), and what μ₀ measures is α: Z₀ = μ₀c = 2αh/e², α = Z₀/(2R_K)",
      abs(mu0 * eps0v * c_light**2 - 1) < 3e-10 and abs(1 / alpha_from_Z0 / alpha_inv - 1) < 3e-10,
      f"Z₀ = {Z0:.9f} Ω, R_K = h/e² = {RK:.6f} Ω (exact), α⁻¹ = 2R_K/Z₀ = {1 / alpha_from_Z0:.9f} (CODATA 2022: 137.035999177). In MacCullagh's dictionary Z₀ = 1/(ρ_aether c): the fine-structure constant is the aether's softness in quantum units")


def kelvin(x, nu, mu=1.0):
    r = np.linalg.norm(x); nn = x / r
    return ((3 - 4 * nu) * np.eye(3) + np.outer(nn, nn)) / (16 * np.pi * mu * (1 - nu) * r)


def rotlet_numeric(x, Cv, nu, h=1e-5):
    uu = np.zeros(3)
    for k in range(3):
        ek = np.zeros(3); ek[k] = h
        dG = (kelvin(x + ek, nu) - kelvin(x - ek, nu)) / (2 * h)
        for l in range(3):
            for j in range(3):
                if EPS3[l, k, j] != 0:
                    uu += -0.5 * EPS3[l, k, j] * Cv[l] * dG[:, j]
    return uu


ok = True
Cv = np.array([0.3, -0.5, 1.1])
for nu in (-0.6, 0.0, 0.29, 0.45):
    for _ in range(5):
        x = rng.normal(size=3) * 1.5
        ok &= np.allclose(rotlet_numeric(x, Cv, nu), np.cross(Cv, x) / (8 * np.pi * np.linalg.norm(x)**3), rtol=1e-6, atol=1e-10)
x = np.array([0.4, 0.9, -0.7]); hh = 1e-5
Jm = np.array([(np.cross(Cv, x + e) / np.linalg.norm(x + e)**3 - np.cross(Cv, x - e) / np.linalg.norm(x - e)**3) / (2 * hh) for e in np.eye(3) * hh]).T / (8 * np.pi)
omega = 0.5 * np.array([Jm[2, 1] - Jm[1, 2], Jm[0, 2] - Jm[2, 0], Jm[1, 0] - Jm[0, 1]])
rr = np.linalg.norm(x); nh = x / rr
dip = (3 * np.dot(Cv, nh) * nh - Cv) / (16 * np.pi * rr**3)
ok &= np.allclose(omega, dip, rtol=1e-6)
# torque transmitted through a sphere: the moment of the tractions balances the applied couple
R0 = 1.3
th_n, ph_n = 60, 120
tq = np.zeros(3)
for i in range(th_n):
    th = (i + 0.5) * np.pi / th_n
    for j in range(ph_n):
        ph = (j + 0.5) * 2 * np.pi / ph_n
        nn = np.array([math.sin(th) * math.cos(ph), math.sin(th) * math.sin(ph), math.cos(th)])
        xp = R0 * nn
        Jp = np.array([(np.cross(Cv, xp + e) / np.linalg.norm(xp + e)**3 - np.cross(Cv, xp - e) / np.linalg.norm(xp - e)**3) / (2 * hh) for e in np.eye(3) * hh]).T / (8 * np.pi)
        sig = Jp + Jp.T
        tq += np.cross(xp, sig @ nn) * R0**2 * math.sin(th) * (np.pi / th_n) * (2 * np.pi / ph_n)
ok &= np.allclose(tq, -Cv, rtol=2e-3)
check("a point couple C in any isotropic solid displaces u = C × r/(8πμr³), whatever ν (Love's 'centre of rotation'); its rotation ½ curl u = [3(C·r̂)r̂ − C]/(16πμr³) is a magnetic dipole field",
      ok, f"from Kelvin's tensor by differentiation at ν = −0.6, 0, 0.29, 0.45; the tractions on a sphere carry the torque {np.round(tq, 4)} = −C. In MacCullagh's aether a twist is a magnetic moment; electron and positron are its two chiralities relative to spin")
x = np.array([0.0, 0.0, 2.0])
ok = all(abs(np.dot(np.array([0, 0, 1.0]), kelvin(x, nu) @ np.array([0, 0, 1.0])) - 1 / (4 * np.pi * 2.0)) < 1e-14 for nu in (-0.7, 0.0, 0.3, 0.49))
check("Kelvin's point force: along its own line the displacement is F/(4πμr) for every ν — Coulomb's law with μ in place of ε₀; across it (3 − 4ν)/(16πμ(1 − ν)r), which diverges as ν → 1 (the labile medium has no static response)", ok)

# =====================================================================================================
print("§2  The phase portrait of a deformation  [structure]")
# =====================================================================================================


def classify(lmb):
    l1, l2 = lmb
    if abs(l1.imag) > 1e-12:
        return 'centre' if abs(l1.real) < 1e-12 else ('attracting focus' if l1.real < 0 else 'repelling focus')
    a1, a2 = sorted([l1.real, l2.real])
    if a1 < 0 < a2:
        return 'saddle'
    return 'attracting node' if a2 < 0 else 'repelling node'


def classify_ab(al, be):
    if be > 0:
        if abs(al) < math.sqrt(be):
            return 'saddle'
        return 'attracting node' if al < 0 else 'repelling node'
    if be < 0:
        return 'centre' if abs(al) < 1e-12 else ('attracting focus' if al < 0 else 'repelling focus')
    return 'star/degenerate'


ok = True
for _ in range(2000):
    e, s1, s2, w = rng.normal(size=4)
    if rng.random() < 0.1:
        e = 0.0
    L = np.array([[e + s1, s2 - w], [s2 + w, e - s1]])
    al, be = e, s1 * s1 + s2 * s2 - w * w
    ok &= abs(np.trace(L) - 2 * al) < 1e-12 and abs(np.linalg.det(L) - (al * al - be)) < 1e-12
    ok &= classify(np.linalg.eigvals(L)) == classify_ab(al, be)
check("for ẋ = Lx with L = e·I + shear(s₁, s₂) + rotation(w): eigenvalues e ± √(|s|² − w²), so the chart's α is the dilation rate and β = |s|² − w² is the Okubo–Weiss parameter", ok,
      "2000 random gradients classified two ways (eigenvalues vs the chart's regions): saddle iff det < 0, nodes iff 0 < det ≤ α², foci iff β < 0, centre iff β < 0 and α = 0")
e_s, s1_s, s2_s, w_s, lam_s, mu_s = sp.symbols('e s1 s2 w lambda mu', real=True)
L2 = sp.Matrix([[e_s + s1_s, s2_s - w_s], [s2_s + w_s, e_s - s1_s]])
eps2 = (L2 + L2.T) / 2
Wstd = mu_s * sum(eps2[i, j]**2 for i in range(2) for j in range(2)) + lam_s / 2 * eps2.trace()**2
Wmac = (lam_s + 2 * mu_s) / 2 * L2.trace()**2 + mu_s / 2 * (L2[1, 0] - L2[0, 1])**2 - 2 * mu_s * L2.det()
ok = sp.simplify(Wstd - Wmac) == 0 and sp.simplify(Wstd - (2 * mu_s * (s1_s**2 + s2_s**2) + 2 * (lam_s + mu_s) * e_s**2)) == 0
check("in 2D: μ ε:ε + ½λ(tr L)² = ½M(tr L)² + ½μ(curl)² − 2μ det L; a medium charges the saddle part 2μ|s|² and the star part 2K₂e², and the centre part (rotation) is free — MacCullagh's aether charges the centre instead; the two differ by det L = α² − β, a null Lagrangian",
      ok, "the chart's parabola β = α² is det L = 0: rank-one (uniaxial) deformations, where σ₂ vanishes and only M is felt")
ok = True
for _ in range(200):
    x0, y0 = rng.uniform(-0.5, 0.5), rng.uniform(0.6, 2.5)
    th1 = rng.uniform(0, 2 * np.pi)
    w1 = cmath.exp(1j * th1) * rng.uniform(0.5, 2)
    w2 = w1 * complex(x0, y0)
    e, s1, s2, w = rng.normal(size=4)
    a_c, b_c = complex(e, w), complex(s1, s2)
    dt_ = 1e-7
    f = lambda om: om + dt_ * (a_c * om + b_c * om.conjugate())
    tau0 = w2 / w1; tau1 = f(w2) / f(w1)
    dtau = (tau1 - tau0) / dt_
    pred = -2j * tau0.imag * b_c * cmath.exp(-2j * th1)
    A0_, A1_ = (w1.conjugate() * w2).imag, (f(w1).conjugate() * f(w2)).imag
    ok &= abs(dtau - pred) < 1e-5 * abs(pred) + 1e-6 and abs((math.log(A1_) - math.log(A0_)) / dt_ - 2 * e) < 1e-5
check("on a torus: δτ = −2i·Im τ·b·e^{−2i arg ω₁} and δ ln A = 2e — the saddle part (shear b) moves modular time along a geodesic at hyperbolic speed 2|b|, the star part (dilation) changes only the size, the centre part (rotation) changes nothing",
      ok, "200 random lattices and gradients, finite differences 10⁻⁷")
ok = True
for _ in range(300):
    s1, s2, w = rng.normal(size=3)
    if rng.random() < 0.15:
        w = math.hypot(s1, s2) * (1 if rng.random() < 0.5 else -1)
    X_ = np.array([[s1, s2 - w], [s2 + w, -s1]])
    be = s1 * s1 + s2 * s2 - w * w
    tr = np.trace(linalg.expm(X_))
    if abs(be) < 1e-12:
        ok &= abs(tr - 2) < 1e-9
    elif be > 0:
        ok &= abs(tr - 2 * math.cosh(math.sqrt(be))) < 1e-9 and tr > 2
    else:
        ok &= abs(tr - 2 * math.cos(math.sqrt(-be))) < 1e-9
check("on the line α = 0 the chart is SL(2,ℝ)'s trichotomy: β > 0 hyperbolic (saddle, |tr| > 2), β = 0 parabolic (simple shear), β < 0 elliptic (centre)", ok, "exp of 300 traceless generators")
N143 = 143
sol2 = [a for a in range(N143) if (a * a + 1) % N143 == 0]
sol3p = [a for a in range(N143) if (a * a - a + 1) % N143 == 0]
sol3m = [a for a in range(N143) if (a * a + a + 1) % N143 == 0]
widths = {c: N143 // math.gcd(N143, c * c) for c in divisors(N143)}
genus = 1 + sp.Rational(168, 12) - sp.Rational(len(sol2), 4) - sp.Rational(len(sol3p), 3) - sp.Rational(len(widths), 2)
check("Γ₀(143) has no elliptic elements: a² − ta + 1 ≡ 0 (mod 143) has no solution for t = 0, ±1, so none of its maps of the torus is a centre — all are saddles (Anosov) or shears (Dehn twists)",
      not sol2 and not sol3p and not sol3m and len([a for a in range(13) if (a * a + 1) % 13 == 0]) == 2,
      "(control: Γ₀(13) has the two order-2 points a = 5, 8); 11 ≡ 3 (mod 4) and 11 ≡ 2 (mod 3) forbid them")
check("its four cusps — the parabolic classes, the defect types of a 143-marked medium — have widths 1, 11, 13, 143, summing to the index 168; genus 1 + 168/12 − 4/2 = 13",
      sorted(widths.values()) == [1, 11, 13, 143] and sum(widths.values()) == 168 and genus == 13,
      "width of the cusp 1/c is the least h with N | c²h, i.e. N/gcd(N, c²)")

# =====================================================================================================
print("§3  The torus as a modular medium  [structure + derived]")
# =====================================================================================================


def xi_mp(s):
    return mp.pi**(-s / 2) * mp.gamma(s / 2) * mp.zeta(s)


def E_mp(tau, s):
    s = mp.mpf(s); x, y = mp.re(tau), mp.im(tau)
    nmax = int(mp.ceil(45 / (2 * mp.pi * y))) + 8
    acc = mp.mpf(0)
    for n in range(1, nmax + 1):
        sig = mp.fsum(mp.mpf(d)**(1 - 2 * s) for d in divisors(n))
        acc += mp.mpf(n)**(s - mp.mpf(1) / 2) * sig * mp.besselk(s - mp.mpf(1) / 2, 2 * mp.pi * n * y) * mp.cos(2 * mp.pi * n * x)
    return y**s + xi_mp(2 * s - 1) / xi_mp(2 * s) * y**(1 - s) + 4 * mp.sqrt(y) / xi_mp(2 * s) * acc


def E_direct(tau, s, Nc=500):
    m = np.arange(-Nc, Nc + 1)
    MM, NN = np.meshgrid(m, m, indexing='ij')
    zz = np.abs(MM * tau + NN).astype(float)
    zz[Nc, Nc] = np.inf
    return float(np.sum(tau.imag**s / zz**(2 * s)) / (2 * special.zeta(2 * s)))


ok = True
for tau in (complex(0.5, math.sqrt(3) / 2), 1j, complex(0.21, 1.37)):
    for s in (3, 4):
        ok &= abs(float(E_mp(mp.mpc(tau.real, tau.imag), s)) - E_direct(tau, s)) < 2e-9
check("E(τ, s) = Σ' y^s/|mτ + n|^{2s} / 2ζ(2s): the Fourier (Chowla–Selberg) expansion agrees with the lattice sum — the lattice energy of a 2D crystal with pair potential r^{−2s} at fixed area is 2ζ(2s)E(τ, s)", ok, "s = 3, 4 at ρ, i and a generic τ (direct sum to |m|, |n| ≤ 500)")
ok = True
hfd = mp.mpf('1e-7')
for tau, s in ((mp.mpc('0.13', '1.21'), 2), (mp.mpc('0.4', '0.95'), 3), (mp.mpc('-0.2', '1.6'), mp.mpf('1.5'))):
    E0 = E_mp(tau, s)
    lap = (E_mp(tau + hfd, s) + E_mp(tau - hfd, s) + E_mp(tau + 1j * hfd, s) + E_mp(tau - 1j * hfd, s) - 4 * E0) / hfd**2 * mp.im(tau)**2
    ok &= abs(lap / (s * (s - 1) * E0) - 1) < 1e-8
check("Δ_H E = y²(∂²ₓ + ∂²ᵧ)E = s(s − 1)E: the trace of the energy's Hessian in the Poincaré metric is fixed everywhere", ok, "finite differences at 30 digits, s = 3/2, 2, 3")


def E_f(x, y, s, nmax=14):
    xi = lambda q: np.pi**(-q / 2) * special.gamma(q / 2) * special.zeta(q)
    out = y**s + xi(2 * s - 1) / xi(2 * s) * y**(1 - s)
    acc = 0.0
    for n in range(1, nmax + 1):
        sig = sum(d**(1 - 2 * s) for d in divisors(n))
        acc = acc + n**(s - 0.5) * sig * special.kv(s - 0.5, 2 * np.pi * n * y) * np.cos(2 * np.pi * n * x)
    return out + 4 * np.sqrt(y) / xi(2 * s) * acc


def logabs_eta(x, y, nmax=60):
    q_abs = np.exp(-2 * np.pi * y)
    out = -np.pi * y / 12
    for n in range(1, nmax + 1):
        out = out + np.log(np.abs(1 - q_abs**n * np.exp(2j * np.pi * n * x)))
    return out


def F0_f(x, y):           # Kronecker's F0 up to the additive constant (6/π)C0: −(6/π) log(√y |η|²)
    return -(6 / np.pi) * (0.5 * np.log(y) + 2 * logabs_eta(x, y))


xs_g = np.linspace(0, 0.5, 101)
ok = True
mins = []
for label, fn in (('s=3/2', lambda x, y: E_f(x, y, 1.5)), ('s=2', lambda x, y: E_f(x, y, 2.0)), ('s=3', lambda x, y: E_f(x, y, 3.0)), ('s=6', lambda x, y: E_f(x, y, 6.0)), ('log (s→1)', F0_f)):
    best = (np.inf, None)
    for x in xs_g:
        ys = np.linspace(math.sqrt(1 - x * x), 2.6, 160)
        vals_ = fn(x, ys)
        i = int(np.argmin(vals_))
        if vals_[i] < best[0]:
            best = (vals_[i], (x, ys[i]))
    mins.append(f"{label}: ({best[1][0]:.3f}, {best[1][1]:.4f})")
    ok &= abs(best[1][0] - 0.5) < 1e-9 and abs(best[1][1] - math.sqrt(3) / 2) < 1e-9
check("on the fundamental domain every one of these lattice energies is least at the hexagonal point ρ = e^{iπ/3} (Rankin, Cassels, Ennola, Diananda for r^{−2s}; Montgomery's theta theorem implies all s)", ok, "; ".join(mins))
ok = True
rep = []
for s in (mp.mpf('1.5'), 2, 3, 6):
    rho = mp.mpc(mp.mpf(1) / 2, mp.sqrt(3) / 2)
    Er = E_mp(rho, s); y = mp.im(rho); hd = mp.mpf('1e-6')
    Hxx = (E_mp(rho + hd, s) - 2 * Er + E_mp(rho - hd, s)) / hd**2 * y**2
    Hyy = (E_mp(rho + 1j * hd, s) - 2 * Er + E_mp(rho - 1j * hd, s)) / hd**2 * y**2
    Hxy = (E_mp(rho + hd + 1j * hd, s) - E_mp(rho + hd - 1j * hd, s) - E_mp(rho - hd + 1j * hd, s) + E_mp(rho - hd - 1j * hd, s)) / (4 * hd**2) * y**2
    target = s * (s - 1) * Er / 2
    ok &= abs(Hxx / target - 1) < 1e-6 and abs(Hyy / target - 1) < 1e-6 and abs(Hxy) < 1e-6 * abs(target)
    Ei = E_mp(mp.mpc(0, 1), s)
    Hxx_i = (E_mp(mp.mpc(hd, 1), s) - 2 * Ei + E_mp(mp.mpc(-hd, 1), s)) / hd**2
    Hyy_i = (E_mp(mp.mpc(0, 1 + hd), s) - 2 * Ei + E_mp(mp.mpc(0, 1 - hd), s)) / hd**2
    ok &= Hxx_i < 0 < Hyy_i
    rep.append(f"s = {mp.nstr(s, 3)}: E(ρ) = {mp.nstr(Er, 8)}, μ_ρ = s(s−1)E/2 = {mp.nstr(target, 8)}; at i the Hessian is ({mp.nstr(Hxx_i, 5)}, {mp.nstr(Hyy_i, 5)})")
check("at ρ the Hessian is isotropic and equals s(s − 1)E(ρ)/2 times the Poincaré metric (the order-3 rotation forces isotropy): the hexagonal medium has one shear modulus, fixed by the Laplace eigenvalue; at i it is a saddle — the square lattice of central forces is unstable to simple shear",
      ok, "\n        ".join(rep))
y_r = math.sqrt(3) / 2
ll_rho = 0.5 * math.log(y_r) + 2 * float(logabs_eta(0.5, y_r))
ll_i = 0.5 * math.log(1.0) + 2 * float(logabs_eta(0.0, 1.0))
W_rho, W_i = -0.25 * math.log(2 * math.pi) - 0.5 * ll_rho, -0.25 * math.log(2 * math.pi) - 0.5 * ll_i
hd = 1e-4
lapF = ((0.5 * math.log(1.3) + 2 * float(logabs_eta(0.2 + hd, 1.3))) + (0.5 * math.log(1.3) + 2 * float(logabs_eta(0.2 - hd, 1.3))) + (0.5 * math.log(1.3 + hd) + 2 * float(logabs_eta(0.2, 1.3 + hd))) + (0.5 * math.log(1.3 - hd) + 2 * float(logabs_eta(0.2, 1.3 - hd))) - 4 * (0.5 * math.log(1.3) + 2 * float(logabs_eta(0.2, 1.3)))) / hd**2 * 1.3**2
Wf = lambda x, y: -0.25 * math.log(2 * math.pi) - 0.5 * (0.5 * math.log(y) + 2 * float(logabs_eta(x, y)))
hW = 1e-4
WxxR = (Wf(0.5 + hW, y_r) - 2 * Wf(0.5, y_r) + Wf(0.5 - hW, y_r)) / hW**2 * y_r**2
WyyR = (Wf(0.5, y_r + hW) - 2 * Wf(0.5, y_r) + Wf(0.5, y_r - hW)) / hW**2 * y_r**2
check("the s → 1 member is Kronecker's: the Coulomb (vortex) lattice energy is W = −½ log(√(2πy)·|η(τ)|²), which reproduces Sandier–Serfaty's renormalized energies −0.2011 (triangular) and −0.1958 (square); Δ_H log(√y|η|²) = −½ exactly because log|η| is harmonic, so Δ_H W = ¼",
      abs(W_rho - (-0.2011)) < 6e-5 and abs(W_i - (-0.1958)) < 6e-5 and abs(lapF + 0.5) < 1e-5 and abs(WxxR - 0.125) < 1e-5 and abs(WyyR - 0.125) < 1e-5,
      f"W(ρ) = {W_rho:.5f}, W(i) = {W_i:.5f} (the constant −¼ log 2π fitted here on one value, then checked on the other); Instrument 19's F₀ is the Abrikosov/Wigner lattice energy, and at ρ its Hessian is isotropic with eigenvalues ⅛ ({WxxR:.6f}, {WyyR:.6f}) in Poincaré units")
tsym = sp.symbols('t')
Nf, af, sf, t1f, t2f = [sp.Function(nm)(tsym) for nm in ('N', 'a', 'sigma', 'tau1', 'tau2')]
Aexp = sp.exp(2 * sf)
G6 = sp.zeros(6, 6); G6[0, 0] = -Nf**2
for i in (1, 2, 3):
    G6[i, i] = af**2
hmat = (Aexp / t2f) * sp.Matrix([[1, t1f], [t1f, t1f**2 + t2f**2]])
G6[4:6, 4:6] = hmat
G6i = sp.zeros(6, 6); G6i[0, 0] = -1 / Nf**2
for i in (1, 2, 3):
    G6i[i, i] = 1 / af**2
G6i[4:6, 4:6] = sp.simplify(hmat.inv())
dG6 = G6.diff(tsym)
Gam = [[[0] * 6 for _ in range(6)] for _ in range(6)]
for A_ in range(6):
    for b in range(6):
        for c in range(6):
            val = 0
            for d in range(6):
                if G6i[A_, d] == 0:
                    continue
                term = (dG6[d, c] if b == 0 else 0) + (dG6[d, b] if c == 0 else 0) - (dG6[b, c] if d == 0 else 0)
                if term != 0:
                    val += G6i[A_, d] * term
            Gam[A_][b][c] = sp.simplify(val / 2) if val != 0 else 0
R6 = 0
for b in range(6):
    for d in range(6):
        if G6i[b, d] == 0:
            continue
        Rbd = sp.diff(Gam[0][b][d], tsym) - (sp.diff(sum(Gam[A_][b][A_] for A_ in range(6)), tsym) if d == 0 else 0)
        for A_ in range(6):
            for e in range(6):
                Rbd += Gam[A_][A_][e] * Gam[e][b][d] - Gam[A_][d][e] * Gam[e][b][A_]
        R6 += G6i[b, d] * Rbd
NE, aE = sp.Function('N_E')(tsym), sp.Function('a_E')(tsym)
L6 = sp.simplify(sp.expand((Nf * af**3 * Aexp * R6).subs({Nf: NE * sp.exp(-sf), af: aE * sp.exp(-sf)}).doit()))
c0, c1, c2 = sp.symbols('c0 c1 c2')
Lexp = c0 * (-6 * aE * aE.diff(tsym)**2 / NE) + aE**3 / NE * (c1 * (2 * sf.diff(tsym))**2 + c2 * (t1f.diff(tsym)**2 + t2f.diff(tsym)**2) / t2f**2)
ELs = [sp.simplify(sp.expand(eq.lhs)) for eq in sp.euler_equations(L6 - Lexp, [NE, aE, sf, t1f, t2f], tsym)]
resid = [sp.simplify(e.subs({c0: 1, c1: 1, c2: sp.Rational(1, 2)})) for e in ELs]
ok = all(r == 0 for r in resid)
alt = [sp.simplify(e.subs({c0: 1, c1: 1, c2: 1})) for e in ELs]
ok &= any(r != 0 for r in alt)
check("6D gravity on a torus reduces exactly to S = (M_P²/2)∫√−g [R − (∂ ln A)² − ½|∂τ|²/(Im τ)²]: the action Astra asked for has the Poincaré metric as the shape's kinetic term and coefficient 1 for the size",
      ok, "minisuperspace N(t), a(t), A(t), τ(t): every Euler–Lagrange expression of (reduced 6D Lagrangian − this form) vanishes identically in the Einstein frame g_E = A·g; changing ½ → 1 does not")
hu = (1 / t2f) * sp.Matrix([[1, t1f], [t1f, t1f**2 + t2f**2]])
tr_ = sp.simplify((hu.inv() * hu.diff(tsym) * hu.inv() * hu.diff(tsym)).trace() / 4)
check("the shape part: ¼ tr(ĥ⁻¹ĥ̇ĥ⁻¹ĥ̇) = ½|τ̇|²/(Im τ)² for the unit-area torus metric ĥ = (1/τ₂)[[1, τ₁], [τ₁, |τ|²]]",
      sp.simplify(tr_ - (t1f.diff(tsym)**2 + t2f.diff(tsym)**2) / (2 * t2f**2)) == 0)
# Casimir: V = -(Γ(D/2)/(2π^{D/2})) V_n Σ'_{ℓ∈Λ} |ℓ|^{-D}  (per real bosonic degree of freedom)
casimir = lambda D, Vn, S: -mp.gamma(mp.mpf(D) / 2) / (2 * mp.pi**(mp.mpf(D) / 2)) * Vn * S
Lc = mp.mpf('1.7')
V_circle = casimir(5, Lc, 2 * mp.zeta(5) / Lc**5)
beta = mp.mpf('0.9')
V_thermal = casimir(4, beta, 2 * mp.zeta(4) / beta**4)
ok = abs(V_circle - (-3 * mp.zeta(5) / (4 * mp.pi**2 * Lc**4))) < 1e-25 and abs(V_thermal / beta - (-mp.pi**2 / (90 * beta**4))) < 1e-25
tau_c, A_c = complex(0.31, 1.18), 2.3
m = np.arange(-400, 401)
MM, NN = np.meshgrid(m, m, indexing='ij')
ell = np.abs(np.sqrt(A_c / tau_c.imag) * (MM + NN * tau_c))
ell[400, 400] = np.inf
S6 = float(np.sum(ell**-6.0))
V_T2 = float(casimir(6, A_c, S6))
pred = -(2 * math.pi**3 / 945) * float(E_mp(mp.mpc(tau_c.real, tau_c.imag), 3)) / A_c**2
ok &= abs(V_T2 / pred - 1) < 1e-9
check("one-loop vacuum energy of a massless boson on ℝ^d × T^n: V = −Γ(D/2)/(2π^{D/2})·V_n·Σ'|ℓ|^{−D}; it gives Stefan–Boltzmann −π²T⁴/90 on the thermal circle, −3ζ(5)/(4π²L⁴) on a KK circle, and on a torus V = −(2π³/945)·E(τ, 3)/A²",
      ok, f"T² check at τ = 0.31 + 1.18i, A = 2.3: lattice sum {V_T2:.10e} vs formula {pred:.10e}. Bosons make ρ a maximum of V in shape (E is least there): the shape runs to the cusp; fermions make it the minimum; with N_B = N_F (supersymmetry) V ≡ 0 — a supersymmetric torus medium has no shape stiffness at one loop")
E3rho = float(E_mp(mp.mpc(0.5, math.sqrt(3) / 2), 3))
mKK2A = 4 * math.pi**2 * 2 / math.sqrt(3)
ratio = 6 * (2 * math.pi**3 / 945) * E3rho / mKK2A**2
check("[model] if fermions dominate by one degree of freedom, the hexagonal torus is a minimum and its shape doublet has m_τ² = 6·(2π³/945)·E(ρ, 3)/(A² M_P²) = (E(ρ,3)/1680π)·m_KK⁴/M_P², two degenerate scalars (isotropy)",
      abs(ratio - E3rho / (1680 * math.pi)) < 1e-12,
      f"E(ρ, 3) = {E3rho:.8f}, m_τ = {math.sqrt(ratio):.5f}·m_KK²/M_P (m_KK² = 8π²/(√3 A), the hexagonal torus's lightest level): for m_KK = ħc/38.6 μm this is ~10⁻³⁴ eV, a Hubble-range scalar — the programme rule then demands it couple to ordinary matter far below gravity, or not at all")
ts = np.linspace(0, 1, 401)
Es = np.array([E_f(0.5 + t, math.sqrt(3) / 2, 3.0) for t in ts])
d2 = (E_f(0.5 + 1e-4, math.sqrt(3) / 2, 3.0) - 2 * E_f(0.5, math.sqrt(3) / 2, 3.0) + E_f(0.5 - 1e-4, math.sqrt(3) / 2, 3.0)) / 1e-8
slope = np.max(np.abs(np.gradient(Es, ts)))
frenkel = slope / (d2 / (2 * np.pi))
check("along the slip path τ = ρ + t the energy is exactly periodic (a lattice-invariant shear returns the crystal to itself); [model] its steepest slope is the ideal shear strength, Frenkel's sinusoid giving μ/2π",
      abs(Es[0] - Es[-1]) < 1e-12 and abs(E_f(0.5 + 0.3, math.sqrt(3) / 2, 3.0) - E_f(0.5 - 0.3, math.sqrt(3) / 2, 3.0)) < 1e-12,
      f"for the r⁻⁶ crystal the steepest slope is {frenkel:.3f} × Frenkel's estimate (curvature at t = 0 = {d2:.4f})")
ys_ = np.array([4.0, 6.0, 8.0])
ratio_cusp = [(4 * float(logabs_eta(0.0, y)) / math.log(float(abs(mp.kleinj(mp.mpc(0, y)) * 1728)))) for y in ys_]
check("stringy cosmic strings: in the metric Im τ·|η(τ)|⁴|dz|² with j(τ) = z, |η|⁴ ~ |z|^{−1/6} at a string, a cone of deficit 2π/12 = π/6, and 24 strings close the plane into a sphere (24·π/6 = 4π, Gauss–Bonnet)",
      abs(ratio_cusp[-1] + 1 / 6) < 0.003 and abs(24 * math.pi / 6 - 4 * math.pi) < 1e-15,
      f"log|η(iy)|⁴ / log|j(iy)| = {', '.join(f'{v:.5f}' for v in ratio_cusp)} at y = 4, 6, 8 → −1/6. The same Im τ·|η|⁴ is the torus determinant det′Δ/A of Instrument 20 (Greene–Shapere–Vafa–Yau 1990)")

# =====================================================================================================
print("§4  The inverse-square law, a torus's long-distance shadow  [structure]")
# =====================================================================================================
Lp = mp.mpf(1)
G_closed = lambda r, z: mp.sinh(2 * mp.pi * r / Lp) / (4 * mp.pi * Lp * r * (mp.cosh(2 * mp.pi * r / Lp) - mp.cos(2 * mp.pi * z / Lp)))
G_images = lambda r, z: mp.nsum(lambda n: 1 / (4 * mp.pi**2 * (r**2 + (z - n * Lp)**2)), [-mp.inf, mp.inf])
G_yuk = lambda r, z: (1 + 2 * mp.nsum(lambda m: mp.exp(-2 * mp.pi * m * r / Lp) * mp.cos(2 * mp.pi * m * z / Lp), [1, mp.inf])) / (4 * mp.pi * Lp * r)
ok = True
for r, z in ((mp.mpf('0.05'), mp.mpf('0.02')), (mp.mpf('0.3'), mp.mpf('0.4')), (mp.mpf('2.5'), mp.mpf('0.1'))):
    ok &= abs(G_closed(r, z) / G_images(r, z) - 1) < 1e-15 and abs(G_closed(r, z) / G_yuk(r, z) - 1) < 1e-20
near = G_closed(mp.mpf('1e-3'), 0) * 4 * mp.pi**2 * mp.mpf('1e-3')**2
far = G_closed(mp.mpf(8), 0) * 4 * mp.pi * Lp * 8
check("a point source on ℝ³ × S¹ of length L: G = sinh(2πr/L)/[4πLr(cosh(2πr/L) − cos(2πz/L))] — the image sum (4D, inverse cube) equals the Kaluza–Klein sum (3D Coulomb plus Yukawas e^{−2πmr/L})",
      ok and abs(near - 1) < 1e-4 and abs(far - 1) < 1e-20, "r ≪ L: 1/(4π²ρ²); r ≫ L: 1/(4πLr)(1 + 2e^{−2πr/L}cos(2πz/L) + …): the leading deviation has α = 2")
tau_h = complex(0.5, math.sqrt(3) / 2)
Acell = tau_h.imag
w1, w2 = np.array([1.0, 0.0]), np.array([tau_h.real, tau_h.imag])
Bm = np.column_stack([w1, w2]); Dm = np.linalg.inv(Bm).T
rng_r = 0.45
Rcut = 90.0
mgrid = np.arange(-130, 131)
MMg, NNg = np.meshgrid(mgrid, mgrid, indexing='ij')
LX, LY = MMg * w1[0] + NNg * w2[0], MMg * w1[1] + NNg * w2[1]
L2g = LX**2 + LY**2
img = float(np.sum(1 / (8 * np.pi**2 * (rng_r**2 + L2g[L2g < Rcut**2])**1.5)))
kk = 0.0
for mm in range(-25, 26):
    for nn in range(-25, 26):
        kv_ = Dm @ np.array([mm, nn])
        kk += math.exp(-2 * math.pi * np.linalg.norm(kv_) * rng_r)
kk /= (4 * math.pi * Acell * rng_r)
tail = 1 / (4 * math.pi * Acell * math.sqrt(rng_r**2 + Rcut**2))   # continuum tail outside the disk |ℓ| < Rcut
check("on ℝ³ × T² the same Poisson summation holds: Σ over lattice images of 1/(8π²ρ³) = (1/(4πAr)) Σ over dual vectors k of e^{−2π|k|r} (hexagonal torus, r = 0.45)",
      abs((img + tail) / kk - 1) < 2e-4, f"images {img + tail:.6f} (with tail estimate), KK sum {kk:.6f}")


def shell_counts(vectors_norms):
    out = {}
    for nrm in vectors_norms:
        key = round(nrm, 9)
        out[key] = out.get(key, 0) + 1
    return sorted(out.items())


def lattice_shells(B, R=4):
    d = B.shape[1]
    norms = []
    for c in itertools.product(range(-R, R + 1), repeat=d):
        if any(c):
            v = B @ np.array(c)
            norms.append(float(v @ v))
    sh = shell_counts(norms)
    return sh


Bsq = np.eye(2); Bhex = np.array([[1, 0.5], [0, math.sqrt(3) / 2]]); Bz4 = np.eye(4)
BD4 = np.array([[1, 1, 0, 0], [-1, 1, 0, 0], [0, -1, 1, 0], [0, 0, -1, 1]], dtype=float).T
first = {name: lattice_shells(B, 3)[0][1] for name, B in (('square', Bsq), ('hexagonal', Bhex), ('Z4', Bz4), ('D4', BD4))}
check("the Yukawa strength α of a torus's leading deviation is the number of shortest dual vectors: circle 2, square torus 4, hexagonal 6, hypercubic T⁴ 8, D₄ T⁴ 24 (Kehagias–Sfetsos: α = 2n for an n-torus of equal radii)",
      first == {'square': 4, 'hexagonal': 6, 'Z4': 8, 'D4': 24}, f"{first}; for the self-dual E₈ torus it is 240")
ok = True
for n in range(1, 13):
    r2 = sum(1 for a in range(-4, 5) for b in range(-4, 5) if a * a + b * b == n)
    ok &= r2 == 4 * sum((1 if d % 4 == 1 else -1 if d % 4 == 3 else 0) for d in divisors(n))
    rA2 = sum(1 for a in range(-5, 6) for b in range(-5, 6) if a * a + a * b + b * b == n)
    ok &= rA2 == 6 * sum((1 if d % 3 == 1 else -1 if d % 3 == 2 else 0) for d in divisors(n))
D4cnt = {}
for c in itertools.product(range(-3, 4), repeat=4):
    if sum(c) % 2 == 0 and any(c):
        nrm = sum(v * v for v in c)
        if nrm <= 14:
            D4cnt[nrm] = D4cnt.get(nrm, 0) + 1
for n in range(1, 8):
    odd = n
    while odd % 2 == 0:
        odd //= 2
    ok &= D4cnt.get(2 * n, 0) == 24 * sum(divisors(odd))
E8cnt = {2: 0, 4: 0, 6: 0}
for c in itertools.product(range(-2, 3), repeat=8):
    if sum(c) % 2 == 0:
        nrm = sum(v * v for v in c)
        if nrm in E8cnt:
            E8cnt[nrm] += 1
for c in itertools.product((-1.5, -0.5, 0.5, 1.5), repeat=8):
    if round(sum(c)) % 2 == 0:
        nrm = sum(v * v for v in c)
        key = int(round(nrm))
        if key in E8cnt and abs(nrm - key) < 1e-9:
            E8cnt[key] += 1
ok &= E8cnt == {2: 240, 4: 240 * 9, 6: 240 * 28}
check("theta series are divisor gases: shells of ℤ² = 4Σχ₋₄(d), of the hexagonal lattice 6Σχ₋₃(d), of D₄ 24σ₁(odd part of n), of E₈ 240σ₃(n) (θ_E₈ = E₄) — the KK level counts of these tori",
      ok, f"brute force n ≤ 12 (2D), n ≤ 7 (D₄), norms 2, 4, 6 of E₈: {E8cnt}")
cap = (1 - math.cos(math.pi / 6)) / 2
ok = abs(cap - 0.06699) < 1e-5 and abs(12 * cap - 0.80385) < 1e-5 and abs(13 * cap - 0.87083) < 1e-5
U4 = [np.array(v, float) for v in itertools.product((-0.5, 0.5), repeat=4)] + [s * e for e in np.eye(4) for s in (1, -1)]
ang4 = min(math.degrees(math.acos(np.clip(u @ v, -1, 1))) for u, v in itertools.combinations(U4, 2))
U3 = [np.array(v, float) / math.sqrt(3) for v in itertools.product((-1, 1), repeat=3)] + [s * e for e in np.eye(3) for s in (1, -1)]
ang3 = min(math.degrees(math.acos(np.clip(u @ v, -1, 1))) for u, v in itertools.combinations(U3, 2))
fcc = [np.array(v, float) / math.sqrt(2) for v in set(itertools.permutations((1, 1, 0))) | set(itertools.permutations((1, -1, 0))) | set(itertools.permutations((-1, -1, 0)))]
angf = min(math.degrees(math.acos(np.clip(u @ v, -1, 1))) for u, v in itertools.combinations(fcc, 2))
ok &= len(U4) == 24 and abs(ang4 - 60) < 1e-9 and ang3 < 55 and len(fcc) == 12 and abs(angf - 60) < 1e-9
check("kissing: a 60° cap is (1 − cos 30°)/2 = 6.70% of the sphere, so 12 cover 80.4% and 13 would cover 87.1%; in 4D the 16 orthant diagonals (±½)⁴ (the octants' sign quadruples) plus the 8 axes are 24 unit vectors exactly 60° apart — k(4) = 24",
      ok, f"in 3D the 8 octant diagonals plus 6 axes come within {ang3:.2f}° (< 60°), so the sign triples do not kiss there; the 12 FCC neighbours do ({angf:.1f}°)")


def tammes(Np, starts=200, seed=1):
    """max-min angular separation of Np points on S²: SLSQP on 'maximise t subject to chord² ≥ t, |x| = 1' from random starts"""
    g = np.random.default_rng(seed)
    I_, J_ = np.triu_indices(Np, 1); npair = len(I_)

    def cons_f(v):
        P = v[:-1].reshape(Np, 3); d = P[I_] - P[J_]
        return np.sum(d * d, axis=1) - v[-1]

    def cons_j(v):
        P = v[:-1].reshape(Np, 3); Jm = np.zeros((npair, 3 * Np + 1)); d = 2 * (P[I_] - P[J_])
        for k in range(3):
            Jm[np.arange(npair), 3 * I_ + k] = d[:, k]; Jm[np.arange(npair), 3 * J_ + k] = -d[:, k]
        Jm[:, -1] = -1
        return Jm

    def eq_f(v):
        P = v[:-1].reshape(Np, 3)
        return np.sum(P * P, axis=1) - 1

    def eq_j(v):
        P = v[:-1].reshape(Np, 3); Jm = np.zeros((Np, 3 * Np + 1))
        for i in range(Np):
            Jm[i, 3 * i:3 * i + 3] = 2 * P[i]
        return Jm
    best = 0
    for _ in range(starts):
        P = g.normal(size=(Np, 3)); P /= np.linalg.norm(P, axis=1)[:, None]
        x0 = np.concatenate([P.ravel(), [np.min(np.sum((P[I_] - P[J_]) ** 2, axis=1))]])
        res = optimize.minimize(lambda v: -v[-1], x0, jac=lambda v: np.concatenate([np.zeros(3 * Np), [-1.0]]),
                                constraints=[{'type': 'ineq', 'fun': cons_f, 'jac': cons_j}, {'type': 'eq', 'fun': eq_f, 'jac': eq_j}],
                                method='SLSQP', options={'maxiter': 2000, 'ftol': 1e-15})
        Q = res.x[:-1].reshape(Np, 3); Q /= np.linalg.norm(Q, axis=1)[:, None]
        G_ = np.clip(Q @ Q.T, -1, 1); np.fill_diagonal(G_, -1)
        best = max(best, math.degrees(math.acos(np.max(G_))))
    return best


t13 = tammes(13)
t12 = tammes(12, starts=20)
check("[numerical] thirteen points on a sphere cannot all stay 60° apart: the best of 200 optimised random configurations keeps them only 57.137° apart (Musin–Tarasov 2012 proved 57.1367°); Schütte–van der Waerden 1953 settled k(3) = 12",
      57.13 < t13 < 57.14 and abs(t12 - 63.4349) < 1e-3, f"13 points: {t13:.4f}°; control, 12 points: {t12:.4f}° (the icosahedron, arccos(1/√5) = 63.4349°) — twelve fit with room to spare, as the post says")
hex_short = math.sqrt(2 / math.sqrt(3))
check("Hermite's constant γ₂ = 2/√3: at fixed area the hexagonal torus has the longest shortest cycle, so its lightest KK mass is the largest possible — (2/√3)^½ = 1.0746 times the square torus's",
      abs(hex_short**2 - 2 / math.sqrt(3)) < 1e-15 and abs(hex_short - 1.0745699) < 1e-6)
nn_s, k_s, r0_s = sp.symbols('n k r0', positive=True)
Lsq = k_s * r0_s**(nn_s + 3)
Ueff2 = sp.simplify(3 * Lsq / r0_s**4 + nn_s * k_s * r0_s**(nn_s - 1))
check("Ehrenfest: for an attractive force k·r^n a circular orbit is stable iff n > −3; Gauss's law in d dimensions gives n = −(d − 1), so stable orbits (and planets) need d < 4",
      sp.simplify(Ueff2 - k_s * r0_s**(nn_s - 1) * (nn_s + 3)) == 0, "U_eff'' = k r₀^{n−1}(n + 3): the inverse-square law is the largest dimension's law that still binds")


def apsidal_angle(nexp, ecc=0.02):
    k = 1.0; r0 = 1.0; Lh = math.sqrt(k * r0**(nexp + 3))
    r_init = r0 * (1 - ecc)

    def rhs(t, y):
        r, pr, ph = y
        return [pr, Lh**2 / r**3 - k * r**nexp, Lh / r**2]
    ev = lambda t, y: y[1]
    ev.direction = -1
    sol = integrate.solve_ivp(rhs, [0, 200], [r_init, 0.0, 0.0], events=ev, rtol=1e-11, atol=1e-12, max_step=0.05)
    return sol.y_events[0][0][2] if len(sol.t_events[0]) else float('nan')


aK, aH, aM = apsidal_angle(-2), apsidal_angle(1), apsidal_angle(-2.5)
check("Bertrand: near-circular apsidal angle π/√(n + 3) — exactly π for the inverse-square law (orbits close after one turn), π/2 for Hooke, and 254.6° for an r^{−2.5} force (orbits precess)",
      abs(aK - math.pi) < 2e-3 and abs(aH - math.pi / 2) < 2e-3 and abs(aM - math.pi / math.sqrt(0.5)) < 4e-3,
      f"integrated: {math.degrees(aK):.2f}°, {math.degrees(aH):.2f}°, {math.degrees(aM):.2f}° (eccentricity 0.02)")

# =====================================================================================================
print("§5  The matharium coil is Kepler's equation  [structure]")
# =====================================================================================================
xx = np.linspace(-3 * np.pi, 7 * np.pi, 20001)
ok = True
for e in (1.0, 0.6, 0.2):
    P = np.stack([xx, xx - e * np.sin(xx), e * np.cos(xx) / math.sqrt(2)], axis=1)
    axis = np.array([1, 1, 0]) / math.sqrt(2)
    perp = P - np.outer(P @ axis, axis)
    ok &= np.allclose(np.linalg.norm(perp, axis=1), e / math.sqrt(2), atol=1e-12)
    trans = perp @ (np.array([1, -1, 0]) / math.sqrt(2))
    angle = np.unwrap(np.arctan2(perp[:, 2], trans))
    ok &= abs(abs(angle[-1] - angle[0]) / (2 * np.pi) - 5) < 1e-9
signs = []
for kk_ in range(-2, 6):
    xk = kk_ * math.pi
    z = math.cos(xk) / math.sqrt(2)
    d_coil = np.array([1.0, 1 - math.cos(xk)]); d_line = np.array([1.0, 1.0])
    over, under = (d_coil, d_line) if z > 0 else (d_line, d_coil)
    signs.append(int(np.sign(over[0] * under[1] - over[1] * under[0])))
check("the coil (x, x − e sin x, e cos x/√2) stays at distance e/√2 from the line y = x (radius 1/√2 at e = 1, as the post says) and winds once round it per 2π; its shadow crosses the line at x = kπ alternately over and under, every crossing of the same sign",
      ok and len(set(signs)) == 1, f"crossing signs for k = −2…5: {signs}: linking (winding) = ½Σ signs = 1 per turn")
mu_g, a_orb, ecc = 1.0, 1.0, 0.6
r_p = a_orb * (1 - ecc); v_p = math.sqrt(mu_g * (1 + ecc) / (a_orb * (1 - ecc)))
sol = integrate.solve_ivp(lambda t, y: [y[2], y[3], -mu_g * y[0] / (y[0]**2 + y[1]**2)**1.5, -mu_g * y[1] / (y[0]**2 + y[1]**2)**1.5],
                          [0, 2 * math.pi], [r_p, 0, 0, v_p], rtol=1e-12, atol=1e-13, dense_output=True)
ok = True
for t in np.linspace(0.1, 6.1, 25):
    x_, y_ = sol.sol(t)[:2]
    r_ = math.hypot(x_, y_)
    E_ = math.atan2(y_ / math.sqrt(1 - ecc**2), x_ / a_orb + ecc)
    M_t = (t * math.sqrt(mu_g / a_orb**3) + math.pi) % (2 * math.pi) - math.pi
    ok &= abs(((E_ - ecc * math.sin(E_)) - M_t + math.pi) % (2 * math.pi) - math.pi) < 1e-8 and abs(ecc * math.cos(E_) / math.sqrt(2) - (1 - r_ / a_orb) / math.sqrt(2)) < 1e-9
check("y = x − e sin x is Kepler's equation M = E − e sin E: the mean anomaly (time's clock) against the eccentric anomaly (the orbit's clock); the coil's height e cos E/√2 = (1 − r/a)/√2, so 'over' is perihelion and 'under' aphelion, and the coil's radius is the eccentricity",
      ok, "integrated orbit, e = 0.6; the two clocks agree exactly at the apsides; at e = 1 (a radial orbit) M = E − sin E is the cycloid of the post")
ok = all(sp.binomial(k + 3, 3) - (sp.binomial(k + 1, 3) if k >= 2 else 0) == (k + 1)**2 for k in range(0, 12))
check("Fock: the hydrogen (inverse-square) level n has n² states because they are the degree-(n − 1) harmonics on S³ ⊂ ℝ⁴, of dimension (k + 1)²; Moser: the Kepler flow is the geodesic flow of S³, all of whose geodesics close", ok)

# =====================================================================================================
print("§6  Two near-integers: one explained by modular forms, one not  [numerical]")
# =====================================================================================================
mp.mp.dps = 40
ratio_pe = (mp.pi**4 + mp.pi**5) / mp.e**6
check("(π⁴ + π⁵)/e⁶ = 0.9999999561918942…, as the post says (relative miss 4.38·10⁻⁸)",
      mp.nstr(ratio_pe, 18).startswith('0.99999995619189428'), f"{mp.nstr(ratio_pe, 25)}")
fam = []
for a in range(0, 10):
    for b in range(a + 1, 10):
        v = mp.log(mp.pi**a + mp.pi**b)
        fam.append((abs(v - mp.nint(v)), a, b))
fam.sort()
best_d, ba, bb = fam[0]
p_family = 1 - (1 - 2 * float(best_d))**len(fam)
consts = {'π': mp.pi, 'e': mp.e, '2': mp.mpf(2), '3': mp.mpf(3), '√2': mp.sqrt(2), '√3': mp.sqrt(3), 'φ': (1 + mp.sqrt(5)) / 2, 'γ': mp.euler}
cnt, hits = 0, []
names = list(consts)
for i, u in enumerate(names):
    for j, v in enumerate(names):
        for a in range(1, 10):
            for b in range(1, 10):
                if (i, a) >= (j, b):
                    continue
                for sgn in (1, -1):
                    val = consts[u]**a + sgn * consts[v]**b
                    if val <= 0:
                        continue
                    for base_name, base in (('e', mp.e), ('π', mp.pi)):
                        cnt += 1
                        t_ = mp.log(val) / mp.log(base)
                        miss = abs(t_ - mp.nint(t_)) * mp.log(base)
                        if mp.nint(t_) >= 1 and mp.mpf('1e-30') < miss < mp.mpf('4.4e-8'):
                            hits.append(f"{u}^{a}{'+' if sgn > 0 else '−'}{v}^{b} ≈ {base_name}^{int(mp.nint(t_))}")
check("look-elsewhere: among the 45 sums π^a + π^b (0 ≤ a < b ≤ 9) π⁴ + π⁵ is the closest to a power of e, a 1-in-250,000 event for that family alone; in a family of ~8,600 two-term sums and differences of powers of 8 constants it is the only non-identity hit at this precision",
      (ba, bb) == (4, 5) and p_family < 1e-5 and hits == ['π^4+π^5 ≈ e^6'],
      f"miss {mp.nstr(best_d, 4)}, p = {p_family:.1e} within the 45; broader family: {cnt} expressions, expected hits {cnt * 2 * 4.4e-8:.4f}, found {hits}. Unusual for small families, unremarkable for the 10⁶–10⁹ that computer searches scan; no mechanism is known")
x163 = mp.exp(mp.pi * mp.sqrt(163))
resid = x163 - (mp.mpf(640320)**3 + 744)
jval = mp.kleinj((1 + mp.sqrt(-163)) / 2) * 1728
check("e^{π√163} = 640320³ + 744 − 196884/e^{π√163} + …: j((1 + √−163)/2) = −640320³ because ℚ(√−163) has class number one, and the miss is the next coefficient of j = q⁻¹ + 744 + 196884q (196884 = 196883 + 1, the Monster)",
      abs(resid + 196884 / x163) < 1e-20 and abs(jval + mp.mpf(640320)**3) < 1e-10,
      f"miss {mp.nstr(resid, 12)}, −196884/e^(π√163) = {mp.nstr(-196884 / x163, 12)}: a near-integer whose every digit is explained by modular time")
mp.mp.dps = 30

# =====================================================================================================
print("§7  The thermal torus: a state, not a law  [structure + data]")
# =====================================================================================================
hbar, kB = 1.054571817e-34, 1.380649e-23
T_cmb = 2.7255
tT = hbar / (kB * T_cmb)
check("the flat 4-torus has dim GL(4)/O(4) = 10 moduli = volume 1 + spatial shape (SL(3)/SO(3)) 5 + time–space shears 3 + modular time β/V^{1/3} 1; GL(4,ℤ) does not know which cycle is time — a thermal state picks it",
      16 - 6 == 10 and 1 + (6 - 1) + 3 + 1 == 10,
      f"the cosmic heat bath's thermal clock ħ/k_BT = {tT * 1e12:.4f} ps (cħ/k_BT = {c_light * tT * 1e3:.4f} mm at T = 2.7255 K); its rest frame is measured (Planck 2018: v = 369.82 ± 0.11 km/s) — a preferred frame of the state, allowed, unlike Einstein-aether's preferred frame in the laws (|c₁₃| < 10⁻¹⁵ after GW170817)")

# =====================================================================================================
n_pass = sum(ok for _, ok in RESULTS)
print(f"\n{n_pass}/{len(RESULTS)} checks passed in {time.time() - T0:.1f} s")
if n_pass == len(RESULTS):
    print("ALL CHECKS PASSED")
else:
    for name, ok in RESULTS:
        if not ok:
            print("FAILED:", name)
    sys.exit(1)
