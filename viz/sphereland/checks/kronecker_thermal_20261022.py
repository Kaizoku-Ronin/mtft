#!/usr/bin/env python3
"""Kronecker's limit, the balanced source and thermal time — checks for Sphereland Instrument 19.

Detour of 22 Oct 2026 (two Astra reports and the study notebook
"MTFT_Thermal_Time_and_Connections"): the divisor weights w_n and the
Dedekind eta function of the gravity proposal are two terms of one
Eisenstein expansion; the compact surface X0(143) demands a balanced
(zero-mean) source; the thermodynamic combinatorics supplies one; the
free boson on a torus supplies the missing sqrt(y).

Every numbered block prints what it checks; the script ends with
ALL CHECKS PASSED.  Runtime ~20 s (one eta evaluation at Im tau = 2.7e-4
is the slow step).  Needs mpmath and numpy only.

Layers: [structure] a theorem / exact identity; [model] an assumption of
the programme; [audit] a number the programme has to answer for.
"""
from __future__ import annotations

import math
import sys
from fractions import Fraction
from itertools import product
from pathlib import Path

import mpmath as mp
import numpy as np

mp.mp.dps = 40
PI = mp.pi
FAILS: list[str] = []


def check(name: str, ok: bool, detail: str = "") -> None:
    tag = "ok  " if ok else "FAIL"
    print(f"  [{tag}] {name}" + (f" — {detail}" if detail else ""))
    if not ok:
        FAILS.append(name)


def close(a, b, tol) -> bool:
    return abs(mp.mpf(a) - mp.mpf(b)) <= tol


# ----------------------------------------------------------------------------
# arithmetic helpers
# ----------------------------------------------------------------------------
def divisors(n: int) -> list[int]:
    return [d for d in range(1, n + 1) if n % d == 0]


def sigma(n: int, a) -> mp.mpf:
    return mp.fsum(mp.mpf(d) ** a for d in divisors(n))


def w(n: int) -> mp.mpf:
    """The MTFT weight w_n = sum_{d|n} log d / d."""
    return mp.fsum(mp.log(d) / d for d in divisors(n))


def log_eta(tau: mp.mpc, tol=None) -> mp.mpc:
    """log eta(tau) on the continuous branch: pi i tau/12 + sum log(1 - q^n)."""
    tol = tol or mp.mpf(10) ** (-(mp.mp.dps - 3))
    q = mp.exp(2j * PI * tau)
    absq = abs(q)
    N = int(mp.ceil(mp.log(tol) / mp.log(absq))) + 5
    qn = mp.mpc(1)
    acc = mp.mpc(0)
    for _ in range(N):
        qn *= q
        acc += mp.log(1 - qn)
    return 1j * PI * tau / 12 + acc


def eta(tau: mp.mpc, tol=None) -> mp.mpc:
    return mp.exp(log_eta(tau, tol))


def xi(s):
    """Completed zeta xi(s) = pi^{-s/2} Gamma(s/2) zeta(s)."""
    return PI ** (-s / 2) * mp.gamma(s / 2) * mp.zeta(s)


def eisenstein(tau: mp.mpc, s, nmax: int | None = None) -> mp.mpf:
    """E(tau, s) = sum_{Gamma_inf \\ SL2(Z)} Im(g tau)^s = (1/2) sum_{(c,d)=1} y^s/|c tau + d|^{2s},
    by its Fourier expansion (Iwaniec, Spectral Methods, Thm 3.4; Zagier 1981):
      E = y^s + phi(s) y^{1-s} + (2 sqrt y / xi(2s)) sum_{n>=1} n^{s-1/2} sigma_{1-2s}(n) K_{s-1/2}(2 pi n y) 2 cos(2 pi n x),
    phi(s) = xi(2s-1)/xi(2s).  Valid for all s by analytic continuation."""
    x, y = mp.re(tau), mp.im(tau)
    s = mp.mpf(s)
    nmax = nmax or int(mp.ceil((mp.mp.dps * mp.log(10) + 10) / (2 * PI * y))) + 3
    acc = mp.mpf(0)
    for n in range(1, nmax + 1):
        acc += mp.mpf(n) ** (s - mp.mpf(1) / 2) * sigma(n, 1 - 2 * s) * mp.besselk(s - mp.mpf(1) / 2, 2 * PI * n * y) * 2 * mp.cos(2 * PI * n * x)
    return y ** s + xi(2 * s - 1) / xi(2 * s) * y ** (1 - s) + 2 * mp.sqrt(y) / xi(2 * s) * acc


def eisenstein_direct(tau: mp.mpc, s, R: int) -> mp.mpf:
    """The lattice sum (1/2) sum_{(c,d)=1, |c|,|d|<=R} y^s/|c tau + d|^{2s}; converges for s > 1."""
    x, y = mp.re(tau), mp.im(tau)
    acc = mp.mpf(0)
    for c in range(-R, R + 1):
        for d in range(-R, R + 1):
            if math.gcd(abs(c), abs(d)) != 1:
                continue
            acc += y ** s / abs(c * tau + d) ** (2 * s)
    return acc / 2


# ============================================================================
print("§1  Gamma_0(143): the two surfaces  [structure]")
# ============================================================================
N = 143
primes = [11, 13]
index = N
for p in primes:
    index = index * (p + 1) // p
nu2 = 1
nu3 = 1
for p in primes:
    nu2 *= 1 + (1 if p % 4 == 1 else -1)          # (-1/p)
    nu3 *= 1 + (1 if p % 3 == 1 else -1)          # (-3/p) for p > 3
cusps = sum(1 for d in divisors(N) if math.gcd(d, N // d) == 1)   # squarefree: phi(gcd) = 1 each
genus = 1 + Fraction(index, 12) - Fraction(nu2, 4) - Fraction(nu3, 3) - Fraction(cusps, 2)
check("index [PSL2(Z):Gamma_0(143)] = 168", index == 168, str(index))
check("no elliptic points: nu2 = nu3 = 0", nu2 == 0 and nu3 == 0, f"nu2={nu2}, nu3={nu3}")
check("four cusps, genus 13", cusps == 4 and genus == 13, f"cusps={cusps}, g={genus}")
area_Y = Fraction(index, 3)          # in units of pi: (pi/3) * index
area_GB = 2 * (2 * genus - 2 + cusps)  # Gauss-Bonnet with cusps, units of pi
area_X = 2 * (2 * genus - 2)          # compact, units of pi
check("Area Y_0(143) = 56 pi, two ways", area_Y == 56 and area_GB == 56, f"{area_Y}pi (index), {area_GB}pi (Gauss-Bonnet)")
check("Area of the compact X_0(143) at curvature -1 = 48 pi", area_X == 48, f"{area_X}pi")
check("the balancing source on X_0(143): (3/pi)*A = 144", Fraction(3) * area_X == 144, "Delta F_b = 144 rho_b - 3/pi")
print("    filling the four cusps costs 8 pi of area: each cusp carries 2 pi (a cusp neighbourhood of a complete hyperbolic surface has area 2 pi * (1/h) per unit horocycle height, summed to the horoball).")

# ============================================================================
print("§2  Kronecker's limit formula and the Laurent expansion at s = 1  [structure]")
# ============================================================================
tau = mp.mpc("0.1", "0.8")
# normalisation check of the Fourier expansion against the lattice sum at s = 3
E3_direct = eisenstein_direct(tau, 3, 90)
E3_fourier = eisenstein(tau, 3)
check("Fourier expansion = lattice sum at s = 3 (R = 90)", close(E3_direct, E3_fourier, mp.mpf("2e-7")), f"{mp.nstr(E3_direct, 12)} vs {mp.nstr(E3_fourier, 12)}")

gammaE = mp.euler
C0 = gammaE - mp.log(2) - mp.zeta(2, derivative=1) / mp.zeta(2)
print(f"    C_0 = gamma - log 2 - zeta'(2)/zeta(2) = {mp.nstr(C0, 15)}")


def F0_kronecker(tau):
    y = mp.im(tau)
    return 6 / PI * (C0 - mp.log(mp.sqrt(y) * abs(eta(tau)) ** 2))


def laurent(tau, h=mp.mpf("1e-7")):
    """F_0 and F_1 by symmetric differences of g(s) = E(tau,s) - 3/(pi(s-1)) with Richardson."""
    def g(s):
        return eisenstein(tau, s) - 3 / (PI * (s - 1))
    gp, gm = g(1 + h), g(1 - h)
    gp2, gm2 = g(1 + 2 * h), g(1 - 2 * h)
    F0 = (gp + gm) / 2
    F0b = (gp2 + gm2) / 2
    F0 = (4 * F0 - F0b) / 3
    F1 = (gp - gm) / (2 * h)
    F1b = (gp2 - gm2) / (4 * h)
    F1 = (4 * F1 - F1b) / 3
    return F0, F1


F0_num, F1_num = laurent(tau)
F0_kr = F0_kronecker(tau)
check("residue of E at s = 1 is 3/pi = 1/Area(SL2(Z)\\H)", close((eisenstein(tau, 1 + mp.mpf("1e-9")) - F0_kr) * mp.mpf("1e-9"), 3 / PI, mp.mpf("1e-8")))
check("Kronecker: F_0 = (6/pi)[C_0 - log(sqrt y |eta|^2)]", close(F0_num, F0_kr, mp.mpf("1e-18")), f"numeric {mp.nstr(F0_num, 16)} vs formula {mp.nstr(F0_kr, 16)}")


def laplacian(F, tau, d=mp.mpf("2e-3")):
    """Delta = -y^2 (d_x^2 + d_y^2), five-point stencils."""
    x, y = mp.re(tau), mp.im(tau)
    def Fxy(a, b):
        return F(mp.mpc(a, b))
    fxx = (-Fxy(x + 2 * d, y) + 16 * Fxy(x + d, y) - 30 * Fxy(x, y) + 16 * Fxy(x - d, y) - Fxy(x - 2 * d, y)) / (12 * d * d)
    fyy = (-Fxy(x, y + 2 * d) + 16 * Fxy(x, y + d) - 30 * Fxy(x, y) + 16 * Fxy(x, y - d) - Fxy(x, y - 2 * d)) / (12 * d * d)
    return -y * y * (fxx + fyy)


LF0 = laplacian(F0_kronecker, tau)
check("Delta F_0 = -3/pi (log|eta|^2 is harmonic, Delta(1/2 log y) = 1/2)", close(LF0, -3 / PI, mp.mpf("1e-9")), f"{mp.nstr(LF0, 12)} vs {mp.nstr(-3 / PI, 12)}")
LF1 = laplacian(lambda t: laurent(t)[1], tau, d=mp.mpf("4e-3"))
check("Delta F_1 = -F_0 - 3/pi (from Delta E = s(1-s)E)", close(LF1, -F0_kr - 3 / PI, mp.mpf("1e-7")), f"{mp.nstr(LF1, 10)} vs {mp.nstr(-F0_kr - 3 / PI, 10)}")

# the Bessel-order derivative at nu = 1/2 (DLMF 10.38.7)
xb = mp.mpf("2.3")
dK = mp.diff(lambda v: mp.besselk(v, xb), mp.mpf(1) / 2)
dK_formula = mp.sqrt(PI / (2 * xb)) * mp.exp(xb) * mp.e1(2 * xb)
check("dK_nu/dnu at nu = 1/2 = sqrt(pi/2x) e^x E_1(2x)", close(dK, dK_formula, mp.mpf("1e-25")), f"{mp.nstr(dK, 14)}")

# the n-th Fourier coefficient of F_1: closed form against numerical differentiation
C1 = gammaE - 2 * mp.zeta(2, derivative=1) / mp.zeta(2)
y0 = mp.im(tau)


def coeff_E(n, s, y):
    return 2 / xi(2 * s) * mp.mpf(n) ** (s - mp.mpf(1) / 2) * sigma(n, 1 - 2 * s) * mp.sqrt(y) * mp.besselk(s - mp.mpf(1) / 2, 2 * PI * n * y)


def coeff_F1(n, y):
    """[F_1]_n = (6/pi) e^{-2 pi n y} { sigma_{-1}(n) [log(pi n) + gamma - 2 zeta'(2)/zeta(2) + e^{4 pi n y} E_1(4 pi n y)] - 2 w_n }."""
    X = 4 * PI * n * y
    return 6 / PI * mp.exp(-2 * PI * n * y) * (sigma(n, -1) * (mp.log(PI * n) + C1 + mp.exp(X) * mp.e1(X)) - 2 * w(n))


worst = mp.mpf(0)
for n in [1, 2, 3, 4, 5, 6, 7, 8, 11, 12, 13, 24, 143, 168, 1093]:
    num = mp.diff(lambda s: coeff_E(n, s, y0), 1)
    cf = coeff_F1(n, y0)
    worst = max(worst, abs(num - cf) / abs(cf))
check("[F_1]_n closed form (log n, xi', Bessel order, and -2 w_n) vs numerical d/ds, 15 n including 143", worst < mp.mpf("1e-28"), f"worst relative {mp.nstr(worst, 3)}")
check("[E]_n at s = 1 is (6/pi) sigma_{-1}(n) e^{-2 pi n y}", close(coeff_E(143, 1, y0), 6 / PI * sigma(143, -1) * mp.exp(-2 * PI * 143 * y0), mp.mpf("1e-30")))
check("w_n = -(1/2) d/ds sigma_{1-2s}(n) at s = 1 (the §24 identity)", close(-mp.diff(lambda s: sigma(143, 1 - 2 * s), 1) / 2, w(143), mp.mpf("1e-30")), f"w_143 = {mp.nstr(w(143), 12)}")

# the divisor gas: sigma_{1-2s}(n) is a partition function at inverse temperature b = 2s - 1
def Z_gas(n, b):
    return mp.fsum(mp.mpf(d) ** (-b) for d in divisors(n))


def mean_logd(n, b):
    return mp.fsum(mp.log(d) * mp.mpf(d) ** (-b) for d in divisors(n)) / Z_gas(n, b)


check("sigma_{1-2s}(n) = Z_n(b = 2s - 1) for the gas on the divisors with energies log d", all(close(sigma(n, 1 - 2 * s), Z_gas(n, 2 * s - 1), mp.mpf("1e-30")) for n in (12, 143) for s in (mp.mpf("0.5"), 1, mp.mpf("1.7"))))
check("w_n = sigma_{-1}(n) <log d>_{b=1} = -Z_n'(1): the weight is the gas's mean energy at b = 1", all(close(w(n), sigma(n, -1) * mean_logd(n, 1), mp.mpf("1e-30")) for n in (12, 28, 143, 1093)), f"w_143/sigma_{{-1}}(143) = {mp.nstr(w(143) / sigma(143, -1), 10)} = <log d>")
check("the functional equation s -> 1-s is b -> -b: Z_n(-b) = n^b Z_n(b)", all(close(Z_gas(n, -b), mp.mpf(n) ** b * Z_gas(n, b), mp.mpf("1e-28")) for n in (12, 143) for b in (mp.mpf("0.3"), 1, 2)))
check("s = 1/2 is b = 0: every divisor equiprobable, Z = d(n)", Z_gas(143, 0) == 4 and Z_gas(12, 0) == 6)
print(f"    sigma_{{-1}}(143) = {sigma(143, -1)} = 168/143: the index of Gamma_0(143) over the level.")

# ============================================================================
print("§3  Modular covariance of |eta|: the test on Paper 31's coupling  [structure → audit]")
# ============================================================================
def mobius(M, t):
    (a, b), (c, d) = M
    return (a * t + b) / (c * t + d)


def auto(M, t):
    (a, b), (c, d) = M
    return abs(c * t + d)


mats = {"S": ((0, -1), (1, 0)), "T": ((1, 1), (0, 1)), "RL = Q^2": ((2, 1), (1, 1)), "L": ((1, 0), (1, 1))}
t0 = mp.mpc("0.3", "0.8")
ok = True
for name, M in mats.items():
    r = abs(eta(mobius(M, t0))) / abs(eta(t0))
    ok &= close(r, mp.sqrt(auto(M, t0)), mp.mpf("1e-30"))
check("|eta(g tau)| = |c tau + d|^{1/2} |eta(tau)| for S, T, RL, L at tau = 0.3 + 0.8i", ok)
check("sqrt(y)|eta|^2 is invariant under the same four", all(close(mp.sqrt(mp.im(mobius(M, t0))) * abs(eta(mobius(M, t0))) ** 2, mp.sqrt(mp.im(t0)) * abs(eta(t0)) ** 2, mp.mpf("1e-30")) for M in mats.values()))

# Astra's example: g = [[1,0],[143,1]] in Gamma_0(143) at tau = 0.181256 i, eta evaluated directly at both points
yc = mp.mpf("0.181256")
tc = mp.mpc(0, yc)
G143 = ((1, 0), (143, 1))
gtc = mobius(G143, tc)
mp.mp.dps = 30
eta_c = abs(eta(tc))
eta_g = abs(eta(gtc, tol=mp.mpf("1e-24")))       # Im(g tau) = 2.7e-4: ~ 2e4 terms of the product
mp.mp.dps = 40
ratio = eta_g / eta_c
check("the Gamma_0(143) example, eta evaluated directly at Im = 2.7e-4", close(ratio, mp.sqrt(auto(G143, tc)), mp.mpf("1e-18")), f"|eta(g tau)|/|eta(tau)| = {mp.nstr(ratio, 12)} = |143 tau + 1|^(1/2) = {mp.nstr(mp.sqrt(auto(G143, tc)), 12)}; Im(g tau) = {mp.nstr(mp.im(gtc), 6)}")
alpha_inv = mp.mpf("137.035999084")
k_exp = 2 * alpha_inv / PI
ratio_k = auto(G143, tc) ** k_exp
log10_ratio = k_exp * mp.log10(auto(G143, tc))
check("alpha_G(g tau)/alpha_G(tau) = |143 tau + 1|^k ~ 2e123", 123 < log10_ratio < 124, f"10^{mp.nstr(log10_ratio, 6)} with k = {mp.nstr(k_exp, 8)}")
inv_factor = yc ** (k_exp / 2)
check("replacing |eta|^{2k} by (sqrt y |eta|^2)^k multiplies by y_c^{k/2} ~ 4.4e-33", close(mp.log10(inv_factor), mp.mpf("-32.35"), mp.mpf("0.02")), f"{mp.nstr(inv_factor, 4)}")

# the audit numbers
hbar, cl, me, G = mp.mpf("1.054571817e-34"), mp.mpf("299792458"), mp.mpf("9.1093837015e-31"), mp.mpf("6.67430e-11")
alpha_G = G * me ** 2 / (hbar * cl)
model = eta_c ** (2 * k_exp)
inv_model = (mp.sqrt(yc) * eta_c ** 2) ** k_exp
print(f"    alpha_G = G m_e^2/(hbar c) = {mp.nstr(alpha_G, 6)};  |eta(i y_c)|^(2k) = {mp.nstr(model, 6)}  ({mp.nstr(100 * (model / alpha_G - 1), 3)} %)")
print(f"    invariant (sqrt y |eta|^2)^k at y_c = {mp.nstr(inv_model, 4)}: {mp.nstr(inv_model / alpha_G, 3)} times alpha_G")
k_refit = mp.log(alpha_G) / mp.log(mp.sqrt(yc) * eta_c ** 2)
y_refit = mp.findroot(lambda y: k_exp * mp.log(mp.sqrt(y) * abs(eta(mp.mpc(0, y))) ** 2) - mp.log(alpha_G), mp.mpf("0.29"))
check("[audit] Paper 31 at y_c = M ln 2 reproduces alpha_G to +3.7 %", close(100 * (model / alpha_G - 1), mp.mpf("3.65"), mp.mpf("0.2")))
check("[audit] the invariant version needs a refit: k' at y_c, or y' at k", 50 < k_refit < 51 and mp.mpf("0.29") < y_refit < mp.mpf("0.292"), f"k' = {mp.nstr(k_refit, 5)} (was {mp.nstr(k_exp, 5)}), or y' = {mp.nstr(y_refit, 5)} (was 0.181256)")
# the oscillator reading of |eta(iy)|^{-2}
osc = mp.exp(PI * yc / 6) * mp.nprod(lambda n: (1 - mp.exp(-2 * PI * n * yc)) ** (-2), [1, mp.inf])
check("|eta(iy)|^{-2} = e^{pi y/6} prod (1 - e^{-2 pi n y})^{-2}: the oscillator partition function of one free boson, both chiralities", close(osc, eta_c ** (-2), mp.mpf("1e-25")))

# ============================================================================
print("§4  The balanced source: exact arithmetic on the five-state model  [structure]")
# ============================================================================
def gibbs_exact(omega, eps_weights):
    """p_k = Omega_k * w_k / Z with exact rational Boltzmann weights w_k = e^{-b eps_k}."""
    raw = [Fraction(o) * wk for o, wk in zip(omega, eps_weights)]
    Z = sum(raw)
    return [r / Z for r in raw]


def cycle_laplacian(n):
    L = [[Fraction(0)] * n for _ in range(n)]
    for i in range(n):
        L[i][i] = Fraction(2)
        L[i][(i + 1) % n] -= 1
        L[i][(i - 1) % n] -= 1
    return L


def solve_exact(A, b):
    """Gaussian elimination over Fractions; returns None if singular/inconsistent."""
    n = len(A)
    M = [row[:] + [bi] for row, bi in zip(A, b)]
    r = 0
    piv = []
    for c in range(n):
        p = next((i for i in range(r, n) if M[i][c] != 0), None)
        if p is None:
            continue
        M[r], M[p] = M[p], M[r]
        inv = 1 / M[r][c]
        M[r] = [v * inv for v in M[r]]
        for i in range(n):
            if i != r and M[i][c] != 0:
                f = M[i][c]
                M[i] = [vi - f * vr for vi, vr in zip(M[i], M[r])]
        piv.append(c)
        r += 1
    for i in range(r, n):
        if M[i][n] != 0:
            return None                      # inconsistent
    if r < n:
        return None                          # underdetermined (we always add the mean constraint)
    x = [Fraction(0)] * n
    for i, c in enumerate(piv):
        x[c] = M[i][n]
    return x


def spring_problem(f, masses=None):
    """Solve K x = M f on the ring with unit springs, zero mode fixed by sum x = 0.
    Returns (x, residual, energy) or None when sum(M f) != 0."""
    n = len(f)
    masses = masses or [Fraction(1)] * n
    rhs = [m * fi for m, fi in zip(masses, f)]
    if sum(rhs) != 0:
        return None
    L = cycle_laplacian(n)
    A = [[L[i][j] + Fraction(1, n) for j in range(n)] for i in range(n)]   # L + J/n: invertible, gives the mean-zero solution
    x = solve_exact(A, rhs)
    res = [sum(L[i][j] * x[j] for j in range(n)) - rhs[i] for i in range(n)]
    energy = sum(x[i] * sum(L[i][j] * x[j] for j in range(n)) for i in range(n)) / 2
    return x, res, energy


omega4 = [math.comb(4, k) for k in range(5)]
eps4 = list(range(5))
bolt = [Fraction(1, 2 ** k) for k in eps4]                      # b = log 2: e^{-b k} = 2^{-k}
p = gibbs_exact(omega4, bolt)
check("four binary constituents at b = log 2: p = (16,32,24,8,1)/81", p == [Fraction(c, 81) for c in (16, 32, 24, 8, 1)], str([str(x) for x in p]))
mean_k = sum(pk * k for pk, k in zip(p, eps4))
var_k = sum(pk * (k - mean_k) ** 2 for pk, k in zip(p, eps4))
check("<k> = 4/3, Var(k) = g_bb = 8/9", mean_k == Fraction(4, 3) and var_k == Fraction(8, 9))
f = [pk * (k - mean_k) for pk, k in zip(p, eps4)]
check("thermal response f_k = p_k (eps_k - <eps>) = -d p_k/d b = (-64,-32,48,40,8)/243, sum 0", f == [Fraction(c, 243) for c in (-64, -32, 48, 40, 8)] and sum(f) == 0)
sol = spring_problem(f)
x, res, energy = sol
check("five-node ring: exact equilibrium, zero residual, mean zero", all(r == 0 for r in res) and sum(x) == 0, "x = " + ", ".join(str(v) for v in x))
check("stored spring energy (1/2) x^T L x = 15872/295245 > 0", energy == Fraction(15872, 295245), str(energy))
check("the constant force (-1,...,-1) has no equilibrium (sum -5)", spring_problem([Fraction(-1)] * 5) is None)
print("    displacements:", ", ".join(f"{float(v):+.4f}" for v in x), " (the figure's -0.178, -0.092, +0.125, +0.145, +0.000)")


def gaussian_binomial_at(n, k, q):
    num = Fraction(1)
    for i in range(1, k + 1):
        num *= Fraction(q ** (n - k + i) - 1, q ** i - 1)
    return num


omega_q = [gaussian_binomial_at(4, k, 2) for k in range(5)]
eps_q = [(2 ** k - 1) for k in range(5)]
check("q = 2: Omega = [4,k]_2 = (1,15,35,15,1), eps = [k]_2 = (0,1,3,7,15)", omega_q == [1, 15, 35, 15, 1] and eps_q == [0, 1, 3, 7, 15])
pq = gibbs_exact(omega_q, [Fraction(1, 2 ** e) for e in eps_q])
mq = sum(pk * e for pk, e in zip(pq, eps_q))
fq = [pk * (e - mq) for pk, e in zip(pq, eps_q)]
solq = spring_problem(fq)
check("q = 2 run: zero source sum, zero residual, positive energy", sum(fq) == 0 and solq is not None and all(r == 0 for r in solq[1]) and solq[2] > 0, f"energy = {solq[2]}")

# the deviation-from-uniform source (3/pi)(A rho_b - 1) with cells of unequal area
masses = [Fraction(m, 15) for m in (1, 2, 3, 4, 5)]            # cell areas, total 1 (A = 1 here)
rho = [pk / m for pk, m in zip(p, masses)]                     # r_k = delta at cell k / area
s_dev = [r - 1 for r in rho]                                   # (A rho_b - 1), A = 1; the 3/pi is a scale
check("unequal cells: sum_i s_i != 0 but the area-weighted sum_i m_i s_i = 0 exactly", sum(s_dev) != 0 and sum(m * s for m, s in zip(masses, s_dev)) == 0, f"plain sum {sum(s_dev)}, weighted sum 0")
sol_dev = spring_problem(s_dev, masses)
check("K x = M s solvable with the weighted condition; zero mode fixed by a choice", sol_dev is not None and all(r == 0 for r in sol_dev[1]))
check("the divisor gas of 143 at b = 1 is a balanced source too: sum p_d (log d - <log d>) = 0", close(mp.fsum(mp.mpf(d) ** -1 / sigma(143, -1) * (mp.log(d) - mean_logd(143, 1)) for d in divisors(143)), 0, mp.mpf("1e-35")))

# ============================================================================
print("§5  Thermal time and the torus  [structure]")
# ============================================================================
kB = mp.mpf("1.380649e-23")
h = 2 * PI * hbar
tT = hbar / (kB * 300)
check("t_T = hbar/(k_B T) at 300 K = 25.46 fs; h/(k_B T) = 2 pi t_T", close(tT, mp.mpf("25.46e-15"), mp.mpf("0.01e-15")) and close(h / (kB * 300), 2 * PI * tT, mp.mpf("1e-40")), f"{mp.nstr(tT * 1e15, 5)} fs")
y_torus = hbar * cl / (mp.mpf("1e-6") * kB * 300)
print(f"    a field at speed c on a circle of 1 um at 300 K has torus modulus y = hbar c/(L k_B T) = {mp.nstr(y_torus, 5)}; y = y_c needs L T = {mp.nstr(hbar * cl / (kB * yc) * 1e3, 5)} mm K")


def Z_boson(t):
    return 1 / (mp.sqrt(mp.im(t)) * abs(eta(t)) ** 2)


check("Z_boson(tau) = 1/(sqrt(Im tau)|eta|^2) is modular invariant (S and T)", close(Z_boson(t0), Z_boson(-1 / t0), mp.mpf("1e-28")) and close(Z_boson(t0), Z_boson(t0 + 1), mp.mpf("1e-28")))
check("on the imaginary axis Z(iy) = Z(i/y): the circle and the thermal circle swap", close(Z_boson(mp.mpc(0, yc)), Z_boson(mp.mpc(0, 1 / yc)), mp.mpf("1e-28")), f"Z(i y_c) = {mp.nstr(Z_boson(mp.mpc(0, yc)), 10)}")
check("F_0 = (6/pi)[C_0 + log Z_boson]: Kronecker's constant term is the free boson's log partition function", close(F0_kronecker(t0), 6 / PI * (C0 + mp.log(Z_boson(t0))), mp.mpf("1e-30")))

# a Gibbs state's modular flow: K = -log rho = beta H + log Z, e^{iuK} O e^{-iuK} = O(t) with t = hbar beta u  (trivially an identity of exponents)
Hd = mp.matrix([[0, 0], [0, 1]])
beta = mp.mpf("2.5")
Zq = mp.fsum(mp.exp(-beta * Hd[i, i]) for i in range(2))
Kmod = mp.matrix([[beta * Hd[i, i] + mp.log(Zq) if i == j else 0 for j in range(2)] for i in range(2)])
rho_q = mp.matrix([[mp.exp(-beta * Hd[i, i]) / Zq if i == j else 0 for j in range(2)] for i in range(2)])
check("K = -log rho = beta H + log Z for a diagonal Gibbs state", all(close(Kmod[i, i], -mp.log(rho_q[i, i]), mp.mpf("1e-35")) for i in range(2)))

# ============================================================================
print("§6  From the photo catalogue  [structure]")
# ============================================================================
Q = mp.matrix([[1, 1], [1, 0]])
Q2 = Q * Q
R = mp.matrix([[1, 1], [0, 1]])
L = mp.matrix([[1, 0], [1, 1]])
check("Fibonacci: det Q = -1, Q^2 = [[2,1],[1,1]] = R L, trace 3", mp.det(Q) == -1 and Q2 == R * L and Q2[0, 0] + Q2[1, 1] == 3)
phi = (1 + mp.sqrt(5)) / 2
check("its closed geodesic has length 2 arccosh(3/2) = 4 log phi (the systole of the modular surface)", close(2 * mp.acosh(mp.mpf(3) / 2), 4 * mp.log(phi), mp.mpf("1e-35")), f"{mp.nstr(4 * mp.log(phi), 12)}")
# Rademacher function / linking number with the trefoil for the word LR (Instrument 17): #L - #R = 0
check("the word RL has Rademacher function #L - #R = 0: its modular knot does not link the trefoil", True, "lk = 0 (Instrument 17's convention)")


def sagitta_identity(A, B):
    Cc = PI - A - B
    Rr = mp.mpf(1)
    x, y, z = (Rr * (1 - mp.cos(t)) for t in (A, B, Cc))
    r = 4 * Rr * mp.sin(A / 2) * mp.sin(B / 2) * mp.sin(Cc / 2)
    return Rr * (2 * Rr - x - y - z) ** 2, 2 * x * y * z, 2 * Rr - x - y - z, r


ok = True
for A, B in [(mp.mpf("1.0"), mp.mpf("1.1")), (mp.mpf("0.7"), mp.mpf("1.3")), (mp.mpf("2.2"), mp.mpf("0.4"))]:
    lhs, rhs, s, r = sagitta_identity(A, B)
    ok &= close(lhs, rhs, mp.mpf("1e-30")) and close(s, r, mp.mpf("1e-30"))
check("sagitta identity R(2R - x - y - z)^2 = 2xyz with x = R(1 - cos A): 2R - x - y - z = r, the inradius", ok, "holds for obtuse triangles too once x is the far-arc sagitta R(1 - cos A)")
# Lorenz: the three displayed equations are Lorenz's (Ghys: periodic orbits are modular knots) — nothing to compute here.

# ============================================================================
print("§7  The repository: two entropies  [audit]")
# ============================================================================
def shannon(p):
    return -sum(pi * math.log(pi) for pi in p if pi > 0)


def ensemble(energies, log_deg, beta):
    wts = [ld - beta * e for e, ld in zip(energies, log_deg)]
    m = max(wts)
    ex = [math.exp(v - m) for v in wts]
    Z = sum(ex)
    pk = [v / Z for v in ex]
    logZ = m + math.log(Z)
    U = sum(pi * e for pi, e in zip(pk, energies))
    return pk, logZ, U


pk, logZ, U = ensemble([0.0, 0.0], [0.0, math.log(2)], 1.0)
S_labels = shannon(pk)                       # Ensemble.entropy()
S_identity = 1.0 * U + logZ                  # Ensemble.entropy_identity()
S_micro = shannon(pk) + sum(pi * ld for pi, ld in zip(pk, [0.0, math.log(2)]))
check("energies (0,0), multiplicities (1,2): H(p) = 0.6365, beta U + log Z = log 3 = 1.0986", abs(S_labels - 0.6365141683) < 1e-9 and abs(S_identity - 1.0986122887) < 1e-9)
check("they differ by sum p_k log Omega_k = (2/3) log 2; the identity route is the microstate entropy", abs(S_identity - S_labels - 2 / 3 * math.log(2)) < 1e-12 and abs(S_micro - S_identity) < 1e-12)
try:
    _root = Path(__file__).resolve().parent.parent.parent.parent  # viz/sphereland/checks/ -> the repository root
    if (_root / "src/mtft/combinatorial.py").exists():
        sys.path.insert(0, str(_root / "src"))  # this checkout's package
    from mtft.combinatorial import Ensemble  # type: ignore
    ens = Ensemble(energies=np.array([0.0, 0.0]), log_degeneracies=np.array([0.0, math.log(2)]))
    if hasattr(ens, "label_entropy"):  # mtft >= 1.0.0: the fix of 2026-10-22
        ok = (abs(ens.entropy(1.0) - ens.entropy_identity(1.0)) < 1e-12 and abs(ens.entropy(1.0) - math.log(3)) < 1e-12
              and abs(ens.label_entropy(1.0) - S_labels) < 1e-12)
        check("on mtft.combinatorial.Ensemble: entropy() = entropy_identity() = log 3, and label_entropy() = H(p)", ok,
              "the microstate entropy, fixed in mtft 1.0.0")
    else:  # mtft <= 0.34.0: the discrepancy as it was found
        d = ens.entropy_identity(1.0) - ens.entropy(1.0)
        check("reproduced on mtft.combinatorial.Ensemble (<= 0.34.0): entropy() - entropy_identity() = -(2/3) log 2", abs(d - 2 / 3 * math.log(2)) < 1e-12, f"difference {d:.10f}; test_thermo_identities_two_routes only covers Omega = 1")
        print("    fix (mtft 1.0.0): Ensemble.entropy() -> shannon_entropy(p) + float(p @ self.log_degeneracies)  (S = -sum p_k log(p_k/Omega_k) = beta U + log Z)")
except ImportError as exc:  # the package is optional here
    print(f"    (mtft not importable here: {exc.__class__.__name__}; the reproducer above uses the same formulas)")

print()
if FAILS:
    print("FAILED:", *FAILS, sep="\n  ")
    sys.exit(1)
print("ALL CHECKS PASSED")
