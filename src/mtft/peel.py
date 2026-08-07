#!/usr/bin/env python3
"""
mtft.peel — Mellin peel engine for the bulk and skeleton stiffness.
===================================================================
MIT License — (c) 2026 Roger Tano — MTFT Research Program

CORRECTED INDEX (supersedes dictionary F(s) = -zeta(s-1)zeta'(s)):
    sum n^k w_n n^{-s} = -zeta(s-k) * zeta'(s-k+1)          [EXACT, convolution]
    stiffness (k=2):  F(s) = -zeta(s-2) * zeta'(s-1)
    rightmost pole s = 3, residue T_inf = -zeta'(2).

CERTIFIED bulk expansion (studies/peel_2026aug/mu_expansion.py):
    mu(y) = T_inf/(4 pi^3 y^3) + (ln y + gamma - 1)/(8 pi^2 y^2)
            - zeta(3) y/(240 pi) - zeta(5) y^3/(252 pi) - ...
    NO 1/y term, NO constant term (skipped orders; falsifiable).
    gamma enters via psi(2) = 1 - gamma; ln(2pi) from zeta'(0) cancels the
    kernel scale exactly; odd zeta values via trivial-zero derivatives.

CERTIFIED skeleton expansion (studies/peel_2026aug/skeleton_zeros.py):
    S(y) = sum n Lambda(n) e^{-2 pi y n}
         = 1/X^2 - sum_rho Gamma(1+rho) X^{-(1+rho)} + (1 - 12 ln A)
           + X [c2 + H1 - gamma - ln X] + (X^2/2) k2
           + (X^3/6)[c4 + H3 - gamma - ln X] + ...     X = 2 pi y
    Constant term = -zeta'(-1)/zeta(-1) = 1 - 12 ln(Glaisher A).
    Residual * X^{3/2} bounded as y -> 0  <=>  RH (abscissa dichotomy;
    mirrors Addendum M L1/L2).
"""
from __future__ import annotations
import mpmath as mp
from mpmath import mpf, mpc, log, exp, pi, zeta, gamma as GAMMA_F, euler

def w_sieve(N):
    lw = [mpf(0)]*(N+1)
    for d in range(2, N+1):
        c = log(mpf(d))/mpf(d)
        for n in range(d, N+1, d):
            lw[n] += c
    return lw

def lambda_sieve(N):
    spf = list(range(N+1)); i = 2
    while i*i <= N:
        if spf[i] == i:
            for j in range(i*i, N+1, i):
                if spf[j] == j: spf[j] = i
        i += 1
    logp = {}; lam = [None]*(N+1)
    for n in range(2, N+1):
        p = spf[n]; m = n
        while m % p == 0: m //= p
        if m == 1:
            if p not in logp: logp[p] = log(mpf(p))
            lam[n] = logp[p]
    return lam

def F_bulk(s):
    """-zeta(s-2) zeta'(s-1): Dirichlet series of n^2 w_n. Pole s=3, res T_inf."""
    return -zeta(s-2)*mp.diff(zeta, s-1)

def mu_bulk_direct(y, N=None, lw=None):
    X = 2*pi*y
    if N is None: N = int(75/float(X)) + 1
    if lw is None: lw = w_sieve(N)
    r = exp(-X); acc = mpf(1); s = mpf(0)
    for n in range(1, min(N, len(lw)-1)+1):
        acc *= r
        s += mpf(n)**2 * lw[n] * acc
    return s

def mu_bulk_expansion(y, odd_terms=2):
    """Certified asymptotic series. odd_terms counts the zeta(3), zeta(5), ... tail."""
    X = 2*pi*y
    Tinf = -mp.diff(zeta, mpf(2))
    out = Tinf/(4*pi**3*y**3) + (log(y) + euler - 1)/(8*pi**2*y**2)
    for k in range(odd_terms):
        s0 = -(2*k+1)
        Fv = -zeta(mpf(s0-2))*mp.diff(zeta, mpf(s0-1))
        out += (mpf(-1)**(2*k+1)/mp.factorial(2*k+1)) * X**(2*k+1) * Fv
    return out

_SKEL = {}
def skeleton_constants():
    """C0 = -zeta'(-1)/zeta(-1) (= 1 - 12 ln A), c2, k2, c4, k4. Cached."""
    if not _SKEL:
        zp  = lambda w: mp.diff(zeta, mpf(w))
        zpp = lambda w: mp.diff(zeta, mpf(w), 2)
        _SKEL.update(C0=-zp(-1)/zeta(mpf(-1)),
                     c2=zpp(-2)/(2*zp(-2)), k2=-zp(-3)/zeta(mpf(-3)),
                     c4=zpp(-4)/(2*zp(-4)), k4=-zp(-5)/zeta(mpf(-5)))
    return dict(_SKEL)

def skeleton_smooth(y):
    K = skeleton_constants()
    X = 2*pi*y; L = log(X)
    return ( X**-2 + K['C0'] + X*(K['c2'] + 1 - euler - L)
             + X**2/2*K['k2'] + X**3/6*(K['c4'] + mpf(11)/6 - euler - L)
             + X**4/24*K['k4'] )

def skeleton_direct(y, N=None, lam=None):
    X = 2*pi*y
    if N is None: N = int(70/float(X)) + 1
    if lam is None: lam = lambda_sieve(N)
    r = exp(-X); acc = mpf(1); s = mpf(0)
    for n in range(1, min(N, len(lam)-1)+1):
        acc *= r
        if lam[n] is not None:
            s += mpf(n)*lam[n]*acc
    return s

def zeta_zero_sum(y, nzeros=30):
    X = 2*pi*y
    return -sum(2*mp.re(GAMMA_F(1+mp.zetazero(j))*X**(-(1+mp.zetazero(j))))
                for j in range(1, nzeros+1))
