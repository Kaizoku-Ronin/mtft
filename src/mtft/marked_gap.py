#!/usr/bin/env python3
"""
mtft.marked_gap — Rung-4 mass gap in two-temperature ensemble language.
=======================================================================
MIT License — (c) 2026 Roger Tano — MTFT Research Program

EXACT factorization (log = Lambda * 1):
    Z_D(beta) = sum w_n n^-beta = zeta(beta) * zeta(beta+1) * sum Lambda(m) m^-(beta+1)
Two gas temperatures one unit apart plus ONE marked excitation.
The n^2-tilted (stiffness) ensemble is Z_D(s-2); condensation at s = 3
(the pole whose residue is T_inf).

IDENTITY [EXACT, 32 digits]: the Rung-4 strong-coupling gap is the spacing of
the marked-primon Boltzmann spectrum at beta_c = 3:
    lambda_m = Lambda(m) m^-3,   eps_m = 3 ln m - ln Lambda(m)
    m_inf = eps_3 - eps_2 = ln(27 ln2/(8 ln3)) = 0.7558345761261298
    R(inf) = m_inf/(1 - e^{-m_inf}) = 1.425077237462  (generator<->transfer)
So T_inf (residue) and m_inf (level spacing) are attributes of ONE spectral
point s = 3.

AUDIT FLAG (open): the three-factor pole bookkeeping places the marked factor
at exponent 2 at condensation; the chain's ledger m_inf sits at exponent 3
(mark carries the full untilted n^-s).  Delta(sigma) = sigma ln(3/2) - ln(log2 3)
is linear: Delta(3) = Delta(2) + ln(3/2).  Derive the sigma = 3 assignment from
the Rung-4 site measure (PR-5.1) rather than matching it.

FALSIFIABLE SPECTRUM PREDICTION at beta_c = 3 (audit vs studies/ chain record):
    ordering  2, 3, 5, 4, 7, 11, 9, 8, 13   (p=11 BELOW the powers 9 and 8)
    second gap eps_5 - eps_3 = ln(125 ln3/(27 ln5)) = 1.1506397
Gap closes at sigma* = ln(log2 3)/ln(3/2) = 1.13588256792 (level crossing;
p=3 becomes the ground mark below it).  [Prox: ln pi is 0.77% away — index only.]
"""
from __future__ import annotations
import mpmath as mp
from mpmath import mpf, log, exp

LN2 = log(mpf(2)); LN3 = log(mpf(3))
M_INF = log(27*LN2/(8*LN3))
R_INF = M_INF/(1 - exp(-M_INF))
SIGMA_STAR = log(LN3/LN2)/log(mpf(3)/2)

def Lambda_of(m):
    p = None; q = m
    for c in range(2, m+1):
        if m % c == 0: p = c; break
    while q % p == 0: q //= p
    return log(mpf(p)) if q == 1 else None

def eps_level(m, sigma=3):
    L = Lambda_of(m)
    if L is None: raise ValueError(f"{m} is not a prime power")
    return sigma*log(mpf(m)) - log(L)

def delta_gap(sigma):
    return sigma*log(mpf(3)/2) - log(LN3/LN2)

def predicted_spectrum(sigma=3, marks=(2,3,4,5,7,8,9,11,13)):
    return sorted((eps_level(m, sigma), m) for m in marks)

def zd_factorization_check(beta, Nw=2000):
    """E2: direct sum w_n n^-beta vs zeta(b) zeta(b+1) (-zeta'/zeta)(b+1)."""
    lw = [mpf(0)]*(Nw+1)
    for d in range(2, Nw+1):
        c = log(mpf(d))/mpf(d)
        for n in range(d, Nw+1, d): lw[n] += c
    direct = sum(lw[n]*mpf(n)**(-beta) for n in range(1, Nw+1))
    closed = -mp.zeta(beta)*mp.diff(mp.zeta, beta+1)
    return direct, closed
