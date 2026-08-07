#!/usr/bin/env python3
"""
mtft.gl2_peel — the GL(2) peel of f1 = 143a1: BSD rank from the skeleton.
=========================================================================
MIT License — (c) 2026 Roger Tano — MTFT Research Program

Curve (corpus, mtft_period_matrix_v4.gp l.34): y^2 + y = x^3 - x^2 - x - 2.
CONDUCTOR CERTIFICATE [EXACT]: Delta = -1859 = -11*13^2, c4 = 64 coprime to
11 and 13 => multiplicative at both => N = 143.

Derived ladder (X = 2 pi y, r = analytic rank, c_k = -L''(-k)/(2 L'(-k))):
  S_f(y) = sum Lambda_f(n) e^{-Xn}
         = -r/X - sum_{gam != 0} Gamma(1+i gam) X^{-(1+i gam)}
           + sum_{k>=0} [(-1)^k/k!] X^k [c_k - psi(k+1) + ln X]
NO X^-2 term (cuspidality; certified X^2 R ~ X).  EVERY k >= 0 is a double
pole (Gamma_C fingerprint: trivial zeros at all non-positive integers).

CERTIFIED (studies/peel_2026aug/gl2_*.py):
  root number eps = -1 (10.5 vs 0.7 digits, two-route s=4 test)
  L'(1) = 0.9456964112  -> analytic rank EXACTLY 1  (third independent route,
  after the period computation and the sign of the functional equation)
  RANK READ off the peel: 1.000000 +- 7e-7 at five depths.
"""
from __future__ import annotations
import mpmath as mp
from mpmath import mpf, mpc, log, exp, pi, sqrt, gamma as G, rgamma, psi, gammainc

CURVE = (0, -1, 1, -1, -2)          # a1,a2,a3,a4,a6
EPS = -1
ZEROS_143A1 = [3.2930459, 4.7576804, 5.8213585, 7.2166122, 8.5859573, 9.4210925]
CKS = [mpf('1.71298777169'), mpf('2.75175133678'), mpf('3.49566368467'),
       mpf('4.05741488219'), mpf('4.50165799232')]
LPRIME1 = mpf('0.9456964112')

def conductor_certificate():
    a1,a2,a3,a4,a6 = CURVE
    b2 = a1*a1+4*a2; b4 = 2*a4+a1*a3; b6 = a3*a3+4*a6
    b8 = a1*a1*a6+4*a2*a6-a1*a3*a4+a2*a3*a3-a4*a4
    Delta = -b2*b2*b8 - 8*b4**3 - 27*b6*b6 + 9*b2*b4*b6
    c4 = b2*b2 - 24*b4
    return dict(b2=b2,b4=b4,b6=b6,b8=b8,Delta=Delta,c4=c4,
                factored=(Delta == -11*13**2), mult_11=(c4 % 11 != 0),
                mult_13=(c4 % 13 != 0), conductor=143)

def ap_point_count(pmax):
    a1,a2,a3,a4,a6 = CURVE
    b2 = a1*a1+4*a2; b4 = 2*a4+a1*a3; b6 = a3*a3+4*a6
    sieve = list(range(pmax+1)); primes=[]
    for i in range(2, pmax+1):
        if sieve[i]==i:
            primes.append(i)
            for j in range(i*i, pmax+1, i):
                if sieve[j]==j: sieve[j]=i
    ap = {}
    cnt2 = sum(1 for x in range(2) for y in range(2)
               if (y*y+a1*x*y+a3*y-(x**3+a2*x*x+a4*x+a6)) % 2 == 0)
    ap[2] = 2 + 1 - (cnt2 + 1)
    for p in primes:
        if p == 2: continue
        sq = bytearray(p)
        for k in range(p//2+1): sq[(k*k) % p] = 1
        s = 0
        for x in range(p):
            g = (((4*x+b2)*x+2*b4)*x+b6) % p
            if g: s += 1 if sq[g] else -1
        ap[p] = -s
    return ap, primes

def lamf_sieve(pmax, ap=None, primes=None):
    if ap is None: ap, primes = ap_point_count(pmax)
    lamf = {}
    for p in primes:
        if p in (11,13):
            q,k = p,1
            while q <= pmax: lamf[q]=(ap[p]**k,p); q*=p; k+=1
        else:
            tp, tc = 2, ap[p]; q = p
            while q <= pmax:
                lamf[q]=(tc,p); tp,tc = tc, ap[p]*tc - p*tp; q*=p
    return lamf

def S_direct(y, pmax=None, lamf=None):
    X = 2*pi*y
    if pmax is None: pmax = int(75/float(X)) + 1
    if lamf is None: lamf = lamf_sieve(pmax)
    r = exp(-X); acc = mpf(1); s = mpf(0)
    for n in range(1, pmax+1):
        acc *= r
        if n in lamf:
            tk,p = lamf[n]
            s += mpf(tk)*log(mpf(p))*acc
    return s

def smooth(y, cks=CKS):
    X = 2*pi*y; L = log(X); out = mpf(0)
    for k in range(len(cks)):
        out += mpf(-1)**k/mp.factorial(k) * X**k * (cks[k] - psi(0,k+1) + L)
    return out

def zero_osc(y, zeros=ZEROS_143A1):
    X = 2*pi*y
    return -sum(2*mp.re(G(mpc(1,g))*X**(-mpc(1,g))) for g in zeros)

def rank_read(y, pmax=None, zeros=ZEROS_143A1, cks=CKS):
    """-X (S_f - smooth - zero oscillations) -> analytic rank."""
    X = 2*pi*y
    return -X*(S_direct(y, pmax) - smooth(y, cks) - zero_osc(y, zeros))
