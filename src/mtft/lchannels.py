#!/usr/bin/env python3
"""
mtft.lchannels — SU(p) gauge filter as Dirichlet L-function channels.
=====================================================================
MIT License — (c) 2026 Roger Tano — MTFT Research Program

SELECTION RULE [EXACT]: the holonomy filter 1 - cos(2 pi n m/N) is even in n,
so the filtered skeleton couples to exactly the EVEN Dirichlet characters
mod N = the characters of the split torus of PSL(2,N) = its principal-series
parameters (see mtft.coset_reps). Quadratic character enters iff N = 1 mod 4.

EXACT SPLIT (prime p, gcd(m,p)=1), verified 34.8 digits (p=5), 35.3 (p=13):
  S_{p,m}(y) = (p/(p-1)) (S - P_p)
             - (1/(p-1)) Re sum_{even chi != chi0} conj(chi(m)) tau(chi) S_conj(chi)
with S_chi(y) = sum n Lambda(n) chi(n) e^{-2 pi y n},  series -L'/L(s-1,chi).
The mode minimum min_m selects the sign/phase of the L-channels.

Even-channel smooth ladder (trivial zeros of even L at 0,-2,-4,...):
  -1/X + C0(chi) + X[c2+1-gamma-lnX] + (X^2/2)k2 + (X^3/6)[c4+11/6-gamma-lnX]+...
Per-channel: residual * X^{3/2} bounded as y->0  <=>  GRH(chi).

Session-recorded zeros (regenerate with scan_channel_zeros; tol ~1e-6):
  chi5:        6.6484533, 9.8314444, 11.958846, 16.033821, ...(11 located)
  chi6 mod13:  3.119341, 7.231591, 8.625427, 10.33642, 12.61701, 15.14833,
               16.27483, 18.75125, 19.54804, 20.95918, 23.59203, 25.3717
  sextic mod13 (j=2): -11.995,-11.1245,-9.04699,-5.42274,-3.66097,
                       4.45485, 6.80983, 7.99513, 11.5303, 13.0356
  cubic  mod13 (j=4): -12.3258,-10.5067,-8.53434,-5.99433,-2.27313,
                       4.93859, 6.72931, 9.41418, 10.4577
Veil |Gamma(3/2+it1)|: zeta 8.07e-9 | chi5 4.87e-4 | chi6 5.90e-2 |
sextic 2.95e-2 | cubic 1.64e-1  (x2.03e7 vs zeta).
"""
from __future__ import annotations
import mpmath as mp
from mpmath import mpf, mpc, log, exp, pi, sqrt, zeta, euler

SHIPPED_ZEROS = {
 (5, 'quadratic'): [6.6484533, 9.8314444, 11.958846, 16.033821],
 (13, 6): [3.119341, 7.231591, 8.625427, 10.33642, 12.61701, 15.14833,
           16.27483, 18.75125, 19.54804, 20.95918, 23.59203, 25.3717],
 (13, 2): [-11.995, -11.1245, -9.04699, -5.42274, -3.66097,
            4.45485, 6.80983, 7.99513, 11.5303, 13.0356],
 (13, 4): [-12.3258, -10.5067, -8.53434, -5.99433, -2.27313,
            4.93859, 6.72931, 9.41418, 10.4577],
}

def primitive_root(p):
    for g in range(2, p):
        seen, x = set(), 1
        for _ in range(p-1):
            x = x*g % p; seen.add(x)
        if len(seen) == p-1: return g
    raise ValueError(p)

def char_table(p):
    """dlog table and chi_j(n) for the cyclic group mod prime p."""
    g = primitive_root(p)
    dlog = {}; pw = 1
    for k in range(p-1):
        dlog[pw] = k; pw = pw*g % p
    def chi(j, n):
        n %= p
        if n == 0: return mpc(0)
        return mp.expjpi(mpf(2*j*dlog[n])/(p-1))
    return chi, dlog

def gauss_sum(p, j, chi=None):
    if chi is None: chi, _ = char_table(p)
    return sum(chi(j, r)*mp.expjpi(mpf(2*r)/p) for r in range(1, p))

def even_js(p):
    return [j for j in range(0, p-1, 2)]

def prime_tower(p, y):
    X = 2*pi*y; s = mpf(0); q = p
    while float(X*q) < 80:
        s += mpf(q)*log(mpf(p))*exp(-X*q); q *= p
    return s

def channel_sums(p, y, lam, js=None):
    """One pass: S, {S_chi_j}, and direct filtered S_{p,m} for m=1..(p-1)//2."""
    chi, _ = char_table(p)
    if js is None: js = [j for j in even_js(p) if j != 0][: (p-1)//2]
    X = 2*pi*y; M = min(len(lam)-1, int(70/float(X))+1)
    r = exp(-X); acc = mpf(1)
    filt = {m: [1 - mp.cos(2*pi*mpf(k*m)/p) for k in range(p)]
            for m in range(1, (p-1)//2 + 1)}
    S = mpf(0); Sx = {j: mpc(0) for j in js}
    Sf = {m: mpf(0) for m in filt}
    for n in range(1, M+1):
        acc *= r
        if lam[n] is not None:
            t = mpf(n)*lam[n]*acc; k = n % p
            S += t
            for j in js: Sx[j] += chi(j, k)*t
            for m in filt: Sf[m] += filt[m][k]*t
    return S, Sx, Sf

def split_formula(p, m, S, Sx, Ppow):
    """(p/(p-1))(S-P_p) - (1/(p-1)) Re sum_{even j!=0} conj(chi_j(m)) tau_j S_{conj chi_j},
    expanding conjugate partners with dedupe (quadratic j = p-1-j counted once)."""
    chi, _ = char_table(p)
    full = {}
    for j, v in Sx.items():
        full[j] = v
        jj = (p - 1 - j) % (p - 1)
        if jj not in full: full[jj] = mp.conj(v)
    acc = mpc(0)
    for j, v in full.items():
        if j == 0 or j % 2: continue
        acc += mp.conj(chi(j, m))*gauss_sum(p, j, chi)*mp.conj(v)
    return (mpf(p)/(p-1))*(S - Ppow) - (mpf(1)/(p-1))*mp.re(acc)

def L_chi(p, j, s, chi=None):
    if chi is None: chi, _ = char_table(p)
    return p**(-s)*sum(chi(j, r)*mp.zeta(s, mpf(r)/p) for r in range(1, p))

def even_channel_smooth(y, Lfunc):
    """Smooth ladder for an even primitive character channel."""
    Ld  = lambda w: mp.re(mp.diff(Lfunc, mpf(w)))
    Ldd = lambda w: mp.re(mp.diff(Lfunc, mpf(w), 2))
    Lr  = lambda w: mp.re(Lfunc(mpf(w)))
    C0 = -Ld(-1)/Lr(-1); c2 = Ldd(-2)/(2*Ld(-2)); k2 = -Ld(-3)/Lr(-3)
    c4 = Ldd(-4)/(2*Ld(-4)); k4 = -Ld(-5)/Lr(-5)
    X = 2*pi*y; L = log(X)
    return ( -1/X + C0 + X*(c2 + 1 - euler - L) + X**2/2*k2
             + X**3/6*(c4 + mpf(11)/6 - euler - L) + X**4/24*k4 )

def scan_channel_zeros(p, j, tlo, thi, step=0.06):
    """Zeros of L(s,chi_j mod p) on Re s = 1/2 via the rotated completed fn."""
    chi, _ = char_table(p)
    eps = gauss_sum(p, j, chi)/sqrt(mpf(p))
    rot = 1/mp.sqrt(eps)
    def Zf(t):
        s = mpc(mpf(1)/2, t)
        return mp.re(rot*(p/pi)**(s/2)*mp.gamma(s/2)*L_chi(p, j, s, chi))
    zs = []; t = mpf(tlo); f0 = Zf(t)
    while t < thi:
        t2 = t + mpf(step); f1 = Zf(t2)
        if f0*f1 < 0:
            lo, hi, flo = t, t2, f0
            for _ in range(45):
                mid = (lo+hi)/2; fm = Zf(mid)
                if flo*fm <= 0: hi = mid
                else: lo, flo = mid, fm
            zs.append((lo+hi)/2)
        t, f0 = t2, f1
    return zs
