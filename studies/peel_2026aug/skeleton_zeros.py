"""
SKELETON STIFFNESS: S(y) = sum_n n*Lambda(n) e^{-2 pi y n}
Dirichlet series: F_L(s) = -zeta'(s-1)/zeta(s-1)   [E2-checked below]

Mellin pole ladder (X = 2 pi y, L = ln X):
  s=2      simple, res +1 of -z'/z at w=1     ->  1/X^2          (PNT term)
  s=1+rho  one per nontrivial zero            ->  -Gamma(1+rho) X^{-(1+rho)}
  s=1      NOTHING (zeta(0) = -1/2 != 0)      ->  skipped order
  s=0      Gamma pole, F_L(0)=-z'(-1)/z(-1)   ->  1 - 12 ln A    (GLAISHER)
  s=-1     DOUBLE (Gamma pole + trivial zero) ->  X[c2 + H1 - gamma - L]
  s=-2     simple Gamma pole                  ->  (X^2/2)(-z'(-3)/z(-3))
  s=-3     DOUBLE                             ->  (X^3/6)[c4 + H3 - gamma - L]
  s=-4     simple                             ->  (X^4/24)(-z'(-5)/z(-5))
with c2 = z''(-2)/(2 z'(-2)), c4 = z''(-4)/(2 z'(-4)), H1=1, H3=11/6.

Then: residual after smooth part  ==  sum over zeros, coefficient Gamma(1+rho).
"""
from mpmath import (mp, mpf, mpc, log, exp, pi, zeta, gamma, euler, diff,
                    zetazero, nstr, re, glaisher)
import time

mp.dps = 40
N = 130000

# ---- von Mangoldt sieve ----
spf = list(range(N+1))
i = 2
while i*i <= N:
    if spf[i] == i:
        for j in range(i*i, N+1, i):
            if spf[j] == j: spf[j] = i
    i += 1
logp = {}
Lam = [None]*(N+1)
for n in range(2, N+1):
    p = spf[n]; m = n
    while m % p == 0: m //= p
    if m == 1:
        if p not in logp: logp[p] = log(mpf(p))
        Lam[n] = logp[p]

def S_direct(y):
    X = 2*pi*y
    M = min(N, int(70/float(X))+1)
    r = exp(-X); acc = mpf(1); s = mpf(0)
    for n in range(1, M+1):
        acc *= r
        if Lam[n] is not None:
            s += mpf(n)*Lam[n]*acc
    return s

# ---- E2: the Dirichlet identity itself, independent route ----
print("="*76)
print("  E2 -- identity check: sum n Lambda(n) n^-s  vs  -zeta'(s-1)/zeta(s-1)")
print("="*76)
for s0 in [6, 8]:
    direct = sum(mpf(n)*Lam[n]*mpf(n)**(-s0) for n in range(2, N+1) if Lam[n] is not None)
    closed = -diff(zeta, mpf(s0-1))/zeta(mpf(s0-1))
    d = float(abs(direct-closed)/abs(closed))
    print(f"  s={s0}:  agree to {float(-mp.log10(d)):.1f} digits")

# ---- smooth-part constants (all computed, none fitted) ----
zp  = lambda w: diff(zeta, mpf(w))
zpp = lambda w: diff(zeta, mpf(w), 2)
C0  = -zp(-1)/zeta(-1)
c2  = zpp(-2)/(2*zp(-2))
c4  = zpp(-4)/(2*zp(-4))
k2  = -zp(-3)/zeta(-3)
k4  = -zp(-5)/zeta(-5)

print("\n  Constant term = -zeta'(-1)/zeta(-1) =", nstr(C0, 20))
print("  1 - 12 ln(Glaisher A)              =", nstr(1-12*log(glaisher), 20))
print("  -> identical: the Glaisher-Kinkelin constant IS the skeleton's constant term.")

def smooth(y):
    X = 2*pi*y; L = log(X)
    return ( X**-2 + C0
             + X*(c2 + 1 - euler - L)
             + (X**2/2)*k2
             + (X**3/6)*(c4 + mpf(11)/6 - euler - L)
             + (X**4/24)*k4 )

# ---- zeros ----
t0 = time.time()
K = 30
rhos = [zetazero(j) for j in range(1, K+1)]
print(f"\n  first {K} nontrivial zeros computed ({time.time()-t0:.1f}s)")
g1 = gamma(1+rhos[0])
print("  |Gamma(1+rho_1)| =", nstr(abs(g1), 8),
      "  (the e^{-pi t/2} veil: zeros are present but exponentially damped)")

def Zsum(y, k):
    X = 2*pi*y
    return -sum(2*re(gamma(1+r)*X**(-(1+r))) for r in rhos[:k])

# ---- the confrontation ----
print("\n" + "="*76)
print("  PEEL AND COMPARE -- residual vs zero sum (nothing fitted)")
print("="*76)
print(f"  {'y':>9} {'R = S - smooth':>16} {'Z_1':>13} {'Z_5':>13} {'Z_30':>13} {'R - Z_30':>12}")
for ystr in ['5e-4','3e-4','2e-4','1.5e-4','1e-4']:
    y = mpf(ystr)
    R  = S_direct(y) - smooth(y)
    z1, z5, z30 = Zsum(y,1), Zsum(y,5), Zsum(y,30)
    print(f"  {ystr:>9} {nstr(R,9):>16} {nstr(z1,7):>13} {nstr(z5,7):>13} "
          f"{nstr(z30,7):>13} {nstr(R-z30,4):>12}")

print("\n  OSCILLATION -- normalized residual R*X^{3/2} changes sign with ln y,")
print("  exactly as -2|G1|cos(t1 ln X - arg G1) demands:")
print(f"  {'y':>10} {'R * X^1.5':>15} {'Z_30 * X^1.5':>15}")
import numpy as _np
for yy in _np.geomspace(1e-4, 6e-4, 10):
    y = mpf(float(yy)); X = 2*pi*y
    R = S_direct(y) - smooth(y)
    print(f"  {float(y):>10.3e} {nstr(R*X**mpf(1.5),7):>15} {nstr(Zsum(y,30)*X**mpf(1.5),7):>15}")
