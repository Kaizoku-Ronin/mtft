"""
SU(13) FILTERED SKELETON -- part A: exact character split + mode structure.

Even characters mod 13 (chi_j(-1) = (-1)^j, generator g=2):
  j=0 principal | j=6 quadratic (Legendre) | j=2,10 sextic pair | j=4,8 cubic pair

Derived split, gcd(m,13)=1:
  S_{13,m}(y) = (13/12)(S - P_13)
              - (1/12) Re SUM_{j in {2,4,6,8,10}} conj(chi_j(m)) tau_j S_{conj(chi_j)}(y)
with S_chi(y) = sum n Lambda(n) chi(n) e^{-2 pi y n},  S_{conj chi} = conj(S_chi).
"""
from mpmath import mp, mpf, mpc, log, exp, pi, sqrt, nstr, re, im, conj, diff, zeta
mp.dps = 35
N = 130000

spf = list(range(N+1)); i = 2
while i*i <= N:
    if spf[i] == i:
        for j in range(i*i, N+1, i):
            if spf[j] == j: spf[j] = i
    i += 1
logp = {}; Lam = [None]*(N+1)
for n in range(2, N+1):
    p = spf[n]; m0 = n
    while m0 % p == 0: m0 //= p
    if m0 == 1:
        if p not in logp: logp[p] = log(mpf(p))
        Lam[n] = logp[p]

# characters mod 13, generator 2
dlog = {}
pw = 1
for k in range(12):
    dlog[pw] = k; pw = (pw*2) % 13
def chi(j, r):
    r %= 13
    if r == 0: return mpc(0)
    return mp.expjpi(mpf(2*j*dlog[r])/12)   # e^{i pi * (2 j k /12)} = e^{2pi i jk/12}

tau = {}
for j in [2,4,6,8,10]:
    tau[j] = sum(chi(j,r)*mp.expjpi(mpf(2*r)/13) for r in range(1,13))
print("  Gauss sums |tau_j| (should all be sqrt(13) =", nstr(sqrt(mpf(13)),10)+"):")
for j in [2,4,6]:
    print(f"    |tau_{j}| = {nstr(abs(tau[j]),10)}   tau_{j} = {nstr(tau[j],8)}")
print(f"    tau_6/sqrt(13) = {nstr(tau[6]/sqrt(mpf(13)),8)}  (epsilon = +1 for the quadratic)")

def P13(y):
    X = 2*pi*y; s = mpf(0); pk = 13
    while float(X*pk) < 80:
        s += mpf(pk)*log(mpf(13))*exp(-X*pk); pk *= 13
    return s

def sums_at(y):
    """One pass: S, S_chi for j=2,4,6, and direct filtered S_{13,m}, m=1..6."""
    X = 2*pi*y
    M = min(N, int(70/float(X))+1)
    r = exp(-X); acc = mpf(1)
    filt = {m: {k: 1 - mp.cos(2*pi*mpf(k*m)/13) for k in range(13)} for m in range(1,7)}
    S = mpf(0); Sx = {2:mpc(0),4:mpc(0),6:mpc(0)}
    Sf = {m: mpf(0) for m in range(1,7)}
    for n in range(1, M+1):
        acc *= r
        if Lam[n] is not None:
            t = mpf(n)*Lam[n]*acc; k = n % 13
            S += t
            for j in (2,4,6): Sx[j] += chi(j,k)*t
            for m in range(1,7): Sf[m] += filt[m][k]*t
    return S, Sx, Sf

print("\n  EXACT SPLIT CHECK at y = 3e-4 (all six modes):")
y0 = mpf('3e-4')
S, Sx, Sf = sums_at(y0)
p13 = P13(y0)
Sx_full = {j: Sx[j] for j in (2,4,6)}
Sx_full[8]  = conj(Sx[4]); Sx_full[10] = conj(Sx[2])
for m in range(1,7):
    pred = (mpf(13)/12)*(S-p13) - (mpf(1)/12)*re(
        sum(conj(chi(j,m))*tau[j]*conj(Sx_full[j]) for j in (2,4,6,8,10)))
    dg = float(-mp.log10(abs(Sf[m]-pred)/abs(Sf[m])))
    print(f"    m={m}:  direct {nstr(Sf[m],10):>15}   split-formula agrees to {dg:.1f} digits")

print("\n  MODE COMPETITION -- which m minimizes (the vacuum's phase choice):")
for ystr in ['1e-2','3e-3','1e-3','3e-4','1e-4']:
    y = mpf(ystr)
    S, Sx, Sf = sums_at(y)
    vals = [(float(Sf[m]), m) for m in range(1,7)]
    vals.sort()
    spread = vals[-1][0]-vals[0][0]
    print(f"    y={ystr:>6}:  argmin m={vals[0][1]}   spread(max-min)={spread:.4g}   "
          f"S_chi6={nstr(re(Sx[6]),6)}")

# E2 on one complex channel: Dirichlet identity for the cubic chi_4
def L13(j, s):
    return 13**(-s)*sum(chi(j,r)*zeta(s, mpf(r)/13) for r in range(1,13))
for j in [6, 4]:
    s0 = 6
    direct = sum(mpf(n)*Lam[n]*chi(j,n%13)*mpf(n)**(-s0)
                 for n in range(2, N+1) if Lam[n] is not None)
    closed = -diff(lambda w: L13(j, w), mpf(s0-1))/L13(j, mpf(s0-1))
    dg = float(-mp.log10(abs(direct-closed)/abs(closed)))
    print(f"\n  E2 chi_{j}: sum n Lam chi n^-6 vs -L'/L(5,chi_{j}):  {dg:.1f} digits"
          f"   (complex channel)" if j!=6 else
          f"\n  E2 chi_{j}: sum n Lam chi n^-6 vs -L'/L(5,chi_{j}):  {dg:.1f} digits")
