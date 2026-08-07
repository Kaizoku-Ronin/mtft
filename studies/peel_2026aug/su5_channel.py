"""
FILTERED SKELETON: S_{N,m}(y) = sum n Lambda(n) e^{-2 pi y n} (1 - cos(2 pi n m / N))

Character decomposition (derived):
  filter is EVEN in n mod N  =>  only EVEN Dirichlet characters enter.
  N=3: quadratic char mod 3 is ODD  => pure zeta:
       S_3(y) = (3/2) [ S(y) - P_3(y) ],  P_3 = 3-power tower (Euler-factor deletion)
  N=5: quadratic char mod 5 is EVEN, tau(chi_5) = sqrt(5):
       S_{5,1}(y) = (5/4)(S - P_5) - (sqrt5/4) S_chi(y)
       S_{5,2}(y) = (5/4)(S - P_5) + (sqrt5/4) S_chi(y)
       with  S_chi(y) = sum n Lambda(n) chi_5(n) e^{-2 pi y n},
             Dirichlet series  -L'/L(s-1, chi_5).

chi_5 channel pole ladder (even primitive chi, trivial zeros of L at 0,-2,-4,...):
  s=2      NO pole  (L entire, L(1,chi_5) = (2/sqrt5) ln phi != 0)
  s=1+rho  zeros of L(s,chi_5)                -> the NEW oscillation family
  s=1      trivial zero of L at 0             -> -1/X   (the skipped order returns!)
  s=0      Gamma pole                          -> C0 = -L'(-1)/L(-1)
  s=-1     DOUBLE                              -> X [c2 + 1 - gamma - ln X]
  s=-2     Gamma pole                          -> (X^2/2)(-L'(-3)/L(-3))
  s=-3     DOUBLE                              -> (X^3/6)[c4 + 11/6 - gamma - ln X]
"""
from mpmath import (mp, mpf, mpc, log, exp, pi, sqrt, zeta, gamma, euler,
                    diff, nstr, re, im, arg, fabs)
import time

mp.dps = 35
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

chi5 = {0:0, 1:1, 2:-1, 3:-1, 4:1}   # Legendre symbol mod 5

def sums_at(y):
    """Direct S, S_chi, S_{5,1}, S_{5,2} in one pass."""
    X = 2*pi*y
    M = min(N, int(70/float(X))+1)
    r = exp(-X); acc = mpf(1)
    c1 = mpf(1) - mp.cos(2*pi*mpf(1)/5)   # filter values by residue, m=1
    c2_ = mpf(1) - mp.cos(2*pi*mpf(2)/5)
    fm1 = {0:mpf(0),1:c1,2:c2_,3:c2_,4:c1}
    fm2 = {0:mpf(0),1:c2_,2:c1,3:c1,4:c2_}   # m=2
    S = Sx = S51 = S52 = mpf(0)
    for n in range(1, M+1):
        acc *= r
        if Lam[n] is not None:
            t = mpf(n)*Lam[n]*acc
            S += t
            k = n % 5
            Sx += chi5[k]*t
            S51 += fm1[k]*t
            S52 += fm2[k]*t
    return S, Sx, S51, S52

def P5(y):
    X = 2*pi*y; s = mpf(0); pk = 5
    while X*pk < 80:
        s += mpf(pk)*log(mpf(5))*exp(-X*pk); pk *= 5
    return s

# ---- L(s, chi_5) via Hurwitz zeta ----
def L5(s):
    return 5**(-s)*( zeta(s, mpf(1)/5) - zeta(s, mpf(2)/5)
                   - zeta(s, mpf(3)/5) + zeta(s, mpf(4)/5) )

print("="*76)
print("  E2 -- chi channel Dirichlet identity, and the exact SU(5) split")
print("="*76)
for s0 in [6, 8]:
    direct = sum(mpf(n)*Lam[n]*chi5[n%5]*mpf(n)**(-s0)
                 for n in range(2, N+1) if Lam[n] is not None)
    closed = -diff(L5, mpf(s0-1))/L5(mpf(s0-1))
    d = float(abs(direct-closed)/abs(closed))
    print(f"  sum n Lam chi n^-{s0}  vs  -L'/L({s0-1},chi5):  {float(-mp.log10(d)):.1f} digits")

y0 = mpf('3e-4')
S, Sx, S51, S52 = sums_at(y0)
p5 = P5(y0)
pred1 = (mpf(5)/4)*(S-p5) - (sqrt(mpf(5))/4)*Sx
pred2 = (mpf(5)/4)*(S-p5) + (sqrt(mpf(5))/4)*Sx
print(f"\n  at y=3e-4:  S_51 direct vs (5/4)(S-P5)-(sqrt5/4)S_chi : "
      f"{float(-mp.log10(abs(S51-pred1)/abs(S51))):.1f} digits")
print(f"              S_52 direct vs (5/4)(S-P5)+(sqrt5/4)S_chi : "
      f"{float(-mp.log10(abs(S52-pred2)/abs(S52))):.1f} digits")
print(f"  min-mode at this y: m={'2' if S52<S51 else '1'} "
      f"(the mode that SUBTRACTS the |chi channel|; S_chi = {nstr(Sx,6)})")
print("  L(1,chi5) =", nstr(L5(mpf(1)),12), "  vs (2/sqrt5) ln phi =",
      nstr(2/sqrt(mpf(5))*log((1+sqrt(mpf(5)))/2),12),
      "\n  -> nonzero (class number formula): certifies NO X^-2 main term in the chi channel.")

# ---- zeros of L(s, chi_5) on the critical line ----
print("\n" + "="*76)
print("  ZEROS of L(s, chi_5) -- located from scratch (completed Lambda real on line)")
print("="*76)
def Zf(t):
    s = mpc(0.5, t)
    val = (5/pi)**(s/2)*gamma(s/2)*L5(s)
    return re(val)

t0 = time.time()
grid = [mpf(2) + mpf('0.05')*k for k in range(561)]   # t in [2, 30]
vals = [Zf(t) for t in grid]
zeros = []
for a, b, fa, fb in zip(grid, grid[1:], vals, vals[1:]):
    if fa*fb < 0:
        lo, hi, flo = a, b, fa
        for _ in range(60):
            mid = (lo+hi)/2; fm = Zf(mid)
            if flo*fm <= 0: hi = mid
            else: lo, flo = mid, fm
        zeros.append((lo+hi)/2)
print(f"  found {len(zeros)} zeros in t in (2,30)  ({time.time()-t0:.0f}s):")
print("  t_j =", ", ".join(nstr(t, 8) for t in zeros[:6]))
print("        " + ", ".join(nstr(t, 8) for t in zeros[6:]))

g1c = gamma(mpc(1.5, zeros[0]))
print(f"\n  veil comparison:  |Gamma(1+rho_1)| for L(chi5) t={nstr(zeros[0],6)}: {nstr(abs(g1c),6)}")
print(f"                    |Gamma(1+rho_1)| for zeta    t=14.135:  8.074e-9")
print(f"  ratio ~ {nstr(abs(g1c)/mpf('8.074e-9'),4)}  -> SU(5) oscillations ~6e4 x less veiled")

# ---- chi channel smooth part (all computed) ----
Ld  = lambda w: diff(L5, mpf(w))
Ldd = lambda w: diff(L5, mpf(w), 2)
C0c = -Ld(-1)/L5(mpf(-1))
c2c = Ldd(-2)/(2*Ld(-2))
k2c = -Ld(-3)/L5(mpf(-3))
c4c = Ldd(-4)/(2*Ld(-4))
k4c = -Ld(-5)/L5(mpf(-5))

def smooth_chi(y):
    X = 2*pi*y; L = log(X)
    return ( -1/X + C0c + X*(c2c + 1 - euler - L)
             + X**2/2*k2c + X**3/6*(c4c + mpf(11)/6 - euler - L) + X**4/24*k4c )

def Zsum_chi(y, k):
    X = 2*pi*y
    return -sum(2*re(gamma(mpc(1.5,t))*X**(-mpc(1.5,t))) for t in zeros[:k])

print("\n" + "="*76)
print("  PEEL AND COMPARE -- chi_5 channel residual vs its OWN zero sum")
print("="*76)
print(f"  {'y':>9} {'S_chi (direct)':>16} {'R = S_chi - smooth':>19} "
      f"{'Z (L-zeros)':>14} {'R - Z':>11}")
for ystr in ['5e-4','3e-4','2e-4','1.5e-4','1e-4']:
    y = mpf(ystr)
    _, Sx, _, _ = sums_at(y)
    R = Sx - smooth_chi(y)
    Z = Zsum_chi(y, len(zeros))
    print(f"  {ystr:>9} {nstr(Sx,9):>16} {nstr(R,9):>19} {nstr(Z,9):>14} {nstr(R-Z,3):>11}")

print("\n  relative visibility at y=1e-4:")
y = mpf('1e-4'); _, Sx, _, _ = sums_at(y)
R = Sx - smooth_chi(y)
print(f"    smooth ~ {nstr(smooth_chi(y),6)},  oscillation ~ {nstr(R,4)}  "
      f"({nstr(100*abs(R)/abs(smooth_chi(y)),3)}% of smooth -- vs 4e-8% in the zeta case)")
