"""
SU(13) part B: zeros of the even L-functions mod 13, peel-and-match on the
quadratic channel, veil comparison across the gauge tower.
"""
from mpmath import (mp, mpf, mpc, log, exp, pi, sqrt, zeta, gamma, euler,
                    diff, nstr, re, im, conj)
import time
mp.dps = 30
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

dlog = {}; pw = 1
for k in range(12):
    dlog[pw] = k; pw = (pw*2) % 13
def chi(j, r):
    r %= 13
    if r == 0: return mpc(0)
    return mp.expjpi(mpf(2*j*dlog[r])/12)

def L13(j, s):
    return 13**(-s)*sum(chi(j,r)*zeta(s, mpf(r)/13) for r in range(1,13))

tau = {j: sum(chi(j,r)*mp.expjpi(mpf(2*r)/13) for r in range(1,13)) for j in (2,4,6)}

# ---- zero location: rotated completed function, real on the line ----
def make_Z(j):
    eps = tau[j]/sqrt(mpf(13))
    rot = 1/mp.sqrt(eps)
    def Zf(t):
        s = mpc(mpf(1)/2, t)
        return re(rot*(13/pi)**(s/2)*gamma(s/2)*L13(j, s))
    # sanity: imaginary part small at a test point
    s = mpc(mpf(1)/2, mpf('1.7'))
    v = rot*(13/pi)**(s/2)*gamma(s/2)*L13(j, s)
    assert abs(im(v)) < mpf('1e-15')*max(abs(v), mpf(1)), "rotation branch wrong"
    return Zf

def scan(Zf, tlo, thi, step):
    zs = []
    t = mpf(tlo); f0 = Zf(t)
    while t < thi:
        t2 = t + step; f1 = Zf(t2)
        if f0*f1 < 0:
            lo, hi, flo = t, t2, f0
            for _ in range(50):
                mid = (lo+hi)/2; fm = Zf(mid)
                if flo*fm <= 0: hi = mid
                else: lo, flo = mid, fm
            zs.append((lo+hi)/2)
        t, f0 = t2, f1
    return zs

t0 = time.time()
Z6 = make_Z(6)
zeros6 = scan(Z6, mpf('0.2'), mpf(26), mpf('0.05'))
print(f"  quadratic chi_6 (real): {len(zeros6)} zeros in (0.2, 26)  ({time.time()-t0:.0f}s)")
print("    t =", ", ".join(nstr(t,7) for t in zeros6[:7]))
print("        ", ", ".join(nstr(t,7) for t in zeros6[7:]))

t0 = time.time()
Z2 = make_Z(2); Z4 = make_Z(4)
zeros2 = scan(Z2, mpf(-13), mpf(13), mpf('0.06'))
zeros4 = scan(Z4, mpf(-13), mpf(13), mpf('0.06'))
print(f"\n  sextic chi_2 (complex): zeros in (-13,13)  ({time.time()-t0:.0f}s)")
print("    t =", ", ".join(nstr(t,6) for t in zeros2))
print("  cubic  chi_4 (complex): zeros in (-13,13)")
print("    t =", ", ".join(nstr(t,6) for t in zeros4))

allfirst = [("zeta", mpf('14.134725')), ("chi_5 (SU5)", mpf('6.6484533')),
            ("chi_6 mod13", min(zeros6, key=abs)),
            ("chi_2 mod13", min(zeros2, key=abs)),
            ("chi_4 mod13", min(zeros4, key=abs))]
print("\n  LOWEST ZERO PER CHANNEL and the Gamma veil |Gamma(3/2+it)|:")
base = None
for name, t in allfirst:
    g = abs(gamma(mpc(mpf(3)/2, t)))
    if base is None: base = g
    print(f"    {name:<14} t1 = {nstr(abs(t),7):>9}   veil = {nstr(g,4):>10}   x{nstr(g/base,4)} vs zeta")

# ---- peel-and-match: quadratic channel ----
Lq  = lambda w: L13(6, w)
Ld  = lambda w: re(diff(Lq, mpf(w)))
Ldd = lambda w: re(diff(Lq, mpf(w), 2))
Lqr = lambda w: re(Lq(mpf(w)))
C0c = -Ld(-1)/Lqr(-1)
c2c = Ldd(-2)/(2*Ld(-2))
k2c = -Ld(-3)/Lqr(-3)
c4c = Ldd(-4)/(2*Ld(-4))
k4c = -Ld(-5)/Lqr(-5)

def smooth6(y):
    X = 2*pi*y; L = log(X)
    return ( -1/X + C0c + X*(c2c + 1 - euler - L)
             + X**2/2*k2c + X**3/6*(c4c + mpf(11)/6 - euler - L) + X**4/24*k4c )

def S6_direct(y):
    X = 2*pi*y; M = min(N, int(70/float(X))+1)
    r = exp(-X); acc = mpf(1); s = mpf(0)
    leg = {r_: re(chi(6,r_)) for r_ in range(13)}
    for n in range(1, M+1):
        acc *= r
        if Lam[n] is not None:
            s += leg[n % 13]*mpf(n)*Lam[n]*acc
    return s

def Zsum6(y):
    X = 2*pi*y
    return -sum(2*re(gamma(mpc(mpf(3)/2,t))*X**(-mpc(mpf(3)/2,t))) for t in zeros6)

print("\n  PEEL AND COMPARE -- quadratic mod-13 channel vs its own zeros:")
print(f"  {'y':>8} {'S_chi6 direct':>15} {'R = S - smooth':>16} {'Z (L-zeros)':>14} {'R - Z':>10}")
for ystr in ['5e-4','3e-4','2e-4','1e-4']:
    y = mpf(ystr)
    Sx = S6_direct(y)
    R  = Sx - smooth6(y)
    Z  = Zsum6(y)
    print(f"  {ystr:>8} {nstr(Sx,9):>15} {nstr(R,9):>16} {nstr(Z,9):>14} {nstr(R-Z,3):>10}")

y = mpf('1e-4')
R = S6_direct(y) - smooth6(y)
print(f"\n  visibility at y=1e-4: oscillation {nstr(abs(R),4)} against smooth "
      f"{nstr(abs(smooth6(y)),5)}  ({nstr(100*abs(R)/abs(smooth6(y)),3)}%)")
