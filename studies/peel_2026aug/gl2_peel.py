"""
GL(2) PEEL, part 2: zeros of L(E,s) for 143a1 located from scratch,
then the confrontation -- read the BSD rank off the skeleton.

Derived ladder (X = 2 pi y, r = analytic rank, c_k = -L''(-k)/2L'(-k)):
  S_f(y) = sum Lambda_f(n) e^{-Xn}
         = -r/X  -  sum_{gam != 0} Gamma(1+i gam) X^{-(1+i gam)}
           + sum_{k>=0} [(-1)^k/k!] X^k [ c_k - psi(k+1) + ln X ]
NO X^{-2} term (L entire: cuspidality). Every k >= 0 is a DOUBLE pole
(trivial zeros at ALL non-positive integers: the Gamma_C(s) fingerprint).
"""
from mpmath import (mp, mpf, mpc, log, exp, pi, sqrt, gamma, rgamma, psi,
                    nstr, re, im, diff, gammainc)
import json, time
mp.dps = 25

D = json.load(open("/home/claude/gl2_data.json"))
ap = {int(p): v for p, v in D["ap"].items()}
cks = [mpf(c) for c in D["cks"]]
primes = sorted(ap.keys())
PMAX = 10000

# rebuild Lambda_f sieve and a_n
LamF = {}
for p in primes:
    if p in (11,13):
        q,k = p,1
        while q <= PMAX: LamF[q]=(ap[p]**k,p); q*=p; k+=1
    else:
        tp, tc = 2, ap[p]; q = p
        while q <= PMAX:
            LamF[q]=(tc,p); tp,tc = tc, ap[p]*tc - p*tp; q*=p
spf = list(range(201))
for i in range(2,201):
    if spf[i]==i:
        for j in range(i*i,201,i):
            if spf[j]==j: spf[j]=i
def a_of(n):
    if n==1: return 1
    v=1; m=n
    while m>1:
        p=spf[m]; k=0
        while m%p==0: m//=p; k+=1
        if p in (11,13): v*=ap[p]**k
        else:
            up,uc=1,ap[p]
            for _ in range(k-1): up,uc=uc,ap[p]*uc-p*up
            v*=uc
    return v
an=[0]+[a_of(n) for n in range(1,201)]

SQN = sqrt(mpf(143)); TP = 2*pi; EPS=-1
def Lam_c(s):
    tot=mpc(0); n=1
    while True:
        x = TP*n/SQN
        if x>34: break
        q = SQN/(TP*n)
        tot += an[n]*( q**s*gammainc(s,x,mp.inf) + EPS*q**(2-s)*gammainc(2-s,x,mp.inf) )
        n+=1
    return tot

# --- zeros: eps=-1  =>  Lam(1+it) purely imaginary; Zf = Im ---
def Zf(t):
    return im(Lam_c(mpc(1,t)))
tv = Lam_c(mpc(1,mpf('1.3')))
print(f"  reality check on the line: |Re Lam(1+1.3i)| = {nstr(abs(re(tv)),3)}"
      f"   (Im = {nstr(im(tv),6)})")

t0=time.time()
zeros=[]; t=mpf('0.30'); f0=Zf(t); step=mpf('0.07')
while t < mpf('10.4'):
    t2=t+step; f1=Zf(t2)
    if f0*f1<0:
        lo,hi,flo=t,t2,f0
        for _ in range(40):
            mid=(lo+hi)/2; fm=Zf(mid)
            if flo*fm<=0: hi=mid
            else: lo,flo=mid,fm
        zeros.append((lo+hi)/2)
    t,f0=t2,f1
print(f"  zeros of L(143a1, s) on Re s = 1, gamma in (0.3, 10.4)  ({time.time()-t0:.0f}s):")
print("    gamma_j =", ", ".join(nstr(z,8) for z in zeros))

# --- smooth ladder and direct sums ---
def smooth(y):
    X=TP*y; L=log(X); s=mpf(0)
    sgn=1; fact=1
    for k in range(5):
        if k>0: fact*=k; 
        sgn = (-1)**k
        s += mpf(sgn)/fact * X**k * (cks[k] - psi(0,k+1) + L)
    return s
def S_direct(y):
    X=TP*y; M=min(PMAX,int(75/float(X))+1)
    rr=exp(-X); acc=mpf(1); s=mpf(0)
    for n in range(1,M+1):
        acc*=rr
        if n in LamF:
            tk,p = LamF[n]
            s += mpf(tk)*log(mpf(p))*acc
    return s
def Zosc(y):
    X=TP*y
    return -sum(2*re(gamma(mpc(1,g))*X**(-mpc(1,g))) for g in zeros)

print("\n  THE CONFRONTATION -- reading BSD rank off the prime sum")
print(f"  {'y':>8} {'S_f (direct)':>14} {'-X (S-smooth)':>14} {'1+X|osc| pred':>14}"
      f" {'RANK READ':>11} {'X^2 R (no pole)':>15}")
for ystr in ['8e-3','5e-3','3e-3','2e-3','1.2e-3']:
    y=mpf(ystr); X=TP*y
    S=S_direct(y); R=S-smooth(y)
    zo=Zosc(y)
    plateau = -X*R
    pred    = 1 - X*zo
    rank_read = -X*(R - zo)
    print(f"  {ystr:>8} {nstr(S,7):>14} {nstr(plateau,8):>14} {nstr(pred,8):>14}"
          f" {nstr(rank_read,8):>11} {nstr(X*X*R,4):>15}")
print("\n  RANK READ column: -X (S_f - smooth - zero-oscillations) -> r = 1")
print("  X^2 R column: shrinking ~X, certifying NO X^-2 pole (cuspidality).")
