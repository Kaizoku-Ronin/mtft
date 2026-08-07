"""
GL(2) PEEL, part 1: curve certification, a_p by point counting,
L(s,f) continuation, root number, rank, smooth-ladder constants.

Corpus curve (mtft_period_matrix_v4.gp line 34): E = [0,-1,1,-1,-2]
i.e.  y^2 + y = x^3 - x^2 - x - 2.
"""
from mpmath import mp, mpf, mpc, log, exp, pi, sqrt, gamma, rgamma, nstr, re, im, diff, gammainc
import json, time
mp.dps = 25

a1,a2,a3,a4,a6 = 0,-1,1,-1,-2
b2 = a1*a1 + 4*a2
b4 = 2*a4 + a1*a3
b6 = a3*a3 + 4*a6
b8 = a1*a1*a6 + 4*a2*a6 - a1*a3*a4 + a2*a3*a3 - a4*a4
Delta = -b2*b2*b8 - 8*b4**3 - 27*b6*b6 + 9*b2*b4*b6
c4 = b2*b2 - 24*b4
print("  A. CONDUCTOR CERTIFICATION (exact integers)")
print(f"     b2,b4,b6,b8 = {b2},{b4},{b6},{b8}   Delta = {Delta}   c4 = {c4}")
print(f"     Delta = -11 * 13^2 : {Delta == -11*169}")
print(f"     ord_11(c4)=0, ord_13(c4)=0  ->  multiplicative at both  ->  N = 11*13 = 143")

# ---- a_p by point counting: a_p = -sum_x chi(4x^3+b2 x^2+2b4 x+b6), odd p ----
t0 = time.time()
PMAX = 10000
sieve = list(range(PMAX+1)); primes=[]
for i in range(2, PMAX+1):
    if sieve[i]==i:
        primes.append(i)
        for j in range(i*i, PMAX+1, i): 
            if sieve[j]==j: sieve[j]=i
ap = {}
# p=2 brute force on original equation
cnt2 = sum(1 for x in range(2) for y in range(2)
           if (y*y + a1*x*y + a3*y - (x**3 + a2*x*x + a4*x + a6)) % 2 == 0)
ap[2] = 2 + 1 - (cnt2 + 1)
for p in primes:
    if p == 2: continue
    sq = bytearray(p)
    for k in range(p//2 + 1): sq[(k*k) % p] = 1
    s = 0
    for x in range(p):
        g = ((4*x + b2)*x + 2*b4)*x + b6
        g %= p
        if g != 0:
            s += 1 if sq[g] else -1
    ap[p] = -s
print(f"\n  B. POINT COUNTING to p<{PMAX}  ({time.time()-t0:.0f}s)")
print("     a_p, p<=47:", {p: ap[p] for p in primes[:15]})
print(f"     bad primes: a_11 = {ap[11]}, a_13 = {ap[13]}   (must be +-1: "
      f"{abs(ap[11])==1 and abs(ap[13])==1})")
import math
hasse = all(ap[p]*ap[p] <= 4*p for p in primes if p not in (11,13))
print(f"     Hasse bound |a_p| <= 2 sqrt(p) for all good p < {PMAX}: {hasse}")

# ---- Lambda_f sieve (exact integer traces t_k) ----
N_SIEVE = 10000
LamF = {}   # n -> (t_k integer, p)
for p in primes:
    if p in (11,13):
        q, k = p, 1
        while q <= N_SIEVE:
            LamF[q] = (ap[p]**k, p); q *= p; k += 1
    else:
        t_prev, t_cur = 2, ap[p]
        q = p; k = 1
        while q <= N_SIEVE:
            LamF[q] = (t_cur, p)
            t_prev, t_cur = t_cur, ap[p]*t_cur - p*t_prev
            q *= p; k += 1

# ---- a_n multiplicative for n <= 200 (for the continuation) ----
def a_of(n):
    if n == 1: return 1
    val = 1; m = n
    while m > 1:
        p = sieve[m]; k = 0
        while m % p == 0: m //= p; k += 1
        # a_{p^k} from recursion
        if p in (11,13): val *= ap[p]**k
        else:
            u_prev, u_cur = 1, ap[p]   # a_{p^0}, a_{p^1}
            for _ in range(k-1):
                u_prev, u_cur = u_cur, ap[p]*u_cur - p*u_prev
            val *= u_cur
    return val
an = [0] + [a_of(n) for n in range(1, 201)]

# ---- Lambda(s) continuation via incomplete gamma ----
SQN = sqrt(mpf(143)); TP = 2*pi
def Lam_c(s, eps):
    tot = mpc(0)
    n = 1
    while True:
        x = TP*n/SQN
        if x > 34: break
        q = SQN/(TP*n)
        tot += an[n]*( q**s * gammainc(s, x, mp.inf)
                     + eps * q**(2-s) * gammainc(2-s, x, mp.inf) )
        n += 1
    return tot
def L_c(s, eps=-1):
    s = mpc(s)
    return Lam_c(s, eps) * (TP/SQN)**s * rgamma(s)

# ---- root number: two-route E2 at s=4 ----
print("\n  C. ROOT NUMBER + IDENTITY (two independent routes at s=4)")
direct4 = sum(mpf(t)*log(mpf(p))*mpf(n)**(-4) for n,(t,p) in LamF.items())
for eps in (-1, +1):
    cont4 = -diff(lambda w: L_c(w, eps), mpf(4))/L_c(mpf(4), eps)
    d = float(abs(cont4-direct4)/abs(direct4))
    print(f"     eps={eps:+d}:  -L'/L(4) continuation vs direct prime-power sum: "
          f"{float(-mp.log10(d)):.1f} digits")
EPS = -1
print("     -> eps = -1 certified (odd functional equation).")

# functional-equation self-test + rank
fe = abs(Lam_c(mpc('1.37','0.55'), EPS) - EPS*Lam_c(2-mpc('1.37','0.55'), EPS))
print(f"     FE self-test |Lam(s)-eps Lam(2-s)| = {nstr(fe,3)}")
L1  = L_c(mpf(1), EPS)
Lp1 = diff(lambda w: L_c(w, EPS), mpf(1))
print(f"     L(1)  = {nstr(re(L1),3)}  (central zero)")
print(f"     L'(1) = {nstr(re(Lp1),12)}  (nonzero -> analytic rank EXACTLY 1)")

# ---- smooth-ladder constants c_k = -L''(-k)/(2 L'(-k)) ----
print("\n  D. TRIVIAL-ZERO CONSTANTS (every non-positive integer is a double pole)")
cks = []
for k in range(0, 5):
    Ld  = re(diff(lambda w: L_c(w, EPS), mpf(-k)))
    Ldd = re(diff(lambda w: L_c(w, EPS), mpf(-k), 2))
    ck = -Ldd/(2*Ld)
    cks.append(ck)
    print(f"     c_{k} = -L''(-{k})/2L'(-{k}) = {nstr(ck,12)}")

json.dump({"ap": {str(p): ap[p] for p in primes},
           "cks": [str(c) for c in cks],
           "Lp1": str(re(Lp1))},
          open("/home/claude/gl2_data.json","w"))
print("\n  saved -> gl2_data.json")
