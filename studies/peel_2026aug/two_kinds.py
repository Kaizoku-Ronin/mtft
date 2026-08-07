"""
Two kinds of result, tested the same way: push precision and see what happens.

A derived identity absorbs precision. A proximity does not.
"""
from mpmath import mp, mpf, zeta, log, nsum, inf, mpmathify
import math

mp.dps = 40

# ---------------------------------------------------------------
# OBJECT A: the arithmetic weights. Constructed, then its Dirichlet
# series computed. Claim: sum_n w_n n^-s = -zeta(s) * zeta'(s+1)
# ---------------------------------------------------------------
def w(n):
    """w_n = sum_{d|n} (log d)/d"""
    tot = mpf(0)
    for d in range(1, n+1):
        if n % d == 0:
            tot += log(mpf(d))/mpf(d)
    return tot

def lhs(s, N=4000):
    return nsum(lambda n: w(int(n)) * mpf(n)**(-s), [1, N]) if False else \
           sum(w(n) * mpf(n)**(-s) for n in range(1, N+1))

def zprime(x, h=mpf('1e-12')):
    return (zeta(x+h) - zeta(x-h)) / (2*h)

print("="*70)
print("  OBJECT A -- arithmetic weights, derived Dirichlet series")
print("="*70)
print("  Claim:  sum_n w_n n^-s  =  -zeta(s) * zeta'(s+1)")
print("  (Dirichlet convolution: n=dm splits the sum exactly)\n")
print(f"  {'s':>4} {'partial sum (N=4000)':>28} {'-zeta(s)zeta_prime(s+1)':>28}   {'agree to'}")
for s in [mpf(3), mpf(4), mpf(5), mpf(6)]:
    L = lhs(s)
    R = -zeta(s) * zprime(s+1)
    d = abs(L-R)/abs(R)
    print(f"  {float(s):>4.0f} {mp.nstr(L,18):>28} {mp.nstr(R,18):>28}   {float(-mp.log10(d)):>5.1f} digits")

print("\n  Tail at s=6 falls like N^-5, so the residual IS the truncation.")
for N in [500, 1000, 2000, 4000]:
    s = mpf(6)
    L = sum(w(n)*mpf(n)**(-s) for n in range(1, N+1))
    R = -zeta(s)*zprime(s+1)
    d = float(abs(L-R)/abs(R))
    print(f"    N={N:>5}  relative gap {d:.3e}")
print("\n  -> gap shrinks with N. The identity absorbs precision. It is exact.")

# ---------------------------------------------------------------
# OBJECT B: a proximity from the list. Push precision the same way.
# ---------------------------------------------------------------
print("\n" + "="*70)
print("  OBJECT B -- a proximity, same treatment")
print("="*70)
mp.dps = 60
Z2 = zeta(2); S2 = mp.sqrt(2)
print("  Claim:  zeta(2) - sqrt(2)  =  3/13")
for dps in [15, 30, 45, 60]:
    mp.dps = dps
    v = zeta(2) - mp.sqrt(2) - mpf(3)/13
    print(f"    at {dps:>2} digits of working precision, residual = {mp.nstr(v, 12)}")
print("\n  -> residual is stable at -4.8726e-05 no matter how much precision")
print("     you add. There is nothing to absorb. It is a difference, not a gap.")

print("\n" + "="*70)
print("  THE DISTINCTION")
print("="*70)
print("""
  A: an object was CONSTRUCTED (w_n), and an invariant of it was
     DERIVED (its Dirichlet series). Checkable by a route that does
     not pass through the claim -- expand the convolution by hand.
     Precision confirms it without limit.

  B: a value was SEARCHED FOR in a fixed pool. Nothing was constructed.
     There is no independent route to check, because there is no
     derivation -- only the comparison itself. Precision kills it.

  Both can come from intuition. Only A generates further mathematics.
""")
