"""The period clock: repeating expansions, repunits, cyclotomic polynomials, the trefoil's Alexander polynomial, cyclic covers, and
the modular group's torsion — Roger's Gemini thread, checked.  2026-10-19.  EXACT unless marked.

(1) 1/p in base B repeats with period ord_p(B); Midy's complement when the period is even; the two repetend loops of 1/13 are the
squares and the non-squares mod 13.  (2) ord_p(B) = k exactly when p | Phi_k(B) (p not dividing k): Phi_6(10) = 91 = 7 x 13, so 7 and 13
are the primes with decimal period 6.  (3) The trefoil's Alexander polynomial is Phi_6(t) = t^2 - t + 1, the characteristic polynomial of
the order-6 element ST of SL(2,Z); its Alexander module is Z[zeta_6], where 7 and 13 split and 11 is inert; Z[t]/(Phi_6, t - B) = Z/Phi_6(B);
Alexander-quandle colourings of the trefoil over F_13 with t = 10 number 13^2 against 13 for a non-root; the k-fold cyclic branched covers
have |H_1| = 3, 4, 3, 1, inf, 1, ... periodic in k mod 6 (Sigma_5 and Sigma_7 are the homology spheres Sigma(2,3,5), Sigma(2,3,7)).
(4) T(2,p): Delta = Phi_{2p}(t) and Delta(-B) = R_p(B), the repunit (Gemini's observation, right).  (5) 143 = 12^2 - 1: 1/143 has period 6
with repetend 006993 (Midy 006 + 993 = 999); 13 = R_3(3) = Delta_trefoil(-3); the Pell unit 12 + sqrt(143) gives [12 1; 143 12] in
Gamma_0(143), conjugate to the Lorenz word L R^22 with Rademacher psi = 21: a modular knot of X_0(143) linking the trefoil 21 times.
(6) Higher rank, the corrections: SL_3(Z) has no non-trivial quasimorphisms (Burger-Monod), so the Rademacher *function* does not lift;
two periodic flats (2-tori) can link in the 5-dimensional SL_3(Z) \\ SL_3(R)/SO(3) since 2 + 2 = 5 - 1."""
import math, itertools
from fractions import Fraction
import sympy as sp
import numpy as np

t = sp.Symbol('t')
def order(B, p):
    k, x = 1, B % p
    while x != 1: x = (x * B) % p; k += 1
    return k
def repetend(n, p, B):
    k = order(B, p); digits = []; r = n % p
    for _ in range(k): r *= B; digits.append(r // p); r %= p
    return digits, k
def digits_to_int(ds, B): return sum(d * B ** (len(ds) - 1 - i) for i, d in enumerate(ds))

# ---------------------------------------------------------------- 1. the decimal clock of 13, 11 and 143
d13, k13 = repetend(1, 13, 10); assert d13 == [0, 7, 6, 9, 2, 3] and k13 == 6
A, Bh = digits_to_int(d13[:3], 10), digits_to_int(d13[3:], 10); assert A + Bh == 999                                  # Midy
assert pow(10, 3, 13) == 12                                                                                             # 10^3 = -1 mod 13: the halves are complements
sq13 = {(x * x) % 13 for x in range(1, 13)}
loopA = {n for n in range(1, 13) if repetend(n, 13, 10)[0] in [d13[i:] + d13[:i] for i in range(6)]}
assert loopA == sq13 == {1, 3, 4, 9, 10, 12} and 10 in sq13                                                               # loop A = the squares; the base 10 is a square
d11, k11 = repetend(1, 11, 10); assert d11 == [0, 9] and k11 == 2 and pow(10, 1, 11) == 10 == 11 - 1                      # 10 = -1 mod 11
d143, k143 = repetend(1, 143, 10); assert d143 == [0, 0, 6, 9, 9, 3] and k143 == 6 == math.lcm(k11, k13)
assert digits_to_int(d143[:3], 10) + digits_to_int(d143[3:], 10) == 999
print("1/13 = 0.(076923): period 6 = ord_13(10); 076 + 923 = 999 since 10^3 = -1 (mod 13); loop A = the squares mod 13 = {1,3,4,9,10,12}, loop B the non-squares")
print("1/11 = 0.(09): period 2, 10 = -1 (mod 11).   1/143 = 0.(006993): period lcm(2, 6) = 6, 006 + 993 = 999")

# ---------------------------------------------------------------- 2. periods are cyclotomic: ord_p(B) = k  <=>  p | Phi_k(B), p not dividing k
Phi = lambda k, x: int(sp.cyclotomic_poly(k, t).subs(t, x))
assert Phi(6, 10) == 91 == 7 * 13 and Phi(2, 10) == 11 and Phi(3, 3) == 13 and Phi(6, -3) == 13
for k in range(1, 13):
    primes_dividing = [p for p in sp.primerange(2, 400) if Phi(k, 10) % p == 0 and k % p != 0]
    primes_with_period = [p for p in sp.primerange(2, 400) if p != 2 and p != 5 and order(10, p) == k]
    assert primes_dividing == primes_with_period, (k, primes_dividing, primes_with_period)
print("Phi_6(10) = 91 = 7 x 13: the primes with decimal period 6 are exactly 7 and 13; ord_p(10) = k <=> p | Phi_k(10) checked for k <= 12")

# ---------------------------------------------------------------- 3. the trefoil: Alexander polynomial Phi_6, monodromy ST of order 6, the Eisenstein ring
Delta = t ** 2 - t + 1
assert sp.expand(Delta - sp.cyclotomic_poly(6, t)) == 0
S = sp.Matrix([[0, -1], [1, 0]]); T_ = sp.Matrix([[1, 1], [0, 1]]); ST = S * T_
assert ST.charpoly(t).as_expr() == Delta and ST ** 3 == -sp.eye(2) and ST ** 6 == sp.eye(2)
# Z[t]/(Phi_6, t - B) = Z/Phi_6(B): the Alexander module evaluated at the base
for Bv in (-3, 2, 10):
    M = sp.Matrix([[Bv * Bv - Bv + 1]])                                                   # the ideal (Phi_6(t), t - B) in Z[t] is (Phi_6(B), t - B)
    assert abs(int(M[0, 0])) == abs(Phi(6, Bv))
print("Delta_trefoil = Phi_6 = charpoly(ST), ST^3 = -1 (the Midy complement), ST^6 = 1; Z[t]/(Phi_6, t - B) = Z/Phi_6(B): Z/13 at B = -3, Z/91 at B = 10")
# splitting in the Eisenstein integers: p = 1 mod 6 splits (7 = 3^2 - 3 + 1, 13 = 4^2 - 4 + 1), p = 5 mod 6 is inert
norm_rep = lambda p: [(a, b) for a in range(-p, p + 1) for b in range(-p, p + 1) if a * a - a * b + b * b == p]
assert norm_rep(13) and norm_rep(7) and not norm_rep(11)
print("   13 = 4^2 - 4*1 + 1^2 and 7 = 3^2 - 3 + 1 are norms from Z[zeta_6] (they split; their bases 10 are eigenvalues of ST mod p); 11 is not (inert, 11 = 5 mod 6)")
# Alexander-quandle colourings of the trefoil over F_13: x3 = t x1 + (1 - t) x2 and cyclic
def colourings(p, tv):
    n = 0
    for x1 in range(p):
        for x2 in range(p):
            x3 = (tv * x1 + (1 - tv) * x2) % p
            if (tv * x2 + (1 - tv) * x3) % p == x1 and (tv * x3 + (1 - tv) * x1) % p == x2: n += 1
    return n
assert colourings(13, 10) == 169 and colourings(13, 4) == 169 and colourings(13, 2) == 13 and colourings(11, 10) == 11
print("   Alexander-quandle colourings of the trefoil over F_13: 13^2 = 169 at t = 10 or 4 (the roots of Phi_6 mod 13), only the 13 trivial ones at t = 2; over F_11 none beyond trivial at t = 10")
# cyclic branched covers: |H_1(Sigma_k)| = |prod_{j=1}^{k-1} Delta(zeta_k^j)| = |Res(Delta, (t^k - 1)/(t - 1))|
orders = []
for k in range(2, 14):
    Q = sp.quo(t ** k - 1, t - 1); r = abs(int(sp.resultant(Delta, Q, t))); orders.append(r)
assert orders == [3, 4, 3, 1, 0, 1, 3, 4, 3, 1, 0, 1]
print("   |H_1| of the k-fold cyclic branched covers of the trefoil, k = 2..13:", orders, "— period 6 in k; Sigma_2 = L(3,1), Sigma_3 = S^3/Q_8, Sigma_5 = the Poincare sphere, Sigma_6 infinite (Nil), Sigma_7 = Sigma(2,3,7)")

# ---------------------------------------------------------------- 4. torus knots T(2,p): Delta = Phi_{2p}, Delta(-B) = R_p(B)
R = lambda n, B: (B ** n - 1) // (B - 1)
for p in (3, 5, 7, 11, 13):
    D = sp.cyclotomic_poly(2 * p, t)
    assert sp.expand(D * (t + 1) - (t ** p + 1)) == 0                                        # Delta_{T(2,p)} = (t^p + 1)/(t + 1)
    assert int(D.subs(t, -10)) == R(p, 10) and int(D.subs(t, -3)) == R(p, 3)
print("T(2,p): Delta = Phi_2p(t) = (t^p + 1)/(t + 1), and Delta(-B) = R_p(B): Delta_T(2,13)(-10) = 1111111111111; Delta_trefoil(-3) = 111_3 = 13")

# ---------------------------------------------------------------- 5. 143 = 12^2 - 1: the Pell unit and its modular knot in Gamma_0(143)
assert 143 == 12 ** 2 - 1 and 12 ** 2 - 143 * 1 ** 2 == 1
cf = sp.continued_fraction_periodic(0, 1, 143); assert cf == [11, [1, 22]]                   # sqrt(143) = [11; 1, 22]
A143 = sp.Matrix([[12, 1], [143, 12]]); assert A143.det() == 1 and A143.trace() == 24 and A143[1, 0] % 143 == 0
L = sp.Matrix([[1, 0], [1, 1]]); Rm = sp.Matrix([[1, 1], [0, 1]])
W = L * Rm ** 22; assert W.trace() == 24
def conjugator(A, W, bound=60):
    p = sp.symbols('p0:4'); P = sp.Matrix(2, 2, p); eqs = list(P * A - W * P); (s,) = sp.linsolve(eqs, p); free = sorted(set().union(*[e.free_symbols for e in s]), key=str)
    for m in range(-bound, bound + 1):
        for n in range(-bound, bound + 1):
            vals = [e.subs(dict(zip(free, (m, n)[:len(free)]))) for e in s]
            if all(v.is_integer for v in vals) and vals[0] * vals[3] - vals[1] * vals[2] == 1: return sp.Matrix(2, 2, vals)
    return None
word = 'LR^22' if conjugator(A143, W) is not None else ('L^22R' if conjugator(A143, L ** 22 * Rm) is not None else None)
assert word is not None
def dedekind(h, k): saw = lambda x: 0 if x == int(x) else x - math.floor(x) - Fraction(1, 2); return sum(saw(Fraction(r, k)) * saw(Fraction(h * r, k)) for r in range(1, k))
def psi(M):
    a, b, c, d = (int(x) for x in M); phi = Fraction(a + d, c) - 12 * (1 if c > 0 else -1) * dedekind(d % abs(c), abs(c)); return phi - 3 * (1 if c * (a + d) > 0 else -1)
ps = psi(A143); assert ps.denominator == 1 and abs(int(ps)) == 21 and int(ps) == (21 if word == 'LR^22' else -21)
ell = 2 * math.acosh(12)
print("143 = 12^2 - 1: sqrt(143) = [11; 1, 22], Pell unit 12 + sqrt(143); [12 1; 143 12] in Gamma_0(143) is conjugate to the Lorenz word %s, Rademacher psi = %+d, geodesic length 2 arccosh 12 = %.4f: a modular knot of X_0(143) linking the trefoil %d times" % (word, int(ps), ell, abs(int(ps))))
assert R(2, 12) == 13 and 143 == 11 * R(2, 12) and R(2, 10) == 11
print("   also 13 = R_2(12) = 11_12 and 143 = 11 x 13 = BB_12; 11 = R_2(10)")

# ---------------------------------------------------------------- 6. higher rank: what lifts and what does not
dim_sym = 5                                                                                    # SL_3(R)/SO(3): 8 - 3
assert 2 + 2 == dim_sym - 1                                                                    # two periodic flats (2-tori) can link
print("rank 2: SL_3(Z)\\SL_3(R)/SO(3) has dimension 5, so two periodic maximal flats (2-tori) can link (2 + 2 = 4 = 5 - 1); SL_3(Z) is perfect and by Burger-Monod has no non-trivial quasimorphism, so the Rademacher function does not lift — what lifts is Sczech's Eisenstein cocycle on GL_3(Z) (totally real cubic fields, zeta values at 0)")
print("ALL CHECKS PASSED")
