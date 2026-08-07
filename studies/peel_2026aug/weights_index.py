"""
Settle the index on F(s), then read it in primon-gas language.
E2: verify each Dirichlet series against a direct partial sum.
"""
from mpmath import mp, mpf, zeta, log, exp, pi as MPI, nstr
mp.dps = 30

N = 200000
lw = [mpf(0)]*(N+1)                      # divisor sieve for w_n
for d in range(2, N+1):
    c = log(mpf(d))/mpf(d)
    for n in range(d, N+1, d):
        lw[n] += c

def zprime(x, h=mpf('1e-10')):
    return (zeta(x+h) - zeta(x-h))/(2*h)

def partial(k, s, M):
    return sum(mpf(n)**k * lw[n] * mpf(n)**(-s) for n in range(1, M+1))

print("="*72)
print("  Which Dirichlet series goes with which weighting?")
print("="*72)
print("  Convolution n = d*m gives  sum_n n^k w_n n^-s  =  -zeta(s-k) zeta'(s-k+1)\n")
cases = [
    (0, "sum w_n n^-s",      lambda s: -zeta(s)*zprime(s+1),     "-zeta(s)zeta'(s+1)"),
    (1, "sum n w_n n^-s",    lambda s: -zeta(s-1)*zprime(s),     "-zeta(s-1)zeta'(s)"),
    (2, "sum n^2 w_n n^-s",  lambda s: -zeta(s-2)*zprime(s-1),   "-zeta(s-2)zeta'(s-1)"),
]
for k, lab, rhs, form in cases:
    s = mpf(6+k)
    L, R = partial(k, s, 60000), rhs(s)
    d = float(abs(L-R)/abs(R))
    print(f"  k={k}  {lab:<20} at s={float(s):.0f}")
    print(f"        partial {nstr(L,16):>24}   {form} {nstr(R,16):>24}")
    print(f"        agree to {float(-mp.log10(d)):.1f} digits\n")

print("  The dictionary states  F(s) = -zeta(s-1)zeta'(s) = sum n^2 w_n n^-s.")
print("  Those are two different objects: that RHS is the k=1 sum, not k=2.")
print("  The stiffness mu_N(y) ~ sum n^2 w_n e^{-2*pi*y*n} needs k=2,")
print("  so the correct pairing is  F(s) = -zeta(s-2)zeta'(s-1).\n")

print("="*72)
print("  Consequence: pole structure, and where T-infinity comes from")
print("="*72)
Tinf = -zprime(mpf(2))
print(f"  T_inf = -zeta'(2) = {nstr(Tinf,12)}\n")
print("  -zeta(s-2)zeta'(s-1):  simple pole s=3 (from zeta),")
print("                          double pole s=2 (from zeta').")
print("  Rightmost pole s=3, residue = -zeta'(2) = T_inf.\n")
print("  Inverse Mellin, mu(y) = (1/2pi i) int Gamma(s)(2*pi*y)^-s F(s) ds:")
print("     leading term = Gamma(3)(2*pi*y)^-3 * T_inf = T_inf/(4*pi^3*y^3)\n")
print("  E2 CHECK -- direct sum vs pole prediction, no shared route:")
print(f"  {'y':>8} {'direct sum':>22} {'T_inf/(4 pi^3 y^3)':>22} {'ratio':>10}")
for y in [mpf('0.05'), mpf('0.02'), mpf('0.01'), mpf('0.005'), mpf('0.002')]:
    a = 2*MPI*y
    M = min(N, int(60/float(a))+10)
    direct = sum(mpf(n)**2*lw[n]*exp(-a*n) for n in range(1, M+1))
    pred = Tinf/(4*MPI**3*y**3)
    print(f"  {float(y):>8.3f} {nstr(direct,12):>22} {nstr(pred,12):>22} {float(direct/pred):>10.5f}")
print("\n  -> ratio -> 1 as y -> 0. The leading asymptotic is fixed by T_inf")
print("     because T_inf IS the residue at the dominant pole. That is the")
print("     derivation of the 'torque normalization', not a fitted scale.\n")

print("="*72)
print("  Primon-gas reading")
print("="*72)
print("""  State n has energy E_n = log n, Boltzmann weight e^{-beta E_n} = n^-beta.
  Bare primon gas:      Z(beta)   = sum n^-beta            = zeta(beta)
  Dirichlet ensemble:   Z_D(beta) = sum w_n n^-beta        = -zeta(beta)zeta'(beta+1)

  Your Z_D(beta) = -zeta(beta)zeta'(beta+1) is exactly the k=0 line above.
  The primon-gas side of the corpus already has the correct form.

  Weighting by n^k is an energy tilt: n^k n^-beta = n^-(beta-k). So

        sum n^k w_n n^-beta  =  Z_D(beta - k)

  The n^2 in the stiffness is not a new object -- it is Z_D evaluated two
  units hotter. mu(y) is the Dirichlet ensemble at beta - 2.
""")
zd_res = -zprime(mpf(2))
print(f"  Z_D(beta) has a simple pole at beta=1, residue -zeta'(2) = {nstr(zd_res,12)} = T_inf.")
print("  Shifting beta -> beta-2 moves that pole to beta=3, which is the s=3")
print("  pole above. Same residue, same T_inf. One object, two temperatures.")
