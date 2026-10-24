"""A gap uniform in the refinement: where it is a theorem, where it is the whole problem.  2026-10-20.  EXACT unless marked.

Roger's prompt: the lab says a complete 4D theory needs a mass gap uniform in the lattice spacing, which no construction has; a search
engine read 'uniform gap in lattice spacing' as crystal planes, d_hkl = a / sqrt(h^2 + k^2 + l^2), and he noticed the three squares.
Two different gaps -- but the three-squares instinct lands on the one place where arithmetic provably supplies a spectral gap uniform
in a refinement.  This script checks:
(1) Refinement families and their spectral gaps: the clock Z/n and the flat torus (Z/L)^3 lose their gap as they grow; the modular
    group's finite quotients PSL(2, F_p), as Cayley graphs on {S, T, T^-1}, keep it (expanders, by Selberg's 3/16 through Lubotzky's
    theorem 4.4.2); the Lubotzky-Phillips-Sarnak graphs X^{5,q}, built from the SIX four-square representations of 5, are Ramanujan:
    every non-trivial eigenvalue has |lambda| <= 2 sqrt(5) (checked in full for q = 13, at the extremes for q = 17, 29).  The Ramanujan
    property is the Riemann hypothesis for the Ihara zeta of Instrument 15.
(2) Sums of squares: r_3(N) = 0 exactly for N = 4^a (8b + 7) (Gauss-Legendre); r_4(N) = 8 sum_{d|N, 4 not | d} d (Jacobi); so the cubic
    crystal has no plane family with d = a/sqrt(7), a/sqrt(15), ..., a/sqrt(143), while Z^4 has r_4(143) = 1344 = 8 * 168 vectors of
    norm 143, 168 = sigma(143) = [PSL(2,Z) : Gamma_0(143)]; r_3 counts class numbers (Gauss), the LPS degree p + 1 is r_4(p)/8.
(3) Two dimensions, where uniformity in a is exact: the character coefficients of a doubled plaquette are the squares of the originals
    (Migdal), so the Wilson loop is u_f^{A/a^2} at every spacing, and beta (-ln u_f) -> 3/2 as beta -> infinity for the Wilson action
    AND for CC-69's repaired action: sigma/g^2 = C_2(f)/2 = 3/8 either way; the tower is a finite shift beta -> beta + Delta(kappa),
    Delta -> (kappa/6) sum a_n n(n+1)(n+2), invisible in the continuum.
(4) Four dimensions: the strong-coupling gap -4 ln u_f against two-loop scaling a Lambda(beta): the ratio is not constant, and cannot be
    (roughening at beta ~ 2); the uniform statement is the open conjecture, with Teper's m(0++)/sqrt(sigma) = 3.74(12) as the number it
    would have to reproduce.
(5) The hyperbolic gap: the bottom of the Laplace spectrum on H^d is ((d-1)/2)^2, and on the congruence quotients lambda_1 >= 3/16
    (Selberg), >= 975/4096 (Kim-Sarnak), conjecturally 1/4 -- uniform in the level; generic surfaces can pinch (Cheeger, Randol).
    This is the gap VE-01 uses; it is geometric and arithmetic, not dynamical."""
import math, itertools
import numpy as np
from scipy.special import ive
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import eigsh
from scipy.optimize import brentq

# ---------------------------------------------------------------- 1. refinement families and their spectral gaps
def psl2_cayley(p):
    """Cayley graph of PSL(2, F_p) on {S, T, T^-1}: vertices are matrices up to sign, edges g -> g s."""
    def canon(m):
        a, b, c, d = [x % p for x in m]
        for x in (a, b, c, d):
            if x: return (a, b, c, d) if x <= (p - 1) // 2 else ((-a) % p, (-b) % p, (-c) % p, (-d) % p)
    elems = {}
    for a in range(p):
        for b in range(p):
            for c in range(p):
                for d in range(p):
                    if (a * d - b * c) % p == 1:
                        k = canon((a, b, c, d)); elems.setdefault(k, len(elems))
    mul = lambda m, n: ((m[0] * n[0] + m[1] * n[2]) % p, (m[0] * n[1] + m[1] * n[3]) % p, (m[2] * n[0] + m[3] * n[2]) % p, (m[2] * n[1] + m[3] * n[3]) % p)
    gens = [(0, p - 1, 1, 0), (1, 1, 0, 1), (1, p - 1, 0, 1)]              # S, T, T^-1
    n = len(elems); rows, cols = [], []
    for g, i in elems.items():
        for s in gens: rows.append(i); cols.append(elems[canon(mul(g, s))])
    A = csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(n, n))
    return A, elems, mul, canon
def spectrum(A):
    return np.sort(np.linalg.eigvalsh(A.toarray()))[::-1]
print("spectral gaps of refinement families (k-regular graphs: gap = k - lambda_2 of the adjacency)")
psl = {}
for p in (5, 7, 11, 13, 17, 19, 23):
    A, elems, mul, canon = psl2_cayley(p); n = A.shape[0]; assert n == p * (p * p - 1) // 2
    if n <= 1200: ev = spectrum(A)
    else: ev = np.concatenate([np.sort(eigsh(A, k=3, which='LA')[0])[::-1], np.sort(eigsh(A, k=2, which='SA')[0])])
    assert abs(ev[0] - 3) < 1e-7 and ev[1] < 3 - 1e-6                                            # connected: S and T generate
    # the period clock sits inside: T has order p, its orbits are p-cycles
    g = (1, 0, 0, 1); k = 0
    while True:
        g = canon(mul(g, (1, 1, 0, 1))); k += 1
        if g == (1, 0, 0, 1): break
    assert k == p
    psl[p] = (n, ev[1], ev[-1]); print("   PSL(2,F_%d) on {S, T, T^-1}: %d vertices, lambda_2 = %.4f, gap = %.4f, lambda_min = %.4f; T-orbits are %d-cycles whose own gap is 2 - 2cos(2pi/%d) = %.4f" % (p, n, ev[1], 3 - ev[1], ev[-1], p, p, 2 - 2 * math.cos(2 * math.pi / p)))
gaps = [3 - psl[p][1] for p in psl]; assert min(gaps) > 0.05, gaps
for n_ in (60, 660, 2448, 6072): print("   the clock Z/%d: gap 2 - 2cos(2pi/n) = %.5f;   the torus (Z/L)^3 with L^3 ~ %d sites: gap 4 sin^2(pi/L) = %.5f" % (n_, 2 - 2 * math.cos(2 * math.pi / n_), n_, 4 * math.sin(math.pi / round(n_ ** (1 / 3))) ** 2))
print("   Alon-Boppana: no family of 3-regular graphs can beat lambda_2 >= 2 sqrt 2 - o(1) = %.4f, i.e. gap <= %.4f; the {S, T} quotients are an expander family by Selberg's theorem (property tau of SL(2,Z) for congruence subgroups: Lubotzky 1994, Thm 4.4.2), but the guaranteed gap is small and the cusp cycles of T make the approach slow -- numerically the gap tracks the T-cycle's 2 - 2cos(2pi/p) out to p = 23, so this family shows the theorem's existence, not its size" % (2 * math.sqrt(2), 3 - 2 * math.sqrt(2)))
# Lubotzky-Phillips-Sarnak X^{5,q}: generators from a^2 + b^2 + c^2 + d^2 = 5, a > 0 odd, b, c, d even
def lps_generators(p):
    out = []
    for a, b, c, d in itertools.product(range(-p, p + 1), repeat=4):
        if a * a + b * b + c * c + d * d == p and a > 0 and a % 2 == 1 and b % 2 == 0 and c % 2 == 0 and d % 2 == 0: out.append((a, b, c, d))
    return out
for p in (5, 13, 17, 29): assert len(lps_generators(p)) == p + 1
print("LPS: the four-square representations of p with a odd positive and b, c, d even number exactly p + 1 (p = 5, 13, 17, 29): for p = 5 they are", lps_generators(5))
def lps_graph(p, q):
    i = next(x for x in range(q) if (x * x + 1) % q == 0)                                        # i^2 = -1 in F_q (q = 1 mod 4)
    def canon(m):                                                                                 # PGL_2(F_q): matrices up to scalars
        a, b, c, d = [x % q for x in m]
        for x in (a, b, c, d):
            if x:
                inv = pow(x, -1, q); return ((a * inv) % q, (b * inv) % q, (c * inv) % q, (d * inv) % q)
    mul = lambda m, n: ((m[0] * n[0] + m[1] * n[2]) % q, (m[0] * n[1] + m[1] * n[3]) % q, (m[2] * n[0] + m[3] * n[2]) % q, (m[2] * n[1] + m[3] * n[3]) % q)
    gens = [canon((a + b * i, c + d * i, -c + d * i, a - b * i)) for a, b, c, d in lps_generators(p)]
    start = canon((1, 0, 0, 1)); idx = {start: 0}; frontier = [start]; rows, cols = [], []
    while frontier:                                                                              # the Cayley graph of the group the generators generate
        nxt = []
        for g in frontier:
            for s in gens:
                h = canon(mul(g, s))
                if h not in idx: idx[h] = len(idx); nxt.append(h)
                rows.append(idx[g]); cols.append(idx[h])
        frontier = nxt
    n = len(idx); return csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(n, n)), n
legendre = lambda a, q: pow(a, (q - 1) // 2, q)
for q in (13, 17, 29):
    A, n = lps_graph(5, q); sym = legendre(5, q) == 1
    assert n == (q * (q * q - 1) // 2 if sym else q * (q * q - 1)), (q, n)
    assert (A - A.T).nnz == 0 and np.all(A.sum(axis=1) == 6)
    if q == 13:
        ev = spectrum(A); nontriv = ev[1:-1] if not sym else ev[1:]
        assert abs(ev[0] - 6) < 1e-9 and (sym or abs(ev[-1] + 6) < 1e-9) and np.max(np.abs(nontriv)) <= 2 * math.sqrt(5) + 1e-9
        lam2, lmin = ev[1], ev[-1]; mx = np.max(np.abs(nontriv))
    else:
        top = np.sort(eigsh(A, k=4, which='LA')[0])[::-1]; bot = np.sort(eigsh(A, k=4, which='SA')[0])
        lam2 = top[1]; lmin = bot[0]; nontriv = np.concatenate([top[1:], bot[1:] if not sym else bot]); mx = np.max(np.abs(nontriv))
        assert abs(top[0] - 6) < 1e-7 and mx <= 2 * math.sqrt(5) + 1e-6, (q, top, bot)
    print("   X^{5,%d} = Cay(%s(2,F_%d)): %d vertices, 6-regular%s; lambda_2 = %.4f, max |nontrivial| = %.4f <= 2 sqrt 5 = %.4f (Ramanujan), gap = %.4f >= 6 - 2 sqrt 5 = %.4f for every q" % (q, 'PSL' if sym else 'PGL', q, n, '' if sym else ', bipartite (lambda_min = -6)', lam2, mx, 2 * math.sqrt(5), 6 - lam2, 6 - 2 * math.sqrt(5)))
print("   the Ramanujan bound is the Riemann hypothesis for the graph's Ihara zeta (Instrument 15): the gap is uniform in q because the arithmetic of the quaternions supplies it (LPS 1988; Selberg's 3/16 behind the PSL(2,F_p) family, Lubotzky 1994 Thm 4.4.2)")

# ---------------------------------------------------------------- 2. sums of squares: the crystal planes and the 4D shells
NMAX = 200
r1 = np.zeros(NMAX + 1, dtype=np.int64); r1[0] = 1
for k in range(1, int(math.isqrt(NMAX)) + 1): r1[k * k] = 2
def conv(a, b):
    out = np.zeros(NMAX + 1, dtype=np.int64)
    for i in range(NMAX + 1):
        if a[i]: out[i:] += a[i] * b[:NMAX + 1 - i]
    return out
r2 = conv(r1, r1); r3 = conv(r2, r1); r4 = conv(r3, r1)
forbidden = lambda N: any(N == 4 ** a * (8 * b + 7) for a in range(8) for b in range(NMAX // 8 + 1))
for N in range(1, NMAX + 1):
    assert (r3[N] == 0) == forbidden(N), N                                                      # Gauss-Legendre
    assert r4[N] == 8 * sum(d for d in range(1, N + 1) if N % d == 0 and d % 4 != 0), N          # Jacobi
sigma = lambda N: sum(d for d in range(1, N + 1) if N % d == 0)
assert r3[143] == 0 and r4[143] == 1344 == 8 * 168 and sigma(143) == 168 and 143 * (1 + 1 / 11) * (1 + 1 / 13) == 168
assert r2[13] == 8 and r2[7] == 0 and r2[11] == 0 and r3[7] == 0 and r3[11] == 24 and r3[13] == 24 and r4[7] == 64
def class_number(D):                                                                               # reduced forms a x^2 + b x y + c y^2, b^2 - 4ac = D < 0
    h = 0
    for a in range(1, int(math.isqrt(-D // 3)) + 2):
        for b in range(-a + 1, a + 1):
            if (b * b - D) % (4 * a): continue
            c = (b * b - D) // (4 * a)
            if c < a or (c == a and b < 0): continue
            h += 1
    return h
assert class_number(-11) == 1 and class_number(-52) == 2 and class_number(-23) == 3 and class_number(-4) == 1
assert r3[11] == 24 * class_number(-11) and r3[13] == 12 * class_number(-52) and r3[19] == 24 * class_number(-19)
forb = [N for N in range(1, 64) if r3[N] == 0]
print("sums of squares: r_3(N) = 0 exactly for N = 4^a(8b+7) (checked to %d): %s ... and 143 -- no cubic-crystal plane family has d = a/sqrt(N) for these; r_4(N) = 8 sum_{4 not| d | N} d (Jacobi, checked to %d), so no 4D shell is empty" % (NMAX, forb, NMAX))
print("   the period-clock primes: 7 = 4+1+1+1 needs four squares (r_3 = 0, r_4 = 64); 11 = 9+1+1 (r_3 = 24 = 24 h(-11)); 13 = 9+4 (r_2 = 8, r_3 = 24 = 12 h(-52)); 143 = 8*17 + 7: r_3 = 0, r_4(143) = 1344 = 8 * 168, and 168 = sigma(143) = 143 (1 + 1/11)(1 + 1/13) = [PSL(2,Z) : Gamma_0(143)]")
t_ = 0.7; lhs = sum(math.exp(-math.pi * n * n * t_) for n in range(-30, 31)); rhs = sum(math.exp(-math.pi * n * n / t_) for n in range(-30, 31)) / math.sqrt(t_)
assert abs(lhs - rhs) < 1e-12
print("   theta(t) = t^{-1/2} theta(1/t) (Poisson, 1e-12): the lattice and its dual lattice are the same object seen at spacing a and 1/a; theta^3 (weight 3/2) counts class numbers, theta^4 (weight 2, Gamma_0(4)) counts divisor sums")

# ---------------------------------------------------------------- 3. two dimensions: uniform in a, exactly
Y = 0.18174; NMAXT = 24
W = np.zeros(NMAXT + 1)
for n in range(1, NMAXT + 1): W[n] = sum(math.log(d) / d for d in range(1, n + 1) if n % d == 0)
A_N = np.array([W[n] * math.exp(-2 * math.pi * Y * n) for n in range(NMAXT + 1)]); A_N[1] = 0
NQ = 40000; th = (np.arange(NQ) + 0.5) * math.pi / NQ; mu = (2 / math.pi) * np.sin(th) ** 2 * (math.pi / NQ)
chi = lambda j2: np.sin((j2 + 1) * th) / np.sin(th)
def action(beta, kappa):
    s = beta * (1 - np.cos(th))
    for n in range(2, NMAXT + 1): s += (kappa / 2) * A_N[n] * ((n + 1) - chi(n))
    return s
def uf(beta, kappa):
    f = np.exp(-(action(beta, kappa) - 0.0)); c0 = np.sum(f * mu); return np.sum(f * chi(1) * mu) / 2 / c0
uW = lambda beta: ive(2, beta) / ive(1, beta)                                                      # Wilson: I_2/I_1, exponentially scaled for large beta
for b in (0.7, 2.3, 10.0): assert abs(uf(b, 0) - uW(b)) < 1e-9
print("two dimensions (Migdal): the doubled plaquette's character coefficients are the squares of the originals, so <W(A)> = u_f^{A/a^2} exactly and sigma a^2 = -ln u_f(beta) at every spacing")
rows = []
for b in (50.0, 100.0, 200.0, 400.0, 800.0):
    sW = b * (-math.log(uW(b))); sR = b * (-math.log(uf(b, 1.0))); rows.append((b, sW, sR))
    print("   beta = %5.0f: beta(-ln u_f) = %.5f (Wilson), %.5f (repaired, kappa = 1); 3/2 + 3/(8 beta) + ... = %.5f" % (b, sW, sR, 1.5 + 3 / (8 * b)))
assert abs(rows[-1][1] - 1.5) < 0.003 and abs(rows[-1][2] - 1.5) < 0.003 and abs(rows[-1][1] - rows[-1][2]) < abs(rows[0][1] - rows[0][2])
print("   -> sigma/g^2 = beta(-ln u_f)/4 -> 3/8 = C_2(fundamental)/2 for both actions with beta = 4/(g^2 a^2): the string tension is uniform in a, and the tower changes nothing in the continuum")
delta_inf = sum(A_N[n] * n * (n + 1) * (n + 2) for n in range(2, NMAXT + 1)) / 6
shifts = []
for b in (2.3, 20.0, 100.0, 400.0, 1600.0, 6400.0):
    beff = brentq(lambda x: uW(x) - uf(b, 1.0), b * 0.5, b * 2.0); shifts.append(beff - b)
print("   the tower as a shift of beta (kappa = 1): beta_eff - beta = %s at beta = 2.3, 20, 100, 400, 1600, 6400 -> the small-angle limit (kappa/6) sum a_n n(n+1)(n+2) = %.4f (section 10's 0.051 is the value at beta = 2.3)" % (", ".join("%.4f" % s for s in shifts), delta_inf))
assert abs(shifts[-1] - delta_inf) < 2e-3 and all(shifts[i] < shifts[i + 1] for i in range(len(shifts) - 1))

# ---------------------------------------------------------------- 4. four dimensions: the strong-coupling gap against asymptotic scaling
def lam_over_mu(g2, N=2):
    b0 = 11 * N / 3; b1 = 34 * N * N / 3; x = b0 * g2 / (16 * math.pi ** 2)
    return x ** (-b1 / (2 * b0 ** 2)) * math.exp(-8 * math.pi ** 2 / (b0 * g2))
print("four dimensions: the leading strong-coupling gap m a = -4 ln u_f(beta) against two-loop scaling a Lambda_L(beta), g^2 = 4/beta (SU(2))")
ratios = []
for b in (1.0, 1.5, 2.0, 2.5, 3.0):
    ma = -4 * math.log(uW(b)); aL = lam_over_mu(4 / b); ratios.append(ma / aL)
    print("   beta = %.1f: m a = %.3f, a Lambda_L = %.3e, m/Lambda_L = %.3e" % (b, ma, aL, ma / aL))
assert max(ratios) / min(ratios) > 10
print("   the ratio moves by a factor %.0f across beta = 1..3: the strong-coupling series cannot deliver a gap uniform in a (its radius of convergence ends at the roughening transition, beta ~ 2 for SU(2)); 'uniform in a' means m(beta(a))/Lambda_L -> a positive constant, which is the open conjecture" % (max(ratios) / min(ratios)))
print("   the number it would have to reproduce (lattice Monte Carlo, continuum-extrapolated): m(0++)/sqrt(sigma) = 3.74(12) for SU(2), 3.64(9) for SU(3) (Teper 1998)")

# ---------------------------------------------------------------- 5. the hyperbolic gap and its arithmetic uniformity
for d in (2, 3, 4): print("   H^%d: the Laplace spectrum starts at ((d-1)/2)^2 = %.4f (curvature units): a massless field in hyperbolic space is gapped, uniformly in any cutoff" % (d, ((d - 1) / 2) ** 2))
print("   congruence quotients: lambda_1(Gamma_0(N)\\H) >= 3/16 = %.4f (Selberg 1965), >= 975/4096 = %.4f (Kim-Sarnak 2003), conjecturally 1/4 -- uniform in the level N; generic surfaces pinch (Cheeger: lambda_1 >= h^2/4, Randol: lambda_1 can be made arbitrarily small)" % (3 / 16, 975 / 4096))
print("the map of theorems: d = 2 -- exact and uniform (Migdal, no propagating modes); d = 3 U(1) -- gap for all beta (Gopfert-Mack 1982) but m ~ e^{-c beta} -> 0 in physical units; d = 4 U(1) -- the gap CLOSES at weak coupling (Guth 1980, Frohlich-Spencer 1982: a Coulomb phase); d = 4 SU(N) -- gap at strong coupling (Osterwalder-Seiler 1978), ultraviolet stability uniform in a (Balaban 1984-89), the uniform gap open; arithmetic surfaces -- gap uniform in the level (Selberg); Ramanujan graphs -- gap uniform in the size (LPS)")
print("ALL CHECKS PASSED")
