"""The Casimir knot: the zeta-regularised vacuum energies are the constant terms of the Eisenstein series, the first of them is the
Dedekind eta's 1/24, and eta's phase round a closed geodesic of the modular surface is Rademacher's function, which Ghys showed is
the linking number of the modular knot with the trefoil.  Then Godel's universe over the modular surface, 4D rotations as torus knots,
radicals as braids, digit dynamics as a base-2 period clock on 143a1, and Siegel zeros as thin three-square shells.  2026-10-20.
EXACT unless marked."""
import math, itertools, bisect, cmath
from fractions import Fraction
import numpy as np
import sympy as sp
import mpmath as mp
mp.mp.dps = 30
rng = np.random.default_rng(5)
t = sp.symbols('t')

# ---------------------------------------------------------------- 1. Casimir energies are the constant terms of Eisenstein series
mp.mp.dps = 50
def regulated(power, eps):                                  # sum_{n>=1} n^power e^{-eps n} = Li_{-power}(e^{-eps}), at 50 digits
    return mp.polylog(-power, mp.exp(-mp.mpf(eps)))
for k, zeta_val in ((1, -1 / 12), (2, 1 / 120), (3, -1 / 252), (4, 1 / 240)):
    p = 2 * k - 1; eps = mp.mpf('0.02')
    finite = regulated(p, eps) - mp.factorial(p) / eps ** (2 * k)                                 # strip the divergent bulk term
    assert abs(finite - zeta_val) < 1e-5, (k, finite, zeta_val)
    assert float(sp.zeta(1 - 2 * k)) == zeta_val and sp.zeta(1 - 2 * k) == -sp.bernoulli(2 * k) / (2 * k)
mp.mp.dps = 30
print("Casimir regularisation: sum n^{2k-1} e^{-eps n} = (2k-1)!/eps^{2k} + zeta(1-2k) + O(eps^2): finite parts -1/12, 1/120, -1/252, 1/240 for k = 1..4 (eps = 0.02, 1e-5)")
print("   Eisenstein series: E_{2k} = (2/zeta(1-2k)) [zeta(1-2k)/2 + sum sigma_{2k-1}(n) q^n], coefficients 2/zeta = %s: the vacuum sum (1/2) sum n^{2k-1}, regularised, IS the constant term" % [int(2 / sp.zeta(1 - 2 * k)) for k in (1, 2, 3, 4)])
print("   the physics: a string of length L, E = (pi/2L) zeta(-1) = -pi/(24 L) (one chirality) or -pi/(6L) (both); plates at distance a, E/A = -(pi^2/6a^3) zeta(-3) = -pi^2/(720 a^3): the numbers 24 and 240 of E_2 and E_4")

# ---------------------------------------------------------------- 2. the knot: log eta round a closed geodesic is Rademacher's function = the linking number with the trefoil (Ghys)
NQ = 4000
SIG = np.array([0] + [sum(d for d in range(1, n + 1) if n % d == 0) for n in range(1, NQ + 1)], dtype=float)
def log_eta(tau):                                           # the continuous branch on H: pi i tau/12 - sum sigma(n)/n q^n
    q = cmath.exp(2j * math.pi * tau); n = np.arange(1, NQ + 1); qn = q ** n
    assert abs(q) ** NQ < 1e-14
    return 1j * math.pi * tau / 12 - complex(np.sum(SIG[1:] / n * qn))
def E2(tau):
    q = cmath.exp(2j * math.pi * tau); n = np.arange(1, NQ + 1); return 1 - 24 * complex(np.sum(SIG[1:] * q ** n))
tau = 0.1 + 0.3j; h = 1e-5
assert abs((log_eta(tau + h) - log_eta(tau - h)) / (2 * h) - (1j * math.pi / 12) * E2(tau)) < 1e-8
assert abs(cmath.exp(log_eta(tau)) - complex(mp.qp(mp.exp(2j * mp.pi * tau)) * mp.exp(2j * mp.pi * tau / 24))) < 1e-14
def dedekind_sum(h, k):
    saw = lambda x: Fraction(0) if x.denominator == 1 else x - (x.numerator // x.denominator) - Fraction(1, 2)
    return sum(saw(Fraction(r, k)) * saw(Fraction(h * r, k)) for r in range(1, k))
def rademacher_Phi(A):                                      # Phi(A) = (a + d)/c - 12 s(d, c) for c > 0
    (a, b), (c, d) = A; assert c > 0; return Fraction(a + d, c) - 12 * dedekind_sum(d % c, c)
def psi_ghys(A):                                            # the lab's psi = Phi - 3 sign(c (a + d)) = lk(k_A, trefoil)
    (a, b), (c, d) = A; return rademacher_Phi(A) - 3 * (1 if c * (a + d) > 0 else -1)
L = np.array([[1, 0], [1, 1]]); R = np.array([[1, 1], [0, 1]])
word = lambda w: [[int(x) for x in row] for row in np.linalg.multi_dot([L if ch == 'L' else R for ch in w]).tolist()] if len(w) > 1 else (L if w == 'L' else R).tolist()
cases = {'LR': word('LR'), 'L²R': word('LLR'), 'LR²': word('LRR'), 'L³R': word('LLLR'), 'L²R²': word('LLRR'), 'L²RLR': word('LLRLR'), '[84 -47; 143 -80]': [[84, -47], [143, -80]], '[12 1; 143 12]': [[12, 1], [143, 12]]}
print("Rademacher from the Casimir term: Phi(A) = (12/pi i)[log eta(A tau) - log eta(tau) - (1/2) log(-i(c tau + d))], tau on the isometric circle |c tau + d| = 1 (Im tau = 1/c)")
for name, A in cases.items():
    (a, b), (c, d) = A; assert a * d - b * c == 1 and c > 0
    tau0 = (-d + 1j) / c; Atau = (a * tau0 + b) / (c * tau0 + d); assert abs(Atau.imag - tau0.imag) < 1e-12
    Phi_num = (12 / (1j * math.pi)) * (log_eta(Atau) - log_eta(tau0) - 0.5 * cmath.log(-1j * (c * tau0 + d)))
    Phi = rademacher_Phi(A); psi = psi_ghys(A)
    assert abs(Phi_num - float(Phi)) < 1e-6 and abs(Phi_num.imag) < 1e-6, (name, Phi_num, Phi)
    assert Phi.denominator == 1
    print("   A = %-18s trace %3d: Phi = %+d (numerically %+.8f), lk(k_A, trefoil) = Phi - 3 = %+d" % (name, a + d, Phi, Phi_num.real, psi))
assert psi_ghys(cases['LR']) == 0 and psi_ghys(cases['L²R']) == -1 and psi_ghys(cases['[84 -47; 143 -80]']) == -1 and psi_ghys(cases['[12 1; 143 12]']) == -21
print("   -> log eta is the integral of (pi i/12) E_2; the Casimir constant 1 of E_2 (i.e. -1/24 in eta) is what makes the phase round the geodesic an integer Phi, and Ghys's theorem reads it as the linking number of the modular knot with the trefoil (Instrument 10's psi: 0, -1, ..., -21 reproduced)")

# ---------------------------------------------------------------- 3. Godel's universe: sl(2,R) + R Killing algebra, the null radius, and the quotient by the braid group
tt, x, y, z, w_, r, phi = sp.symbols('t x y z omega r phi', real=True)
coords = [tt, x, y, z]
g = sp.Matrix(4, 4, lambda i, j: 0)
# ds^2 = (1/2 omega^2) [-(dt + e^x dy)^2 + dx^2 + e^{2x} dy^2 / 2 + dz^2]
g[0, 0] = -1; g[0, 2] = g[2, 0] = -sp.exp(x); g[2, 2] = -sp.exp(2 * x) + sp.exp(2 * x) / 2; g[1, 1] = 1; g[3, 3] = 1
g = g / (2 * w_ ** 2)
K = {'K1': [1, 0, 0, 0], 'K2': [0, 0, 1, 0], 'K3': [0, 0, 0, 1], 'K4': [0, 1, -y, 0], 'K5': [-2 * sp.exp(-x), y, sp.exp(-2 * x) - y ** 2 / 2, 0]}
def lie_derivative_metric(X):
    Lg = sp.zeros(4, 4)
    for i in range(4):
        for j in range(4):
            Lg[i, j] = sum(X[k] * sp.diff(g[i, j], coords[k]) + g[k, j] * sp.diff(X[k], coords[i]) + g[i, k] * sp.diff(X[k], coords[j]) for k in range(4))
    return sp.simplify(Lg)
for name, X in K.items(): assert lie_derivative_metric(X) == sp.zeros(4, 4), name
def bracket(X, Y): return [sp.simplify(sum(X[k] * sp.diff(Y[i], coords[k]) - Y[k] * sp.diff(X[i], coords[k]) for k in range(4))) for i in range(4)]
eq = lambda X, Y: all(sp.simplify(X[i] - Y[i]) == 0 for i in range(4))
assert eq(bracket(K['K4'], K['K2']), K['K2']) and eq(bracket(K['K4'], K['K5']), [-c for c in K['K5']]) and eq(bracket(K['K2'], K['K5']), K['K4'])
assert all(eq(bracket(K['K1'], K[n]), [0, 0, 0, 0]) for n in K)
assert g[0, 0].subs(w_, 1) < 0                                 # the dust worldlines d/dt are timelike everywhere
print("Godel: the five Killing fields of Wikipedia's chart are Killing (Lie derivative of g vanishes, sympy); K2, K4, K5 close on [K4,K2] = K2, [K4,K5] = -K5, [K2,K5] = K4 -- sl(2,R), Bianchi VIII -- and d/dt commutes with all: isometry algebra sl(2,R) + R on the (t,x,y) part, d/dt timelike (g_tt = -1/(2 omega^2))")
# the cylindrical chart: g(d_phi, d_phi) = (sinh^2 / omega^2) [1/2 - (cosh - 1)/(cosh + 1)] vanishes at cosh(sqrt 2 omega r) = 3
u = sp.symbols('u', positive=True)                            # u = sqrt(2) omega r
gphiphi = sp.sinh(u) ** 2 * (sp.Rational(1, 2) - (sp.cosh(u) - 1) / (sp.cosh(u) + 1))
assert sp.simplify(gphiphi.subs(u, sp.acosh(3))) == 0 and gphiphi.subs(u, 1.0) > 0 and gphiphi.subs(u, 2.0) < 0
print("   the circle d/dphi turns null at cosh(sqrt2 omega r) = 3, r_c = arccosh(3)/(sqrt2 omega) = %.4f/omega, timelike beyond: Godel's closed timelike curves" % (math.acosh(3) / math.sqrt(2)))
print("   Rooman-Spindel 1998: the (t,x,y) part is squashed AdS_3 = SL(2,R) with a left-invariant metric, isometries SL(2,R)_L x U(1)_R, U(1)_R = d/dt; Milnor: SL(2,Z)\\SL(2,R) = S^3 minus the trefoil with the U(1)_R orbits as Seifert fibres (regular fibre links K six times)")
print("   so the quotient of Godel's universe by the braid group B_3 (the preimage of SL(2,Z) in the universal cover, = pi_1 of the trefoil complement) is (S^3 - trefoil) x R_z with a locally-Godel metric in which every dust worldline is a closed timelike geodesic, a Seifert fibre linking the trefoil six times (three and two for the fibres over i and rho)")

# ---------------------------------------------------------------- 4. 4D rotations: orbits of double rotations are torus knots; quaternion pairs
def alexander(components, tries=6):                        # Alexander polynomial of a polygonal link from a generic projection (Wirtinger + Fox), as in bessel_knots_20261019.py
    for _ in range(tries):
        Q, _ = np.linalg.qr(rng.normal(size=(3, 3))); comps = [np.asarray(c)[:, :3] @ Q.T for c in components]; cross = []; ok = True
        for ci, P in enumerate(comps):
            for cj, Qp in enumerate(comps):
                if cj < ci: continue
                N, M = len(P), len(Qp); p = P[:, :2]; rr = np.roll(P, -1, 0)[:, :2] - p; q = Qp[:, :2]; s = np.roll(Qp, -1, 0)[:, :2] - q
                dp = P[:, 2]; dq = Qp[:, 2]; dpn = np.roll(dp, -1); dqn = np.roll(dq, -1)
                for i in range(N):
                    js = np.arange(M)
                    if ci == cj: js = js[(js > i + 1) & ~((i == 0) & (js == N - 1))]
                    if len(js) == 0: continue
                    den = rr[i, 0] * s[js, 1] - rr[i, 1] * s[js, 0]; good = np.abs(den) > 1e-14; qp = q[js] - p[i]
                    ti = np.where(good, (qp[:, 0] * s[js, 1] - qp[:, 1] * s[js, 0]) / np.where(good, den, 1), -1); tj = np.where(good, (qp[:, 0] * rr[i, 1] - qp[:, 1] * rr[i, 0]) / np.where(good, den, 1), -1)
                    hit = good & (ti >= 0) & (ti < 1) & (tj >= 0) & (tj < 1)
                    for k in np.nonzero(hit)[0]:
                        j = js[k]
                        if min(ti[k], 1 - ti[k], tj[k], 1 - tj[k]) < 1e-6: ok = False
                        di = dp[i] + ti[k] * (dpn[i] - dp[i]); dj = dq[j] + tj[k] * (dqn[j] - dq[j]); over_i = di > dj
                        d_over, d_under = (rr[i], s[j]) if over_i else (s[j], rr[i]); sign = np.sign(d_over[0] * d_under[1] - d_over[1] * d_under[0])
                        cross.append((ci, i + ti[k], cj, j + tj[k], over_i, sign))
        if not ok: continue
        unders = [[] for _ in comps]
        for ci, si, cj, sj, over_i, sign in cross:
            if over_i: unders[cj].append(sj)
            else: unders[ci].append(si)
        for uu in unders: uu.sort()
        offset = []; narcs = 0
        for uu in unders: offset.append(narcs); narcs += max(1, len(uu))
        def arc_of(c, s_):
            uu = unders[c]
            if not uu: return offset[c]
            return offset[c] + (bisect.bisect_right(uu, s_) - 1) % len(uu)
        n = len(cross)
        if n == 0: return sp.Integer(1), 0
        A = sp.zeros(n, narcs)
        for row, (ci, si, cj, sj, over_i, sign) in enumerate(cross):
            if over_i: co, so, cu, su = ci, si, cj, sj
            else: co, so, cu, su = cj, sj, ci, si
            kk = arc_of(co, so); uu = unders[cu]; idx = bisect.bisect_left(uu, su); out_arc = offset[cu] + idx; in_arc = offset[cu] + (idx - 1) % len(uu)
            A[row, kk] += 1 - t
            if sign > 0: A[row, in_arc] += t; A[row, out_arc] += -1
            else: A[row, in_arc] += -1; A[row, out_arc] += t
        D = sp.expand(A[:-1, :-1].det(method='berkowitz'))
        if D == 0: return sp.Integer(0), n
        P_ = sp.Poly(D, t); low = min(mon[0] for mon in P_.monoms()); D = sp.expand(sp.cancel(D / t ** low))
        if sp.Poly(D, t).LC() < 0: D = -D
        return D, n
    raise RuntimeError('no generic projection')
def same_poly(D, E): return sp.expand(D - E) == 0 or sp.expand(D - sp.expand(sp.cancel(E.subs(t, 1 / t) * t ** sp.degree(E, t)))) == 0
def stereo(v): return np.array([v[0], v[1], v[2]]) / (1 - v[3])
def orbit(p, q, base=(1 / math.sqrt(2), 0, 1 / math.sqrt(2), 0), n=720, shift=0.0):
    pts = []
    for i in range(n):
        th = 2 * math.pi * i / n; z0 = complex(base[0], base[1]) * cmath.exp(1j * (p * th + shift)); z1 = complex(base[2], base[3]) * cmath.exp(1j * q * th)
        pts.append(stereo(np.array([z0.real, z0.imag, z1.real, z1.imag]) @ ROT.T))
    return pts
ROT, _ = np.linalg.qr(rng.normal(size=(4, 4)))                                                   # a generic rotation of S^3 before projecting, to keep the pole off the orbits
torus_poly = lambda p, q: sp.expand(sp.cancel((t ** (p * q) - 1) * (t - 1) / ((t ** p - 1) * (t ** q - 1))))
for p, q in ((2, 3), (2, 5), (3, 4), (3, 5)):
    D, nc = alexander([orbit(p, q)]); assert same_poly(D, torus_poly(p, q)), (p, q, D)
    print("   the orbit of a point of the Clifford torus under the double rotation R(p th, q th), (p,q) = (%d,%d): Alexander polynomial %s = T(%d,%d) (%d crossings in the projection)" % (p, q, D, p, q, nc))
D, nc = alexander([orbit(1, 1), orbit(1, 1, base=(0.6, 0, 0.8, 0), shift=1.0)]); assert same_poly(D, t - 1)
print("   isoclinic (1,1): two orbits form the Hopf link (Alexander t - 1); the Hopf fibration is the orbit structure of a left-isoclinic one-parameter subgroup of SO(4)")
# quaternion pairs: x -> a x conj(b) with a = e^{alpha i}, b = e^{beta i} is the double rotation by alpha - beta in the (1, i) plane and alpha + beta in the (j, k) plane
def qmul(p, q):
    a1, b1, c1, d1 = p; a2, b2, c2, d2 = q
    return np.array([a1 * a2 - b1 * b2 - c1 * c2 - d1 * d2, a1 * b2 + b1 * a2 + c1 * d2 - d1 * c2, a1 * c2 - b1 * d2 + c1 * a2 + d1 * b2, a1 * d2 + b1 * c2 - c1 * b2 + d1 * a2])
def rot_matrix(a, b):                                                                            # the 4x4 matrix of x -> a x conj(b)
    bc = np.array([b[0], -b[1], -b[2], -b[3]]); return np.array([qmul(qmul(a, e), bc) for e in np.eye(4)]).T
for alpha, beta in ((0.7, 0.2), (1.1, -0.4)):
    a = np.array([math.cos(alpha), math.sin(alpha), 0, 0]); b = np.array([math.cos(beta), math.sin(beta), 0, 0]); M = rot_matrix(a, b)
    th1, th2 = alpha - beta, alpha + beta
    expected = np.block([[np.array([[math.cos(th1), -math.sin(th1)], [math.sin(th1), math.cos(th1)]]), np.zeros((2, 2))], [np.zeros((2, 2)), np.array([[math.cos(th2), -math.sin(th2)], [math.sin(th2), math.cos(th2)]])]])
    assert np.allclose(M, expected) and np.allclose(M @ M.T, np.eye(4)) and abs(np.linalg.det(M) - 1) < 1e-12
for _ in range(5):                                                                                # random unit quaternion pairs give SO(4); so(4) = so(3) + so(3)
    a = rng.normal(size=4); a /= np.linalg.norm(a); b = rng.normal(size=4); b /= np.linalg.norm(b); M = rot_matrix(a, b); assert np.allclose(M @ M.T, np.eye(4)) and abs(np.linalg.det(M) - 1) < 1e-12
print("4D rotations: x -> a x conj(b) with unit quaternions a = e^{alpha i}, b = e^{beta i} is the double rotation by (alpha - beta, alpha + beta): so a (p,q) double rotation is the pair ((p+q)/2, (q-p)/2), isoclinic iff p = q (Hopf); SO(4) = SU(2) x SU(2) / +-1, and the three matrices of the picture are the diagonal SO(3), x -> a x conj(a)")
# the lattice action (t^4 g2, t^6 g3): a (4,6) double rotation whose orbits are (2,3) curves, the trefoil g2^3 = 27 g3^2 among them
rr_ = 0.5722153; traj = [(3 * rr_ ** 2 * cmath.exp(4j * th), rr_ ** 3 * cmath.exp(6j * th)) for th in np.linspace(0, math.pi, 400, endpoint=False)]
assert all(abs(g2 ** 3 - 27 * g3 ** 2) < 1e-12 for g2, g3 in traj) and all(abs(abs(g2) ** 2 + abs(g3) ** 2 - 1) < 1e-6 for g2, g3 in traj)
print("   the lattice rotation (e^{4i th} g2, e^{6i th} g3) is the (4,6) = 2(2,3) double rotation: its orbits are the Seifert fibres of Instrument 10, and the trefoil g2^3 = 27 g3^2 is the orbit through (3r^2, r^3), r = 0.57222, closing at th = pi")

# ---------------------------------------------------------------- 5. radicals are braids: (a^{1/n})^m = (a^m)^{1/n} as sets iff gcd(m,n) = 1 iff T(m,n) is a knot
a0 = 2.0 * cmath.exp(0.7j)
for m, n in itertools.product(range(1, 7), repeat=2):
    key = lambda z: (round(z.real, 8), round(z.imag, 8))
    left = sorted({key(abs(a0) ** (m / n) * cmath.exp(1j * (m * (cmath.phase(a0) + 2 * math.pi * k) / n))) for k in range(n)})
    right = sorted({key(abs(a0) ** (m / n) * cmath.exp(1j * ((m * cmath.phase(a0) + 2 * math.pi * j) / n))) for j in range(n)})
    d = math.gcd(m, n); assert len(left) == n // d and len(right) == n and ((left == right) == (d == 1)), (m, n)
    perm = {k: (k + m) % n for k in range(n)}; seen = set(); cycles = []
    for k in range(n):
        if k in seen: continue
        c_ = 0; j = k
        while j not in seen: seen.add(j); j = perm[j]; c_ += 1
        cycles.append(c_)
    assert len(cycles) == d and all(c_ == n // d for c_ in cycles)
print("radicals: (a^{1/n})^m has n/gcd(m,n) values, (a^m)^{1/n} has n, and the two sets agree exactly when gcd(m,n) = 1 (checked m, n <= 6); carrying a once round 0 permutes the n branches by k -> k + m: gcd(m,n) cycles of length n/gcd -- the torus braid (sigma_1 ... sigma_{n-1})^m whose closure is T(m,n), a knot iff the identity holds")
print("   the real-line rules of the picture (n odd: a; n even: |a|) are the one real branch; the full n-valued root is the braid, and the link of w^n = a^m is the torus link with gcd(m,n) components")

# ---------------------------------------------------------------- 6. digit dynamics (arithmetic dynamics, base by base) and the doubling clock on 143a1
def kaprekar_cycles(B, k):
    N = B ** k; f = np.zeros(N, dtype=np.int64)
    for v in range(N):
        digs = []; x = v
        for _ in range(k): digs.append(x % B); x //= B
        digs.sort(); asc = 0; desc = 0
        for dg in digs: asc = asc * B + dg
        for dg in reversed(digs): desc = desc * B + dg
        f[v] = desc - asc
    cycles = []; state = np.zeros(N, dtype=np.int8)                                                 # 0 unvisited, 1 on path, 2 done
    for v in range(N):
        if state[v]: continue
        path = []; x = v
        while state[x] == 0: state[x] = 1; path.append(x); x = f[x]
        if state[x] == 1:
            cyc = path[path.index(x):]; cycles.append(tuple(int(v) for v in cyc))
        for y in path: state[y] = 2
    return sorted(cycles, key=lambda c: (len(c), c))
def dyn_zeta(cycles):
    zz = sp.symbols('z'); return sp.factor(sp.prod(1 / (1 - zz ** len(c)) for c in cycles))
exp10 = {3: [(0,), (495,)], 4: [(0,), (6174,)]}
for k in (3, 4, 5, 6):
    cyc = kaprekar_cycles(10, k); lens = sorted(len(c) for c in cyc)
    if k in exp10: assert cyc == sorted(exp10[k], key=lambda c: (len(c), c)), cyc
    if k == 5: assert lens == [1, 2, 4, 4] and (53955, 59994) in [tuple(sorted(c)) for c in cyc]
    if k == 6: assert lens == [1, 1, 1, 7] and {549945, 631764} <= {c[0] for c in cyc if len(c) == 1}
    print("   Kaprekar, base 10, %d digits: cycles of lengths %s (%s); dynamical zeta exp(sum #Fix(f^n) z^n/n) = %s" % (k, lens, ", ".join(str(c[0]) if len(c) == 1 else "(%s ...)" % c[0] for c in cyc), dyn_zeta(cyc)))
narc = [v for v in range(1, 10 ** 7) if v == sum(int(ch) ** len(str(v)) for ch in str(v))]
assert narc[:16] == [1, 2, 3, 4, 5, 6, 7, 8, 9, 153, 370, 371, 407, 1634, 8208, 9474] and narc[16:] == [54748, 92727, 93084, 548834, 1741725, 4210818, 9800817, 9926315]
fact = [v for v in range(1, 2000000) if v == sum(math.factorial(int(ch)) for ch in str(v))]; assert fact == [1, 2, 145, 40585]
print("   narcissistic numbers below 10^7 (fixed points of the digit-power map): %s; factorions: %s; 407 = 11 x 37 and 6174 = 2 x 3^2 x 7^3, 495 = 5 x 9 x 11 -- the period-clock moduli, but these are base-10 accidents (the maps depend on the base)" % (narc[9:], fact))
for B in (2, 3, 7, 12, 16):
    print("      base %2d: Kaprekar 3-digit cycles %s, 4-digit cycle lengths %s" % (B, [c if len(c) > 1 else c[0] for c in kaprekar_cycles(B, 3)], sorted(len(c) for c in kaprekar_cycles(B, 4))))
# the doubling map on 143a1 over F_p: y^2 + y = x^3 - x^2 - x - 2; its cycles are the period clock with base 2 and modulus the order of the point
def curve_points(p):
    pts = [None]                                                                                   # None = the point at infinity
    for X in range(p):
        for Y in range(p):
            if (Y * Y + Y - (X ** 3 - X * X - X - 2)) % p == 0: pts.append((X, Y))
    return pts
def add(P, Q, p):                                                                                  # a1 = 0, a2 = -1, a3 = 1, a4 = -1, a6 = -2
    if P is None: return Q
    if Q is None: return P
    x1, y1 = P; x2, y2 = Q
    if x1 == x2 and (y1 + y2 + 1) % p == 0: return None
    if P == Q: lam = ((3 * x1 * x1 - 2 * x1 - 1) * pow(2 * y1 + 1, -1, p)) % p
    else: lam = ((y2 - y1) * pow(x2 - x1, -1, p)) % p
    x3 = (lam * lam + 1 - x1 - x2) % p; y3 = (lam * (x1 - x3) - y1 - 1) % p
    return (x3, y3)
def order(P, p):
    k = 1; Q = P
    while Q is not None: Q = add(Q, P, p); k += 1
    return k
def ord_mod(B, m):
    if m == 1: return 1
    k, x = 1, B % m
    while x != 1: x = (x * B) % m; k += 1
    return k
print("   the doubling map P -> 2P on 143a1(F_p): every cycle length is ord_m(2) for m the odd part of the point's order -- the period clock of Instrument 14 with base 2:")
for p in (7, 17, 19, 23):                                                                        # good primes (11 and 13 are bad: the fibre is singular)
    pts = curve_points(p); n = len(pts); idx = {P: i for i, P in enumerate(pts)}
    f = [idx[add(P, P, p)] for P in pts]
    assert add(None, None, p) is None
    cyc_len = {}
    for i in range(n):
        seen = []; j = i
        while j not in seen: seen.append(j); j = f[j]
        if j == i: cyc_len[i] = len(seen)
    report = {}
    for i, l in cyc_len.items():
        m = order(pts[i], p) if pts[i] is not None else 1
        while m % 2 == 0: m //= 2
        assert l == ord_mod(2, m), (p, pts[i], l, m); report[m] = l
    print("      p = %2d: #E = %2d; points of odd order m -> cycle length ord_m(2): %s" % (p, n, ", ".join("m = %d: %d" % (m, l) for m, l in sorted(report.items()))))

# ---------------------------------------------------------------- 7. Siegel zeros are thin three-square shells
def class_number(D):
    hcount = 0
    for a_ in range(1, int(math.isqrt(-D // 3)) + 2):
        for b_ in range(-a_ + 1, a_ + 1):
            if (b_ * b_ - D) % (4 * a_): continue
            c_ = (b_ * b_ - D) // (4 * a_)
            if c_ < a_ or (c_ == a_ and b_ < 0): continue
            hcount += 1
    return hcount
def kronecker(D, n):                                                                              # the Kronecker symbol (D/n) for fundamental D
    res = 1
    for p_, e in sp.factorint(n).items():
        if p_ == 2: s = 0 if D % 2 == 0 else (1 if D % 8 in (1, 7) else -1)
        else: s = 0 if D % p_ == 0 else (1 if pow(D % p_, (p_ - 1) // 2, p_) == 1 else -1)
        res *= s ** e
    return res
def L1(D, N=200000):                                                                              # L(1, chi_D) by a smoothed partial sum
    s = 0.0
    for n in range(1, N + 1): s += kronecker(D, n) / n * (1 - n / N)
    return s
r3 = {5: 24, 13: 24, 21: 48, 11: 24, 19: 24}
for n in (5, 13, 21):                                                                             # n = 1, 2 mod 4 squarefree: r_3(n) = 12 h(-4n) = (24 sqrt n / pi) L(1, chi_{-4n})
    D = -4 * n; h = class_number(D); L = L1(D)
    assert r3[n] == 12 * h and abs(L - math.pi * h / (2 * math.sqrt(n))) < 2e-3, (n, h, L)
for n in (11, 19):                                                                                # n = 3 mod 8 squarefree: r_3(n) = 24 h(-n), L(1, chi_{-n}) = pi h / sqrt n (w = 2)
    D = -n; h = class_number(D); L = L1(D); assert r3[n] == 24 * h and abs(L - math.pi * h / math.sqrt(n)) < 2e-3, (n, h, L)
h143 = class_number(-143); L143 = L1(-143); assert abs(L143 - math.pi * h143 / math.sqrt(143)) < 3e-3
print("Siegel zeros as thin shells: r_3(n) = 12 h(-4n) = (24 sqrt n/pi) L(1, chi_{-4n}) for n = 5, 13, 21 and r_3(n) = 24 h(-n) for n = 11, 19 (Gauss + Dirichlet, L(1,chi) summed directly, 2e-3); a Siegel zero beta would make L(1,chi) ~ (1 - beta) log|d| small, hence a thin shell")
print("   143: h(-143) = %d, L(1, chi_{-143}) = %.4f, and r_3(143) = 0 anyway -- that empty shell is Gauss's (143 = 7 mod 8), not Siegel's" % (h143, L143))
print("   status: Yitang Zhang's 2022 preprint claims L(sigma, chi) != 0 for sigma > 1 - c/(log D)^2024 (a weak Landau-Siegel statement); as of my knowledge it is unverified and not accepted, so nothing here depends on it; Siegel's ineffective bound and Goldfeld-Gross-Zagier's effective h(d) >> log|d| are the theorems")
print("ALL CHECKS PASSED")
