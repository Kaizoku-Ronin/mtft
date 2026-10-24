"""The trefoil complement as the modular group's unit tangent bundle; modular knots and their Rademacher linking.  2026-10-18.

EXACT unless marked.  Milnor: the space of unimodular lattices SL(2,R)/SL(2,Z), mapped by (g2, g3) to C^2 minus the discriminant
{g2^3 = 27 g3^2} and normalised along the weighted scaling (t^4 g2, t^6 g3) to S^3, is S^3 minus a trefoil.  Ghys: the periodic orbits of
the geodesic flow (closed geodesics of the modular surface, lifted) are Lorenz knots whose linking number with the trefoil is Rademacher's
function psi(A).  This script builds the picture numerically and checks: the fibration (regular fibres are (2,3) curves linking the trefoil
6 times, the exceptional fibres over i and rho link it 3 and 2 times), the modular knots for several words in L = [[1,0],[1,1]] and
R = [[1,1],[0,1]] with their Gauss linking integrals against psi(A) from Dedekind sums (psi = #R - #L on words), and the arithmetic of
Gamma_0(143): torsion-free, genus 13, four cusps, index 168, its trace-4 element and that element's Lorenz word."""
import math, itertools
from fractions import Fraction
import numpy as np
import sympy as sp

L = np.array([[1, 0], [1, 1]]); R = np.array([[1, 1], [0, 1]])
word_matrix = lambda w: __import__('functools').reduce(lambda M, ch: M @ (L if ch == 'L' else R), w, np.eye(2, dtype=int))

# ---------------------------------------------------------------- Eisenstein series and the lattice map
def E4E6(tau, N=40):
    q = np.exp(2j * np.pi * tau); e4 = 1; e6 = 1
    for n in range(1, N):
        s3 = sum(d ** 3 for d in range(1, n + 1) if n % d == 0); s5 = sum(d ** 5 for d in range(1, n + 1) if n % d == 0)
        e4 += 240 * s3 * q ** n; e6 -= 504 * s5 * q ** n
    return e4, e6
def reduce_tau(tau):                       # returns tau' in F and (c, d) of the matrix used, so that omega1' = omega1 (c tau + d)
    a, b, c, d = 1, 0, 0, 1                # track M = [[a,b],[c,d]] with tau' = M tau
    for _ in range(200):
        n = round(tau.real); tau = tau - n; a, b = a - n * c, b - n * d
        if abs(tau) < 1 - 1e-13: tau = -1 / tau; a, b, c, d = -c, -d, a, b
        else: break
    return tau, c, d
def g2g3_lattice(omega1, tau):             # g2, g3 of the lattice omega1 (Z + tau Z)
    t, c, d = reduce_tau(tau); w1 = omega1 * (c * tau + d)
    e4, e6 = E4E6(t)
    return w1 ** -4 * (4 * np.pi ** 4 / 3) * e4, w1 ** -6 * (8 * np.pi ** 6 / 27) * e6
def weighted_normalize(g2, g3):            # the point of the orbit (t^4 g2, t^6 g3) on |u|^2 + |v|^2 = 1
    lo, hi = 0.0, 1.0
    f = lambda t: t ** 8 * abs(g2) ** 2 + t ** 12 * abs(g3) ** 2 - 1
    while f(hi) < 0: hi *= 2
    for _ in range(200):
        mid = (lo + hi) / 2
        if f(mid) < 0: lo = mid
        else: hi = mid
    t = (lo + hi) / 2; return t ** 4 * g2, t ** 6 * g3
KAPPA = None
def display(u, v):                         # (u, v) on S^3 -> stretched so the trefoil sits on the Clifford torus, then back to S^3
    w = np.array([u.real, u.imag, KAPPA * v.real, KAPPA * v.imag]); return w / np.linalg.norm(w)
POLE = None
def stereo(x):                              # S^3 subset R^4 -> R^3 from the pole POLE (rotate POLE to e4 first)
    e4 = POLE; basis = [e4]
    for e in np.eye(4):
        v = e - sum(np.dot(e, b) * b for b in basis)
        if np.linalg.norm(v) > 1e-6: basis.append(v / np.linalg.norm(v))
        if len(basis) == 4: break
    Q = np.array(basis[1:] + [basis[0]])       # rows: new coordinates, last = pole direction
    y = Q @ x; return y[:3] / (1 - y[3])

r_tref = math.sqrt(sp.nsolve(sp.Symbol('x') ** 3 + 9 * sp.Symbol('x') ** 2 - 1, 0.33))     # 9 r^4 + r^6 = 1
KAPPA = 3 / r_tref                            # puts |u| = |v| for the trefoil
POLE = np.array([0.25 * math.cos(0.4), 0.25 * math.sin(0.4), 0.968246 * math.cos(2.1), 0.968246 * math.sin(2.1)]); POLE /= np.linalg.norm(POLE)   # display coordinates, near the rho-core
assert abs(9 * r_tref ** 4 + r_tref ** 6 - 1) < 1e-12
print("trefoil radius r = %.6f (9 r^4 + r^6 = 1): |u| = %.4f, |v| = %.4f on S^3; display stretch kappa = %.4f" % (r_tref, 3 * r_tref ** 2, r_tref ** 3, KAPPA))

def curve_trefoil(n=600):
    return np.array([stereo(display(3 * r_tref ** 2 * np.exp(2j * p), r_tref ** 3 * np.exp(3j * p))) for p in np.linspace(0, 2 * np.pi, n, endpoint=False)])
def curve_fiber(tau, n=600, mult=1):          # the fibre over tau: lattices e^{i theta}(Z + tau Z)/sqrt(y); over i (rho) the lattice repeats after pi/2 (pi/3): multiplicity 2 (3)
    y = tau.imag; pts = []
    for th in np.linspace(0, np.pi / mult, n, endpoint=False):
        g2, g3 = g2g3_lattice(np.exp(1j * th) / math.sqrt(y), tau); pts.append(stereo(display(*weighted_normalize(g2, g3))))
    return np.array(pts)
def curve_modular_knot(A, n=800):             # the periodic orbit of the geodesic flow a(t) = diag(e^{t/2}, e^{-t/2}) on lattices g Z^2 with g^{-1} a(l) g = A
    tr = A[0, 0] + A[1, 1]; lam = (tr + math.sqrt(tr * tr - 4)) / 2; ell = 2 * math.log(lam)
    a, b, c, d = A.ravel()
    vp = np.array([lam - d, c], float) if abs(c) > 1e-12 else np.array([1.0, 0.0])            # A v = lam v
    vm = np.array([1 / lam - d, c], float) if abs(c) > 1e-12 else np.array([b / (1 / lam - a), 1.0])
    ginv = np.column_stack([vp, vm]); ginv = ginv / math.sqrt(abs(np.linalg.det(ginv)))
    if np.linalg.det(ginv) < 0: ginv[:, 1] *= -1
    g = np.linalg.inv(ginv); pts = []; taus = []
    for t in np.linspace(0, ell, n, endpoint=False):
        M = np.diag([math.exp(t / 2), math.exp(-t / 2)]) @ g
        w1 = complex(M[0, 0], M[1, 0]); w2 = complex(M[0, 1], M[1, 1]); tau = w2 / w1
        if tau.imag < 0: w1, w2 = w2, w1; tau = w2 / w1
        g2, g3 = g2g3_lattice(w1, tau); pts.append(stereo(display(*weighted_normalize(g2, g3)))); taus.append(reduce_tau(tau)[0])
    return np.array(pts), ell, np.array(taus)
def linking(P, Q, chunk=1500):                # Gauss linking integral of two closed polylines (midpoint rule, chunked)
    dP = np.roll(P, -1, 0) - P; dQ = np.roll(Q, -1, 0) - Q; mP = P + dP / 2; mQ = Q + dQ / 2; total = 0.0
    for i in range(0, len(P), chunk):
        r = mP[i:i + chunk, None, :] - mQ[None, :, :]; n3 = np.linalg.norm(r, axis=2) ** 3
        cr = np.cross(dP[i:i + chunk, None, :], dQ[None, :, :]); total += float(np.sum(np.einsum('ijk,ijk->ij', r, cr) / n3))
    return total / (4 * np.pi)

K = curve_trefoil(2400)
assert min(np.linalg.norm(K, axis=1)) > 0.05 and max(np.linalg.norm(K, axis=1)) < 50
# the Seifert fibration: regular fibre is a (2,3) curve linking the trefoil 6 times; the cores over i and rho link it 3 and 2 times
F_gen = curve_fiber(0.25 + 1.2j); F_i = curve_fiber(1j, mult=2); F_rho = curve_fiber(np.exp(2j * np.pi / 3), mult=3)
lk_gen, lk_i, lk_rho = linking(F_gen, K), linking(F_i, K), linking(F_rho, K)
assert abs(abs(lk_gen) - 6) < 0.02 and abs(abs(lk_i) - 3) < 0.02 and abs(abs(lk_rho) - 2) < 0.02, (lk_gen, lk_i, lk_rho)
assert abs(abs(linking(F_i, F_rho)) - 1) < 0.02                               # the two cores form a Hopf link
# the theta-parametrised fibre over i (rho) runs over its core twice (three times): the Seifert multiplicities
assert abs(abs(linking(curve_fiber(1j), K)) - 6) < 0.02 and abs(abs(linking(curve_fiber(np.exp(2j * np.pi / 3)), K)) - 6) < 0.02
print("Seifert fibration of the trefoil complement: |lk(regular fibre, K)| = %.3f (6 = 2 x 3), |lk(core over i, K)| = %.3f (3), |lk(core over rho, K)| = %.3f (2); the cores form a Hopf link (%.3f); multiplicities 2 and 3 confirmed"
      % (abs(lk_gen), abs(lk_i), abs(lk_rho), abs(linking(F_i, F_rho))))

# ---------------------------------------------------------------- Rademacher's function by Dedekind sums, and on words
def dedekind_sum(h, k):
    saw = lambda x: 0 if x == int(x) else x - math.floor(x) - Fraction(1, 2)
    return sum(saw(Fraction(r, k)) * saw(Fraction(h * r, k)) for r in range(1, k))
def rademacher_psi(A):                        # psi(A) = Phi(A) - 3 sign(c (a + d)), Phi(A) = (a + d)/c - 12 sign(c) s(d, |c|)   (Rademacher-Grosswald; Ghys)
    a, b, c, d = (int(x) for x in A.ravel())
    if c == 0: return Fraction(b, d)
    phi = Fraction(a + d, c) - 12 * (1 if c > 0 else -1) * dedekind_sum(d, abs(c))
    return phi - 3 * (1 if c * (a + d) > 0 else -1)
WORDS = ['LR', 'LLR', 'LRR', 'LLLR', 'LLRR', 'LLLRR', 'LRLRR', 'LLRLR']
for w in WORDS:
    A = word_matrix(w); psi = rademacher_psi(A)
    assert psi.denominator == 1 and int(psi) == w.count('R') - w.count('L'), (w, psi)
print("Rademacher psi by Dedekind sums equals #R - #L on the words", WORDS)

# ---------------------------------------------------------------- Ghys: lk(modular knot, trefoil) = psi(A)
results = []
for w in WORDS:
    A = word_matrix(w); knot, ell, taus = curve_modular_knot(A, n=6000); lk = linking(knot, K)
    results.append((w, int(A[0, 0] + A[1, 1]), ell, int(rademacher_psi(A)), lk, max(taus.imag)))
    assert min(np.linalg.norm(knot, axis=1)) > 0.02 and max(np.linalg.norm(knot, axis=1)) < 80, w
signs = {round(lk / psi) for w, tr, ell, psi, lk, _ in results if psi != 0}
assert len(signs) == 1 and all(abs(lk - list(signs)[0] * psi) < 0.03 for w, tr, ell, psi, lk, _ in results), results
SIGN = list(signs)[0]
print("modular knots: lk(k_A, K) = %+d x psi(A) for every word (orientation convention of this construction); the knots stay at Im tau' <= %.2f in F" % (SIGN, max(r[5] for r in results)))
for w, tr, ell, psi, lk, _ in results:
    print("   %-6s trace %2d  length 2 arccosh(tr/2) = %.4f  psi = %+d  lk = %+.3f" % (w, tr, ell, psi, lk))

# ---------------------------------------------------------------- Gamma_0(143): torsion-free, genus 13, 4 cusps, index 168; its trace-4 element
legendre = lambda a, p: pow(a % p, (p - 1) // 2, p)
N = 143; mu = N * (1 + Fraction(1, 11)) * (1 + Fraction(1, 13))
nu2 = math.prod(1 + (1 if legendre(-1, p) == 1 else -1) for p in (11, 13))
nu3 = math.prod(1 + (1 if legendre(-3, p) == 1 else -1) for p in (11, 13))
cusps = 4; genus = 1 + mu / 12 - Fraction(nu2, 4) - Fraction(nu3, 3) - Fraction(cusps, 2)
assert mu == 168 and nu2 == 0 and nu3 == 0 and genus == 13
print("Gamma_0(143): index 168, no elliptic points (nu2 = nu3 = 0: torsion-free), 4 cusps, genus 13 -> its unit tangent bundle is a circle bundle over Y_0(143) and a 168-sheeted cover of the trefoil complement")
A143 = np.array([[84, -47], [143, -80]])
assert round(np.linalg.det(A143)) == 1 and A143[0, 0] + A143[1, 1] == 4 and A143[1, 0] % 143 == 0
# its SL(2,Z)-conjugacy class among the trace-4 words: solve P A = W P for integer P with det 1
def conjugator(A, W, bound=300):
    p = sp.symbols('p0:4'); P = sp.Matrix(2, 2, p); eqs = list(P * sp.Matrix(A.tolist()) - sp.Matrix(W.tolist()) * P)
    sol = sp.linsolve(eqs, p); (s,) = sol; free = sorted(set().union(*[e.free_symbols for e in s]), key=str)
    if len(free) != 2: return None
    for m in range(-bound, bound + 1):
        for n in range(-bound, bound + 1):
            vals = [e.subs({free[0]: m, free[1]: n}) for e in s]
            if all(v.is_integer for v in vals) and vals[0] * vals[3] - vals[1] * vals[2] == 1: return sp.Matrix(2, 2, vals)
    return None
word143 = None
for w in ('LLR', 'LRR'):
    P = conjugator(A143, word_matrix(w), bound=200)
    if P is not None: word143 = w; break
assert word143 is not None
ell143 = 2 * math.acosh(2)
print("the trace-4 element [[84,-47],[143,-80]] of Gamma_0(143) is conjugate in SL(2,Z) to the word %s: closed geodesic of length 2 arccosh 2 = %.4f (VE-01's shortest), Rademacher psi = %+d, the same modular knot as %s" % (word143, ell143, int(rademacher_psi(A143)), word143))
assert abs(curve_modular_knot(A143)[1] - ell143) < 1e-9
print("ALL CHECKS PASSED")
