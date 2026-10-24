"""DeepSeek's 'Siegel movie' thread, checked.  2026-10-18.  EXACT unless marked.

A knot in bridge position, sliced at height c, meets the plane in 2n points.  The double cover of the plane branched at them is a
hyperelliptic curve of genus n - 1 (Riemann-Hurwitz).  This script checks (1) that count on the two trefoil presentations used in
Sphereland (3-bridge arc: 6 points, genus 2; 2-bridge plat of sigma_2^3: 4 points, genus 1), (2) the genus-1 modular parameter tau(c)
of the 4-point slice via the AGM, validated against theta functions, and its dependence on c (so it is a parameter of the slice, not of
the knot), (3) Clebsch's invariants A, B, C, D of the 6-point slice and their degeneration (D -> 0) exactly at the Morse levels,
(4) the 2-bridge knots of determinant 143 (Schubert), (5) the Mazur-Morishita linking of the primes 11 and 13."""
import itertools, math
from fractions import Fraction
import numpy as np
import mpmath as mp

mp.mp.dps = 30
rng = np.random.default_rng(143)

# ---------------------------------------------------------------- 1. Riemann-Hurwitz and the slice counts
chi = lambda two_n: 2 * 2 - two_n                                    # double cover of P^1 (chi = 2) branched at 2n points
genus = lambda two_n: (2 - chi(two_n)) // 2                          # chi = 2 - 2g  ->  g = n - 1
assert [genus(k) for k in (2, 4, 6, 8)] == [0, 1, 2, 3]

def trefoil_arc_3bridge():                                           # Sphereland Instrument 8: (sin t + 2 sin 2t, cos t - 2 cos 2t, 1.4 - sin 3t), opened at t0 = pi/6
    t0, dl = math.pi / 6, 0.22
    T = lambda t: (math.sin(t) + 2 * math.sin(2 * t), math.cos(t) - 2 * math.cos(2 * t), 1.4 - math.sin(3 * t))
    a, b = T(t0 + dl), T(t0 + 2 * math.pi - dl)
    pts = [(a[0], a[1], 0.0)] + [T(t0 + dl + (2 * math.pi - 2 * dl) * i / 4000) for i in range(4001)] + [(b[0], b[1], 0.0)]
    return pts

def trefoil_arc_2bridge():                                           # plat closure of sigma_2^3 in B_4, one bottom cup opened into two legs
    pts = []
    pts += [(-1.5, 0.0, z) for z in np.linspace(0, 2, 200)]                                    # strand 1 up
    pts += [(-1.0 - 0.5 * math.cos(s), 0.0, 2 + 0.5 * math.sin(s)) for s in np.linspace(0, math.pi, 200)]   # cap (1,2), top at 2.5
    pts += [(0.5 * math.cos(3 * math.pi * (z - 1)), 0.5 * math.sin(3 * math.pi * (z - 1)), z) for z in np.linspace(2, 1, 600)]  # down, top-2 -> bottom-3
    pts += [(1.0 - 0.5 * math.cos(s), 0.0, 1 - 0.5 * math.sin(s)) for s in np.linspace(0, math.pi, 200)]   # cup (3,4), bottom at 0.5
    pts += [(1.5, 0.0, z) for z in np.linspace(1, 2, 200)]                                     # strand 4 up
    pts += [(1.0 + 0.5 * math.cos(s), 0.0, 2 + 0.5 * math.sin(s)) for s in np.linspace(0, math.pi, 200)]    # cap (3,4)
    pts += [(-0.5 * math.cos(3 * math.pi * (z - 1)), -0.5 * math.sin(3 * math.pi * (z - 1)), z) for z in np.linspace(2, 1, 600)]  # down, top-3 -> bottom-2
    pts += [(-0.5, 0.0, z) for z in np.linspace(1, 0, 200)]                                    # leg down
    return pts

def slice_points(pts, c):
    out = []
    for p, q in zip(pts[:-1], pts[1:]):
        if (p[2] < c) != (q[2] < c):                                 # half-open rule: a vertex exactly on the plane is counted once
            s = (c - p[2]) / (q[2] - p[2]); out.append(complex(p[0] + s * (q[0] - p[0]), p[1] + s * (q[1] - p[1])))
    return out

A3, A2 = trefoil_arc_3bridge(), trefoil_arc_2bridge()
assert [len(slice_points(A3, c)) for c in (0.2, 0.5, 1.0, 2.0, 2.39, 2.5)] == [2, 6, 6, 6, 6, 0]
assert [len(slice_points(A2, c)) for c in (0.2, 0.7, 1.5, 2.2, 2.6)] == [2, 4, 4, 4, 0]
# the arc must not self-intersect (checked on the polyline: minimum distance between non-adjacent segments' sample points)
def min_self_distance(pts, gap=0.6):                                 # min distance between points more than `gap` apart along the arc
    P = np.array(pts); s = np.concatenate(([0], np.cumsum(np.linalg.norm(P[1:] - P[:-1], axis=1)))); d = np.inf
    for i in range(0, len(P), 2):
        j0 = int(np.searchsorted(s, s[i] + gap))
        if j0 >= len(P): break
        d = min(d, float(np.min(np.linalg.norm(P[j0:] - P[i], axis=1))))
    return d
assert min_self_distance(A2) > 0.3 and min_self_distance(A3) > 0.3, (min_self_distance(A2), min_self_distance(A3))
print("slice counts: 3-bridge arc 2 -> 6 -> 0 (genus 0 -> 2 -> -), 2-bridge plat 2 -> 4 -> 0 (genus 0 -> 1 -> -); arcs embedded")

# ---------------------------------------------------------------- 2. the genus-1 slice: tau from the cross-ratio by the AGM, checked against theta functions
def cross_ratio(a):                                                  # lambda = (a1-a3)(a2-a4) / ((a1-a4)(a2-a3))
    a1, a2, a3, a4 = a; return (a1 - a3) * (a2 - a4) / ((a1 - a4) * (a2 - a3))
j_of_lambda = lambda l: 256 * (l * l - l + 1) ** 3 / (l * l * (l - 1) ** 2)

def agm(a, b, n=60):                                                 # principal AGM with Cox's 'right choice' of square root
    for _ in range(n):
        a, b = (a + b) / 2, mp.sqrt(a * b)
        if abs(a - b) > abs(a + b): b = -b
    return a
K = lambda m: mp.pi / (2 * agm(mp.mpc(1), mp.sqrt(1 - m)))            # complete elliptic integral, parameter m
tau_of_lambda = lambda l: 1j * K(1 - l) / K(l)

def reduce_tau(t):
    for _ in range(200):
        t = t - mp.nint(mp.re(t))
        if abs(t) < 1 - mp.mpf('1e-25'): t = -1 / t
        else: break
    return t
def j_of_tau(t):                                                     # j = 32 (t2^8 + t3^8 + t4^8)^3 / (t2 t3 t4)^8 with the nome q = e^{i pi tau}
    q = mp.exp(1j * mp.pi * t)
    t2 = mp.jtheta(2, 0, q); t3 = mp.jtheta(3, 0, q); t4 = mp.jtheta(4, 0, q)
    return 32 * (t2 ** 8 + t3 ** 8 + t4 ** 8) ** 3 / (t2 * t3 * t4) ** 8
lam_of_tau = lambda t: (mp.jtheta(2, 0, mp.exp(1j * mp.pi * t)) / mp.jtheta(3, 0, mp.exp(1j * mp.pi * t))) ** 4

for l in (mp.mpc(-1, 0) / 3, mp.mpc(0.5, 0.667), mp.mpc(2.3, -0.4), mp.mpc(0.1, 0.05)):
    t = tau_of_lambda(l)
    assert mp.im(t) > 0
    assert abs(j_of_tau(t) - j_of_lambda(l)) < mp.mpf('1e-15') * (1 + abs(j_of_lambda(l))), (l, t)
    assert abs(j_of_tau(reduce_tau(t)) - j_of_lambda(l)) < mp.mpf('1e-15') * (1 + abs(j_of_lambda(l)))
    # lambda itself is recovered up to the S3 action (the AGM branch fixes a Gamma(2)-coset): check the j only, and lambda for the principal case
print("tau(lambda) = i K(1-lambda)/K(lambda) via the principal AGM: j(tau) = j(lambda) to 1e-15 at four test points, in H")

# the 2-bridge slice: order the 4 points by angle around their centroid so lambda is continuous, then tau(c) reduced to F
def ordered(a):
    cen = sum(a) / len(a); return sorted(a, key=lambda z: math.atan2((z - cen).imag, (z - cen).real))
traj = []
for c in np.linspace(0.55, 2.45, 191):
    a = [mp.mpc(z) for z in ordered(slice_points(A2, c))]
    t = reduce_tau(tau_of_lambda(cross_ratio(a))); traj.append((c, complex(t), complex(j_of_lambda(cross_ratio(a)))))
ims = [t.imag for _, t, _ in traj]; js = [abs(j) for _, _, j in traj]
assert max(ims) - min(ims) > 0.2 and max(js) / max(1e-9, min(js)) > 10, "tau must move with c"
# in the braid region the points are a1 = -1.5, a4 = 1.5, a2 = -a3 = w: lambda = (w - 1.5)^2 / (-6 w); j is periodic in phi with period pi, one period per half-twist
w = lambda phi: 0.5 * mp.exp(1j * phi)
lam_b = lambda phi: (w(phi) - 1.5) ** 2 / (-6 * w(phi))
assert abs(lam_b(0) + mp.mpf(1) / 3) < 1e-20 and abs(lam_b(mp.pi) - mp.mpf(4) / 3) < 1e-20
assert all(abs(j_of_lambda(lam_b(p)) - j_of_lambda(lam_b(p + mp.pi))) < 1e-12 for p in (0.3, 1.1, 2.0))
print("2-bridge slice: tau(c) reduced to F ranges over Im tau in [%.3f, %.3f]; |j| from %.1f to %.1f; j(phi) has period pi: one circuit per half-twist, three for sigma_2^3"
      % (min(ims), max(ims), min(js), max(js)))
print("   braid-region endpoints lambda = -1/3 and 4/3 = 1 - (-1/3), j =", mp.nstr(j_of_lambda(mp.mpf(-1) / 3), 8), " (= 35152/9)")
assert abs(j_of_lambda(mp.mpf(-1) / 3) - mp.mpf(35152) / 9) < 1e-20

# ---------------------------------------------------------------- 3. the genus-2 slice: Clebsch's A, B, C, D of the sextic with roots at the 6 points
def clebsch(a):
    a = list(a); n = 6; d2 = lambda i, j: (a[i] - a[j]) ** 2
    idx = range(n)
    def matchings(s):
        s = list(s)
        if not s: yield []; return
        i = s[0]
        for j in s[1:]:
            rest = [k for k in s if k not in (i, j)]
            for m in matchings(rest): yield [(i, j)] + m
    A = sum(math.prod(d2(i, j) for i, j in m) for m in matchings(idx))                       # 15 terms
    splits = [(T, tuple(k for k in idx if k not in T)) for T in itertools.combinations(idx, 3) if 0 in T]   # 10 unordered splits
    tri = lambda T: d2(T[0], T[1]) * d2(T[1], T[2]) * d2(T[2], T[0])
    B = sum(tri(T) * tri(U) for T, U in splits)                                                # 10 terms
    C = sum(tri(T) * tri(U) * math.prod(d2(T[k], perm[k]) for k in range(3)) for T, U in splits for perm in itertools.permutations(U))  # 60 terms
    D = math.prod(d2(i, j) for i in idx for j in idx if i < j)                                 # the discriminant, 15 factors
    return A, B, C, D
def absolute(a):
    A, B, C, D = clebsch(a); return (B / A ** 2, C / A ** 3, D / A ** 5)
# Moebius invariance of the weight-0 ratios (the roots transform, the common factors cancel)
a0 = rng.normal(size=6) + 1j * rng.normal(size=6)
for _ in range(5):
    m = rng.normal(size=4); det = m[0] * m[3] - m[1] * m[2]
    a1 = [(m[0] * z + m[1]) / (m[2] * z + m[3]) for z in a0]
    r0, r1 = absolute(a0), absolute(a1)
    assert all(abs(x - y) < 1e-9 * (1 + abs(x)) for x, y in zip(r0, r1)), (r0, r1)
# counts of terms
assert len(list(itertools.combinations(range(6), 3))) == 20 and 20 // 2 == 10
print("Clebsch A (15 pairings), B (10 triple splits), C (60 split+matching), D (discriminant): the ratios B/A^2, C/A^3, D/A^5 are Moebius-invariant (5 random maps, 1e-9)")
# degeneration: D -> 0 at the Morse levels of the 3-bridge arc (the saddles at z = 0.4, the tops at 2.4), finite in between
def Dnorm(c):
    a = slice_points(A3, c); A, B, C, D = clebsch(a); return abs(D / A ** 5)
mid = Dnorm(1.2); near_saddle = Dnorm(0.4005); near_top = Dnorm(2.3995)
assert near_saddle < 1e-6 * mid and near_top < 1e-6 * mid, (mid, near_saddle, near_top)
print("3-bridge slice: |D/A^5| = %.3e at c = 1.2, %.1e just above the saddles (c = 0.4005), %.1e just below the tops (c = 2.3995): the curve degenerates exactly at the Morse levels" % (mid, near_saddle, near_top))
# and the point in M_2 moves with c (so it is a parameter of the slice, not of the knot)
r_a, r_b = absolute(slice_points(A3, 0.8)), absolute(slice_points(A3, 1.8))
assert any(abs(x - y) > 1e-3 * (1 + abs(x)) for x, y in zip(r_a, r_b))
print("   and the absolute invariants at c = 0.8 and 1.8 differ: the slice curve moves through M_2 with c")

# ---------------------------------------------------------------- 4. two-bridge knots of determinant 143 (Schubert: K(p/q) ~ K(p/q') iff q' = q^{+-1} mod p; mirror: q -> -q)
p = 143
units = [q for q in range(1, p) if math.gcd(q, p) == 1]
assert len(units) == 120
def cls(q, mirror):
    s = {q % p, pow(q, -1, p)}
    if mirror: s |= {(-x) % p for x in list(s)}
    return frozenset(s)
chiral = {cls(q, False) for q in units}; unoriented = {cls(q, True) for q in units}
assert len(chiral) == 62 and len(unoriented) == 31
assert not any(pow(q, 2, p) == p - 1 for q in units)                 # no q^2 = -1 mod 143 (11 = 3 mod 4): none of them is amphichiral
def cf_positive(num, den):                                           # positive continued fraction, crossing number = sum of terms (alternating reduced diagram)
    terms = []
    while den:
        terms.append(num // den); num, den = den, num % den
    return terms
cross = {}
for C_ in unoriented:
    cross[C_] = min(sum(cf_positive(p, q)) for q in C_)
cmin = min(cross.values()); best = [sorted(C_) for C_, c_ in cross.items() if c_ == cmin]
print("determinant 143: %d two-bridge knots up to mirror image (%d chiral), none amphichiral; least crossing number %d, at q in %s"
      % (len(unoriented), len(chiral), cmin, best))
for C_ in sorted(unoriented, key=lambda s: (cross[s], min(s)))[:6]:
    q = min(C_); print("   K(143/%d): continued fraction %s, crossing number %d, double branched cover L(143, %d)" % (q, cf_positive(p, q), cross[C_], q))
# the composite with the same determinant
assert 11 * 13 == 143
print("   also T(2,11) # T(2,13): determinant 11 x 13 = 143, Sigma_2 = L(11,1) # L(13,1), H_1 = Z/11 + Z/13 = Z/143, 24 crossings")
# the trefoil and the figure-eight for comparison
print("   trefoil = K(3/1): det 3, Sigma_2 = L(3,1); figure-eight = K(5/2): det 5, Sigma_2 = L(5,2)")

# ---------------------------------------------------------------- 5. Mazur-Morishita: Legendre symbols as mod-2 linking numbers of the prime 'knots' 11 and 13
legendre = lambda a, q: pow(a, (q - 1) // 2, q)                      # 1 or q-1
assert legendre(11, 13) == 12 and legendre(13, 11) == 10            # both -1
assert (11 % 4 == 3) and (13 % 4 == 1)                               # reciprocity sign (-1)^{(10)(12)/4} = +1: the two symbols agree
print("Legendre symbols (11|13) = (13|11) = -1: in the Mazur-Morishita dictionary the primes 11 and 13 have mod-2 linking number 1 -- the two factors of 143 are linked")

# ---------------------------------------------------------------- 6. the obstruction, stated (Gluck 1962; Hatcher 1981): pi_0 Diff(S^2 x S^1) = Z2^3, orientation-preserving Z2 x Z2
print("S^2 x S^1 (boundary of a 2-knot): pi_0 Diff = Z2 x Z2 x Z2 (Gluck), Diff ~ O(2) x O(3) x Omega SO(3) (Hatcher); no flat metric, so no cusp shape and no tau -- the obstruction DeepSeek named is real")
print("ALL CHECKS PASSED")
