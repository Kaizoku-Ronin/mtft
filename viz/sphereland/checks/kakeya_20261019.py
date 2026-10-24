"""Kakeya: Besicovitch's Perron tree in the plane, and the finite-field Kakeya sets of F_q^2.  2026-10-19.  EXACT unless marked.

(1) Perron tree: the triangle T (base [0,1], apex (1/2, 1)) cut into 2^k elementary triangles with the common apex; at level j the right
half of every block of 2^j is translated left by (1 - alpha) alpha^{j-1} 2^j b0 and block m by m (1 - alpha) alpha^{j-1} 2^j b0, so that the
'hearts' of each level are the elementary triangles of alpha^j T (Falconer, The Geometry of Fractal Sets, Lemma 7.1 and Theorem 7.2).
Area of the union <= [alpha^{2k} + 2 (1 - alpha)^2 (1 - alpha^{2k}) / (1 - alpha^2)] area(T), checked by rasterisation, and every
direction of T's range still has a unit segment in the union (the elementary triangles are only translated).  (2) F_q^2, q odd: K =
{(x, y) : y + x^2 is a square} together with the line x = 0 contains a line in every direction, has q(q+1)/2 + (q-1)/2 points — the
minimum (Blokhuis-Mazzocca 2008) — and meets Dvir's bound |K| >= C(q + n - 1, n) = q(q+1)/2 for n = 2 up to the (q-1)/2."""
import math, itertools
import numpy as np

# ---------------------------------------------------------------- 1. the Perron tree
def perron_offsets(k, alpha):
    n = 2 ** k; b0 = 1.0 / n; off = np.zeros(n)
    for j in range(1, k + 1):
        bj = 2 ** j * b0; scale = alpha ** (j - 1)
        for i in range(n):
            block, pos = divmod(i, 2 ** j)
            if pos >= 2 ** (j - 1): off[i] -= (1 - alpha) * scale * bj           # the right half of the block, as in the lemma
            off[i] -= block * (1 - alpha) * scale * bj                              # re-align the hearts to the elementary triangles of alpha^j T
    return off
def perron_triangles(k, alpha):
    n = 2 ** k; b0 = 1.0 / n; off = perron_offsets(k, alpha); tris = []
    for i in range(n):
        tris.append([(i * b0 + off[i], 0.0), ((i + 1) * b0 + off[i], 0.0), (0.5 + off[i], 1.0)])
    return tris
from shapely.geometry import Polygon
from shapely.ops import unary_union
def union_area(tris):                                                                    # exact polygon union
    return unary_union([Polygon(t) for t in tris]).area
def bound(k, alpha): return (alpha ** (2 * k) + 2 * (1 - alpha) ** 2 * (1 - alpha ** (2 * k)) / (1 - alpha ** 2)) * 0.5
assert abs(union_area(perron_triangles(0, 0.8)) - 0.5) < 2e-3
rows = []
for alpha in (0.7, 0.8, 0.9):
    prev = 1.0
    for k in (1, 2, 3, 4, 5, 6, 7):
        a = union_area(perron_triangles(k, alpha)); b = bound(k, alpha)
        assert a <= b * (1 + 1e-9) + 1e-12, (alpha, k, a, b)
        rows.append((alpha, k, a, b)); prev = a
print("Perron tree, area of the union of the 2^k translated elementary triangles against the bound (area(T) = 0.5):")
for alpha, k, a, b in rows:
    if k in (1, 3, 5, 7): print("   alpha = %.1f  k = %d  area = %.4f  bound = %.4f" % (alpha, k, a, b))
# the areas decrease in k for fixed alpha beyond the first levels, and the bound tends to 2(1-alpha)/(1+alpha) * 0.5
for alpha in (0.7, 0.8, 0.9):
    seq = [a for al, k, a, b in rows if al == alpha]; assert seq[-1] < seq[0] and seq[-1] < 0.4 * 0.5, seq
# directions: every direction between the edges of T has a segment of length >= 1 in the union (apex to base of its elementary triangle)
def point_in_tri(p, t):
    (x1, y1), (x2, y2), (x3, y3) = t; x, y = p
    d1 = (x - x2) * (y1 - y2) - (x1 - x2) * (y - y2); d2 = (x - x3) * (y2 - y3) - (x2 - x3) * (y - y3); d3 = (x - x1) * (y3 - y1) - (x3 - x1) * (y - y1)
    return not ((d1 < -1e-12 or d2 < -1e-12 or d3 < -1e-12) and (d1 > 1e-12 or d2 > 1e-12 or d3 > 1e-12))
rng = np.random.default_rng(1); tris = perron_triangles(6, 0.8); n = len(tris); b0 = 1 / n
for _ in range(300):
    s = rng.uniform(0, 1); i = min(n - 1, int(s * n))                                   # the base point s of T, in elementary triangle i
    apex = tris[i][2]; foot = (tris[i][0][0] + (s - i * b0), 0.0)                       # same direction as in T: apex -> (s, 0)
    seg = [(apex[0] + (foot[0] - apex[0]) * u, apex[1] + (foot[1] - apex[1]) * u) for u in np.linspace(0, 1, 50)]
    assert all(point_in_tri(p, tris[i]) for p in seg)
    assert math.hypot(foot[0] - apex[0], foot[1] - apex[1]) >= 1.0
print("   every direction in T's range keeps a segment of length >= 1 inside its translated elementary triangle (300 random directions, k = 6)")

# ---------------------------------------------------------------- 2. finite-field Kakeya sets in F_q^2
def kakeya_fq2(q):
    sq = {(x * x) % q for x in range(q)}
    K = {(x, y) for x in range(q) for y in range(q) if (y + x * x) % q in sq} | {(0, y) for y in range(q)}
    return K
def contains_line_every_direction(K, q):
    inv2 = pow(2, -1, q)
    for m in range(q):                                                                   # direction (1, m): the line y = m x + m^2/4 lies in K
        c = (m * m * inv2 * inv2) % q
        if not all(((x, (m * x + c) % q) in K) for x in range(q)): return False
    return all(((0, y) in K) for y in range(q))                                          # direction (0, 1): the line x = 0
for q in (5, 7, 11, 13):
    K = kakeya_fq2(q); size = len(K); formula = q * (q + 1) // 2 + (q - 1) // 2; dvir = math.comb(q + 1, 2)
    assert contains_line_every_direction(K, q) and size == formula and size == dvir + (q - 1) // 2
    print("F_%d^2: the parabola set has %d points = q(q+1)/2 + (q-1)/2 (the minimum, Blokhuis-Mazzocca) against Dvir's bound C(q+1,2) = %d and q^2 = %d; a line in each of the q + 1 directions" % (q, size, dvir, q * q))
# exhaustive check of minimality for q = 3 (small enough to enumerate): no Kakeya set with fewer than 3*4/2 + 1 = 7 points
q = 3; pts = [(x, y) for x in range(q) for y in range(q)]
lines = [[(x, (m * x + c) % q) for x in range(q)] for m in range(q) for c in range(q)] + [[(c, y) for y in range(q)] for c in range(q)]
dirs = [(m, c) for m in range(q) for c in range(q)]
def is_kakeya(S):
    S = set(S)
    for m in range(q):
        if not any(all(((x, (m * x + c) % q) in S) for x in range(q)) for c in range(q)): return False
    return any(all(((c, y) in S) for y in range(q)) for c in range(q))
best = min(r for r in range(1, 10) if any(is_kakeya(S) for S in itertools.combinations(pts, r)))
assert best == 7 == q * (q + 1) // 2 + (q - 1) // 2
print("F_3^2: exhaustive search gives the minimum Kakeya set size 7 = 3*4/2 + 1, matching the formula")
print("Dvir (2008): every Kakeya set in F_q^n has at least C(q + n - 1, n) >= q^n / n! points; in R^n the measure can be zero (Besicovitch) and the dimension conjecture is proven for n = 2 (Davies 1971) and n = 3 (Wang-Zahl 2025), open for n >= 4")
print("ALL CHECKS PASSED")
