"""SQ-01 (2026-10-08): the star quotients X0(pq)* for supersingular p < q — is X0(143) special?

X0(N)* = X0(N)/<W_p, W_q> is how X0(143) folds onto 143a1.  This wave computes its genus at every level N = pq with p < q supersingular (105 levels)
by two independent routes, and records what does and does not single out 143.

1. Two routes (EXACT).  Route 1 (PARI/GP, frozen in ROUTE1): the level-N newforms whose Atkin–Lehner signs are all +1; old forms add nothing, since
   every newform of supersingular prime level has w_p = -1 (Ogg).  Route 2 (pure Python, computed live): each involution's fixed points from class
   numbers of imaginary quadratic orders with Eichler symbols (Ogg's count, plus the Z[i] term for W_2), then the holomorphic Lefschetz formula
   tr(W | H^0(K)) = 1 - nu/2 and g* = (1/|W|) sum tr.  Route 2 covers all 105 levels and agrees with route 1 at all 99 that route 1 reached, and at the
   75 levels with three or four supersingular primes up to 700.
2. The census (EXACT).  X0(N)* is a sphere at 21 levels and an elliptic curve at 13 — 57, 58, 65, 77, 82, 91, 118, 123, 141, 142, 143, 145, 155 — and
   has genus >= 2 at the other 71.  Every elliptic quotient has analytic rank 1 (NUMERICAL, PARI lfunorderzero); odd rank is forced, since all signs
   +1 give Fricke sign +1 and root number -1.  X0(143) is not alone.
3. What singles out 143 (EXACT).  The twin-prime levels 15, 35, 143, 323, 899 have star genus 0, 0, 1, 5, 12: 143 is the only twin level with an
   elliptic quotient.  77 and 143 are the only consecutive-prime levels with one.  dim S2^new = p, the smaller prime, at 11 levels, and with p, q
   consecutive only at 143 (read literally, 'the larger prime' with consecutive primes picks 221).  The Galois-orbit pattern (1, 4, 6) occurs only at 143.
4. The nearest sibling (EXACT).  123 = 3 x 41 is the only other level of index 168; with 143 it is the only torsion-free level of genus 13; it has four
   cusps and a rank-1 elliptic quotient.  Its new space has dimension 7 (orbits 1, 1, 2, 3) against 143's 11 (1, 4, 6).
5. Moonshine cross-check (EXACT).  Every level in the scan whose star quotient is a sphere — 21 with two primes and 7 with three — is the order of an
   element of the Monster (GAP's character table library: 73 element orders); 143 is not, and among semiprime orders the converse fails at 57 and 93.
6. Denser factorizations (EXACT).  With three or four supersingular primes and N <= 700 the star quotient is an elliptic curve at 102, 114, 130, 138,
   174, 182, 190, 195, 210, 231, 238; each extra prime doubles the cusps and the Atkin–Lehner group, so these change the arena rather than extend it.
7. Flavour (STRUCTURAL).  The exact obstructions behind pieces IX and X do not use the level: DS-1's stack pattern, the odd-mode count (an odd B-stack
   mode costs 3(2g - 2) chiral families, at least 24 at every torus level and 72 at 143, OM-3), QL-1's dressing bound and CI-2's universality.  DS-4's
   tail cap and BL-01's Wilson-line numbers were computed on X0(143) x 143a1 only.
8. Pre-registration (STRUCTURAL).  A rule that selects the level must be stated before it is evaluated on the 13 torus levels, and passes only if it
   selects 143 and rejects the other twelve, 123 and 77 included.  The distinctions of item 3 were found after the fact and are not such a rule.
Provenance: scripts/star_quotient_scan.py (needs cypari2) reproduces ROUTE1, DENSE_ROUTE1 and RANKS; MONSTER_ELEMENT_ORDERS is
Set(OrdersClassRepresentatives(CharacterTable("M"))) in GAP 4.12 with CTblLib 1.3.7."""
from __future__ import annotations

import functools
import itertools
import math
from fractions import Fraction

SUPERSINGULAR = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 41, 47, 59, 71)
DENSE_BOUND = 700

# route 1, PARI/GP 2.17.2 via cypari2 (mfinit([N,2],0), mfsplit, mfatkineigenvalues): level-N Galois orbits as (dimension, w_p, w_q), p < q
ROUTE1 = {6: (), 10: (), 14: ((1,1,-1),), 15: ((1,1,-1),), 21: ((1,-1,1),), 22: (), 26: ((1,1,-1), (1,-1,1)), 33: ((1,1,-1),), 34: ((1,-1,1),),
          35: ((1,1,-1), (2,-1,1)), 38: ((1,1,-1), (1,-1,1)), 39: ((1,1,-1), (2,-1,1)), 46: ((1,1,-1),), 51: ((1,-1,1), (2,1,-1)),
          55: ((1,-1,1), (2,1,-1)), 57: ((1,1,1), (1,-1,1), (1,-1,1)), 58: ((1,1,1), (1,-1,1)), 62: ((1,-1,1), (2,1,-1)),
          65: ((1,1,1), (2,1,-1), (2,-1,1)), 69: ((1,-1,1), (2,1,-1)), 77: ((1,1,1), (1,1,-1), (1,-1,1), (2,-1,1)), 82: ((1,1,1), (2,-1,1)),
          85: ((1,1,-1), (2,1,1), (2,-1,1)), 87: ((2,-1,1), (3,1,-1)), 91: ((1,1,1), (1,-1,-1), (2,-1,1), (3,1,-1)), 93: ((2,1,1), (3,-1,1)),
          94: ((1,-1,1), (2,1,-1)), 95: ((3,-1,1), (4,1,-1)), 115: ((1,1,-1), (2,1,1), (4,-1,1)), 118: ((1,1,1), (1,1,-1), (1,-1,1), (1,-1,1)),
          119: ((4,-1,1), (5,1,-1)), 123: ((1,1,1), (1,-1,-1), (2,-1,1), (3,1,-1)), 133: ((2,1,1), (2,-1,1), (2,-1,-1), (3,1,-1)),
          141: ((1,1,1), (1,1,-1), (1,-1,1), (1,-1,1), (1,-1,-1), (2,1,-1)), 142: ((1,1,1), (1,1,-1), (1,1,-1), (1,-1,1), (1,-1,-1)),
          143: ((1,1,1), (4,-1,1), (6,1,-1)), 145: ((1,1,1), (2,-1,-1), (3,1,-1), (3,-1,1)), 155: ((1,1,1), (1,1,-1), (1,-1,-1), (4,1,-1), (4,-1,1)),
          161: ((1,-1,1), (2,1,1), (3,1,-1), (5,-1,1)), 177: ((2,1,1), (2,-1,1), (2,-1,-1), (3,1,-1)),
          187: ((1,1,-1), (1,-1,1), (2,-1,1), (2,-1,-1), (3,1,1), (4,1,-1)),
          203: ((1,-1,1), (1,-1,1), (1,-1,-1), (2,1,-1), (2,1,-1), (3,1,1), (5,-1,1)),
          205: ((1,1,-1), (1,-1,1), (1,-1,-1), (2,1,1), (2,-1,-1), (3,1,-1), (3,-1,1)), 209: ((1,-1,-1), (2,1,1), (5,-1,1), (7,1,-1)),
          213: ((1,-1,1), (2,1,1), (2,-1,1), (2,-1,-1), (4,1,-1)), 217: ((3,1,1), (3,-1,-1), (4,-1,1), (5,1,-1)),
          221: ((1,1,-1), (1,1,-1), (2,1,1), (2,1,-1), (2,1,-1), (3,-1,-1), (6,-1,1)), 235: ((1,1,-1), (1,1,-1), (1,-1,-1), (5,1,1), (7,-1,1)),
          247: ((2,-1,1), (3,-1,-1), (4,1,1), (5,1,-1), (5,-1,1)), 253: ((3,1,-1), (3,-1,-1), (5,1,1), (6,-1,1)),
          287: ((2,1,1), (2,-1,-1), (3,1,-1), (3,-1,1), (5,-1,1), (6,1,-1)), 295: ((3,1,1), (3,-1,-1), (6,-1,1), (7,1,-1)),
          299: ((2,1,1), (2,-1,1), (2,-1,1), (2,-1,1), (2,-1,-1), (3,-1,1), (10,1,-1)), 319: ((1,1,-1), (3,-1,-1), (4,1,1), (7,-1,1), (8,1,-1)),
          323: ((1,1,-1), (2,1,-1), (4,-1,-1), (5,1,1), (6,1,-1), (7,-1,1)),
          329: ((1,1,-1), (2,-1,1), (3,1,1), (3,1,-1), (3,-1,-1), (5,1,-1), (6,-1,1)), 341: ((2,-1,-1), (4,1,1), (8,1,-1), (11,-1,1)),
          355: ((1,-1,1), (4,1,1), (4,-1,-1), (6,-1,1), (8,1,-1)), 377: ((1,-1,-1), (2,1,-1), (5,1,1), (5,-1,-1), (7,1,-1), (9,-1,1)),
          391: ((2,1,1), (3,1,1), (3,-1,-1), (9,1,-1), (12,-1,1)), 403: ((2,-1,1), (6,-1,-1), (7,1,-1), (8,1,1), (8,-1,1)),
          413: ((2,1,-1), (3,1,-1), (5,1,1), (5,1,-1), (5,-1,-1), (9,-1,1)),
          437: ((1,1,-1), (1,-1,-1), (2,1,1), (2,-1,-1), (2,-1,-1), (5,1,1), (8,1,-1), (12,-1,1)),
          451: ((1,1,1), (5,1,1), (5,-1,-1), (10,-1,1), (12,1,-1)),
          493: ((1,-1,1), (1,-1,-1), (2,1,-1), (4,-1,1), (5,-1,-1), (6,-1,1), (8,1,-1), (10,1,1)),
          497: ((1,-1,-1), (2,-1,-1), (8,1,-1), (9,1,1), (15,-1,1)),
          517: ((1,1,-1), (1,-1,1), (1,-1,-1), (2,1,-1), (2,1,-1), (2,-1,1), (2,-1,-1), (3,1,-1), (3,-1,-1), (10,1,1), (10,-1,1)),
          527: ((2,1,-1), (5,-1,-1), (7,1,1), (11,1,-1), (16,-1,1)), 533: ((2,1,-1), (7,-1,-1), (8,1,-1), (11,1,1), (13,-1,1)),
          551: ((1,-1,-1), (1,-1,-1), (1,-1,-1), (1,-1,-1), (2,1,1), (3,1,1), (16,-1,1), (18,1,-1)),
          589: ((3,-1,1), (9,-1,1), (9,-1,-1), (10,1,1), (14,1,-1)), 611: ((1,-1,1), (5,-1,-1), (9,1,1), (14,1,-1), (18,-1,1)),
          649: ((1,-1,-1), (2,-1,-1), (4,-1,-1), (11,1,-1), (12,1,1), (17,-1,1)), 667: ((10,-1,-1), (12,1,1), (13,1,-1), (16,-1,1)),
          697: ((3,1,-1), (3,-1,-1), (7,-1,-1), (10,1,-1), (15,1,1), (15,-1,1)),
          713: ((1,-1,-1), (2,1,1), (5,1,1), (5,1,1), (9,-1,-1), (15,1,-1), (18,-1,1)),
          767: ((2,1,1), (3,-1,-1), (4,-1,-1), (10,1,1), (17,1,-1), (23,-1,1)), 779: ((2,-1,-1), (8,-1,-1), (11,1,1), (19,-1,1), (21,1,-1)),
          781: ((1,1,-1), (1,-1,-1), (2,-1,-1), (6,-1,-1), (12,1,-1), (15,1,1), (20,-1,1)),
          799: ((1,-1,1), (1,-1,-1), (2,-1,-1), (8,-1,-1), (12,1,1), (17,-1,1), (20,1,-1)), 893: ((12,-1,-1), (16,1,1), (18,1,-1), (23,-1,1)),
          899: ((1,1,1), (1,-1,1), (10,-1,-1), (11,1,1), (23,1,-1), (25,-1,1)), 923: ((1,1,-1), (10,-1,-1), (16,1,1), (18,1,-1), (26,-1,1)),
          943: ((1,1,-1), (13,-1,-1), (16,1,1), (19,1,-1), (24,-1,1)),
          1003: ((1,1,1), (1,-1,1), (1,-1,1), (1,-1,-1), (3,1,1), (4,1,1), (10,1,1), (16,-1,-1), (18,-1,1), (22,1,-1)),
          1081: ((17,-1,-1), (20,1,-1), (21,1,1), (25,-1,1)), 1121: ((16,-1,-1), (17,1,1), (26,-1,1), (28,1,-1)),
          1189: ((1,1,-1), (2,1,1), (3,1,1), (3,-1,-1), (15,-1,-1), (19,1,1), (24,1,-1), (26,-1,1)),
          1207: ((16,-1,-1), (22,1,1), (24,1,-1), (31,-1,1)), 1271: ((1,-1,-1), (13,-1,-1), (17,1,1), (33,1,-1), (37,-1,1)),
          1349: ((1,1,1), (1,1,-1), (19,1,1), (19,-1,-1), (32,-1,1), (33,1,-1)),
          1357: ((1,1,1), (2,1,-1), (5,1,-1), (18,1,-1), (22,-1,-1), (26,1,1), (31,-1,1)),
          1363: ((2,1,1), (2,1,-1), (4,-1,-1), (18,-1,-1), (24,1,1), (25,1,-1), (32,-1,1)),
          1457: ((1,-1,1), (2,1,1), (2,1,-1), (24,-1,-1), (26,1,1), (27,1,-1), (33,-1,1)),
          1633: ((3,1,-1), (6,-1,-1), (21,-1,-1), (27,1,-1), (33,1,1), (37,-1,1)), 1711: ((26,-1,-1), (28,1,1), (39,-1,1), (42,1,-1)),
          1829: ((1,-1,-1), (28,-1,-1), (34,1,1), (38,1,-1), (44,-1,1)), 1927: ((2,1,1), (32,-1,-1), (34,1,1), (40,1,-1), (45,-1,1)),
          2059: ((1,1,-1), (36,-1,-1), (38,1,1), (43,-1,1), (45,1,-1))}

# route 1 for three or four supersingular primes, N <= 700: sum over divisors M with >= 2 primes of the level-M newforms with every sign +1
DENSE_ROUTE1 = {30: 0, 42: 0, 66: 0, 70: 0, 78: 0, 102: 1, 105: 0, 110: 0, 114: 1, 130: 1, 138: 1, 154: 2, 165: 2, 170: 2, 174: 1, 182: 1, 186: 2,
                190: 1, 195: 1, 210: 1, 230: 2, 231: 1, 238: 1, 246: 3, 255: 2, 266: 2, 273: 4, 282: 3, 285: 2, 286: 2, 290: 3, 310: 3, 322: 4, 330: 2,
                345: 4, 354: 4, 357: 2, 374: 4, 385: 4, 390: 2, 399: 4, 406: 5, 410: 5, 418: 5, 426: 4, 429: 3, 434: 4, 435: 5, 442: 5, 455: 3, 462: 3,
                465: 5, 470: 6, 483: 4, 494: 5, 506: 6, 510: 3, 546: 4, 561: 6, 570: 4, 574: 5, 590: 6, 595: 5, 598: 6, 609: 6, 615: 6, 627: 6, 638: 7,
                646: 7, 651: 8, 658: 7, 663: 5, 665: 8, 682: 8, 690: 6}

# analytic rank of the elliptic curve X0(N)* (the rational newform with signs (+, +)), PARI lfunorderzero
RANKS = {57: 1, 58: 1, 65: 1, 77: 1, 82: 1, 91: 1, 118: 1, 123: 1, 141: 1, 142: 1, 143: 1, 145: 1, 155: 1}

# GAP 4.12, CTblLib 1.3.7: Set(OrdersClassRepresentatives(CharacterTable("M"))) — 73 orders over 194 classes
MONSTER_ELEMENT_ORDERS = (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34,
                          35, 36, 38, 39, 40, 41, 42, 44, 45, 46, 47, 48, 50, 51, 52, 54, 55, 56, 57, 59, 60, 62, 66, 68, 69, 70, 71, 78, 84, 87, 88,
                          92, 93, 94, 95, 104, 105, 110, 119)

ODD_MODE_TORUS_UNITS = 3                    # an odd mode needs torus degree 3 (OM-1); each unit costs 2g - 2 chiral families (TP-1)


# ------------------------------------------------------------------------------------------------------------------------------- arithmetic
def _is_prime(n):
    return n > 1 and all(n % k for k in range(2, int(n ** 0.5) + 1))

def _check_primes(primes):
    ps = tuple(primes)
    if not ps or len(set(ps)) != len(ps) or not all(isinstance(p, int) and _is_prime(p) for p in ps):
        raise ValueError(f"expected distinct primes (a squarefree level), got {primes!r}")
    return tuple(sorted(ps))

def kronecker(a, l):
    """Kronecker symbol (a / l) for a prime l."""
    if l == 2:
        return 0 if a % 2 == 0 else (1 if a % 8 in (1, 7) else -1)
    a %= l
    return 0 if a == 0 else (1 if pow(a, (l - 1) // 2, l) == 1 else -1)

def conductor(D):
    """f with D = D0 f^2 and D0 a fundamental discriminant."""
    f = 1
    for k in range(2, int(abs(D) ** 0.5) + 1):
        while D % (k * k) == 0 and (D // (k * k)) % 4 in (0, 1):
            D //= k * k; f *= k
    return f

def eichler(D, l):
    """Eichler symbol: 1 when l divides the conductor of D, otherwise the Kronecker symbol (D / l)."""
    return 1 if conductor(D) % l == 0 else kronecker(D, l)

@functools.lru_cache(maxsize=None)
def class_number(D):
    """Class number of the imaginary quadratic order of discriminant D: primitive reduced forms (a, b, c), b^2 - 4ac = D."""
    if D >= 0 or D % 4 not in (0, 1):
        raise ValueError(f"not a negative discriminant: {D}")
    n, a = 0, 1
    while 3 * a * a <= -D:
        for b in range(-a + 1, a + 1):
            if (b * b - D) % (4 * a):
                continue
            c = (b * b - D) // (4 * a)
            if c < a or (b < 0 and a == c) or math.gcd(math.gcd(a, abs(b)), c) != 1:
                continue
            n += 1
        a += 1
    return n

def index(primes):
    return math.prod(p + 1 for p in _check_primes(primes))

def elliptic_points(primes):
    ps = _check_primes(primes)
    return math.prod(1 + kronecker(-4, p) for p in ps), math.prod(1 + kronecker(-3, p) for p in ps)

def genus(primes):
    """Genus of X0(N) for squarefree N, from the index, the elliptic points and the 2^k cusps."""
    ps = _check_primes(primes); nu2, nu3 = elliptic_points(ps)
    g = 1 + Fraction(index(ps), 12) - Fraction(nu2, 4) - Fraction(nu3, 3) - Fraction(2 ** len(ps), 2)
    return int(g)

def new_dim(p, q):
    """dim S2^new(Gamma0(pq)) = g(pq) - 2 g(p) - 2 g(q)."""
    return genus((p, q)) - 2 * genus((p,)) - 2 * genus((q,))

def fixed_points(Q_primes, other_primes=()):
    """Fixed points of W_Q on X0(N), N = Q * prod(other_primes): Ogg's class-number count with Eichler symbols, plus the Z[i] term for Q = 2."""
    qs, others = _check_primes(Q_primes), tuple(other_primes)
    if others:
        _check_primes(qs + others)
    Q = math.prod(qs)
    a = class_number(-4 * Q) * math.prod(1 + eichler(-4 * Q, l) for l in others)
    b = class_number(-Q) * math.prod(1 + eichler(-Q, l) for l in others) if Q % 4 == 3 else 0
    c = math.prod(1 + eichler(-4, l) for l in others) if Q == 2 else 0
    return a + b + c

def star_genus(primes):
    """Route 2: genus of X0(N)/W(N) by the holomorphic Lefschetz formula, tr(W_Q | H^0(K)) = 1 - nu(Q)/2."""
    ps = _check_primes(primes)
    tot = Fraction(genus(ps))
    for k in range(1, len(ps) + 1):
        for sub in itertools.combinations(ps, k):
            tot += 1 - Fraction(fixed_points(sub, [l for l in ps if l not in sub]), 2)
    s = tot / 2 ** len(ps)
    return int(s) if s.denominator == 1 else s

def star_genus_route1(N):
    """Route 1, frozen from PARI: level-N newforms with every Atkin–Lehner sign +1 (two primes), or the frozen divisor sum (three or four)."""
    if N in ROUTE1:
        return sum(d for d, wp, wq in ROUTE1[N] if wp == 1 and wq == 1)
    if N in DENSE_ROUTE1:
        return DENSE_ROUTE1[N]
    raise ValueError(f"level {N} is not in the frozen PARI tables")

def root_number(N):
    """Root number of the elliptic curve X0(N)* (two primes, genus 1): -w_N with w_N = w_p w_q = +1, so -1 — odd rank is forced."""
    if star_genus_route1(N) != 1 or N not in ROUTE1:
        raise ValueError(f"X0({N})* is not an elliptic curve in the two-prime table")
    (d, wp, wq), = [o for o in ROUTE1[N] if o[1] == 1 and o[2] == 1]
    return -(wp * wq)


# ------------------------------------------------------------------------------------------------------------------------------- the census
def levels():
    """The 105 pairs p < q of supersingular primes."""
    return list(itertools.combinations(SUPERSINGULAR, 2))

def census():
    """One row per level N = pq: genus, index, torsion-freeness, new dimension, the star genus by both routes, and the level's shape."""
    return [dict(r) for r in _census()]

@functools.lru_cache(maxsize=None)
def _census():
    rows = []
    for p, q in levels():
        N = p * q; i, j = SUPERSINGULAR.index(p), SUPERSINGULAR.index(q)
        g1 = star_genus_route1(N) if N in ROUTE1 else None
        nu2, nu3 = elliptic_points((p, q))
        rows.append({"N": N, "p": p, "q": q, "genus": genus((p, q)), "index": index((p, q)), "cusps": 4, "torsion_free": nu2 == 0 and nu3 == 0,
                     "new_dim": new_dim(p, q), "orbits": tuple(sorted(d for d, _, _ in ROUTE1[N])) if N in ROUTE1 else None,
                     "gstar": star_genus((p, q)), "gstar_route1": g1, "routes": 1 + (g1 is not None),
                     "twin": q - p == 2, "consecutive": j == i + 1})
    return tuple(rows)

def _row(N):
    for r in _census():
        if r["N"] == N:
            return dict(r)
    raise ValueError(f"{N} is not a product of two supersingular primes")

def spheres():
    return sorted(r["N"] for r in census() if r["gstar"] == 0)

def tori():
    return sorted(r["N"] for r in census() if r["gstar"] == 1)

def twin_ladder():
    return [(r["N"], r["gstar"]) for r in sorted(census(), key=lambda r: r["N"]) if r["twin"]]

def consecutive_ladder():
    return [(r["N"], r["gstar"]) for r in sorted(census(), key=lambda r: r["N"]) if r["consecutive"]]

def new_dim_levels(which="smaller"):
    if which not in ("smaller", "larger"):
        raise ValueError("which must be 'smaller' or 'larger'")
    return sorted(r["N"] for r in census() if r["new_dim"] == (r["p"] if which == "smaller" else r["q"]))

def distinctions(N=143):
    """Which arithmetic properties single out the level N among the 105 (each True means: N has it and no other level does)."""
    rows, me = census(), _row(N)
    def only(pred):
        return pred(me) and all(not pred(r) for r in rows if r["N"] != N)
    return {"only_twin_torus": only(lambda r: r["twin"] and r["gstar"] == 1),
            "only_new_dim_smaller_and_consecutive": only(lambda r: r["new_dim"] == r["p"] and r["consecutive"]),
            "only_orbit_pattern": me["orbits"] is not None and only(lambda r: r["orbits"] == me["orbits"]),
            "only_consecutive_torus": only(lambda r: r["consecutive"] and r["gstar"] == 1),
            "only_index": only(lambda r: r["index"] == me["index"]),
            "only_torus": only(lambda r: r["gstar"] == 1)}

def siblings(N=143):
    """Other levels sharing N's index, genus, cusps, torsion-freeness and star genus."""
    me = _row(N); keys = ("index", "genus", "cusps", "torsion_free", "gstar")
    return sorted(r["N"] for r in census() if r["N"] != N and all(r[k] == me[k] for k in keys))

def ranks():
    return dict(RANKS)

def dense_census(bound=DENSE_BOUND):
    """Levels with three or four supersingular primes up to the bound: both routes where route 1 was run."""
    rows = []
    for k in (3, 4):
        for ps in itertools.combinations(SUPERSINGULAR, k):
            N = math.prod(ps)
            if N <= bound:
                rows.append({"N": N, "primes": ps, "cusps": 2 ** k, "gstar": star_genus(ps), "gstar_route1": DENSE_ROUTE1.get(N)})
    return sorted(rows, key=lambda r: r["N"])

def monster_check():
    """Star-quotient spheres against the Monster's element orders (Ogg's observation, carried to products)."""
    M = set(MONSTER_ELEMENT_ORDERS)
    sph = spheres() + [r["N"] for r in dense_census() if r["gstar"] == 0]
    semiprime_orders = sorted(r["N"] for r in census() if r["N"] in M)
    return {"spheres": sph, "spheres_are_orders": all(N in M for N in sph), "orders_not_spheres": [N for N in semiprime_orders if N not in sph],
            "143_is_an_order": 143 in M, "n_orders": len(M)}

def flavour_transfer():
    """Cost of an odd B-stack mode, 3 units of torus flux at 2g - 2 chiral families each, at every torus level."""
    return {N: ODD_MODE_TORUS_UNITS * (2 * _row(N)["genus"] - 2) for N in tori()}


# ------------------------------------------------------------------------------------------------------------------------------- views
def report():
    rows = census(); d = dense_census()
    both = [r for r in rows if r["routes"] == 2]
    return {"levels": len(rows), "spheres": spheres(), "tori": tori(), "higher": sum(1 for r in rows if r["gstar"] >= 2),
            "two_routes": len(both), "routes_agree": all(r["gstar"] == r["gstar_route1"] for r in both),
            "dense_agree": all(r["gstar"] == r["gstar_route1"] for r in d if r["gstar_route1"] is not None), "dense_levels": len(d),
            "dense_tori": [r["N"] for r in d if r["gstar"] == 1], "dense_spheres": [r["N"] for r in d if r["gstar"] == 0],
            "twin_ladder": twin_ladder(), "consecutive_ladder": consecutive_ladder(), "new_dim_smaller": new_dim_levels("smaller"),
            "distinctions_143": distinctions(143), "siblings_143": siblings(143), "ranks": ranks(), "monster": monster_check(),
            "odd_mode_cost": flavour_transfer()}

def text(level=None):
    if level is not None:
        r = _row(level); dist = distinctions(level)
        shape = "torsion-free" if r["torsion_free"] else "elliptic points"
        g1 = "" if r["gstar_route1"] is None else f", {r['gstar_route1']} (route 1)"
        lines = [f"X0({level}) = X0({r['p']} x {r['q']}): genus {r['genus']}, index {r['index']}, {shape}, "
                 f"new dimension {r['new_dim']}, Galois orbits {r['orbits']}",
                 f"  star quotient genus {r['gstar']} (route 2){g1}",
                 "  unique among the 105: " + (", ".join(k for k, v in dist.items() if v) or "nothing"),
                 "  siblings (same index, genus, cusps, torsion, star genus): " + (", ".join(map(str, siblings(level))) or "none")]
        return "\n".join(lines)
    R = report()
    return "\n".join([
        f"{R['levels']} levels N = pq (p < q supersingular); two routes at {R['two_routes']}, agreeing: {R['routes_agree']}",
        f"star quotient a sphere at {len(R['spheres'])}: {', '.join(map(str, R['spheres']))}",
        f"an elliptic curve at {len(R['tori'])} (all rank 1): {', '.join(map(str, R['tori']))}",
        f"genus >= 2 at {R['higher']}",
        "twin ladder: " + ", ".join(f"{N} -> {g}" for N, g in R["twin_ladder"]),
        "143 unique in: " + ", ".join(k for k, v in R["distinctions_143"].items() if v),
        "siblings of 143: " + ", ".join(map(str, R["siblings_143"])),
        f"three or four primes (N <= {DENSE_BOUND}): tori {', '.join(map(str, R['dense_tori']))}; spheres {', '.join(map(str, R['dense_spheres']))}",
        f"every sphere level is a Monster element order: {R['monster']['spheres_are_orders']}; 143 is an order: {R['monster']['143_is_an_order']}"])
