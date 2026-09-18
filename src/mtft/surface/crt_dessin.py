"""INT-04a (HOPF-02 §§1–3): P^1(Z/143) = P^1(F_11) x P^1(F_13) by CRT and the canonical dessin of X0(143).

CRT idempotents 78, 66 (78 = 1 mod 11, 0 mod 13; 66 = 0 mod 11, 1 mod 13).  S = [[0,-1],[1,0]] and T = [[1,1],[0,1]] act on
the 168 = 12 * 14 projective points; T has cycles of lengths 1, 11, 13, 143 (the cusp widths), S has 84 orbits (edges),
ST has 56 orbits (triangles): 4 - 84 + 56 = 2 - 2g gives g = 13.  Level 143 is a modular level, not a sphere dimension
(H-10)."""
from math import gcd
import itertools

E11, E13 = 78, 66

def crt(a, b, N=143): return (E11 * a + E13 * b) % N

def canon_prime(c, d, p):
    c %= p; d %= p
    if not (c or d): raise ValueError("zero vector")
    return (0, 1) if c == 0 else (1, (d * pow(c, -1, p)) % p)

def crt_projective_line():
    """Ordered pairs (P^1(F_11) x P^1(F_13)) with their CRT representatives (c, d) mod 143."""
    P11 = [(0, 1)] + [(1, t) for t in range(11)]; P13 = [(0, 1)] + [(1, t) for t in range(13)]
    pairs = list(itertools.product(P11, P13)); reps = [(crt(u[0], v[0]), crt(u[1], v[1])) for u, v in pairs]
    return {"pairs": pairs, "reps": reps, "count": len(pairs), "all_primitive": all(gcd(gcd(c, d), 143) == 1 for c, d in reps)}

def act(pair, M, p):
    c, d = pair; a, b, cc, dd = M; return canon_prime(c * a + d * cc, c * b + d * dd, p)

def permutation(M, pairs=None):
    pairs = pairs or crt_projective_line()["pairs"]; index = {pair: i for i, pair in enumerate(pairs)}
    return [index[(act(u, M, 11), act(v, M, 13))] for u, v in pairs]

def cycles(perm):
    seen = set(); out = []
    for j in range(len(perm)):
        if j in seen: continue
        cyc = []; k = j
        while k not in seen: seen.add(k); cyc.append(k); k = perm[k]
        out.append(cyc)
    return out

def dessin():
    S, T = (0, -1, 1, 0), (1, 1, 0, 1); pairs = crt_projective_line()["pairs"]; PS, PT = permutation(S, pairs), permutation(T, pairs)
    ST = [PT[PS[j]] for j in range(len(pairs))]
    widths = sorted(len(c) for c in cycles(PT)); edges = len(cycles(PS)); triangles = len(cycles(ST))
    return {"darts": len(pairs), "cusp_widths": widths, "edge_orbits": edges, "triangle_orbits": triangles,
            "genus": (2 - (len(widths) - edges + triangles)) // 2, "S_involution": all(PS[PS[j]] == j for j in range(len(pairs))),
            "ST_order_three": all(ST[ST[ST[j]]] == j for j in range(len(pairs)))}

def genus_squarefree(N):
    import sympy as sp
    primes = list(sp.factorint(N)); mu = sp.Integer(N) * sp.prod(1 + sp.Rational(1, p) for p in primes)
    e2 = sp.prod(1 + sp.legendre_symbol(-1, p) for p in primes); e3 = sp.prod(1 + sp.legendre_symbol(-3, p) for p in primes); cusps = 2 ** len(primes)
    return {"level": N, "index": mu, "e2": e2, "e3": e3, "cusps": cusps, "genus": 1 + mu / 12 - e2 / 4 - e3 / 3 - sp.Rational(cusps, 2)}
