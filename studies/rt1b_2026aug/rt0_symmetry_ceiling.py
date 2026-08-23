"""RT-0: how large a continuous symmetry can the canonical ideal carry?

Sol proposes pre-registering an "SM representation gate suite" to see
whether su(3) + su(2) + u(1) emerges from v0.17.0.  Before spending a
session on that, compute the ceiling.

The object: the Lie algebra of the stabiliser of I_2 inside gl_13,

    stab(I_2) = { M in gl_13 : D_M(Q) in I_2 for every Q in I_2 }

where D_M is the derivation induced by x -> x + eps M x.  Any algebraic
group acting linearly on P^12 and preserving the canonical curve has its
Lie algebra inside this.  Scalars always qualify (D_(lam I) Q = 2 lam Q),
so dim >= 1 unconditionally; dim = 1 means nothing continuous acts.

Positive control: the rational normal curve of degree 12 in P^12, whose
ideal is C(12,2) = 66 quadrics and whose automorphisms are PGL_2.  There
the answer must be 3 + 1 = 4.  If the code returns 4 there and 1 here,
the 1 is a real no-go and not a bug.
"""

from itertools import combinations_with_replacement

from mtft.canonical import MONOMIALS, ideal_basis

PRIMES = (2147483647, 1000003)
G = 13
MIDX = {m: k for k, m in enumerate(MONOMIALS)}


def rref_mod_p(rows, ncols, p):
    rows = [[v % p for v in r] for r in rows]
    piv, r = [], 0
    for c in range(ncols):
        k = next((i for i in range(r, len(rows)) if rows[i][c]), None)
        if k is None:
            continue
        rows[r], rows[k] = rows[k], rows[r]
        inv = pow(rows[r][c], p - 2, p)
        rows[r] = [v * inv % p for v in rows[r]]
        for i in range(len(rows)):
            if i != r and rows[i][c]:
                f = rows[i][c]
                rows[i] = [(a - f * b) % p for a, b in zip(rows[i], rows[r])]
        piv.append(c)
        r += 1
        if r == len(rows):
            break
    return rows[:r], piv


def rank_mod_p(rows, ncols, p):
    return len(rref_mod_p(rows, ncols, p)[0])


def derivation(quad, a, b):
    """D_{E_ab} applied to a quadric given as {monomial index: coeff}."""
    out = {}
    for m, c in quad.items():
        i, j = MONOMIALS[m]
        if i == a:
            k = MIDX[tuple(sorted((b, j)))]
            out[k] = out.get(k, 0) + c
        if j == a:
            k = MIDX[tuple(sorted((i, b)))]
            out[k] = out.get(k, 0) + c
    return out


def stabiliser_dim(quadrics, p):
    """dim of { M in gl_13 : D_M(I) subset I }, computed over F_p."""
    rows = [[q.get(m, 0) for m in range(91)] for q in quadrics]
    ech, piv = rref_mod_p(rows, 91, p)
    pivset = set(piv)
    free = [c for c in range(91) if c not in pivset]

    def reduce_mod_I(vec):
        v = [x % p for x in vec]
        for r, c in zip(ech, piv):
            if v[c]:
                f = v[c]
                v = [(a - f * b) % p for a, b in zip(v, r)]
        return [v[c] for c in free]

    big = []
    for a in range(G):
        for b in range(G):
            row = []
            for q in quadrics:
                d = derivation(q, a, b)
                vec = [0] * 91
                for m, c in d.items():
                    vec[m] = c
                row.extend(reduce_mod_I(vec))
            big.append(row)
    r = rank_mod_p(big, len(big[0]), p)
    return G * G - r


def rational_normal_curve_quadrics():
    """I_2 of the degree-12 rational normal curve.

    Monomials x_i x_j with i + j = s all restrict to t^s on the curve, so
    consecutive differences within each block vanish.  Taking differences
    of consecutive monomials block by block gives 91 - 25 = 66 quadrics
    with entries +-1 that are manifestly independent over Z — no modular
    reduction anywhere, so the same integer basis is valid at every prime.
    """
    blocks = {}
    for i, j in MONOMIALS:
        blocks.setdefault(i + j, []).append((i, j))
    quads = []
    for s in sorted(blocks):
        mons = sorted(blocks[s])
        for a, b in zip(mons, mons[1:]):
            quads.append({MIDX[a]: 1, MIDX[b]: -1})
    return quads


if __name__ == "__main__":
    I2 = ideal_basis()
    quads = [{m: I2[m][c] for m in range(91) if I2[m][c]} for c in range(55)]

    print("CONTROL — rational normal curve of degree 12 in P^12")
    rnc = rational_normal_curve_quadrics()
    print(f"  dim I_2 = {len(rnc)}   [expected C(12,2) = 66]")
    for p in PRIMES:
        print(f"  dim stab(I_2) mod {p} = {stabiliser_dim(rnc, p)}"
              f"   [expected 4 = dim sl_2 + scalars]")

    print()
    print("X0(143) canonical ideal")
    print(f"  dim I_2 = {len(quads)}")
    for p in PRIMES:
        print(f"  dim stab(I_2) mod {p} = {stabiliser_dim(quads, p)}"
              f"   [1 = scalars only means nothing continuous acts]")

    print()
    print("For reference: dim(su(3) + su(2) + u(1)) = 8 + 3 + 1 = 12")
