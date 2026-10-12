"""Provenance script for SQ-01: route 1 of the star-quotient census, from PARI/GP's modular-forms package.

Requires cypari2 (pip install cypari2).  For each level N = pq with p < q supersingular, splits S2^new(Gamma0(N)) into Galois orbits (mfsplit) and
reads each orbit's Atkin–Lehner signs (mfatkineigenvalues); for three or four supersingular primes (N <= 700) it sums, over divisors M with at least
two primes, the level-M newforms whose signs are all +1; at the genus-1 levels it computes the analytic rank of the (+, +) newform (lfunorderzero).
These are frozen as ROUTE1, DENSE_ROUTE1 and RANKS in mtft.matter.star_quotients; route 2 there is pure Python and independent of this script.

    python scripts/star_quotient_scan.py               # two-prime levels up to --max-level (default 2100), dense levels, ranks; JSON on stdout
    python scripts/star_quotient_scan.py --check       # also compare with the frozen tables; exit status 1 on any difference

Levels 2201, 2419, 2773, 2911, 3337 and 4189 were not reached within the time allowed (route 2 covers them)."""
import argparse
import itertools
import json
import math
import sys

import cypari2

SS = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 41, 47, 59, 71)


def orbits(pari, N, primes):
    pari(f"mf = mfinit([{N},2],0)")
    if int(pari("mfdim(mf)")) == 0:
        return []
    pari("S = mfsplit(mf)")
    degs = [int(x) for x in pari("apply(poldegree, S[2])")]
    signs = [[int(v[0]) for v in pari(f"mfatkineigenvalues(mf,{p})")] for p in primes]
    return [(d, *(s[i] for s in signs)) for i, d in enumerate(degs)]


def all_plus(pari, M, primes, cache):
    if M not in cache:
        cache[M] = sum(o[0] for o in orbits(pari, M, primes) if all(s == 1 for s in o[1:]))
    return cache[M]


def rank(pari, N, p, q):
    o = orbits(pari, N, (p, q))
    k = [i for i, t in enumerate(o) if t == (1, 1, 1)]
    assert len(k) == 1, (N, o)
    pari(f"F = mfeigenbasis(mf)[{k[0] + 1}]; L = lfunmf(mf, F)")
    return int(pari("lfunorderzero(L)"))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--max-level", type=int, default=2100)
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args(argv)
    pari = cypari2.Pari(); pari.allocatemem(2 * 10**9)
    route1 = {p * q: sorted(orbits(pari, p * q, (p, q)), key=lambda o: (o[0], -o[1], -o[2]))
              for p, q in itertools.combinations(SS, 2) if p * q <= a.max_level}
    cache, dense = {}, {}
    for k in (3, 4):
        for ps in itertools.combinations(SS, k):
            N = math.prod(ps)
            if N <= 700:
                dense[N] = sum(all_plus(pari, math.prod(sub), sub, cache) for r in range(2, k + 1) for sub in itertools.combinations(ps, r))
    tori = sorted(N for N, o in route1.items() if sum(t[0] for t in o if t[1] == 1 and t[2] == 1) == 1)
    ranks = {N: rank(pari, N, *[p for p in SS if N % p == 0]) for N in tori}
    out = {"pari": ".".join(map(str, pari.version())), "route1": {N: [list(t) for t in o] for N, o in route1.items()}, "dense": dense, "ranks": ranks}
    print(json.dumps(out))
    if a.check:
        from mtft.matter import star_quotients as S
        bad = [N for N, o in route1.items() if N in S.ROUTE1 and tuple(map(tuple, o)) != S.ROUTE1[N]]
        bad += [N for N, g in dense.items() if S.DENSE_ROUTE1.get(N) != g] + [N for N, r in ranks.items() if S.RANKS.get(N) != r]
        print("differences from the frozen tables:", bad or "none", file=sys.stderr)
        return 1 if bad else 0
    return 0


if __name__ == "__main__":
    sys.exit(main())
