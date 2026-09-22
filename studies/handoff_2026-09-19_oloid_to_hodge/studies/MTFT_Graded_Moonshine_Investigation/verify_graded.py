#!/usr/bin/env python3
"""Exact graded-ring investigation from frozen MTFT 0.32.0 coefficients.

Requires Python 3.10+ and NumPy. No mtft installation or network is needed.
The geometric interpretation assumes the shipped, provenance-checked modular
forms are the claimed basis. Finite-field ranks are lower bounds over Q;
Riemann--Roch and holomorphic Lefschetz provide the matching upper bounds.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from fractions import Fraction as F
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
SECTORS = ("(+,+)", "(+,-)", "(-,+)", "(-,-)")
# Sector indices use XOR for multiplication of characters.
COORD_SECTORS = [0] + [1] * 6 + [2] * 5 + [3]
PRIMES = (1000003, 1000033)
HASHES = {
    "X0_143_AL_adapted_basis.txt": "b24d315e1f4d3a67651b057266352630ac9dce517609abf7a7b1c38c35e8721f",
    "X0_143_AL_adapted_qexpansions.txt": "cdb9f3caf6fb6eb529931f41a2b47f0657489328737eb5abf0ab82150b56eb52",
    "X0_143_S2_qexpansions.txt": "e7c2aaed96af74c4f21b1aa6d811a427b1a2dc7597abab9d5d5a6308888bb69f",
}


def read_rows(name):
    lines = [s for s in (ROOT / "data" / name).read_text().splitlines()
             if s.strip() and not s.startswith("#")]
    return [[int(x) for x in s.split(",")[1:]] for s in lines[1:]]


def canonical_dim(n):
    return 1 if n == 0 else 13 if n == 1 else 12 * (2 * n - 1)


def canonical_sectors(n):
    if n == 0:
        return [1, 0, 0, 0]
    if n == 1:
        return [1, 6, 5, 1]
    eps = (-1) ** n
    return [6*n-3+3*eps, 6*n-3-3*eps, 6*n-3-2*eps, 6*n-3+2*eps]


def log_sectors(n):
    if n == 0:
        return [1, 0, 0, 0]
    eps = (-1) ** n
    return [7*n-3+3*eps, 7*n-3-3*eps, 7*n-3-2*eps, 7*n-3+2*eps]


def trace_from_sectors(d):
    return {"identity": sum(d), "W11": d[0]+d[1]-d[2]-d[3],
            "W13": d[0]-d[1]+d[2]-d[3],
            "W143": d[0]-d[1]-d[2]+d[3]}


def is_prime(n):
    return n >= 2 and all(n % k for k in range(2, math.isqrt(n) + 1))


def insert_basis(vector, basis, p):
    """Incremental row echelon basis; retain all 141 coefficients."""
    v = vector.copy() % p
    for pivot in sorted(basis):
        if v[pivot]:
            v = (v - int(v[pivot]) * basis[pivot]) % p
    nz = np.flatnonzero(v)
    if len(nz) == 0:
        return False
    pivot = int(nz[0])
    v = (v * pow(int(v[pivot]), -1, p)) % p
    basis[pivot] = v
    return True


def ring_ranks(E, p, max_degree=5):
    assert is_prime(p)
    size = len(E)
    # Products and elimination use int64. These explicit bounds exclude overflow.
    assert size * (p-1)**2 < np.iinfo(np.int64).max
    generators = [np.array([row[i] % p for row in E], dtype=np.int64)
                  for i in range(13)]
    one = np.zeros(size, dtype=np.int64)
    one[0] = 1
    previous = [{0: one}, {}, {}, {}]
    out = []
    for n in range(1, max_degree + 1):
        current = [{}, {}, {}, {}]
        attempted = 0
        for sector, basis in enumerate(previous):
            for v in basis.values():
                for i, f in enumerate(generators):
                    product = np.convolve(v, f)[:size] % p
                    insert_basis(product, current[sector ^ COORD_SECTORS[i]], p)
                    attempted += 1
        ranks = [len(b) for b in current]
        assert ranks == canonical_sectors(n), (p, n, ranks)
        # Pivots within the weight-2n Sturm bound certify that the supplied
        # precision is sufficient, rather than just counting long vectors.
        pivots = [[int(k) for k in sorted(b)] for b in current]
        assert max(k for ks in pivots for k in ks) <= 28*n <= 140
        out.append({"degree": n, "weight": 2*n, "sturm_bound": 28*n,
                    "rank_by_sector": ranks, "total_rank": sum(ranks),
                    "products_processed": attempted, "pivots_by_sector": pivots})
        previous = current
    return out


def convolve(a, b, length):
    return [sum(a[k]*b[n-k] for k in range(n+1)
                if k < len(a) and n-k < len(b)) for n in range(length)]


def sigma(n, power):
    return sum(d**power for d in range(1, n+1) if n % d == 0)


def normalized_moonshine_dimensions(length=11):
    """Compute q(j-744) using E4^3 and prod(1-q^n)^(-24), exact integers."""
    e4 = [1] + [240*sigma(n, 3) for n in range(1, length)]
    oscillators = [1]
    for n in range(1, length):
        numerator = 24 * sum(sigma(k, 1)*oscillators[n-k] for k in range(1, n+1))
        assert numerator % n == 0
        oscillators.append(numerator // n)
    j_times_q = convolve(convolve(convolve(e4, e4, length), e4, length),
                         oscillators, length)
    j_times_q[1] -= 744
    assert j_times_q[:6] == [1, 0, 196884, 21493760, 864299970, 20245856256]
    return j_times_q


def hilbert_coeffs(numerator, count):
    return [sum(a*(n-k+1) for k, a in enumerate(numerator) if k <= n)
            for n in range(count)]


def partition_at_half():
    t = F(1, 2)
    p = 1 + 11*t + 11*t*t + t**3
    dp = 11 + 22*t + 3*t*t
    ddp = 22 + 6*t
    z = p / (1-t)**2
    mean = t*dp/p + 2*t/(1-t)
    variance = t*dp/p + t*t*(ddp*p-dp*dp)/p**2 + 2*t/(1-t)**2
    assert z == F(75, 2) and mean == F(241, 75) and variance == F(24794, 5625)
    return {"t": str(t), "Z": str(z), "mean_degree": str(mean),
            "variance_degree": str(variance),
            "status": "EXACT for the defined grading ensemble; no physical energy scale fixed"}


def investigate():
    for name, wanted in HASHES.items():
        assert hashlib.sha256((ROOT / "data" / name).read_bytes()).hexdigest() == wanted
    E = read_rows("X0_143_AL_adapted_qexpansions.txt")
    S = read_rows("X0_143_S2_qexpansions.txt")
    B = read_rows("X0_143_AL_adapted_basis.txt")
    assert len(E) == 141 and all(len(row) == 13 for row in E)
    assert E == [[sum(S[n][i]*B[i][j] for i in range(13)) for j in range(13)]
                 for n in range(141)]
    mu, cusps, elliptic2, elliptic3 = 168, 4, 0, 0
    assert F(1) + F(mu, 12) - F(cusps, 2) - F(elliptic2, 4) - F(elliptic3, 3) == 13
    assert sum([143, 13, 11, 1]) == mu
    dims = [canonical_dim(n) for n in range(51)]
    assert hilbert_coeffs([1, 11, 11, 1], 51) == dims
    logdims = [1] + [28*n-12 for n in range(1, 51)]
    assert hilbert_coeffs([1, 14, 13], 51) == logdims
    for n in range(51):
        assert sum(canonical_sectors(n)) == dims[n]
        assert sum(log_sectors(n)) == logdims[n]
        traces = trace_from_sectors(canonical_sectors(n))
        if n >= 2:
            assert [traces[k] for k in ("W11", "W13", "W143")] == [0, 2*(-1)**n, 10*(-1)**n]
            assert [b-a for a,b in zip(canonical_sectors(n),log_sectors(n))] == [n]*4
    ranks = {str(p): ring_ranks(E, p) for p in PRIMES}
    monster_order = 808017424794512875886459904961710757005754368000000000
    assert monster_order > 84*12
    assert max(13, 16, 26) < 196883
    moonshine = normalized_moonshine_dimensions()
    rows = [{"degree": n, "canonical": canonical_dim(n),
             "canonical_AL_sectors": canonical_sectors(n),
             "all_even_modular_forms": 1 if n == 0 else 28*n-12,
             "logcanonical_AL_sectors_theoretical": log_sectors(n),
             "moonshine_conformal_weight_dim": moonshine[n]}
            for n in range(len(moonshine))]
    return {
        "investigation": "MTFT graded geometry and Monster compatibility",
        "upstream_version": "0.32.0", "input_sha256": HASHES,
        "sector_order": SECTORS, "genus": 13, "cusps": 4, "index": 168,
        "canonical_hilbert_series": "(1+11t+11t^2+t^3)/(1-t)^2",
        "logcanonical_hilbert_series": "(1+14t+13t^2)/(1-t)^2",
        "twined_canonical_series": {"W11": "1+t", "W13": "(1+t^2)/(1+t)",
                                    "W143": "(1-8t+t^2)/(1+t)"},
        "finite_field_certificates": ranks, "dimension_comparison": rows,
        "partition_at_t_half": partition_at_half(),
        "obstructions": {"Hurwitz_bound": 1008, "Monster_order": str(monster_order),
                         "least_nontrivial_complex_Monster_rep_dimension": 196883,
                         "scope": "grading-preserving algebra actions on these rings; not all conceivable extensions"},
        "verification_status": "PASS",
        "limits": ["Frozen q-expansion modularity and AL labels are upstream mathematical inputs.",
                   "Matching Q-rank uses geometric upper bounds, not agreement of primes alone.",
                   "Modular/logcanonical dimensions use theorems; no full M_2n basis was constructed.",
                   "No VOA, physical Hamiltonian, RG fixed point, Feigenbaum value, or coupling constant was derived."]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "graded_results.json")
    args = parser.parse_args()
    result = investigate()
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print("PASS: two-prime product ranks through degree five, exact Hilbert series,")
    print("AL characters, input basis transformation, moonshine coefficients, and grading ensemble.")
    print("Sector order:", ", ".join(SECTORS))
    for row in result["dimension_comparison"][:6]:
        print(f"n={row['degree']}: canonical={row['canonical']}, sectors={row['canonical_AL_sectors']}")
    print("Results:", args.output)
