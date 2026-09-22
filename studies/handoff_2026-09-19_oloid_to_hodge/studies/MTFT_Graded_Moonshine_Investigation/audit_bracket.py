#!/usr/bin/env python3
"""Canonical first Rankin--Cohen / Wahl bracket on frozen X0(143) data.

Run: python audit_bracket.py --output bracket_results.json
Requires Python 3.10+ and NumPy. Does not import mtft or modify its source.

For f in R_n=H0(K^n), g in R_m, use
    {f,g} = n*f*theta(g) - m*g*theta(f), theta=q*d/dq.
This is half the conventional first Rankin--Cohen bracket for weights
2n,2m, and lands in R_(n+m+1). At a fixed h in R1, {h,-} raises
canonical degree by TWO. Canonical degree is distinct from Fourier order.

The cancellation of coordinate-change terms and the cusp vanishing bound
give the geometric global closure theorem; finite-field membership checks
below are consistency checks of the supplied Fourier data, not a substitute
for that theorem. Jacobi and Leibniz tests use exact integer coefficients.
No real physical Hamiltonian, unitary evolution, or Feigenbaum RG is inferred.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
PRECISION = 141
PRIMES = (1000003, 1000033)
SECTOR_NAMES = ("(+,+)", "(+,-)", "(-,+)", "(-,-)")
SECTORS = [0] + [1]*6 + [2]*5 + [3]
R3_DIMS = (12, 18, 17, 13)
FIXTURE = "X0_143_AL_adapted_qexpansions.txt"
FIXTURE_HASH = "cdb9f3caf6fb6eb529931f41a2b47f0657489328737eb5abf0ab82150b56eb52"


def read_forms(data_dir):
    raw = (data_dir / FIXTURE).read_bytes()
    assert hashlib.sha256(raw).hexdigest() == FIXTURE_HASH
    rows = list(csv.reader(line for line in raw.decode().splitlines()
                           if line and not line.startswith("#")))
    coefficients = [[int(x) for x in row[1:]] for row in rows[1:]]
    assert len(coefficients) == PRECISION
    assert [int(row[0]) for row in rows[1:]] == list(range(PRECISION))
    return [[row[j] for row in coefficients] for j in range(13)]


def product(f, g):
    return [sum(f[i] * g[k-i] for i in range(k+1)) for k in range(PRECISION)]


def bracket(f, g, n=1, m=1):
    return [sum((n*(k-i) - m*i) * f[i] * g[k-i] for i in range(k+1))
            for k in range(PRECISION)]


def add(*vectors):
    return [sum(xs) for xs in zip(*vectors)]


def minus(f):
    return [-x for x in f]


def reduce_mod(vector, basis, p):
    v = np.array([int(x) % p for x in vector], dtype=np.int64)
    for pivot in sorted(basis):
        if v[pivot]:
            v = (v - int(v[pivot])*basis[pivot]) % p
    return v


def basis_mod(vectors, p):
    # No modular convolution: each exact Python product is reduced first.
    # Elimination multiplies residues at most (p-1)^2, far below int64 max.
    assert (p-1)**2 + (p-1) < np.iinfo(np.int64).max
    basis = {}
    independent_indices = []
    for j, vector in enumerate(vectors):
        v = reduce_mod(vector, basis, p)
        nz = np.flatnonzero(v)
        if len(nz):
            pivot = int(nz[0])
            basis[pivot] = (v * pow(int(v[pivot]), -1, p)) % p
            independent_indices.append(j)
    return basis, independent_indices


def nonzero_terms(v, count=8):
    return [{"q_power": n, "coefficient": a}
            for n, a in enumerate(v) if a][:count]


def audit(data_dir):
    E = read_forms(data_dir)
    # Independent construction: all 455 degree-three generator monomials,
    # not the recursive basis reduction used in verify_graded.py.
    pair_products = {(i,j): product(E[i], E[j])
                     for i,j in itertools.combinations_with_replacement(range(13), 2)}
    cubics = [[] for _ in SECTOR_NAMES]
    cubic_labels = [[] for _ in SECTOR_NAMES]
    for i,j,k in itertools.combinations_with_replacement(range(13), 3):
        sector = SECTORS[i] ^ SECTORS[j] ^ SECTORS[k]
        cubics[sector].append(product(pair_products[i,j], E[k]))
        cubic_labels[sector].append([i+1,j+1,k+1])
    assert sum(map(len,cubics)) == 455

    brackets = [[] for _ in SECTOR_NAMES]
    bracket_labels = [[] for _ in SECTOR_NAMES]
    for i,j in itertools.combinations(range(13), 2):
        sector = SECTORS[i] ^ SECTORS[j]
        brackets[sector].append(bracket(E[i], E[j]))
        bracket_labels[sector].append([i+1,j+1])
    assert [len(v) for v in brackets] == [25,11,11,31]

    fields = {}
    for p in PRIMES:
        records = []
        for sector in range(4):
            cbasis, cpivots = basis_mod(cubics[sector], p)
            bbasis, bpivots = basis_mod(brackets[sector], p)
            assert len(cbasis) == R3_DIMS[sector]
            assert max(cbasis) <= 84
            failures = [j+1 for j,v in enumerate(brackets[sector])
                        if np.any(reduce_mod(v, cbasis, p))]
            assert not failures
            records.append({
                "sector": SECTOR_NAMES[sector],
                "cubic_monomials_tested": len(cubics[sector]),
                "R3_basis_rank": len(cbasis),
                "R3_basis_pivot_q_orders": sorted(int(x) for x in cbasis),
                "R3_independent_monomials_one_based": [cubic_labels[sector][j] for j in cpivots],
                "wedge_domain_dimension": len(brackets[sector]),
                "bracket_image_rank": len(bbasis),
                "independent_bracket_pairs_one_based": [bracket_labels[sector][j] for j in bpivots],
                "brackets_outside_cubic_span_mod_p": failures,
            })
        fields[str(p)] = records

    exact_ranks = []
    exact_rank_method = []
    for sector in range(4):
        a, b = (fields[str(p)][sector]["bracket_image_rank"] for p in PRIMES)
        assert a == b
        upper = min(len(brackets[sector]), R3_DIMS[sector])
        assert a == upper
        rank_q = a
        exact_rank_method.append("finite-field lower bound saturates minimum of domain dimension and geometric target dimension")
        exact_ranks.append(rank_q)
    assert exact_ranks == [12,11,11,13]

    assert all(x % 72 == 0 for x in E[0])
    assert all((x+y) % 2 == 0 for x,y in zip(E[7], E[12]))
    h = [x//72 for x in E[0]]
    g = [(x+y)//2 for x,y in zip(E[7], E[12])]
    gh = bracket(g, h)
    components = [bracket(E[7], E[0]), bracket(E[12], E[0])]
    assert all(144*x == a+b for x,a,b in zip(gh,*components))
    assert gh[3] == 2 and all(x == 0 for x in gh[:3])

    # Exact formal-series checks, with the shifted canonical degrees included.
    f1,f2,f8 = E[0],E[1],E[7]
    jacobi = add(bracket(f1,bracket(f2,f8),1,3),
                 bracket(f2,bracket(f8,f1),1,3),
                 bracket(f8,bracket(f1,f2),1,3))
    leibniz = add(bracket(f1,product(f2,f8),1,2),
                  minus(product(bracket(f1,f2),f8)),
                  minus(product(f2,bracket(f1,f8))))
    skew = [add(bracket(f,g),bracket(g,f)) for f,g in [(f1,f2),(f2,f8),(f8,f1)]]
    assert all(x == 0 for x in jacobi + leibniz)
    assert all(x == 0 for v in skew for x in v)

    return {
        "source_release": "MTFT 0.32.0",
        "fixture_sha256": {FIXTURE: FIXTURE_HASH},
        "convention": "{f,g}=n*f*theta(g)-m*g*theta(f) for canonical degrees n,m; theta=q*d/dq",
        "output_degree": "n+m+1; fixed h of degree 1 gives degree +2 derivation",
        "normalization": "half the usual first Rankin-Cohen bracket for modular weights 2n,2m",
        "sector_order": list(SECTOR_NAMES),
        "q_coefficients_checked_inclusive": [0,140],
        "weight6_Sturm_bound": 84,
        "finite_field_cubic_membership_checks": fields,
        "Wahl_map": {
            "domain": "wedge^2 R1, dimension 78",
            "codomain": "R3=H0(K^3), dimension 60",
            "domain_dimensions_by_sector": [25,11,11,31],
            "target_dimensions_by_sector": list(R3_DIMS),
            "rational_image_ranks_by_sector": exact_ranks,
            "rank_certification_methods": exact_rank_method,
            "total_rank": sum(exact_ranks),
            "kernel_dimension": 78-sum(exact_ranks),
            "cokernel_dimensions_by_sector": [b-a for a,b in zip(exact_ranks,R3_DIMS)],
            "total_cokernel_dimension": sum(R3_DIMS)-sum(exact_ranks),
            "interpretation": "The bracket saturates the maximum permitted in every AL sector, but the target has 7 more (+,-) and 6 more (-,+) dimensions than the source. These 13 missing directions are an explicit representation-theoretic deficit of brackets of degree-one generators.",
        },
        "normalized_example": {
            "g": "(e8+e13)/2 = eta(tau)^2 eta(11*tau)^2",
            "h": "e1/72 = normalized 143a1 form",
            "bracket_first_nonzero_terms": nonzero_terms(gh),
            "bracket_coefficient_q3": 2,
            "decomposition": "{g,h}=({e8,e1}+{e13,e1})/144",
            "AL_sectors": ["(-,+)","(-,-)"],
            "all_78_generator_brackets_checked_in_R3_mod_both_primes": True,
        },
        "Poisson_identity_controls": {
            "forms": "e1,e2,e8 (unnormalized exact integral frozen forms)",
            "Jacobi_max_abs_residual": max(map(abs,jacobi)),
            "Jacobi_canonical_degree": 5,
            "Jacobi_modular_weight": 10,
            "Jacobi_Sturm_bound": 140,
            "Leibniz_max_abs_residual": max(map(abs,leibniz)),
            "Leibniz_canonical_degree": 4,
            "Leibniz_modular_weight": 8,
            "Leibniz_Sturm_bound": 112,
            "skew_symmetry_max_abs_residual": max(abs(x) for v in skew for x in v),
            "scope": "Exact data checks for these forms through q^140. The full identities are algebraic consequences of the bracket formula; global closure relies on geometric coordinate cancellation and cusp behavior.",
        },
        "physical_status": "This is a canonical graded Poisson algebra and its algebraic Hamiltonian derivations. No Hilbert adjoint, physical energy scale, self-adjoint Hamiltonian, unitary flow, or renormalization fixed point has been derived.",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir",type=Path,default=HERE/"data")
    parser.add_argument("--output",type=Path,default=HERE/"bracket_results.json")
    args = parser.parse_args()
    result = audit(args.data_dir)
    args.output.write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps({"output":str(args.output),
                      "Wahl_rank":result["Wahl_map"]["total_rank"],
                      "Wahl_ranks_by_sector":result["Wahl_map"]["rational_image_ranks_by_sector"],
                      "Wahl_cokernel_by_sector":result["Wahl_map"]["cokernel_dimensions_by_sector"],
                      "normalized_bracket_q3":2,
                      "Jacobi_and_Leibniz_residuals":0},indent=2))


if __name__ == "__main__":
    main()
