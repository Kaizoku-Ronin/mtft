#!/usr/bin/env python3
"""Replay MTFT's numerical D4 fingerprint and an exact triality control.

Run with the package's dependencies available, for example:
  OPENBLAS_NUM_THREADS=1 PYTHONPATH=/path/to/deps python d4_check.py \
      --source /path/to/mtft-0.26.2 --output d4_results.json

The exact rational control is standard D4 representation theory. It is
not a claim that an arithmetic triality operator has been constructed.
"""
import argparse
from fractions import Fraction as F
import hashlib
import itertools
import json
import os
from pathlib import Path
import sys
import time

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")


def matrix_product(a, b):
    return [[sum(x * y for x, y in zip(row, col))
             for col in zip(*b)] for row in a]


def transpose(a):
    return [list(row) for row in zip(*a)]


def matrix_vector(a, v):
    return tuple(sum(x * y for x, y in zip(row, v)) for row in a)


def triality_control():
    identity = [[F(i == j) for j in range(4)] for i in range(4)]
    h = [[F(x, 2) for x in row] for row in
         [[1, 1, 1, 1], [1, 1, -1, -1],
          [1, -1, 1, -1], [1, -1, -1, 1]]]
    reflection = [[F((1 if i < 3 else -1) if i == j else 0)
                   for j in range(4)] for i in range(4)]
    triality = matrix_product(h, reflection)
    roots = set()
    for i, j in itertools.combinations(range(4), 2):
        for si, sj in itertools.product((-1, 1), repeat=2):
            v = [F(0)] * 4
            v[i], v[j] = F(si), F(sj)
            roots.add(tuple(v))
    vector = set()
    for i in range(4):
        for sign in (-1, 1):
            v = [F(0)] * 4
            v[i] = F(sign)
            vector.add(tuple(v))
    halfspin_even, halfspin_odd = set(), set()
    for signs in itertools.product((-1, 1), repeat=4):
        target = halfspin_even if signs.count(-1) % 2 == 0 else halfspin_odd
        target.add(tuple(F(x, 2) for x in signs))
    weight_sets = {"vector": vector, "halfspin_even": halfspin_even,
                   "halfspin_odd": halfspin_odd}
    mapping = {}
    for name, weights in weight_sets.items():
        image = {matrix_vector(triality, v) for v in weights}
        mapping[name] = next(n for n, w in weight_sets.items() if image == w)
    checks = {
        "orthogonal": matrix_product(transpose(triality), triality) == identity,
        "order_three": matrix_product(matrix_product(triality, triality), triality) == identity,
        "preserves_all_24_D4_roots": {matrix_vector(triality, v) for v in roots} == roots,
        "three_distinct_8_weight_sets": len({frozenset(w) for w in weight_sets.values()}) == 3
            and all(len(w) == 8 for w in weight_sets.values()),
        "cycles_three_representations": all(mapping[n] != n for n in weight_sets),
    }
    assert all(checks.values()), checks
    return {
        "status": "EXACT rational arithmetic; standard D4 control, not an MTFT-derived intertwiner",
        "matrix": [[str(x) for x in row] for row in triality],
        "checks": checks,
        "weight_set_permutation": mapping,
        "implication": "The same positive-definite D4 root metric admits vector and both half-spin labels. The dimension 8 does not distinguish them or select a spacetime signature.",
        "even_Clifford_dimension_check": {
            "Clifford_even_on_4_generators_real_dimension": 2 ** (4 - 1),
            "Mat_2_H_real_dimension": 2 * 2 * 4,
            "claimed_isomorphism_impossible_by_dimension": 2 ** (4 - 1) != 2 * 2 * 4,
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=Path(__file__).with_name("d4_results.json"))
    args = parser.parse_args()
    source = args.source.resolve()
    sys.path.insert(0, str(source / "src"))
    import mtft
    from mtft import liealg
    import numpy as np

    started = time.perf_counter()
    report = liealg.d4_report(screen=True, dps=50)
    elapsed = time.perf_counter() - started
    assert report["dim"] == 28
    assert report["structure"]["killing_signature"] == [28, 0, 0]
    assert report["representation"]["active_dim"] == 8
    assert report["representation"]["common_fixed_dim"] == 5
    assert all(report["symmetry_screen"][w]["normalization_residual"] > .1
               for w in ("W11", "W13", "W143"))
    files = ["src/mtft/liealg.py", "src/mtft/periods/hamiltonian.py",
             "src/mtft/periods/physics.py", "studies/TH2_PREREGISTRATION.md",
             "tests/test_v0230.py", "CHANGELOG_v0230.md"]
    result = {
        "package_version": getattr(mtft, "__version__", "unreported"),
        "source_path": str(source),
        "numpy_version": np.__version__,
        "elapsed_d4_seconds": elapsed,
        "requested_period_precision_dps": 50,
        "matrix_computation": "numpy double precision, despite 50-digit period inputs",
        "status": "DIAGNOSTIC replay of package CERT(tol, E2) fingerprint; no exact reconstruction or interval proof performed",
        "source_sha256": {f: hashlib.sha256((source / f).read_bytes()).hexdigest() for f in files},
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "d4_report": report,
        "exact_standard_triality_control": triality_control(),
        "conditional_fixed_subspace_lemma": {
            "hypothesis": "A nonzero common fixed space F is annihilated by every X in a Lie algebra g.",
            "lie_algebra_consequence": "Every X in g has X^2|F = 0, so X^2 = c I with c nonzero is impossible on the full module.",
            "unital_associative_consequence": "Every element of the unital associative algebra generated by g restricts to a scalar on F. If A and B are invertible then these scalars a and b are nonzero; {A,B}|F = 2ab I cannot vanish in characteristic zero.",
            "scope": "Exact linear-algebra implication conditional on the decomposition. The numerical fixed-five fingerprint is evidence for that hypothesis, not an exact proof of it. Restricting to the active block or adjoining other operators changes the hypothesis.",
        },
        "primary_sources": {
            "triality_three_representations": "https://math.ucr.edu/home/baez/octonions/node7.html",
            "Clifford_algebra_tables": "https://math.ucr.edu/home/baez/octonions/node6.html",
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"output": str(args.output), "dimension": report["dim"],
                      "elapsed_seconds": elapsed,
                      "exact_triality_checks": result["exact_standard_triality_control"]["checks"]}, indent=2))


if __name__ == "__main__":
    main()
