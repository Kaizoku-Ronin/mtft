#!/usr/bin/env python3
"""Exact rational CRT projectors for the four packaged X0(143) Hecke blocks.

The projectors are derived from T2 and pairwise coprime polynomials, without
using the package's block bases. Those bases provide a second construction.
No floating-point arithmetic or active-sector data enter this calculation.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys
import time

import sympy as sp


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def rational_matrix(rows):
    return sp.Matrix([[sp.Rational(v.numerator, v.denominator)
                       if hasattr(v, "numerator") else sp.Rational(v)
                       for v in row] for row in rows])


def polynomial_matrix(poly, matrix):
    answer = sp.zeros(matrix.rows)
    identity = sp.eye(matrix.rows)
    for coefficient in poly.all_coeffs():
        answer = answer * matrix + coefficient * identity
    return answer


def serialize_matrix(matrix):
    return [[str(v) for v in row] for row in matrix.tolist()]


def run(output_dir):
    from mtft import hecke as H
    from mtft.periods import involutions

    started = time.perf_counter()
    x = sp.Symbol("x")
    factors = {
        "ell": sp.Poly(x, x, domain=sp.QQ),
        "old": sp.Poly(x + 2, x, domain=sp.QQ),
        "q4": sp.Poly(x**4 - 3*x**3 - x**2 + 5*x + 1, x, domain=sp.QQ),
        "q6": sp.Poly(x**6 - 10*x**4 + 2*x**3 + 24*x**2 - 7*x - 12,
                       x, domain=sp.QQ),
    }
    dims = {"ell": 2, "old": 4, "q4": 8, "q6": 12}
    minimal_poly = sp.prod(factors.values())
    t2 = rational_matrix(H.cuspidal_hecke(2))
    zero, identity = sp.zeros(26), sp.eye(26)
    assert minimal_poly.degree() == 12
    assert sp.gcd(minimal_poly, minimal_poly.diff()).degree() == 0
    assert polynomial_matrix(minimal_poly, t2) == zero
    assert all(p.is_irreducible for p in factors.values())
    expected_charpoly = factors["ell"]**2 * factors["old"]**4
    expected_charpoly *= factors["q4"]**2 * factors["q6"]**2
    assert sp.Poly(t2.charpoly(x).as_expr(), x, domain=sp.QQ) == expected_charpoly

    projectors, crt_polys, checks = {}, {}, {}
    for name, factor in factors.items():
        complement = minimal_poly.exquo(factor)
        inverse = sp.invert(complement, factor)
        crt = (complement * inverse).rem(minimal_poly)
        projector = polynomial_matrix(crt, t2)
        residues = {other: str(crt.rem(other_factor).as_expr())
                    for other, other_factor in factors.items()}
        assert residues == {other: "1" if other == name else "0"
                            for other in factors}
        assert projector * projector == projector
        assert projector.rank() == dims[name]
        assert projector.trace() == dims[name]
        projectors[name], crt_polys[name] = projector, crt
        checks[name] = {"dimension_real": dims[name],
                        "dimension_complex": dims[name] // 2,
                        "factor": str(factor.as_expr()),
                        "crt_polynomial": str(crt.as_expr()),
                        "crt_residues": residues,
                        "idempotent_exact": True,
                        "rank_exact": dims[name],
                        "trace_exact": str(projector.trace())}
    assert sum(projectors.values(), zero) == identity
    assert all(a * b == zero for i, a in enumerate(projectors.values())
               for j, b in enumerate(projectors.values()) if i != j)

    # Independent nullspace-based package construction, with column convention.
    blocks = H.blocks()
    basis = sp.Matrix.hstack(*[rational_matrix(blocks[name]).T
                              for name in factors])
    assert basis.rank() == 26
    inverse_basis = basis.inv()
    offset = 0
    for name in factors:
        selector = sp.zeros(26)
        for i in range(offset, offset + dims[name]):
            selector[i, i] = 1
        assert basis * selector * inverse_basis == projectors[name]
        checks[name]["matches_package_nullspace_projector_exact"] = True
        offset += dims[name]

    operators = {f"T{p}": rational_matrix(H.cuspidal_hecke(p))
                 for p in (2, 3, 5, 7, 17, 19)}
    operators.update({f"W{q}": rational_matrix(involutions.al_matrix(q))
                      for q in (11, 13, 143)})
    operators["STAR"] = rational_matrix(H.star_involution())
    operator_checks = {}
    for name, operator in operators.items():
        commutations = {block: operator * projector == projector * operator
                        for block, projector in projectors.items()}
        assert all(commutations.values())
        coordinate_operator = inverse_basis * operator * basis
        cursor, block_rows = 0, {}
        reconstructed = sp.zeros(26)
        for block, dimension in dims.items():
            block_matrix = coordinate_operator[cursor:cursor+dimension,
                                                cursor:cursor+dimension]
            reconstructed[cursor:cursor+dimension, cursor:cursor+dimension] = block_matrix
            block_rows[block] = {"trace_exact": str(block_matrix.trace())}
            if name.startswith("W") or name == "STAR":
                assert block_matrix * block_matrix == sp.eye(dimension)
                plus = (dimension + block_matrix.trace()) / 2
                block_rows[block].update({"plus_eigenspace_dimension": int(plus),
                                         "minus_eigenspace_dimension": int(dimension-plus)})
            cursor += dimension
        assert reconstructed == coordinate_operator
        operator_checks[name] = {"commutes_with_each_projector_exact": commutations,
                                 "all_cross_hecke_blocks_zero_exact": True,
                                 "block_data": block_rows}

    # Explanatory follow-up to the registered observed W143 coupling rank.
    # This uses no numerical active projector and proves an upper bound for
    # every complementary pair P,Q, not just the current local decomposition.
    w143 = operators["W143"]
    plus143 = (identity + w143) / 2
    assert plus143 * plus143 == plus143
    assert w143 == 2 * plus143 - identity
    assert plus143.rank() == 4
    plus_support = {name: (p * plus143).rank()
                    for name, p in projectors.items()}
    assert plus_support == {"ell": 2, "old": 2, "q4": 0, "q6": 0}
    assert projectors["ell"] * plus143 == projectors["ell"]
    assert projectors["q4"] * plus143 == zero
    assert projectors["q6"] * plus143 == zero
    assert (projectors["ell"] + projectors["old"]) * plus143 == plus143
    assert projectors["q4"] * w143 == -projectors["q4"]
    assert projectors["q6"] * w143 == -projectors["q6"]
    w143_lemma = {
        "status": "EXACT explanatory follow-up, not a preregistered endpoint",
        "plus_projector_definition": "R_plus=(I+W143)/2",
        "plus_projector_idempotent_exact": True,
        "reflection_identity_exact": "W143=2 R_plus-I",
        "plus_projector_rank_exact": 4,
        "plus_projector_rank_by_hecke_block_exact": plus_support,
        "plus_projector_supported_in_ell_plus_old_exact": True,
        "q4_and_q6_W143_scalar_exact": -1,
        "universal_coupling_identity": "For any P^2=P and Q=I-P: Q W143 P = 2 Q R_plus P.",
        "universal_coupling_rank_bound": 4,
        "proof": "QP=0 cancels the identity part of W143. Rank of a product cannot exceed rank(R_plus)=4. Orthogonality of P is not required.",
        "cross_q4_q6_commutator_identity": "For every matrix P, E_q4 [P,W143] E_q6 = E_q6 [P,W143] E_q4 = 0, since W143 acts as -I on both blocks.",
        "limitation": "This certifies the rank upper bound and exact zero commutator blocks for packaged W143; saturation at rank four for the selected active projector remains a numerical result.",
    }

    result = {
        "epistemic_status": "EXACT rational arithmetic on the packaged matrices",
        "mtft_version": __import__("mtft").__version__,
        "matrix_basis": "mtft.hecke cuspidal homology, column action, dimension 26 over Q",
        "minimal_polynomial": str(minimal_poly.as_expr()),
        "minimal_polynomial_annihilates_exact": True,
        "minimal_polynomial_squarefree_exact": True,
        "charpoly_factorization_exact": str(expected_charpoly.as_expr()),
        "projectors_sum_to_identity_exact": True,
        "projectors_pairwise_algebraically_orthogonal_exact": True,
        "orthogonality_scope": "E_i E_j = 0; this check alone does not claim Euclidean orthogonal projectors in raw coordinates.",
        "blocks": checks,
        "operators": operator_checks,
        "W143_rank_explanation": w143_lemma,
        "interpretation": "Every tested arithmetic operator preserves the four Hecke blocks. Active/fixed mixing under these operators is caused by a different, incompatible decomposition; it is not transition between distinct Hecke blocks.",
        "source_caveat": "The packaged W matrices were originally integer-recognized from high-precision periods. Here their algebraic block-preservation properties are checked exactly; geometric identification is inherited from the package.",
        "source_sha256": {"hecke.py": sha256(H.__file__),
                          "periods/involutions.py": sha256(involutions.__file__),
                          "X0_143_atkin_lehner_v022.json": sha256(involutions.data_path("X0_143_atkin_lehner_v022.json")),
                          "exact_blocks.py": sha256(__file__)},
        "elapsed_seconds": time.perf_counter() - started,
    }
    matrices = {"matrix_basis": result["matrix_basis"],
                "encoding": "Each entry is an exact rational string accepted by fractions.Fraction or sympy.Rational.",
                "projectors": {name: serialize_matrix(p) for name, p in projectors.items()},
                "crt_polynomials_coefficients_ascending": {
                    name: [str(v) for v in reversed(poly.all_coeffs())]
                    for name, poly in crt_polys.items()}}
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "exact_blocks_results.json").write_text(json.dumps(result, indent=2)+"\n")
    (output_dir / "exact_projectors.json").write_text(json.dumps(matrices, indent=2)+"\n")
    print(json.dumps({"all_exact_checks_passed": True, "block_ranks": dims,
                      "operators_checked": list(operators),
                      "elapsed_seconds": result["elapsed_seconds"]}, indent=2))
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path)
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).resolve().parent)
    args = parser.parse_args()
    if args.source:
        sys.path.insert(0, str(args.source.resolve() / "src"))
    run(args.output_dir.resolve())
