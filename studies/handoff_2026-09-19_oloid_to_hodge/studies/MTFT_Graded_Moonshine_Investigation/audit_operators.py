#!/usr/bin/env python3
"""Exact checks of candidate operators on the X0(143) canonical ring.

Run: python audit_operators.py --output operator_results.json
Uses Python's standard library and the two frozen files in ./data/.

Pinned frame: AL-adapted e1,...,e13 from MTFT 0.32.0.  Variables and
monomial pairs in reported relations are ONE-based.  Quadrics within a
sector are also numbered from ONE.  Sector order is (++),(+-),(-+),(--),
with degree-one dimensions (1,6,5,1).  Degree is canonical-ring degree,
not q-exponent.  q is the Fourier parameter at the cusp infinity.

For p=2 away from level 143, the weight-2 Hecke formula is
  (T2 f)[n] = f[2n] + 2 f[n/2], second term zero for n odd.
Thus q-expansions through q^140 give exact T2 coefficients through q^70.
This is enough for exact nonzero witnesses, with no numerical tolerance.
Zero residuals are finite controls; the shipped quadrics' global identity
is supported by MTFT's separate Sturm-bound certificate.

The substitution e_i -> T2(e_i) tests Sym^2(T2) preservation of I2.
The Leibniz extension D_T(e_i)=T2(e_i) tests derivation preservation of I2.
These are distinct candidate lifts; neither is the weight-4 Hecke action.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import json
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).resolve().parent
NVAR = 13
NMAX = 70
SECTOR_ORDER = ("(+,+)", "(+,-)", "(-,+)", "(-,-)")
SECTORS = ((1, 1),) + ((1, -1),) * 6 + ((-1, 1),) * 5 + ((-1, -1),)
MONOMIALS = tuple(itertools.combinations_with_replacement(range(NVAR), 2))
EXPECTED_HASHES = {
    "X0_143_AL_adapted_qexpansions.txt": "cdb9f3caf6fb6eb529931f41a2b47f0657489328737eb5abf0ab82150b56eb52",
    "X0_143_I2_by_AL_sector.txt": "13a6b37c188bc0f4f49f2d2b95f432f15adaa00fc9535276695f2c32df14335d",
}


def load_data(data_dir):
    hashes = {}
    for name, expected in EXPECTED_HASHES.items():
        got = hashlib.sha256((data_dir / name).read_bytes()).hexdigest()
        if got != expected:
            raise ValueError(f"Frozen fixture hash mismatch: {name}")
        hashes[name] = got
    lines = (data_dir / "X0_143_AL_adapted_qexpansions.txt").read_text().splitlines()
    rows = list(csv.reader(line for line in lines if line and not line.startswith("#")))
    E = [[int(v) for v in row[1:]] for row in rows[1:]]
    assert [int(row[0]) for row in rows[1:]] == list(range(141))
    assert len(E) == 141 and all(len(row) == NVAR for row in E)

    quadrics, sector, columns = [], None, None
    for line in (data_dir / "X0_143_I2_by_AL_sector.txt").read_text().splitlines():
        line = line.strip()
        if line.startswith("## class"):
            sector, columns = line.split()[2], None
            continue
        if not line or line.startswith("#"):
            continue
        parts = line.split(",")
        if parts[0] == "monomial":
            columns = [dict() for _ in parts[1:]]
            quadrics.extend((sector, i + 1, col) for i, col in enumerate(columns))
            continue
        label = parts[0].replace("^2", "")
        indices = [int(v) - 1 for v in label.split("y") if v]
        if len(indices) == 1:
            indices *= 2
        m = MONOMIALS.index(tuple(sorted(indices)))
        for col, value in zip(columns, parts[1:]):
            if int(value):
                col[m] = int(value)
    assert len(quadrics) == 55
    assert [sum(sec == s for sec, _, _ in quadrics) for s in SECTOR_ORDER] == [26, 5, 4, 20]
    return E, quadrics, hashes


def product(A, B, i, j):
    return [sum(A[k][i] * B[n - k][j] for k in range(n + 1))
            for n in range(NMAX + 1)]


def residual_report(products, quadrics):
    nonzero = []
    for sector, number, poly in quadrics:
        residual = [sum(a * products[m][n] for m, a in poly.items())
                    for n in range(NMAX + 1)]
        nz = [(n, value) for n, value in enumerate(residual) if value]
        if nz:
            nonzero.append({
                "sector": sector, "quadric_number_one_based": number,
                "first_nonzero_coefficients": [{"q_power": n, "coefficient": v} for n, v in nz[:4]],
                "relation_terms": [{"variables_one_based": [i + 1, j + 1], "coefficient": a}
                                   for m, a in poly.items() for i, j in [MONOMIALS[m]]],
            })
    return {
        "quadrics_checked": len(quadrics),
        "zero_residuals_through_q70": len(quadrics) - len(nonzero),
        "nonzero_residuals": len(nonzero),
        "first_witness": nonzero[0] if nonzero else None,
        "shortest_witness": min(nonzero, key=lambda x: len(x["relation_terms"])) if nonzero else None,
    }


def eta11(N=140):
    """Exact q eta-product: q prod_(n>=1)(1-q^n)^2(1-q^(11n))^2."""
    a = [1] + [0] * N
    for shift in (1, 11):
        for _ in range(2):
            for d in range(shift, N + 1, shift):
                for k in range(N, d - 1, -1):
                    a[k] -= a[k - d]
    return [0] + a[:N]


def audit(data_dir):
    E, quadrics, hashes = load_data(data_dir)
    T = [[E[2*n][j] + (2*E[n//2][j] if n % 2 == 0 else 0)
          for j in range(NVAR)] for n in range(NMAX + 1)]
    original = [product(E, E, i, j) for i, j in MONOMIALS]
    multiplication = [product(T, T, i, j) for i, j in MONOMIALS]
    derivation = [[a + b for a, b in zip(product(T, E, i, j), product(E, T, i, j))]
                  for i, j in MONOMIALS]
    report = {
        "source_release": "MTFT 0.32.0",
        "source_path_within_release": "src/mtft/canonical/_data",
        "fixture_sha256": hashes,
        "arithmetic": "exact Python integers; no floating point",
        "degree": "canonical-ring degree; distinct from Fourier q-exponent",
        "sector_order": list(SECTOR_ORDER),
        "degree_one_sector_dimensions": [1, 6, 5, 1],
        "quadrics_per_sector": [26, 5, 4, 20],
        "q_coefficients_checked_inclusive": [0, NMAX],
        "identity_control": residual_report(original, quadrics),
        "Sym2_T2_substitution": residual_report(multiplication, quadrics),
        "Leibniz_D_T2": residual_report(derivation, quadrics),
        "Atkin_Lehner_controls": {},
    }
    for label, powers in [("W11", (1, 0)), ("W13", (0, 1)), ("W143", (1, 1))]:
        signs = [a**powers[0] * b**powers[1] for a, b in SECTORS]
        transformed = [[signs[i] * signs[j] * a for a in original[m]]
                       for m, (i, j) in enumerate(MONOMIALS)]
        report["Atkin_Lehner_controls"][label] = residual_report(transformed, quadrics)

    # Euler grading D_N(e_i)=e_i gives D_N(Q)=2Q, hence a genuine derivation.
    report["Euler_grading_derivation_control"] = residual_report(
        [[2*a for a in row] for row in original], quadrics)

    # Independent eta-product check fixes the normalized oldform and convention.
    g = eta11()
    old = [Fraction(E[n][7] + E[n][12], 2) for n in range(len(E))]
    assert old == g
    tg = [g[2*n] + (2*g[n//2] if n % 2 == 0 else 0) for n in range(NMAX + 1)]
    assert all(tg[n] == -2*g[n] for n in range(NMAX + 1))
    g2 = [sum(g[k] * g[n-k] for k in range(n + 1)) for n in range(len(g))]
    # Weight-4 Hecke action uses factor 2^(4-1)=8.
    t4g2_q1 = g2[2]
    t2g_squared_q1 = sum(tg[k] * tg[1-k] for k in range(2))
    report["simple_weighted_Hecke_counterexample"] = {
        "form": "g=eta(tau)^2 eta(11 tau)^2=(e8+e13)/2",
        "eta_product_matches_shipped_form_through_q140": True,
        "T2_weight2_g_equals_minus2g_through_q70": True,
        "g_first_coefficients_q1_through_q6": g[1:7],
        "coefficient_q1_of_T2_weight4_g_squared": t4g2_q1,
        "coefficient_q1_of_T2_weight2_g_then_squared": t2g_squared_q1,
        "conclusion": "T2 across weights is not multiplicative. T2(weight4)(g^2) has cusp order 1, so does not represent a holomorphic section of 2K at this cusp.",
    }

    # Serre D_k=q*d/dq - k E2/12 is modular, but changes cusp vanishing.
    report["Serre_derivative_cusp_gate"] = {
        "operator": "D2(g)=q dg/dq - E2*g/6",
        "input_canonical_degree": 1,
        "output_modular_weight": 4,
        "coefficient_q1": "5/6",
        "minimum_q_order_for_holomorphic_2K_section": 2,
        "preserves_canonical_ring": False,
        "reason": "At cusp infinity d tau=dq/(2*pi*i*q); coefficient g=O(q) gives a pole for g*(d tau)^2.",
    }
    assert report["identity_control"]["nonzero_residuals"] == 0
    assert report["Sym2_T2_substitution"]["nonzero_residuals"] == 55
    assert report["Leibniz_D_T2"]["nonzero_residuals"] == 55
    assert all(r["nonzero_residuals"] == 0 for r in report["Atkin_Lehner_controls"].values())
    assert report["Euler_grading_derivation_control"]["nonzero_residuals"] == 0
    assert (t4g2_q1, t2g_squared_q1) == (1, 0)
    report["interpretation"] = (
        "These failures obstruct these two specific T2 lifts to the canonical ring. "
        "They do not obstruct Hecke correspondences on H1 or weightwise modular-form spaces, "
        "and do not rule out a different state-space construction. A scalar shift T2-lambda*Id "
        "cannot repair the derivation obstruction: D_(T2-lambda*Id)(Q)=D_T2(Q)-2lambda*Q."
    )
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=HERE / "data")
    parser.add_argument("--output", type=Path, default=HERE / "operator_results.json")
    args = parser.parse_args()
    result = audit(args.data_dir)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"output": str(args.output), "identity_zero": 55,
                      "Sym2_T2_failures": 55, "Leibniz_D_T2_failures": 55,
                      "AL_controls": "W11,W13,W143 each preserve all 55 relations",
                      "Euler_grading_control": "preserves all 55 relations",
                      "weighted_T2_counterexample_q1": [1, 0]}, indent=2))


if __name__ == "__main__":
    main()
