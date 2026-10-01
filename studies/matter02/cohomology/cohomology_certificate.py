#!/usr/bin/env python3
"""Exact certificate for the specified Matter-02 collinear holomorphic bundle.

Only Python's standard library is used. All row reduction is over Fraction.
No period matrix, metric, or particle-spectrum identification is assumed.
"""
from __future__ import annotations
import argparse
import csv
from fractions import Fraction
import hashlib
from itertools import combinations
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
GENUS = 13
EXPECTED_SHA1 = {
    "X0_143_S2_qexpansions.txt": "479316a0952d72dbd4731f50389003263cbfc4e1",
    "X0_143_AL_adapted_qexpansions.txt": "e71eaabe5cd23f3a3a35c204373ad6dfd5115358",
    "X0_143_AL_adapted_basis.txt": "5e4537c805a48e76027d4d12bdf0cea300cd1472",
}
STACK_POWERS = {"c": 0, "L": 1, "A1": 2, "A2": -3, "B": -1}


def load_matrix(path):
    rows = list(csv.reader(line for line in path.read_text().splitlines()
                           if line.strip() and not line.lstrip().startswith("#")))
    return [[int(value) for value in row[1:]] for row in rows[1:]]


def rref(matrix):
    if not matrix:
        return [], []
    a = [[Fraction(x) for x in row] for row in matrix]
    ncols = len(a[0])
    pivots = []
    rr = 0
    for col in range(ncols):
        pivot = next((i for i in range(rr, len(a)) if a[i][col]), None)
        if pivot is None:
            continue
        a[rr], a[pivot] = a[pivot], a[rr]
        scale = a[rr][col]
        a[rr] = [value / scale for value in a[rr]]
        for i in range(len(a)):
            if i != rr and a[i][col]:
                factor = a[i][col]
                a[i] = [u - factor * v for u, v in zip(a[i], a[rr])]
        pivots.append(col)
        rr += 1
        if rr == len(a):
            break
    return a, pivots


def determinant(matrix):
    a = [[Fraction(x) for x in row] for row in matrix]
    det = Fraction(1)
    for col in range(len(a)):
        pivot = next((i for i in range(col, len(a)) if a[i][col]), None)
        if pivot is None:
            return Fraction(0)
        if pivot != col:
            a[col], a[pivot] = a[pivot], a[col]
            det *= -1
        scale = a[col][col]
        det *= scale
        for i in range(col + 1, len(a)):
            factor = a[i][col] / scale
            for j in range(col, len(a)):
                a[i][j] -= factor * a[col][j]
    return det


def kernel(matrix, ncols=GENUS):
    a, pivots = rref(matrix)
    free = [j for j in range(ncols) if j not in pivots]
    basis = []
    for col in free:
        vector = [Fraction(0)] * ncols
        vector[col] = Fraction(1)
        for row, pivot in enumerate(pivots):
            vector[pivot] = -a[row][col]
        basis.append(vector)
    return basis


def chi(q):
    return (-9 * q * q - 3 * q) // 2


def cohomology(q):
    if not isinstance(q, int) or isinstance(q, bool):
        raise TypeError("q must be an integer")
    if q == 0:
        return (1, 14, 13)
    h2 = max(GENUS - q, 0) if q > 0 else max(GENUS + 1 + 2 * q, 0)
    return (0, h2 - chi(q), h2)


def curve_h0(degree):
    """h0 O_X(degree P) at the certified ordinary cusp."""
    if degree < 0:
        return 0
    return max(1, degree + 1 - GENUS)


def curve_k_h0(minimum_vanishing):
    return max(GENUS - minimum_vanishing, 0)


def run(data_dir=HERE / "data"):
    data_dir = Path(data_dir)
    gates = []
    def check(name, truth, detail=None):
        row = {"name": name, "status": "PASS" if truth else "FAIL"}
        if detail is not None:
            row["detail"] = detail
        gates.append(row)
        if not truth:
            raise AssertionError(name)

    hashes = {}
    for filename, expected in EXPECTED_SHA1.items():
        raw = (data_dir / filename).read_bytes()
        actual_sha1 = hashlib.sha1(raw).hexdigest()
        hashes[filename] = {"sha1": actual_sha1, "sha256": hashlib.sha256(raw).hexdigest()}
        check("archived provenance " + filename, actual_sha1 == expected)
    f = load_matrix(data_dir / "X0_143_S2_qexpansions.txt")
    g = load_matrix(data_dir / "X0_143_AL_adapted_qexpansions.txt")
    b = load_matrix(data_dir / "X0_143_AL_adapted_basis.txt")
    check("original basis shape", len(f) == 141 and all(len(row) == 13 for row in f))
    check("adapted basis shape", len(g) == 141 and all(len(row) == 13 for row in g))
    check("basis change shape", len(b) == 13 and all(len(row) == 13 for row in b))
    check("both q0 rows vanish", f[0] == g[0] == [0] * 13)
    transformed = [[sum(row[k] * b[k][j] for k in range(13))
                    for j in range(13)] for row in f]
    check("exact basis transformation across all141 rows", transformed == g)
    det_b = determinant(b)
    check("archived change of basis determinant", det_b == -1078272)

    jets = []
    for m in range(1, 14):
        original = f[1:m+1]
        adapted = g[1:m+1]
        _, pivot_cols = rref(original)
        _, adapted_pivots = rref(adapted)
        minor = [[row[c] for c in pivot_cols] for row in original]
        det = determinant(minor)
        check(f"cusp jet rank m={m} both bases", len(pivot_cols) == len(adapted_pivots) == m)
        check(f"cusp jet nonzero rational minor m={m}", det != 0)
        basis = kernel(original)
        check(f"cusp jet nullity m={m}", len(basis) == 13-m)
        check(f"cusp jet kernel annihilation m={m}",
              all(sum(row[i] * vec[i] for i in range(13)) == 0
                  for row in original for vec in basis))
        jets.append({"minimum_canonical_vanishing": m, "rank": m,
                     "nullity": 13-m, "pivot_columns_0based": pivot_cols,
                     "minor_determinant": str(det),
                     "rational_kernel_basis_original": [[str(x) for x in vec] for vec in basis]})
    det13 = determinant(f[1:14])
    check("full canonical jet determinant", det13 == -18196056815164784640)

    # Independent local residue check of the Yukawa rank-three witness.
    # Take a rational canonical form of exact order3, scaled by a4=1.
    order3_candidates = kernel(f[1:4])
    phi_vector = next(vec for vec in order3_candidates
                      if sum(f[4][i]*vec[i] for i in range(13)))
    a4 = sum(f[4][i]*phi_vector[i] for i in range(13))
    phi_vector = [value/a4 for value in phi_vector]
    canonical_coefficients = [sum(row[i]*phi_vector[i] for i in range(13)) for row in f]
    check("canonical exact order3 witness", canonical_coefficients[1:4] == [0,0,0]
          and canonical_coefficients[4] == 1)
    # Relative to K_C tensor O(2P), f(z)dz/z contributes coefficient z^2*f(z)/z.
    phi_coefficient = {n+1: value for n,value in enumerate(canonical_coefficients) if value}
    positive_jets = [(0,1),(0,2),(1,1)]
    negative_principal_parts = [(6,2),(7,2),(6,3)]
    residue_matrix = [[phi_coefficient.get(ii-i-1,Fraction(0)) if jj-j == 1 else Fraction(0)
                       for ii,jj in negative_principal_parts] for i,j in positive_jets]
    check("rank3 local residue matrix", residue_matrix == [[1,phi_coefficient.get(6,0),0],
                                                           [0,0,1],[0,1,0]])
    check("rank3 local residue determinant", determinant(residue_matrix) == -1)
    residue_witness = {"canonical_basis_coefficients_original": [str(x) for x in phi_vector],
                       "canonical_a4": "1", "canonical_a5":str(canonical_coefficients[5]),
                       "positive_jets_z_w_powers":positive_jets,
                       "negative_principal_parts_minus_z_w_powers":negative_principal_parts,
                       "residue_matrix":[[str(x) for x in row] for row in residue_matrix],
                       "determinant":"-1", "meaning":"algebraic coordinate normalization; no L2 metric supplied"}

    table = []
    for q in range(-5, 6):
        h0, h1, h2 = cohomology(q)
        check(f"RR chi q={q}", h0-h1+h2 == chi(q))
        # Geometry calculations do not use the production cohomology formula.
        # Positive q uses a length-k fat point with k=3q. Evaluation is injective,
        # because a section of degree q cannot vanish to order 3q at P.
        if q > 0:
            curve0 = curve_h0(q)
            curve1 = curve_k_h0(q)
            fat_length = sum(i+1 for i in range(3*q))
            route_h1 = curve1 + curve0 + fat_length - curve0
            route_h2 = curve1
            route = {"method": "negative exceptional Leray and fat-point sequence",
                     "fat_point_order": 3*q, "fat_point_length": fat_length,
                     "evaluation_rank": curve0, "base_h1": curve1+curve0,
                     "h1": route_h1, "h2": route_h2}
        elif q < 0:
            r = -q
            curve1 = r+12
            r1_length = sum(k-1 for k in range(1,3*r+1))
            dual_vanishing = 2*r-1
            dual_jet_rank = min(dual_vanishing,13)
            route_h2 = 13-dual_jet_rank
            connecting_rank = curve1-route_h2
            route_h1 = curve1+r1_length-connecting_rank
            route = {"method": "positive exceptional Leray with dual jet connecting map",
                     "R1_direct_image_length": r1_length,
                     "base_h1": curve1, "base_h2": curve1,
                     "connecting_map_rank": connecting_rank,
                     "dual_minimum_canonical_vanishing": dual_vanishing,
                     "h1": route_h1, "h2": route_h2}
        else:
            route = {"method": "birational invariance and Kunneth", "h1": 14, "h2": 13}
        check(f"Leray/jet route agrees q={q}", (h1,h2) == (route["h1"],route["h2"]))
        table.append({"q":q, "class_(a,b,e)":[q,0,-3*q], "chi":chi(q),
                      "h0":h0,"h1":h1,"h2":h2,
                      "h0_K_tensor_Lq": cohomology(-q)[2], "route": route})

    hom_rows = []
    for i, j in combinations(STACK_POWERS, 2):
        q = STACK_POWERS[i]-STACK_POWERS[j]
        forward = cohomology(q)
        reverse = cohomology(-q)
        check(f"net index orientation {i}/{j}", chi(q)-chi(-q) == -3*q)
        hom_rows.append({"i":i,"j":j,"q_i_minus_q_j":q,
                         "forward_cohomology":list(forward),"reverse_cohomology":list(reverse),
                         "chi_forward_minus_reverse":chi(q)-chi(-q)})
    # The mathematical theorem covers all powers because ordinary-cusp canonical
    # orders are exactly 0..12. This wider grid checks the arithmetic API edges.
    for q in range(-25,26):
        hh = cohomology(q)
        check(f"all-power RR/nonnegativity q={q}", min(hh)>=0 and hh[0]-hh[1]+hh[2]==chi(q))
    return {"milestone":"Matter-02", "status":"EXACT CONDITIONAL GEOMETRY",
            "geometry": {"surface":"Bl_(P,Q)(X0(143) x 143a1)",
                         "P":"rational cusp infinity, width1",
                         "Q":"identity of elliptic curve 143a1",
                         "line_bundle":"pi*pr_X*O_X(P) tensor O(-3E_exc)",
                         "flat_twists":"none", "coefficient_field":"Q, then base change to C",
                         "stack_powers":STACK_POWERS},
            "formula": {"chi":"(-9q^2-3q)/2", "q0":[1,14,13],
                        "h0_nonzero_q":0,
                        "h2_positive_q":"max(13-q,0)",
                        "h2_negative_q":"max(14+2q,0)",
                        "h1_nonzero_q":"h2+(9q^2+3q)/2"},
            "input_hashes":hashes, "cusp_gap_sequence":list(range(1,14)),
            "canonical_vanishing_sequence":list(range(13)),
            "basis_change_determinant":str(det_b), "full13_jet_determinant":str(det13),
            "jet_certificates":jets, "rank3_yukawa_witness":residue_witness,
            "cohomology_table":table,"off_diagonal_hom_table":hom_rows,
            "unweighted_bookkeeping":{
                "h1_all20_Hom_orientations":sum(row["forward_cohomology"][1]+row["reverse_cohomology"][1] for row in hom_rows),
                "h2_all20_Hom_orientations":sum(row["forward_cohomology"][2]+row["reverse_cohomology"][2] for row in hom_rows),
                "note":"Unweighted vector-space dimensions, not particle multiplicities."},
            "gates":gates,"gate_count":len(gates),
            "open":["8D twisted field dictionary and its representation orientation",
                    "normalized Hodge/fermion bases and Yukawa integrals",
                    "Higgs vev and vector-like lifting",
                    "physical massless photon and coupling matching"]}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=HERE/"data")
    parser.add_argument("--output", type=Path, default=HERE/"cohomology_exact.json")
    args=parser.parse_args()
    result=run(args.data_dir)
    args.output.write_text(json.dumps(result,indent=2)+"\n")
    print(f"PASS: {result['gate_count']} exact arithmetic gates; ordinary cusp gap sequence1..13")
    print("q     h0    h1    h2    chi")
    for row in result["cohomology_table"]:
        print(f"{row['q']:2d} {row['h0']:6d} {row['h1']:5d} {row['h2']:5d} {row['chi']:6d}")
    print("Written:",args.output)

if __name__=="__main__":
    main()
