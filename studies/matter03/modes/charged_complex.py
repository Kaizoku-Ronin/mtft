"""Matter-03 charged holomorphic deformation-complex certificate.

The output counts sheaf/hypercohomology classes, not observed particles.
The paired two-form Higgs representative is an input.  Existence of a
Hermitian BPS solution and physical L2 normalization are separate gates.
"""
from __future__ import annotations
import argparse
import csv
from fractions import Fraction as F
import json
from pathlib import Path
import sys


def rref(matrix):
    a = [[F(x) for x in row] for row in matrix]
    pivots = []
    rr = 0
    if not a:
        return a, pivots
    for col in range(len(a[0])):
        pivot = next((i for i in range(rr, len(a)) if a[i][col]), None)
        if pivot is None:
            continue
        a[rr], a[pivot] = a[pivot], a[rr]
        div = a[rr][col]
        a[rr] = [x / div for x in a[rr]]
        for i in range(len(a)):
            if i != rr and a[i][col]:
                k = a[i][col]
                a[i] = [x - k * y for x, y in zip(a[i], a[rr])]
        pivots.append(col)
        rr += 1
        if rr == len(a):
            break
    return a, pivots


def monomials(k):
    return [(i, j) for degree in range(k) for i in range(degree + 1)
            for j in [degree - i]]


def matrix_rank(matrix):
    return len(rref(matrix)[1])


def local_map(a, shift, coefficients):
    """Multiply fat jets and project to the exact target cokernel.

    a is positive except for the special a=0 global map. Source O(aP)
    has one global section with coefficient z**a. Target K_C((a+shift)P)
    has all one-variable jets except the forbidden simple-pole residue.
    These ranks use a<=4 and the ordinary-cusp certificate.
    """
    if a == 0:
        return {"source_basis": [], "target_basis": [], "matrix": [], "rank": 0}
    ks = 3 * a
    kt = 3 * (a + shift) - 1
    if kt <= 0:
        return {"source_basis": [], "target_basis": [], "matrix": [], "rank": 0}
    source = [x for x in monomials(ks) if x != (a, 0)]
    residue_power = a + shift - 1
    target = [x for x in monomials(kt)
              if x[1] > 0 or x == (residue_power, 0)]
    matrix = [[coefficients.get(i - ii, F(0)) if j == jj else F(0)
               for ii, jj in source] for i, j in target]
    return {"source_basis": source, "target_basis": target,
            "matrix": matrix, "rank": matrix_rank(matrix)}


def canonical_up(repo):
    path = repo / "src/mtft/canonical/_data/X0_143_S2_qexpansions.txt"
    lines = [x for x in path.read_text().splitlines()
             if x.strip() and not x.lstrip().startswith("#")]
    rows = list(csv.reader(lines))[1:]
    coeffs = [[int(x) for x in row[1:]] for row in rows]
    equations = [row + [int(n == 3)]
                 for n, row in enumerate(coeffs[1:5], 1)]
    reduced, pivots = rref(equations)
    if len(pivots) != 4 or 13 in pivots:
        raise ValueError("Order-two canonical representative is not certified")
    vector = [F(0)] * 13
    for i, pivot in enumerate(pivots):
        vector[pivot] = reduced[i][-1]
    series = [sum(row[j] * vector[j] for j in range(13)) for row in coeffs]
    if series[1:5] != [0, 0, 1, 0]:
        raise ValueError("Wrong canonical order-two series")
    # canonical differential f(z)dz/z divided by the frame z**2 dz.
    phi = {n - 3: value for n, value in enumerate(series) if value}
    return {"canonical_basis_coefficients": vector,
            "canonical_q_coefficients": series, "phi_coefficients": phi}


def report(repo):
    sys.path.insert(0, str(repo / "src"))
    from mtft.matter.collinear import cohomology, cusp_jet_certificate, yukawa_witness
    cusp = cusp_jet_certificate()
    down = dict(yukawa_witness().phi_coefficients)
    up = canonical_up(repo)
    checks = []

    def check(name, actual, expected):
        passed = actual == expected
        checks.append({"name": name, "actual": actual, "expected": expected,
                       "passed": passed})
        if not passed:
            raise AssertionError(name)

    check("down local leading order", min(down), 5)
    check("down leading coefficient", down[5], F(1))
    check("up local leading order", min(up["phi_coefficients"]), 0)
    check("up leading coefficient", up["phi_coefficients"][0], F(1))
    check("correct N(-3)", cohomology(-3).holo_chiral_dim, 54)
    check("rank-three original minor", yukawa_witness().determinant, F(-1))
    check("ordinary-cusp certificate", cusp["passed"], True)

    plus = {a: local_map(a, 2, down) for a in (0, 1, 2)}
    minus = {a: local_map(a, -2, up["phi_coefficients"]) for a in (2, 3, 4)}
    # The curve residue functional adds one for plus maps a=0,1,2.
    # For minus a=2 its target is H1(K_C), giving one curve contribution.
    # Every map has an independent rank-one elliptic-H1 contribution.
    plus_ranks = {a: plus[a]["rank"] + 2 for a in plus}
    minus_ranks = {a: minus[a]["rank"] + 1 + int(a == 2) for a in minus}
    for a, want in {0: 0, 1: 3, 2: 15}.items():
        check(f"down a={a} exact local rank", plus[a]["rank"], want)
    for a, want in {2: 0, 3: 2, 4: 11}.items():
        check(f"up a={a} exact local rank", minus[a]["rank"], want)
    check("down full ranks", plus_ranks, {0: 2, 1: 5, 2: 17})
    check("up full ranks", minus_ranks, {2: 2, 3: 3, 4: 12})

    # U8 component index, stack exponent, electric charge.
    components = [(f"c{i}", 0, F(1, 6)) for i in range(1, 4)] + [
        ("L+", 1, F(1, 2)), ("L-", 1, F(-1, 2)),
        ("A1", 2, F(-1, 2)), ("A2", -3, F(-1, 2)),
        ("B", -1, F(1, 2))]
    charged = {"-1": [], "+1": []}
    for i, ni, qi in components:
        for j, nj, qj in components:
            if qi - qj not in (1, -1):
                continue
            key = "+1" if qi - qj == 1 else "-1"
            cc = cohomology(ni - nj)
            charged[key].append({"root": f"{i} <- {j}", "q": ni - nj,
                                 "a": cc.h1, "Phi": cc.phi_dim,
                                 "gauge": cc.h0})
    totals = {key: {name: sum(row[name] for row in table)
                    for name in ("a", "Phi", "gauge")}
              for key, table in charged.items()}
    check("Qminus totals", totals["-1"], {"a": 216, "Phi": 63, "gauge": 1})
    check("Qplus totals", totals["+1"], {"a": 217, "Phi": 62, "gauge": 1})
    check("equal total split chirals", totals["-1"]["a"] + totals["-1"]["Phi"],
          totals["+1"]["a"] + totals["+1"]["Phi"])

    blocks = [
        {"spectator": "L-", "Higgs": "down", "a_q": 0, "partner_q": -2,
         "shape": [14, 25], "rank": plus_ranks[0]},
        {"spectator": "L-", "Higgs": "up", "a_q": 2, "partner_q": 0,
         "shape": [32, 14], "rank": minus_ranks[2]},
        {"spectator": "A1", "Higgs": "down", "a_q": 1, "partner_q": -3,
         "shape": [18, 44], "rank": plus_ranks[1]},
        {"spectator": "A1", "Higgs": "up", "a_q": 3, "partner_q": -1,
         "shape": [55, 15], "rank": minus_ranks[3]},
        {"spectator": "A2", "Higgs": "down", "a_q": -4, "partner_q": 2,
         "shape": [72, 32], "rank": plus_ranks[2]},
        {"spectator": "A2", "Higgs": "up", "a_q": -2, "partner_q": 4,
         "shape": [25, 87], "rank": minus_ranks[4]}]
    for block in blocks:
        check(f"{block['spectator']} {block['Higgs']} exponent closure",
              block["a_q"] + block["partner_q"] + (2 if block["Higgs"] == "down" else -2), 0)
        check(f"{block['spectator']} {block['Higgs']} rank bound",
              block["rank"] <= min(block["shape"]), True)
    pair_rank = sum(block["rank"] for block in blocks)
    max_rank = sum(min(block["shape"]) for block in blocks)
    check("paired a-map exact rank", pair_rank, 41)
    check("paired a-map universal size bound", max_rank, 118)
    hyper1 = {key: totals[key]["a"] - pair_rank + totals[key]["Phi"] - 1
              for key in charged}
    check("Qminus paired hyperH1", hyper1["-1"], 237)
    check("Qplus paired hyperH1", hyper1["+1"], 237)

    return {"status": "EXACT CONDITIONAL HOLOMORPHIC deformation complex",
            "claim_type": "Sheaf/hypercohomology dimensions; no physical masses",
            "assumptions": ["Specified Matter-02 cusp representative and bulk BHV field dictionary",
                            "Nonzero neutral paired Phi with canonical orders down=3, up=2",
                            "a background expectation value zero",
                            "Two-term complex End(V) --adPhi--> K_S tensor End(V)",
                            "Hermitian BPS existence is established separately for the stable paired branch; the actual metric is not computed"],
            "charged_roots": charged, "split_totals": totals,
            "down_maps": plus, "up_maps": minus,
            "down_full_ranks": plus_ranks, "up_full_ranks": minus_ranks,
            "six_a_mass_blocks": blocks, "paired_a_map_rank": pair_rank,
            "rank_size_bound": max_rank, "paired_hyperH1": hyper1,
            "up_representative": up, "cusp_certificate": cusp,
            "checks": checks, "gates_passed": len(checks),
            "proof_review": "Root and vacuum investigators independently reviewed the Leray extension argument; analytic contributions are proved, not numerically estimated.",
            "open": ["explicit Hermitian metric", "L2 Gram matrices and singular values",
                     "Additional interactions to remove excess charged modes",
                     "Observed light-family selection"]}


def safe(obj):
    if isinstance(obj, F):
        return str(obj)
    if isinstance(obj, dict):
        return {str(k): safe(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [safe(v) for v in obj]
    return obj


def default_repo():
    for ancestor in Path(__file__).resolve().parents:
        if (ancestor / "src/mtft/matter/collinear.py").is_file():
            return ancestor
        candidate = ancestor / "matter03/repo_candidate"
        if (candidate / "src/mtft/matter/collinear.py").is_file():
            return candidate
    raise RuntimeError("Pass --repo to locate the MTFT candidate source")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=default_repo())
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = report(args.repo.resolve())
    output = json.dumps(safe(result), indent=2) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(output)
        print(f"PASS: {result['gates_passed']} exact gates; paired charged hyperH1 =237 per sign")
    else:
        print(output, end="")
