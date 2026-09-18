#!/usr/bin/env python3
"""SC7-01: exact tests of a declared smooth, commuting, twisted 7D extension.

This script reads the supplied MTFT source without editing it. All matrix
ranks and geometric identities use exact arithmetic. Its negative matter
verdict applies to the model specified in PREREGISTRATION.md.
"""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
from functools import lru_cache
import hashlib
import itertools
import json
from math import gcd
from pathlib import Path
import platform
import sys

import sympy as sp
from sympy.matrices.normalforms import smith_normal_form


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def clean(x):
    return sp.simplify(sp.expand_complex(x))


def matrix_strings(m):
    return [[str(x) for x in row] for row in m.tolist()]


def say(message):
    print(message, flush=True)


def commutator(a, b):
    # Generator indices are one based; a negative index means inverse.
    return [a, b, -a, -b]


def seifert_relators(genus, euler):
    fiber = 2 * genus + 1
    central = [commutator(j, fiber) for j in range(1, fiber)]
    surface = []
    for i in range(genus):
        surface += commutator(2 * i + 1, 2 * i + 2)
    surface += [-fiber] * euler
    return central + [surface]


def fox_row(word, character):
    """Left Fox derivative evaluated in a rank-one complex character."""
    row = [sp.S.Zero] * len(character)
    prefix = sp.S.One
    for letter in word:
        j = abs(letter) - 1
        if letter > 0:
            row[j] += prefix
            prefix = clean(prefix * character[j])
        else:
            prefix = clean(prefix / character[j])
            row[j] -= prefix
    return [clean(x) for x in row], clean(prefix)


@lru_cache(maxsize=None)
def low_cohomology(genus, euler, character):
    n = 2 * genus + 1
    d0 = sp.Matrix([clean(x - 1) for x in character])
    rows = []
    for relator in seifert_relators(genus, euler):
        row, value = fox_row(relator, character)
        assert value == 1, "Input does not define a representation of pi_1."
        rows.append(row)
    d1 = sp.Matrix(rows)
    assert (d1 * d0).applyfunc(clean) == sp.zeros(n, 1)
    r0, r1 = d0.rank(), d1.rank()
    return {"rank_d0": r0, "rank_d1": r1,
            "h0": 1 - r0, "h1": n - r0 - r1}


def all_cohomology(genus, euler, character):
    low = low_cohomology(genus, euler, tuple(character))
    dual_character = tuple(clean(1 / x) for x in character)
    dual = low_cohomology(genus, euler, dual_character)
    # Poincare duality on the closed oriented 3-manifold, with dual local
    # coefficients, supplies the dimensions in degrees 2 and 3.
    dims = [low["h0"], low["h1"], dual["h1"], dual["h0"]]
    assert sum((-1) ** i * v for i, v in enumerate(dims)) == 0
    return {**low, "dual_rank_d0": dual["rank_d0"],
            "dual_rank_d1": dual["rank_d1"], "cohomology_dimensions": dims}


def genus_and_topology(source):
    N = 143
    primes = [11, 13]
    index = sp.Integer(N)
    for p in primes:
        index *= 1 + sp.Rational(1, p)
    elliptic2 = sp.prod(1 + sp.kronecker_symbol(-1, p) for p in primes)
    elliptic3 = sp.prod(1 + sp.kronecker_symbol(-3, p) for p in primes)
    cusps = 4
    genus = int(1 + index / 12 - elliptic2 / 4 - elliptic3 / 3 - sp.Rational(cusps, 2))
    basis = json.loads((source / "src/mtft/surface/_data/x0143_weight2_basis.json").read_text())
    assert genus == 13 == len(basis["coefficients"])
    euler = genus - 1
    n = 2 * genus + 1
    abelian_relations = []
    for word in seifert_relators(genus, euler):
        row = [0] * n
        for letter in word:
            row[abs(letter) - 1] += 1 if letter > 0 else -1
        abelian_relations.append(row)
    snf = smith_normal_form(sp.Matrix(abelian_relations), domain=sp.ZZ)
    nonzero = [abs(int(snf[i, i])) for i in range(n) if snf[i, i]]
    assert nonzero == [12]
    # Gysin: H^0(X,Z) -- cup e --> H^2(X,Z) is multiplication by 12.
    gysin_snf = smith_normal_form(sp.Matrix([[euler]]), domain=sp.ZZ)
    assert int(gysin_snf[0, 0]) == nonzero[0]
    return {
        "genus": genus, "euler_class": euler, "modular_index": int(index),
        "cusps": cusps, "elliptic_points": [int(elliptic2), int(elliptic3)],
        "abelianized_relations_smith_nonzero": nonzero,
        "gysin_cup_e_matrix": [[euler]],
        "homology_Z": ["Z", "Z^26 + Z/12", "Z^26", "Z"],
        "cohomology_Z": ["Z", "Z^26", "Z^26 + Z/12", "Z"],
        "base_pullback_degree_map": "m -> m mod 12 in the torsion summand of H^2",
        "fiber_in_pi1": "infinite central cyclic subgroup",
        "fiber_in_H1": "order 12"
    }, basis


def spectrum_scan(genus, euler, stacks, rep_table):
    # All roots live exactly in Q(i,sqrt(3)); no floating-point tolerance.
    zeta = (sp.sqrt(3) + sp.I) / 2
    roots = [clean(zeta ** (-m)) for m in range(euler)]
    assert len(set(roots)) == euler and all(clean(x ** euler) == 1 for x in roots)
    scenarios = {
        "trivial_base": [sp.S.One] * (2 * genus),
        "unitary_nontrivial_base": [-sp.S.One] + [sp.S.One] * (2 * genus - 1),
        "complex_nonunitary_base": [sp.Integer(2), sp.Integer(3)] + [sp.S.One] * (2 * genus - 2)
    }
    scans = {}
    for name, base in scenarios.items():
        rows = []
        for residue, lam in enumerate(roots):
            data = all_cohomology(genus, euler, [*base, lam])
            expected = ([1, 2 * genus, 2 * genus, 1] if name == "trivial_base"
                        else [0, 2 * genus - 2, 2 * genus - 2, 0])
            if residue:
                expected = [0, 0, 0, 0]
            assert data["cohomology_dimensions"] == expected
            rows.append({"residue": residue, "fiber_holonomy": str(lam), **data})
        scans[name] = rows
        say(f"Exact Fox scan: {name}, all {euler} residues match the independent fiber/Gysin calculation.")

    ordered = list(stacks["m"])
    pairs = []
    for x, y in itertools.combinations(ordered, 2):
        degree = stacks["m"][x] - stacks["m"][y]
        residue = degree % euler
        data = scans["trivial_base"][residue]
        pairs.append({
            "ordered_pair": [x, y], "signed_degree": degree, "residue": residue,
            "torsion_order": euler // gcd(euler, degree),
            "fiber_holonomy": str(roots[residue]),
            "zero_modes_for_trivial_base_holonomy": data["cohomology_dimensions"],
            "zero_modes_for_nontrivial_base_holonomy":
                scans["unitary_nontrivial_base"][residue]["cohomology_dimensions"]
        })
    model = {"degrees": stacks["m"],
             "y": {x: Fraction(v) for x, v in stacks["Y"].items()}}
    oriented = []
    for row in rep_table(model):
        # rep_table conjugates negative-degree sectors; use the displayed
        # representation's positive degree, not its retained original sign.
        degree = sum(row["charges"][x] * stacks["m"][x] for x in ordered)
        assert degree == row["multiplicity"] > 0
        oriented.append({
            "sector": row["sector"], "representation": [str(x) for x in row["rep"]],
            "six_dimensional_multiplicity": row["multiplicity"],
            "oriented_bundle_degree": degree, "fiber_holonomy": str(roots[degree % euler]),
            "seven_dimensional_cohomology_dimensions":
                scans["trivial_base"][degree % euler]["cohomology_dimensions"]
        })
    assert len(oriented) == 8
    assert all(x["seven_dimensional_cohomology_dimensions"] == [0, 0, 0, 0] for x in oriented)
    return {"exact_field": "Q(i,sqrt(3))", "characters_scanned": 36,
            "scans": scans, "all_ten_stack_pairs": pairs, "oriented_M1_sectors": oriented,
            "statement": "All nonzero-degree M1 sectors are acyclic for every smooth commuting flat/Higgs background on the fixed line bundles."}


def symmetric_square(m):
    a, b, c, d = list(m)
    return sp.Matrix([[a*a, a*b, b*b], [2*a*c, a*d+b*c, 2*b*d],
                      [c*c, c*d, d*d]]).applyfunc(sp.simplify)


def group_summary(generators):
    ident = sp.eye(2)
    seen = {tuple(ident): ident}
    pending = [ident]
    while pending:
        left = pending.pop()
        for right in generators:
            product = (left * right).applyfunc(sp.simplify)
            if tuple(product) not in seen:
                seen[tuple(product)] = product
                pending.append(product)
                assert len(seen) <= 8
    orders = Counter(next(n for n in range(1, 9) if (m**n).applyfunc(sp.simplify) == ident)
                     for m in seen.values())
    return {"order": len(seen), "element_order_counts": dict(sorted(orders.items()))}


def symmetry_check(basis):
    A = sp.diag(sp.I, -sp.I)
    B = sp.Matrix([[0, -1/sp.sqrt(13)], [sp.sqrt(13), 0]])
    ident = sp.eye(2)
    assert A*A == B*B == (A*B)**2 == -ident
    assert (A*B*A.inv()*B.inv()).applyfunc(sp.simplify) == -ident
    QA = sp.diag(-1, 1, -1)
    QB = sp.Matrix([[0, 0, sp.Rational(1, 13)], [0, -1, 0], [13, 0, 0]])
    assert symmetric_square(A) == QA
    assert symmetric_square(B) == QB
    group = group_summary([A, B])
    assert group["element_order_counts"] == {1: 1, 2: 1, 4: 6}
    W11 = sp.Matrix(basis["W11"]).applyfunc(sp.Rational)
    W13 = sp.Matrix(basis["W13"]).applyfunc(sp.Rational)
    assert W11**2 == W13**2 == sp.eye(13) and W11*W13 == W13*W11
    traces = {name: 2*int(W.trace())
              for name, W in [("W11", W11), ("W13", W13), ("W143", W11*W13)]}
    assert traces == {"W11": 2, "W13": -2, "W143": -18}
    return {"theta_compatible_A": matrix_strings(A), "theta_compatible_B": matrix_strings(B),
            "group": group, "square_map_W11": matrix_strings(QA),
            "square_map_W13": matrix_strings(QB),
            "traces_on_H1_P_R": traces,
            "central_half_turn_on_H1": "identity (fiber rotation is isotopic to identity)",
            "original_M1_W11_invariance": False,
            "full_Q8_equivariance_of_flat_7D_background": "not established",
            "M1_note": "W11 is broken by the original odd-degree CM divisor. Choosing a flat 7D background does not automatically retain that divisor or establish a restored symmetry.",
            "Hecke_correspondence_lift": "not established"}


def curvature_check(euler):
    # Local metric for a circle bundle over the compact uniformized surface:
    # R^2(dx^2+dy^2)/y^2 + r^2(dtheta+k dx/y)^2, k=e/24.
    x, y, theta = sp.symbols("x y theta", real=True)
    R, r = sp.symbols("R r", positive=True)
    k = sp.symbols("k", real=True)
    coords = [x, y, theta]
    metric = sp.Matrix([[(R**2+r**2*k**2)/y**2, 0, r**2*k/y],
                        [0, R**2/y**2, 0], [r**2*k/y, 0, r**2]])
    inv = metric.inv().applyfunc(sp.simplify)
    gamma = [[[sp.simplify(sum(
        inv[a, d] * (sp.diff(metric[d, b], coords[c]) +
                     sp.diff(metric[d, c], coords[b]) -
                     sp.diff(metric[b, c], coords[d])) / 2
        for d in range(3))) for c in range(3)] for b in range(3)] for a in range(3)]
    ricci = sp.zeros(3)
    for a in range(3):
        for b in range(3):
            ricci[a, b] = sp.simplify(sum(
                sp.diff(gamma[c][a][b], coords[c]) - sp.diff(gamma[c][a][c], coords[b]) +
                sum(gamma[c][a][b] * gamma[d][c][d] -
                    gamma[d][a][c] * gamma[c][b][d] for d in range(3))
                for c in range(3)))
    scalar = sp.simplify(sum(inv[a, b] * ricci[a, b] for a in range(3) for b in range(3)))
    predicted = -2/R**2 - r**2*k**2/(2*R**4)
    assert sp.simplify(scalar - predicted) == 0
    assert sp.simplify(metric.det() - R**4*r**2/y**4) == 0
    # Dual orthonormal vector frame to (R dx/y, R dy/y, r(dtheta+k dx/y)).
    frame = sp.Matrix([[y/R, 0, 0], [0, y/R, 0], [-k/R, 0, 1/r]])
    assert (frame.T*metric*frame).applyfunc(sp.simplify) == sp.eye(3)
    ricci_frame = (frame.T*ricci*frame).applyfunc(sp.simplify)
    expected_ricci = sp.diag(-1/R**2-r**2*k**2/(2*R**4),
                            -1/R**2-r**2*k**2/(2*R**4),
                            r**2*k**2/(2*R**4))
    assert (ricci_frame - expected_ricci).applyfunc(sp.simplify) == sp.zeros(3)
    value = sp.simplify(scalar.subs(k, sp.Rational(euler, 24)))
    assert sp.simplify(value + 2/R**2 + r**2/(8*R**4)) == 0

    rho, a = sp.symbols("rho a", positive=True)
    b, c = sp.symbols("b c", nonnegative=True)
    # Internal dimension 3; rescaling to 4D Einstein frame gives powers
    # -(n+2), -(n+4), -n for curvature, internal two-form flux, bulk energy.
    potential = a/rho**5 + b/rho**7 + c/rho**3
    derivative = sp.diff(potential, rho)
    assert derivative.is_negative is True
    return {"metric": matrix_strings(metric),
            "connection_k": str(sp.Rational(euler, 24)),
            "metric_determinant": str(sp.factor(metric.det())),
            "Ricci_in_orthonormal_frame": matrix_strings(ricci_frame),
            "scalar_curvature_general_k": str(scalar),
            "scalar_curvature_e12": str(value),
            "volume": "96*pi^2*R^2*r",
            "M4_squared": "M7^5 * (96*pi^2*R^2*r)",
            "inverse_g4_squared": "(96*pi^2*R^2*r) / g7_squared",
            "Einstein_frame_uniform_scale_potential": str(potential),
            "uniform_scale_derivative": str(derivative),
            "derivative_strictly_negative": True,
            "actual_flat_BPS_background_flux_coefficient_b": 0,
            "status": "CONDITIONAL: additional unwarped 7D Einstein action and nonnegative classical energy assumptions; not a derived supergravity solution"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--archive", type=Path)
    parser.add_argument("--output", type=Path, default=Path("results.json"))
    args = parser.parse_args()
    source = args.source.resolve()
    sys.path.insert(0, str(source / "src"))
    import mtft
    from mtft.surface.smflux import M1_STACKS, rep_table
    assert Path(mtft.__file__).resolve().is_relative_to(source)
    tracked = ["src/mtft/surface/smflux.py", "src/mtft/surface/arithspin.py",
               "src/mtft/surface/_data/x0143_weight2_basis.json"]
    before = {p: digest(source / p) for p in tracked}
    result = {"experiment": "SC7-01", "recorded_date": "2026-09-16",
              "source": str(source), "python": platform.python_version(),
              "sympy": sp.__version__, "input_sha256": before,
              "preregistration_sha256": digest(Path(__file__).with_name("PREREGISTRATION.md")),
              "script_sha256": digest(__file__),
              "scope": "closed spin-circle bundle, standard twisted 7D SYM, smooth commuting independent stack backgrounds",
              "M1_stack_data": M1_STACKS}
    if args.archive:
        result["archive_sha256"] = digest(args.archive)
        assert result["archive_sha256"] == "2e6fcd0fde3b1601e3cf58aaec338790aa760034e548ecded3086bfa2a110db5"
    result["topology"], basis = genus_and_topology(source)
    say("Topology: genus 13, e=12; Gysin and presentation/SNF agree.")
    result["spectrum"] = spectrum_scan(13, 12, M1_STACKS, rep_table)
    result["symmetry"] = symmetry_check(basis)
    say("Symmetry: Q8 square map retained; original M1 breaks W11; no full 7D vacuum symmetry is inferred.")
    result["gravity"] = curvature_check(12)
    say("Gravity: direct Ricci calculation agrees; restricted Einstein-frame scale potential is strictly decreasing.")
    assert before == {p: digest(source / p) for p in tracked}
    result["source_files_unchanged"] = True
    result["computation_status"] = "ALL_IMPLEMENTED_EXACT_CHECKS_PASSED"
    result["physical_gate"] = "FAIL: none of the eight nonzero-degree M1 matter sectors has a massless twisted 7D bulk mode"
    result["radius_gate"] = "FAIL within the separately declared nonnegative-energy unwarped Einstein ansatz"
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    say(result["physical_gate"])
    say(f"Results written to {args.output}")


if __name__ == "__main__":
    main()
