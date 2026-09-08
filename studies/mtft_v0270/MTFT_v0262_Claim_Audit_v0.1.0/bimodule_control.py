#!/usr/bin/env python3
"""Exact finite doubled-algebra control for the v0.26.2 review.

This is a synthetic diagonal R^3 construction, NOT an X0(143) construction.
It tests the algebraic conditions measured by surface.bimodule plus selected
reality/grading conditions. It does not certify every spectral-triple axiom.
On complexification the real structure is swap followed by conjugation.

All control matrix arithmetic uses Python integers, without floating point.
Optional --verify-theta rechecks the finite characteristic census in both
packaged frames; it does not evaluate theta functions or run TH-2.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys


def eye(n):
    return [[int(i == j) for j in range(n)] for i in range(n)]


def zero(n):
    return [[0] * n for _ in range(n)]


def transpose(a):
    return [list(row) for row in zip(*a)]


def add(a, b):
    return [[x + y for x, y in zip(ar, br)] for ar, br in zip(a, b)]


def neg(a):
    return [[-x for x in row] for row in a]


def mul(a, b):
    bt = transpose(b)
    return [[sum(x * y for x, y in zip(ar, bc)) for bc in bt] for ar in a]


def comm(a, b):
    return add(mul(a, b), neg(mul(b, a)))


def block(a, b, c, d):
    return [ar + br for ar, br in zip(a, b)] + [cr + dr for cr, dr in zip(c, d)]


def norm_squared(a):
    return sum(x * x for row in a for x in row)


def control(name, w):
    n = 3
    i, z = eye(n), zero(n)
    z6 = zero(2 * n)
    j = block(z, i, i, z)
    grading = block(i, z, z, neg(i))
    wt = transpose(w)
    m = add(w, wt)
    d = block(z, m, transpose(m), z)
    alphabet = [[[int(r == c == k) for c in range(n)] for r in range(n)] for k in range(n)]
    left = [block(a, z, z, mul(mul(w, a), wt)) for a in alphabet]
    opposite = [block(mul(mul(w, transpose(a)), wt), z, z, transpose(a)) for a in alphabet]
    forms = [comm(d, a) for a in left]
    checks = {
        "orthogonal_twist": mul(w, wt) == i,
        "faithful_idempotent_algebra_representation": all(
            mul(left[r], left[c]) == (left[r] if r == c else z6)
            for r in range(n) for c in range(n)
        ),
        "order_zero_all_generator_pairs": all(comm(a, b) == z6 for a in left for b in opposite),
        "first_order_all_generator_pairs": all(comm(f, b) == z6 for f in forms for b in opposite),
        "self_adjoint_D": d == transpose(d),
        "J_squared_identity": mul(j, j) == eye(2 * n),
        "JD_equals_DJ": mul(j, d) == mul(d, j),
        "grading_anticommutes_D": mul(grading, d) == neg(mul(d, grading)),
        "grading_anticommutes_J": mul(grading, j) == neg(mul(j, grading)),
    }
    assert all(checks.values()), (name, checks)
    form_norms = [norm_squared(f) for f in forms]
    expected = [4, 4, 4] if name == "three_cycle" else [0, 0, 0]
    assert form_norms == expected, (name, form_norms)
    return {
        "arithmetic": "EXACT: Python integer matrix operations",
        "W": w,
        "M_equals_W_plus_transpose_W": m,
        "checks": checks,
        "one_form_Frobenius_norms_squared": form_norms,
        "has_nonzero_represented_one_forms": any(form_norms),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True, help="Extracted mtft-0.26.2 source directory")
    parser.add_argument("--output", type=Path, default=Path(__file__).with_name("bimodule_control_results.json"))
    parser.add_argument("--verify-theta", action="store_true")
    args = parser.parse_args()
    source = args.source.resolve()
    tracked = [
        "src/mtft/surface/bimodule.py",
        "studies/TH2_PREREGISTRATION.md",
        "src/mtft/thetachar/__init__.py",
        "src/mtft/thetachar/_data/x0143_theta_census.json",
        "src/mtft/liealg.py",
        "CHANGELOG_v0230.md",
        "CHANGELOG_v0260.md",
    ]
    hashes = {rel: hashlib.sha256((source / rel).read_bytes()).hexdigest() for rel in tracked}
    prereg = (source / "studies/TH2_PREREGISTRATION.md").read_text()
    th2_paths = sorted(str(p.relative_to(source)) for p in source.rglob("*")
                       if p.is_file() and ("th2" in p.name.lower() or "th-2" in p.name.lower()))
    out = {
        "scope": "Synthetic R^3 doubled-algebra control; NOT an X0(143) construction or a signature derivation.",
        "construction": {
            "algebra": "R^3 represented by all diagonal 3x3 matrices; complexification available",
            "left_action": "pi(a)=diag(a, W a W^T)",
            "opposite_action": "b^o=diag(W b^T W^T, b^T)",
            "Dirac": "D=[[0,M],[M^T,0]], M=W+W^T",
            "real_structure": "J=swap on real matrices; swap composed with conjugation after complexification",
            "grading": "Gamma=diag(I,-I)",
            "commutator": "ordinary [D,pi(a)], not a sigma-twisted commutator",
        },
        "cases": {
            "identity": control("identity", eye(3)),
            "transposition": control("transposition", [[0, 1, 0], [1, 0, 0], [0, 0, 1]]),
            "three_cycle": control("three_cycle", [[0, 0, 1], [1, 0, 0], [0, 1, 0]]),
        },
        "interpretation": [
            "A noninvolutive permutation can produce nonzero ordinary one-forms while satisfying the tested order-zero, first-order, reality and grading conditions.",
            "This changes the algebra and does not provide a noninvolutive rational field automorphism of the quartic or sextic Hecke fields.",
            "Aut_Q(K4)=Aut_Q(K6)={1} constrains field-automorphism twists on those factors; it is not a theorem excluding all bimodules, larger algebras, multiplicity-space transformations or scalar extensions.",
            "All represented commutators vanishing implies no represented inner one-forms in that construction. It does not identify a physical flat gauge connection.",
        ],
        "TH2_status": {
            "preregistration_present": True,
            "preregistration_states_no_primary_quantities_computed_at_filing": "No quantity defined below has been computed at filing time" in prereg,
            "filenames_matching_TH2_or_TH_minus_2": th2_paths,
            "review_of_supplied_release": "Only preregistration found; no packaged primary gradient ledger or H-A/H-B/H-C verdict found. This does not exclude unpublished external results.",
            "what_H_A_would_establish": "All 32 gradient support ratios satisfy the preregistered tolerance for the active 8, in a correctly matched frame and metric.",
            "what_H_A_would_not_establish": "D4 equivariance of the gradient construction, vector-versus-spinor labeling under triality, or a real Clifford signature.",
            "roundoff_note": "Preregistration certifies truncation tails but labels floating-point roundoff estimate NOT certified; final intervals inherit this limitation.",
            "theta_function_evaluations_performed_in_this_control": 0,
        },
        "source_sha256": hashes,
    }
    if args.verify_theta:
        sys.path.insert(0, str(source / "src"))
        from mtft import thetachar
        counts = {}
        for name, make in [("periods", thetachar.x0143_periods_frame),
                           ("independent_packaged_GP_frame", thetachar.x0143_gp_frame)]:
            action = make()
            points, parities = action.invariant_characteristics()
            _, basis = action.joint_fixed_locus()
            thetachar.verify_action(action, samples=50, seed=143)
            counts[name] = {"affine_dimension": int(basis.shape[1]),
                            "invariants": len(points), "even": int((parities == 0).sum()),
                            "odd": int((parities == 1).sum()),
                            "affine_action_50_samples_per_generator": "PASS"}
            assert (counts[name]["affine_dimension"], counts[name]["even"], counts[name]["odd"]) == (7, 96, 32)
        out["TH2_status"]["finite_input_census_rechecked"] = counts
    else:
        out["TH2_status"]["finite_input_census_rechecked"] = "Not run; use --verify-theta"
    args.output.write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps({"output": str(args.output.resolve()), "control_cases": 3,
                      "three_cycle_one_form_norms_squared": out["cases"]["three_cycle"]["one_form_Frobenius_norms_squared"],
                      "theta_functions_evaluated": 0}, indent=2))


if __name__ == "__main__":
    main()
