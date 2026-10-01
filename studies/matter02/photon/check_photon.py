#!/usr/bin/env python3
"""Matter-02 exact photon-kernel audit and explicitly conditional completions.

No Green--Schwarz/Stueckelberg action is inferred from flux projections or
anomaly rows.  Every mass kernel below belongs to a named INPUT probe action.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import sympy as sp


def serial(value):
    if isinstance(value, sp.MatrixBase):
        return [[serial(x) for x in value.row(i)] for i in range(value.rows)]
    if isinstance(value, sp.Basic):
        return str(value)
    if isinstance(value, dict):
        return {str(k): serial(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [serial(v) for v in value]
    return value


def audit(repo: Path):
    sys.path.insert(0, str(repo / "src"))
    from mtft.research import anomaly8d, anomalies
    from mtft.surface import blowup

    names = ["c", "L", "A1", "A2", "B"]
    ranks = sp.Matrix([3, 2, 1, 1, 1])
    powers = sp.Matrix([0, 1, 2, -3, -1])
    Y = sp.Matrix([sp.Rational(1, 6), 0, -sp.Rational(1, 2),
                   -sp.Rational(1, 2), sp.Rational(1, 2)])
    center = sp.ones(5, 1)
    v = (1, 0, -3)
    v2 = sp.Integer(blowup.dot(v, v))
    s = sp.Symbol("s", real=True)
    Ys = Y + s * center
    gates = []

    def check(name, condition, classification="EXACT conditional arithmetic"):
        passed = bool(condition)
        gates.append({"name": name, "pass": passed, "classification": classification})
        if not passed:
            raise AssertionError(name)

    def moment(q, p):
        return sp.expand(sum(ranks[i] * q[i] * powers[i] ** p for i in range(5)))

    check("collinear_self_intersection", v2 == -9)
    check("traceless_flux", moment(center, 1) == 0)
    check("traceless_hypercharge", moment(Y, 0) == 0)
    check("linear_flux_projection_vanishes", moment(Y, 1) == 0)
    check("quadratic_flux_projection_is_54", moment(Y, 2) * v2 == 54)
    check("center_shift_linear_projection_unchanged", moment(Ys, 1) == 0)
    check("center_shift_quadratic_projection", moment(Ys, 2) * v2 == 54 - 144 * s)
    check("center_shift_trace_projection", moment(Ys, 0) == 8 * s)
    check("root_charges_center_invariant", all(sp.expand(Ys[i] - Ys[j] - Y[i] + Y[j]) == 0
                                               for i in range(5) for j in range(5)))
    check("hypercharge_kinetic_trace", sum(ranks[i] * Ys[i] ** 2 for i in range(5)).expand()
          == sp.Rational(5, 6) + 8 * s ** 2)

    linear_row = sp.Matrix([[ranks[i] * powers[i] for i in range(5)]])
    quadratic_row = sp.Matrix([[ranks[i] * powers[i] ** 2 for i in range(5)]])
    trace_row = ranks.T
    K12 = linear_row.col_join(quadratic_row)
    K012 = trace_row.col_join(K12)
    check("linear_quadratic_probe_rank", K12.rank() == 2)
    check("shifted_hypercharge_in_linear_quadratic_kernel", K12 * Ys.subs(s, sp.Rational(3, 8)) == sp.zeros(2, 1))
    check("trace_plus_quadratic_has_no_SM_hypercharge_lift", sp.solve(list(K012 * Ys), [s]) == [])

    # Independent oriented adjoint net-index ledger.  This is not a supplied GS
    # matrix.  The positive chiral index is phi_i - phi_j = -3(n_i-n_j).
    q = sp.Matrix(sp.symbols("q_c q_L q_A1 q_A2 q_B"))
    phi = -3 * powers
    P = sp.expand(sum((phi[i] - phi[j]) * ranks[i] * ranks[j] * (q[i] - q[j]) ** 3
                      for i in range(5) for j in range(i + 1, 5)))
    Agr = sp.expand(sum((phi[i] - phi[j]) * ranks[i] * ranks[j] * (q[i] - q[j])
                        for i in range(5) for j in range(i + 1, 5)))
    mixed = [sp.expand(sum(sp.Rational(1, 2) * ranks[j] * (phi[k] - phi[j]) * (q[k] - q[j])
                           for j in range(5) if j != k)) for k in (0, 1)]
    A = sp.Matrix([[sp.diff(p, z) for z in q] for p in mixed + [Agr]])
    model = {"order": names, "N": list(ranks), "m": list(phi)}
    published = anomalies.anomaly_polynomials(model, list(q))
    check("independent_cubic_matches_repository", sp.expand(P - published["cubic"]) == 0)
    check("independent_gravity_matches_repository", sp.expand(Agr - published["gravity"]) == 0)
    check("independent_mixed_matches_repository", all(sp.expand(mixed[i] - published["mixed"][k]) == 0
                                                      for i, k in enumerate(("c", "L"))))
    check("mixed_anomaly_matrix_rank_two", A.rank() == 2)
    check("mixed_anomaly_hypercharge_zero", A * Y == sp.zeros(3, 1))
    check("mixed_anomaly_center_zero", A * center == sp.zeros(3, 1))
    check("cubic_hypercharge_zero", P.subs(dict(zip(q, Y))) == 0)
    T = sp.Matrix.hstack(*A.nullspace())
    t = sp.Matrix(sp.symbols("t0 t1 t2"))
    restricted = sp.factor(P.subs(dict(zip(q, T * t))))
    check("mixed_kernel_not_cubic_anomaly_free", restricted != 0)
    check("restricted_cubic_polynomial", restricted == -72 * (t[0] - t[1]) ** 2 * (4 * t[0] - t[1] - 3 * t[2]))
    tr_center_Y = sp.Matrix.hstack(center, Y)
    a, b = sp.symbols("a b")
    check("center_hypercharge_plane_cubic_zero", sp.expand(P.subs(dict(zip(q, tr_center_Y * sp.Matrix([a, b]))))) == 0)

    # An example INPUT axion-shift matrix: two mixed-anomaly rows, an extra
    # cubic-compatible row, then a trace row to remove the decoupled center.
    # This demonstrates what a completion must achieve, not what B3 derives.
    Lrow = sp.Matrix([[3, -8, -1, 4, 2]])
    cubic_row = sp.Matrix([[0, 0, 1, -1, 0]])
    K3 = linear_row.col_join(Lrow).col_join(cubic_row)
    K4 = K3.col_join(trace_row)
    check("anomaly_row_spaces_match", linear_row.col_join(Lrow).col_join(A).rank() == 2)
    check("input_K3_rank_and_kernel", K3.rank() == 3 and K3 * tr_center_Y == sp.zeros(3, 2))
    check("input_K4_rank_and_hypercharge_kernel", K4.rank() == 4 and K4 * Y == sp.zeros(4, 1))
    Kneutral = sp.zeros(K4.rows, 1).row_join(K4)
    Hd = sp.Matrix([[sp.Rational(1, 2), 0, 1, 0, 0, -1]])
    Hu = sp.Matrix([[-sp.Rational(1, 2), 0, 1, -1, 0, 0]])
    Q = sp.Matrix([1] + list(Y))
    combined = Kneutral.col_join(Hd)
    check("neutral_Hd_is_Q_neutral", Hd * Q == sp.zeros(1, 1))
    check("neutral_Hu_is_Q_neutral", Hu * Q == sp.zeros(1, 1))
    check("input_neutral_combined_exactly_one_kernel", combined.rank() == 5 and combined * Q == sp.zeros(5, 1))
    check("second_neutral_Higgs_preserves_same_Q", combined.col_join(Hu).rank() == 5)

    # Positive definite supplied kinetic matrices cannot repair a failed
    # kernel; an exact nonorthogonal canonical change makes that explicit.
    G = sp.diag(1, 2, 3, 4, 5)
    Hhalf = sp.eye(6)
    Hhalf[0, 1] = sp.Rational(1, 3)
    Hhalf[2, 4] = -sp.Rational(2, 5)
    H = Hhalf.T * Hhalf
    C = Hhalf.inv()
    M = combined.T * G * combined
    canonical_M = C.T * M * C
    check("supplied_gauge_kinetic_is_positive", all(sp.det(H[:i, :i]) > 0 for i in range(1, 7)))
    check("canonical_kinetic_identity", C.T * H * C == sp.eye(6))
    check("canonical_mass_rank_unchanged", canonical_M.rank() == M.rank() == 5)
    check("canonical_mass_kernel_transformed", canonical_M * Hhalf * Q == sp.zeros(6, 1))

    c1_baseline = [(0, 0, 0), (-3, 1, 21), (-1, 0, 8), (4, -1, -29), (3, -1, -21)]
    baseline_linear = sp.Matrix([sum(ranks[i] * Y[i] * c1_baseline[i][j] for i in range(5)) for j in range(3)])
    baseline_quadratic = sum(ranks[i] * Y[i] * blowup.dot(c1_baseline[i], c1_baseline[i]) for i in range(5))
    check("Matter01_baseline_linear_projection", baseline_linear == sp.zeros(3, 1))
    check("Matter01_baseline_quadratic_projection", baseline_quadratic == 233)
    I10 = anomaly8d.I10([(anomaly8d.adjoint("A", 8), 1)])
    check("ordinary_8D_adjoint_local_I10_zero", I10 == 0)
    check("ordinary_8D_adjoint_local_GS_count_zero", anomaly8d.gs_decomposition(I10)["green_schwarz_fields_needed"] == 0)

    traces = {str(p): moment(Y, p) for p in range(5)}
    mixed_moments = {f"Y^{u}_n^{p}": sum(ranks[i] * Y[i] ** u * powers[i] ** p for i in range(5))
                     for u, p in [(1, 1), (1, 2), (1, 3), (1, 4), (2, 1), (2, 2), (3, 1), (4, 0)]}
    source_files = ["src/mtft/surface/blowup.py", "src/mtft/research/anomaly8d.py",
                    "src/mtft/research/anomalies.py", "src/mtft/research/charge_lattices.py",
                    "src/mtft/parent.py", "docs/SM/PV02_BLOWUP_SUSY_REGISTER.md"]
    hashes = {rel: hashlib.sha256((repo / rel).read_bytes()).hexdigest() for rel in source_files}
    return serial({
        "study": "Matter-02 collinear photon-kernel audit", "all_pass": all(g["pass"] for g in gates),
        "exact_gates_passed": len(gates), "gates": gates, "source_hashes": hashes,
        "frozen": {"stack_order": names, "ranks": ranks, "line_powers": powers, "v_(a,b,e)": v,
                   "v_squared": v2, "Y_traceless": Y, "Y_center_lift": Ys},
        "flux_projections": {"Tr_Y_c1_class": [0, 0, 0], "integral_Tr_Y_c1_squared": 54,
                             "integral_Tr_Ys_c1_squared": 54 - 144 * s, "Tr_Ys": 8 * s,
                             "Tr_Ys_squared": sp.Rational(5, 6) + 8 * s ** 2,
                             "matrix_moments_Tr_Y_n_power": traces, "mixed_matrix_moments": mixed_moments,
                             "pure_internal_c1_cubed_and_fourth": "zero in real dimension 4; matrix moments 9 and -48 are not those integrals",
                             "formal_Tr_F_cubed_fY_barF_squared_coefficient_integrated": 162,
                             "formal_Tr_F_fourth_fY_squared_barF_squared_coefficient_integrated": -189,
                             "baseline_Matter01_quadratic": baseline_quadratic},
        "ordinary_parent_local_anomaly": {"I10_adjoint": I10, "GS_fields_for_this_polynomial": 0,
                         "status": "EXACT local polynomial only; spin^c/R-twist/global and 4D anomaly completion OPEN"},
        "4D_anomaly": {"mixed_row_matrix": A, "rank": A.rank(), "nullspace": A.nullspace(),
                       "cubic_polynomial": P, "cubic_on_mixed_kernel": restricted,
                       "nullspace_parameter_basis": T, "mixed_kernel_is_full_anomaly_free_space": False,
                       "hypercharge_and_center_plane_cubic_zero": True},
        "INPUT_probe_kernels": {
            "linear_flux_only": {"K": linear_row, "rank": linear_row.rank(), "nullspace": linear_row.nullspace(), "Y_kernel": True},
            "linear_plus_quadratic_flux": {"K": K12, "rank": K12.rank(), "nullspace": K12.nullspace(), "Y_kernel": False,
                                          "same_SM_charges_shifted_kernel": Ys.subs(s, sp.Rational(3, 8))},
            "trace_linear_quadratic": {"K": K012, "rank": K012.rank(), "nullspace": K012.nullspace(), "any_Y_plus_center_lift": False},
            "anomaly_compatible_K3": {"K": K3, "rank": 3, "kernel_basis": [Y, center]},
            "anomaly_compatible_K4_trace_center_lifted": {"K": K4, "rank": 4, "kernel_basis": [Y]},
            "K4_plus_neutral_Hd": {"basis": ["T3", *names], "stacked_constraint_matrix": combined,
                                   "rank": 5, "kernel_basis": [Q], "gauge_mass_rank_unchanged_by_positive_kinetics": True}},
        "physical_photon_status": "OPEN: no B3-specific axion/GS kinetic, shift and counterterm action has been supplied and reduced",
        "requirements": ["Specify parent form fields and their kinetic/shift/counterterm action, including spin^c twist data",
                         "Reduce on independent cohomology classes, respecting form degrees and field dualities",
                         "Compute full 4D anomaly polynomial and cancel anomalies, including cubic anomalies",
                         "Derive K in the same neutral connection basis as every Higgs charge covector",
                         "Show a positive kinetic mass operator with ker equal to the desired electromagnetic direction",
                         "Normalize gauge and matter mode kinetic forms and match measured coupling separately"],
    })


def main():
    here = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=here.parents[2])
    parser.add_argument("--output", type=Path, default=here / "photon_checks.json")
    args = parser.parse_args()
    result = audit(args.repo.resolve())
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(f"Matter-02 photon audit: {result['exact_gates_passed']} exact gates PASS; physical photon OPEN")
    print(args.output)


if __name__ == "__main__":
    main()
