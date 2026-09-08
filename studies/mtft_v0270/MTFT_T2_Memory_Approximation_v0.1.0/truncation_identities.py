#!/usr/bin/env python3
"""Exact synthetic controls for structure-preserving auxiliary truncation.

No MTFT measurements, prior study matrices, or numerical selection results
are read. The rational examples check algebra; the accompanying note proves
the arbitrary-dimensional norm bounds.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import sympy as sp


def run() -> dict:
    z = sp.symbols("z")
    x, t = sp.symbols("x t", real=True)
    eta = sp.symbols("eta", positive=True, real=True)
    checks = {}

    def zero(name, expr):
        entries = list(expr) if isinstance(expr, sp.MatrixBase) else [expr]
        checks[name] = all(sp.simplify(sp.expand(q)) == 0 for q in entries)
        if not checks[name]:
            raise AssertionError(name)

    def truth(name, expr):
        checks[name] = bool(expr)
        if not checks[name]:
            raise AssertionError(name)

    rat = sp.Rational
    A = sp.Matrix([[1, rat(1, 3)], [rat(1, 3), rat(-1, 2)]])
    B = sp.Matrix([[rat(1, 2), rat(1, 4)], [rat(1, 5), rat(2, 3)]])
    D = sp.Matrix([[2, rat(3, 7)], [rat(3, 7), -1]])
    Qr = sp.Matrix([rat(3, 5), rat(4, 5)])
    Qs = sp.Matrix([rat(-4, 5), rat(3, 5)])
    Q = Qr.row_join(Qs)
    V = sp.diag(sp.eye(2), Q)
    H = A.row_join(B).col_join(B.T.row_join(D))
    Br, Bs = B * Qr, B * Qs
    Dr, Ds, M = Qr.T * D * Qr, Qs.T * D * Qs, Qr.T * D * Qs
    Hr = A.row_join(Br).col_join(Br.T.row_join(Dr))
    C = Bs.col_join(M)
    Hrot = V.T * H * V
    H0 = sp.diag(Hr, Ds)
    E = Hrot - H0
    zero("orthonormal_retained_basis", Qr.T * Qr - sp.eye(1))
    zero("orthonormal_discarded_basis", Qs.T * Qs - sp.eye(1))
    zero("retained_discarded_orthogonality", Qr.T * Qs)
    zero("complete_complement_basis", Q * Q.T - sp.eye(2))
    zero("original_generator_Hermitian", H - H.T)
    zero("retained_generator_Hermitian", Hr - Hr.T)
    zero("embedded_generator_Hermitian", H0 - H0.T)
    zero("rotated_full_block_decomposition", Hrot - Hr.row_join(C).col_join(C.T.row_join(Ds)))
    E_expected = sp.zeros(3, 3).row_join(C).col_join(C.T.row_join(sp.zeros(1, 1)))
    zero("error_contains_only_C_off_diagonal", E - E_expected)
    zero("squared_error_blocks", E * E - sp.diag(C * C.T, C.T * C))
    c2 = (C.T * C)[0, 0]
    zero("error_characteristic_polynomial", E.charpoly(z).as_expr() - z**2 * (z**2 - c2))
    truth("direct_discarded_coupling_nonzero", Bs != sp.zeros(2, 1))
    truth("complement_mixing_nonzero", M[0, 0] != 0)
    truth("exact_norm_squared_positive", c2 > 0)

    R = (z * sp.eye(4) - Hrot).inv()
    R0 = (z * sp.eye(4) - H0).inv()
    Rr = (z * sp.eye(3) - Hr).inv()
    zero("embedded_projected_resolvent_equals_reduced", R0[:2, :2] - Rr[:2, :2])
    # Multiplication form avoids a needlessly expanded four-by-four rational
    # inverse product while checking precisely the resolvent difference identity.
    zero("resolvent_difference_identity", (z * sp.eye(4) - Hrot) * (R - R0) - E * R0)
    dr = Dr[0, 0]
    sigma = Br * Br.T / (z - dr)
    zero("reduced_Schur_identity", (z * sp.eye(2) - A - sigma) * Rr[:2, :2] - sp.eye(2))
    zero("self_energy_residue", sigma.applyfunc(lambda q: sp.cancel((z-dr)*q)) - Br * Br.T)
    residue = Br * Br.T
    truth("self_energy_residue_PSD_rank_one", residue.det() == 0 and residue.trace() > 0)
    sig_up = sigma.subs(z, x + sp.I * eta)
    neg_im_sig = -(sig_up - sig_up.conjugate().T) / (2 * sp.I)
    zero("upper_half_plane_self_energy_sign", neg_im_sig - eta * Br * Br.T / ((x-dr)**2 + eta**2))
    z0, eta0 = rat(1, 2) + sp.I * rat(3, 5), rat(3, 5)
    Rrat = Rr.subs(z, z0)
    zero("upper_half_plane_full_resolvent_sign", -(Rrat-Rrat.conjugate().T)/(2*sp.I) - eta0 * Rrat.conjugate().T * Rrat)
    zero("active_second_moment_error", (Hrot**2)[:2, :2] - (Hr**2)[:2, :2] - Bs * Bs.T)
    zero("active_first_moment_error_zero", Hrot[:2, :2] - Hr[:2, :2])
    zero("memory_initial_error", B * B.T - Br * Br.T - Bs * Bs.T)

    # Exact zero-direct / nonzero-indirect counterexample. The auxiliary
    # truncation keeps coordinate 2 and discards coordinate 3.
    K = sp.Matrix([[0, 1, 0], [1, 0, 1], [0, 1, 0]])
    Kr = K[:2, :2]
    truth("indirect_example_discarded_B_exactly_zero", K[0, 2] == 0)
    truth("indirect_example_complement_mixing_nonzero", K[1, 2] != 0)
    for power in range(4):
        zero(f"indirect_example_active_moment_{power}_matches", (K**power)[0, 0] - (Kr**power)[0, 0])
    zero("indirect_example_fourth_moment_difference", (K**4)[0, 0] - (Kr**4)[0, 0] - 1)
    c = sp.cos(sp.sqrt(2) * t)
    s = sp.sin(sp.sqrt(2) * t)
    U = sp.eye(3) - sp.I * s * K / sp.sqrt(2) + (c-1) * K**2 / 2
    Ur = sp.cos(t) * sp.eye(2) - sp.I * sp.sin(t) * Kr
    zero("indirect_example_full_unitarity", U.conjugate().T * U - sp.eye(3))
    zero("indirect_example_reduced_unitarity", Ur.conjugate().T * Ur - sp.eye(2))
    zero("indirect_example_full_evolution_equation", sp.I * U.diff(t) - K * U)
    zero("indirect_example_reduced_evolution_equation", sp.I * Ur.diff(t) - Kr * Ur)
    zero("indirect_example_full_norm_balance", (U[:, 0:1].conjugate().T * U[:, 0:1])[0, 0] - 1)
    zero("indirect_example_reduced_norm_balance", (Ur[:, 0:1].conjugate().T * Ur[:, 0:1])[0, 0] - 1)
    zero("indirect_example_active_error_t4", sp.diff(U[0, 0] - Ur[0, 0], t, 4).subs(t, 0)/24 - rat(1, 24))
    zero("indirect_example_discarded_amplitude_t2", sp.diff(U[2, 0], t, 2).subs(t, 0)/2 + rat(1, 2))
    zero("indirect_example_discarded_norm_t4", sp.diff(sp.conjugate(U[2, 0])*U[2, 0], t, 4).subs(t, 0)/24 - rat(1, 4))
    Gk = (z * sp.eye(3) - K).inv()[0, 0]
    Gkr = (z * sp.eye(2) - Kr).inv()[0, 0]
    zero("indirect_example_nonzero_resolvent_difference", Gk - Gkr - 1/(z*(z**2-2)*(z**2-1)))
    zero("indirect_example_resolvent_moment_tied_to_t4", sp.limit(z**5*(Gk-Gkr), z, sp.oo)-1)
    # Complementary memory is cos(t), whereas the retained complement is
    # scalar zero and has constant memory 1.
    zero("indirect_example_memory_initial_agreement", sp.cos(t).subs(t, 0)-1)
    zero("indirect_example_memory_second_derivative_difference", sp.diff(sp.cos(t)-1, t, 2).subs(t, 0)+1)

    return {
        "status": "EXACT_SYNTHETIC_CONTROLS",
        "scope": "No T2 or other MTFT measured matrices read; generic proofs in TRUNCATION_IDENTITIES.md.",
        "checks": checks,
        "checks_passed": sum(checks.values()),
        "checks_total": len(checks),
        "all_checks_passed": all(checks.values()),
        "sympy_version": sp.__version__,
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "first_control": {
            "active_dimension": 2, "complement_dimension": 2, "retained_dimension": 1,
            "A": str(A), "B": str(B), "D": str(D), "Qr": str(Qr), "Qs": str(Qs),
            "C": str(C), "error_operator_norm_squared": str(c2),
            "norm_equality_route": "E^2=diag(CC*,C*C), with exact characteristic polynomial control"
        },
        "indirect_control": {
            "H": str(K), "Hr": str(Kr), "discarded_direct_coupling": 0,
            "complement_mixing": 1, "active_propagator_error_leading_term": "t^4/24",
            "discarded_norm_leading_term": "t^4/4",
            "active_resolvent_difference": "1/(z*(z^2-2)*(z^2-1))"
        },
        "bounds_proved_in_note": {
            "norm_convention": "absolute spectral operator norm",
            "domain": "H Hermitian; Qr orthonormal; Im(z)=eta>0; t real",
            "resolvent": "||G(z)-G_r(z)||_2 <= ||C||_2/eta^2",
            "propagator": "||U_PP(t)-U_r,PP(t)||_2 <= min(2,abs(t)*||C||_2)",
            "relative_budget_warning": "These bounds alone do not certify a 1% relative Frobenius error."
        },
        "interpretation": "Direct-sum amplitude coordinates and dimensionless Hermitian evolution; no density-matrix partial trace or physical Hamiltonian identification."
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=Path(__file__).with_name("truncation_identity_results.json"))
    args = parser.parse_args()
    results = run()
    args.output.write_text(json.dumps(results, indent=2) + "\n")
    print(json.dumps({"passed": results["checks_passed"], "total": results["checks_total"], "output": str(args.output)}))


if __name__ == "__main__":
    main()
