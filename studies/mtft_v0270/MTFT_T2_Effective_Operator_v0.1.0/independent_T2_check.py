#!/usr/bin/env python3
"""Independent T2 response from the real 26-dimensional Hodge eigensystem.

No primary block-inverse or reduction helpers are used. The localization
projector and Hodge frame are shared upstream data, so this is a distinct
computation route rather than independent geometric or physical evidence.
All Hodge-dependent results are float64/complex128 DIAGNOSTIC. Only the
fresh arithmetic matrix's characteristic polynomial is exact here.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import scipy.linalg as sla
import sympy as sp


def fro(a):
    return float(np.linalg.norm(a, "fro"))


def linear_branch(a):
    """Complex-linear branch, including for complex-valued real-frame arrays.

    At complex z, an ambient real-frame resolvent is complex valued and is
    not the literal realification of the desired 13-dimensional response.
    This algebraic branch formula is applied without taking real/imag parts.
    """
    n = a.shape[-1] // 2
    return (a[..., :n, :n] + a[..., n:, n:]
            + 1j*(a[..., n:, :n] - a[..., :n, n:])) / 2


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def exact_spectral_data(h):
    x = sp.Symbol("x")
    cp = h.charpoly(x).as_poly()
    factorization = sp.factor_list(cp.as_expr(), x)
    factors = []
    all_roots = []
    max_imaginary = 0.0
    for expression, multiplicity in factorization[1]:
        poly = sp.Poly(expression, x)
        roots = poly.nroots(n=50, maxsteps=200)
        for value in roots:
            re, im = value.as_real_imag()
            max_imaginary = max(max_imaginary, abs(float(im)))
            all_roots.extend([float(re)] * multiplicity)
        factors.append({"polynomial": str(poly.as_expr()),
                        "coefficients_descending": [str(c) for c in poly.all_coeffs()],
                        "multiplicity": int(multiplicity),
                        "numerical_roots_50_digits": [str(v) for v in roots]})
    square = all(multiplicity % 2 == 0 for _, multiplicity in factorization[1])
    return {
        "epistemic": "EXACT rational characteristic polynomial and factorization; 50-digit factor roots are numerical diagnostics.",
        "characteristic_polynomial": str(cp.as_expr()),
        "coefficients_descending": [str(c) for c in cp.all_coeffs()],
        "leading_factor": str(factorization[0]),
        "factorization": factors,
        "all_factor_multiplicities_even": square,
        "square_root_monic_polynomial": str(sp.prod(expr**(mult//2) for expr, mult in factorization[1]).expand()) if square else None,
        "arithmetic_trace": str(sp.trace(h)),
        "arithmetic_determinant": str(h.det()),
        "maximum_numerical_root_imaginary_part": max_imaginary,
    }, np.sort(all_roots)


def run(source, previous, output):
    sys.path.insert(0, str(source / "src"))
    import mtft
    from mtft import hecke

    tolerance = 1e-9
    path = previous / "kernel_check_projectors.npz"
    data = np.load(path)
    r, ri = data["R"], data["Ri"]
    p_raw = data["localization_active_real"]
    h_exact = sp.Matrix(hecke.cuspidal_hecke(2))
    h_raw = r @ np.asarray(h_exact, dtype=float) @ ri
    eye = np.eye(26)
    j = np.block([[np.zeros((13, 13)), -np.eye(13)],
                  [np.eye(13), np.zeros((13, 13))]])
    p_scale, h_scale = max(1, fro(p_raw)), max(1, fro(h_raw))
    gate_residuals = {
        "frame_inverse": fro(r @ ri-eye)/fro(eye),
        "H_selfadjoint_relative": fro(h_raw-h_raw.T)/h_scale,
        "H_complex_linear_relative": fro(h_raw @ j-j @ h_raw)/h_scale,
        "P_selfadjoint_relative": fro(p_raw-p_raw.T)/p_scale,
        "P_idempotent_relative": fro(p_raw @ p_raw-p_raw)/p_scale,
        "P_complex_linear_relative": fro(p_raw @ j-j @ p_raw)/p_scale,
    }
    if max(gate_residuals.values()) > tolerance:
        raise ValueError("Raw real-frame Hermitian/projector gate failed: "+str(gate_residuals))

    h = (h_raw+h_raw.T)/2
    pvals, pu = np.linalg.eigh((p_raw+p_raw.T)/2)
    e, f = pu[:, pvals > 0.5], pu[:, pvals <= 0.5]
    if (e.shape[1], f.shape[1]) != (16, 10):
        raise ValueError("Expected active16/complement10 real dimensions")
    p = e @ e.T
    q = f @ f.T
    cleanup = {"H_symmetrization_change_relative": fro(h-h_raw)/h_scale,
               "P_eigenspace_projector_change_relative": fro(p-p_raw)/p_scale}
    if max(cleanup.values()) > tolerance:
        raise ValueError("Roundoff cleanup exceeded preregistered tolerance")

    eigenvalues, eigenvectors = np.linalg.eigh(h)
    exact, exact_roots = exact_spectral_data(h_exact)
    exact["real26_eigenspectrum_vs_factor_roots_maximum_difference"] = float(np.max(np.abs(eigenvalues-exact_roots)))
    a, b, d = e.T @ h @ e, e.T @ h @ f, f.T @ h @ f
    bsingular = np.linalg.svd(b, compute_uv=False)
    deigen, dvectors = np.linalg.eigh(d)
    normalized = bsingular / (fro(h)/np.sqrt(2))
    rank_real = int(np.count_nonzero(normalized >= 1e-7))
    ambiguous = bool(np.any((normalized > 1e-9) & (normalized < 1e-7)))
    pairing_errors = {
        "full_spectrum_pairing_maximum_difference": float(np.max(np.abs(eigenvalues[::2]-eigenvalues[1::2]))),
        "compression_spectrum_pairing_maximum_difference": float(np.max(np.abs(deigen[::2]-deigen[1::2]))),
        "coupling_singular_pairing_maximum_difference": float(np.max(np.abs(bsingular[::2]-bsingular[1::2]))),
    }
    grid = np.concatenate([np.linspace(-3, 3, 601)+1j*eta for eta in (0.05, 0.2, 0.7)]
                          + [np.asarray([-4, 4, 2j, 0.17+0.37j])])
    projected_eigenvectors = p @ eigenvectors
    responses = np.empty((len(grid), 13, 13), dtype=complex)
    row_errors, condition_numbers = [], []
    skipped = []
    for i, z in enumerate(grid):
        distances = np.abs(z-eigenvalues)
        if z.imag == 0 and np.min(distances) <= 1e-12*max(1, np.max(np.abs(eigenvalues))):
            responses[i] = np.nan + 1j*np.nan
            skipped.append(i)
            row_errors.append(float("nan"))
            condition_numbers.append(float("inf"))
            continue
        # Full real26 eigensystem and ambient real P; no sector Schur inverse.
        spectral = (projected_eigenvectors / (z-eigenvalues)[None, :]) @ projected_eigenvectors.T
        responses[i] = linear_branch(spectral)
        direct = p @ np.linalg.solve(z*eye-h, p)
        direct_branch = linear_branch(direct)
        row_errors.append(fro(direct_branch-responses[i])/max(fro(responses[i]), np.finfo(float).tiny))
        condition_numbers.append(float(np.max(distances)/np.min(distances)))

    h13, p13 = linear_branch(h), linear_branch(p)
    q13 = np.eye(13)-p13
    time_controls = []
    for t in (0, 0.17, np.pi/4, np.pi/2, np.pi, 2*np.pi):
        spectral26 = (eigenvectors*np.exp(-1j*t*eigenvalues)[None,:]) @ eigenvectors.T
        exponential26 = sla.expm(-1j*t*h)
        spectral13, exponential13 = linear_branch(spectral26), linear_branch(exponential26)
        pp = p13 @ exponential13 @ p13
        transfer = q13 @ exponential13 @ p13
        # D is diagonalized independently as a real symmetric 10x10 block.
        memory_spectral = (b @ dvectors*np.exp(-1j*t*deigen)[None,:]) @ (b @ dvectors).T
        memory_exponential = b @ sla.expm(-1j*t*d) @ b.T
        time_controls.append({
            "t": float(t),
            "real26_expm_vs_spectral_relative": fro(exponential26-spectral26)/fro(exponential26),
            "complex13_expm_vs_spectral_relative": fro(exponential13-spectral13)/fro(exponential13),
            "ambient_projected_expm_vs_spectral_relative": fro(pp-p13 @ spectral13 @ p13)/max(fro(pp), np.finfo(float).tiny),
            "complex_branch_unitarity_relative": fro(exponential13.conj().T @ exponential13-np.eye(13))/np.sqrt(13),
            "projected_norm_feedback_identity_relative": fro(p13-pp.conj().T @ pp-transfer.conj().T @ transfer)/max(1, fro(p13)),
            "memory_expm_vs_spectral_relative": fro(memory_exponential-memory_spectral)/max(fro(memory_exponential), np.finfo(float).tiny),
            "complementary_average_squared_norm": fro(transfer)**2/8,
            "complementary_worst_squared_norm": float(np.linalg.norm(transfer, 2)**2),
        })

    response_error = max(x for x in row_errors if np.isfinite(x))
    control_errors = [value for row in time_controls for key, value in row.items() if key.endswith("relative")]
    checks = {
        "raw_real_frame_gates": max(gate_residuals.values()) <= tolerance,
        "roundoff_cleanup_within_gate": max(cleanup.values()) <= tolerance,
        "active16_and_complement10_real": e.shape[1] == 16 and f.shape[1] == 10,
        "fresh_charpoly_has_degree26": len(exact["coefficients_descending"])-1 == 26,
        "exact_factor_roots_numerically_real": exact["maximum_numerical_root_imaginary_part"] < 1e-30,
        "eigenspectrum_agrees_with_exact_polynomial": exact["real26_eigenspectrum_vs_factor_roots_maximum_difference"] < tolerance,
        "all_real_spectra_pair_to_complex_multiplicities": max(pairing_errors.values()) < tolerance,
        "regular_grid_has1807_points": len(grid) == 1807,
        "all_admissible_grid_direct_and_spectral_responses_agree": response_error <= tolerance,
        "six_independent_expm_controls_agree": max(control_errors) <= tolerance,
    }
    inputs = {str(path): digest(path), str(Path(hecke.__file__)): digest(Path(hecke.__file__)),
              str(Path(__file__)): digest(Path(__file__))}
    protocol = output / "PROTOCOL.md"
    if protocol.exists():
        inputs[str(protocol)] = digest(protocol)
    result = {
        "epistemic": "EXACT arithmetic characteristic polynomial and factorization; all Hodge geometry and response measurements are DIAGNOSTIC.",
        "route": "Fresh T2 arithmetic matrix; real26 Hodge transport; localization projector; full real symmetric eigensystem; algebraic J-linear branch extraction at complex z. No primary code, outputs, direct projector, or Schur helper read.",
        "shared_upstream_data": "Prior frozen Hodge frame R,Ri and localization geometry; not independent physical data.",
        "version": mtft.__version__, "tolerance": tolerance,
        "input_sha256": inputs, "raw_gate_residuals": gate_residuals,
        "cleanup": cleanup, "projector_raw_eigenvalues": pvals.tolist(),
        "exact_arithmetic_spectrum": exact,
        "full_real26_eigenvalues": eigenvalues.tolist(),
        "full_complex13_eigenvalues_from_real_pair_means": ((eigenvalues[::2]+eigenvalues[1::2])/2).tolist(),
        "compression_real10_eigenvalues": deigen.tolist(),
        "compression_complex5_eigenvalues_from_real_pair_means": ((deigen[::2]+deigen[1::2])/2).tolist(),
        "coupling_real_singular_values": bsingular.tolist(),
        "coupling_complex_singular_values_from_real_pair_means": ((bsingular[::2]+bsingular[1::2])/2).tolist(),
        "coupling_real_singular_values_over_complex_operator_Frobenius_norm": normalized.tolist(),
        "coupling_rank_real_at_frozen_bands": rank_real,
        "coupling_rank_ambiguous": ambiguous,
        "pairing_residuals": pairing_errors,
        "regular_grid": {"point_count": len(grid), "skipped_indices": skipped,
                         "maximum_independent_direct_vs_spectral_relative_error": response_error,
                         "maximum_full_operator_condition_number": max(condition_numbers),
                         "complex_z_branch_note": "Responses are extracted from a complex-valued 26x26 projected resolvent via (aa+dd+i(cc-bb))/2. No literal realification assumption is made at complex z."},
        "time_controls": time_controls,
        "checks": {k: bool(v) for k,v in checks.items()},
        "all_checks_passed": all(checks.values()),
    }
    output.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output / "independent_T2_matrices.npz", grid=grid,
                        ambient_projected_responses=responses,
                        direct_spectral_relative_errors=np.asarray(row_errors),
                        full_condition_numbers=np.asarray(condition_numbers),
                        H_real=h, P_real=p, H_complex=h13, P_complex=p13,
                        active_basis_real=e, complement_basis_real=f,
                        A_real=a, B_real=b, D_real=d,
                        H_real_eigenvalues=eigenvalues, H_real_eigenvectors=eigenvectors,
                        D_real_eigenvalues=deigen, D_real_eigenvectors=dvectors)
    (output / "independent_T2_results.json").write_text(json.dumps(result, indent=2)+"\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--previous", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=Path(__file__).parent)
    args = parser.parse_args()
    result = run(args.source, args.previous, args.output)
    print(json.dumps({"characteristic_polynomial": result["exact_arithmetic_spectrum"]["characteristic_polynomial"],
                      "factorization": result["exact_arithmetic_spectrum"]["factorization"],
                      "coupling_singular_values": result["coupling_complex_singular_values_from_real_pair_means"],
                      "compression_eigenvalues": result["compression_complex5_eigenvalues_from_real_pair_means"],
                      "regular_grid": result["regular_grid"],
                      "checks": result["checks"]}, indent=2))
    if not result["all_checks_passed"]:
        raise SystemExit(1)
