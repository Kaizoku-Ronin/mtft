#!/usr/bin/env python3
"""Independent W143 reduction from its positive eigenspace, without QWP SVD.

Let W = 2FF* - I be a Hermitian involution, and P an orthogonal projector.
Diagonalize F*PF with eigenvalues rho. For 0 < rho < 1, the normalized
projections x=PFv/sqrt(rho), y=(I-P)Fv/sqrt(1-rho) give W on (x,y) as
[[2rho-1, 2sqrt(rho(1-rho))], [2sqrt(rho(1-rho)), 1-2rho]].
The formula is exact algebra. The Hodge-frame realization below is a
double-precision DIAGNOSTIC using previously frozen geometric inputs.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import sympy as sp


def fro(a):
    return float(np.linalg.norm(a, "fro"))


def herm(a):
    return (a + a.conj().T) / 2


def complex_linear(a):
    n = a.shape[0] // 2
    return (a[:n, :n] + a[n:, n:] + 1j * (a[n:, :n] - a[:n, n:])) / 2


def realify(a):
    return np.block([[a.real, -a.imag], [a.imag, a.real]])


def hash_file(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(previous, source, output, tolerance):
    sys.path.insert(0, str(source / "src"))
    import mtft
    from mtft.periods.involutions import al_matrix
    from mtft.periods.core import data_path

    input_path = previous / "kernel_check_projectors.npz"
    data = np.load(input_path)
    r, ri = data["R"], data["Ri"]
    preal_raw = data["localization_active_real"]
    jreal = np.block([[np.zeros((13, 13)), -np.eye(13)], [np.eye(13), np.zeros((13, 13))]])
    if max(fro(preal_raw-preal_raw.T), fro(preal_raw @ preal_raw-preal_raw),
           fro(preal_raw @ jreal-jreal @ preal_raw)) >= tolerance:
        raise ValueError("Localization input failed the real orthogonal-projector / complex-linearity gate")
    p_raw = complex_linear(preal_raw)
    pvals, pvecs = np.linalg.eigh(herm(p_raw))
    active_basis = pvecs[:, pvals > 0.5]
    p = active_basis @ active_basis.conj().T
    q = np.eye(13) - p

    # Fresh exact arithmetic operator; no previous coupling matrices.
    w_exact = sp.Matrix(al_matrix(143))
    rplus_exact = (w_exact + sp.eye(26)) / 2
    exact_checks = {
        "involution": w_exact * w_exact == sp.eye(26),
        "positive_projector_idempotent": rplus_exact * rplus_exact == rplus_exact,
        "positive_projector_real_rank": int(rplus_exact.rank()),
        "operator_real_trace": int(sp.trace(w_exact)),
    }
    wreal = r @ np.asarray(w_exact, dtype=float) @ ri
    if max(fro(wreal-wreal.T), fro(wreal @ wreal-np.eye(26)),
           fro(wreal @ jreal-jreal @ wreal)) >= tolerance:
        raise ValueError("W143 failed the real selfadjoint-involution / complex-linearity gate")
    w_raw = complex_linear(wreal)
    w = herm(w_raw)
    wvals, wvecs = np.linalg.eigh(w)
    f = wvecs[:, wvals > 0]
    rplus = f @ f.conj().T
    overlap = herm(f.conj().T @ p @ f)
    rho, v = np.linalg.eigh(overlap)
    # Ordered by increasing rho, not by coupling singular value.
    plus_vectors = f @ v
    if not (np.all(rho > tolerance) and np.all(rho < 1-tolerance)):
        raise ValueError("Endpoint or ambiguous principal angle: generic coupled-block formula not applicable")
    x = p @ plus_vectors / np.sqrt(rho)[None, :]
    y = q @ plus_vectors / np.sqrt(1-rho)[None, :]
    a = 2*rho - 1
    b = 2*np.sqrt(rho*(1-rho))
    p_decoupled = herm(p - x @ x.conj().T)
    q_decoupled = herm(q - y @ y.conj().T)
    ep, up = np.linalg.eigh(p_decoupled)
    eq, uq = np.linalg.eigh(q_decoupled)
    xp = up[:, ep > 0.5]
    yq = uq[:, eq > 0.5]
    coupled_basis = np.column_stack([vec for j in range(len(rho)) for vec in (x[:, j], y[:, j])])
    u = np.column_stack((coupled_basis, xp, yq))
    expected = -np.eye(13, dtype=complex)
    for j in range(len(rho)):
        expected[2*j:2*j+2, 2*j:2*j+2] = [[a[j], b[j]], [b[j], -a[j]]]
    actual = u.conj().T @ w_raw @ u

    rows = []
    for j in range(len(rho)):
        theta = np.arccos(np.sqrt(rho[j]))
        e = np.column_stack((x[:, j], y[:, j]))
        block = np.array([[a[j], b[j]], [b[j], -a[j]]])
        rows.append({
            "index_in_increasing_rho": j,
            "rho_positive_eigenvector_active_overlap": float(rho[j]),
            "positive_eigenvector_fixed_overlap": float(1-rho[j]),
            "principal_angle_radians": float(theta),
            "principal_angle_degrees": float(np.degrees(theta)),
            "a": float(a[j]), "b": float(b[j]),
            "canonical_block": block.tolist(),
            "block_eigenvalues": np.linalg.eigvalsh(block).tolist(),
            "a_squared_plus_b_squared_minus_one": float(a[j]**2+b[j]**2-1),
            "eigenvector_projection_reconstruction": float(np.linalg.norm(
                plus_vectors[:, j] - np.sqrt(rho[j])*x[:, j] - np.sqrt(1-rho[j])*y[:, j])),
            "invariant_block_residual": fro(w_raw @ e - e @ block),
            "compression_active_eigenvalue": float(a[j]),
            "compression_fixed_eigenvalue": float(-a[j]),
        })

    errors = {
        "frame_inverse": fro(r @ ri-np.eye(26)),
        "localization_projector_complex_linearity": fro(preal_raw @ jreal-jreal @ preal_raw),
        "localization_complex_extraction_realification": fro(realify(p_raw)-preal_raw),
        "localization_hermitian_projection_change": fro(p-p_raw),
        "W_complex_linearity": fro(wreal @ jreal-jreal @ wreal),
        "W_complex_extraction_realification": fro(realify(w_raw)-wreal),
        "W_hermiticity": fro(w_raw-w_raw.conj().T),
        "W_involution": fro(w_raw @ w_raw-np.eye(13)),
        "positive_projector_eigenspace_reconstruction": fro(rplus-(w_raw+np.eye(13))/2),
        "x_orthonormality": fro(x.conj().T @ x-np.eye(len(rho))),
        "y_orthonormality": fro(y.conj().T @ y-np.eye(len(rho))),
        "cross_orthogonality": fro(x.conj().T @ y),
        "canonical_basis_unitarity": fro(u.conj().T @ u-np.eye(13)),
        "canonical_matrix_reconstruction": fro(actual-expected),
        "full_operator_reconstruction": fro(w_raw-u @ expected @ u.conj().T),
        "active_decoupled_eigenvalue_minus_one": fro(w_raw @ xp+xp),
        "fixed_decoupled_eigenvalue_minus_one": fro(w_raw @ yq+yq),
        "coupling_reconstruction_without_SVD": fro(q @ w_raw @ p-y @ np.diag(b) @ x.conj().T),
        "positive_eigenvectors_reconstruction": fro(plus_vectors-x*np.sqrt(rho)[None,:]-y*np.sqrt(1-rho)[None,:]),
    }
    checks = {
        "exact_W_is_involution": exact_checks["involution"],
        "exact_Rplus_is_projector": exact_checks["positive_projector_idempotent"],
        "exact_positive_rank_four_real": exact_checks["positive_projector_real_rank"] == 4,
        "positive_eigenspace_two_complex": f.shape[1] == 2,
        "active_rank_eight_complex": active_basis.shape[1] == 8,
        "two_strictly_mixed_principal_angles": len(rho) == 2 and np.all(rho > tolerance) and np.all(rho < 1-tolerance),
        "active_decoupled_six_complex": xp.shape[1] == 6,
        "fixed_decoupled_three_complex": yq.shape[1] == 3,
        "all_geometry_residuals_within_tolerance": max(errors.values()) < tolerance,
        "all_block_residuals_within_tolerance": max(row["invariant_block_residual"] for row in rows) < tolerance,
    }
    checks = {key: bool(value) for key, value in checks.items()}
    inputs = {str(input_path): hash_file(input_path),
              str(data_path("X0_143_atkin_lehner_v022.json")): hash_file(data_path("X0_143_atkin_lehner_v022.json")),
              str(Path(__file__)): hash_file(Path(__file__))}
    protocol = output / "PROTOCOL.md"
    if protocol.exists():
        inputs[str(protocol)] = hash_file(protocol)
    result = {
        "epistemic": "EXACT: integer W143 involution and positive-projector rank; conditional canonical algebra. DIAGNOSTIC: every Hodge-frame coefficient, principal angle, and residual.",
        "mtft_version": mtft.__version__,
        "route": "Fresh W143 arithmetic matrix, prior localization projector, positive-eigenspace Gram eigenproblem F*PF. No QWP SVD or primary reduction helper is used.",
        "shared_inputs": "Frozen Hodge frame R,Ri and localization projector from the preceding study; this is an independent computation route, not independent geometric source data.",
        "tolerance": tolerance,
        "input_sha256": inputs,
        "exact_checks": exact_checks,
        "raw_localization_projector_spectrum": pvals.tolist(),
        "raw_hermitian_W_spectrum": wvals.tolist(),
        "canonical_modes": rows,
        "coupling_singular_values_predicted_by_canonical_algebra_descending": sorted(b.tolist(), reverse=True),
        "decoupled_complex_dimensions": {"active": xp.shape[1], "fixed": yq.shape[1]},
        "decoupled_eigenvalue": -1,
        "total_complex_dimension": u.shape[1],
        "residuals": errors,
        "checks": checks,
        "all_checks_passed": all(checks.values()),
        "scope": "The dimensionless W143 spectral parameter is not physical energy, mass, or time. This canonical involution decomposition does not establish a new Clifford or gauge representation.",
    }
    output.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output / "canonical_check_matrices.npz", P=p, Q=q, W_raw=w_raw,
                        Rplus=rplus, positive_basis=f, positive_vectors=plus_vectors,
                        rho=rho, a=a, b=b, x=x, y=y, active_decoupled=xp,
                        fixed_decoupled=yq, canonical_basis=u,
                        canonical_matrix_expected=expected, canonical_matrix_actual=actual)
    (output / "canonical_check_results.json").write_text(json.dumps(result, indent=2)+"\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--previous", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=Path(__file__).parent)
    parser.add_argument("--tolerance", type=float, default=1e-9)
    args = parser.parse_args()
    result = run(args.previous, args.source, args.output, args.tolerance)
    print(json.dumps({"canonical_modes": result["canonical_modes"],
                      "maximum_geometry_residual": max(result["residuals"].values()),
                      "checks": result["checks"]}, indent=2))
    if not result["all_checks_passed"]:
        raise SystemExit(1)
