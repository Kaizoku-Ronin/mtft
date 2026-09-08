#!/usr/bin/env python3
"""Independent reconstruction of the MTFT active module.

There are two different sets of three operators: B_i represent J-antilinear
Hamiltonians z -> B_i conjugate(z), while C_ij = B_i conjugate(B_j) -
B_j conjugate(B_i) are the complex-linear Lie generators. Their kernels
are compared, not identified by definition.

The strongest independent route uses no channel matrices or Lie closure:
let W be the exact rational harmonic matrix, H its four independent columns
on the eight selected edges, Gg = W W^T, and G the Hodge metric. Every
localized V_i has image in im(Gg^-1 H), and its G-adjoint has image in
im(G^-1 H). Thus every A_i^- = (A_i + J A_i J)/2 has image in

  E = span(Gg^-1 H, G^-1 H, J Gg^-1 H, J G^-1 H), dim_R(E) <= 16.

For an exact compatible Hodge pair (G,J), A_i^- is G-self-adjoint and
J-antilinear. Consequently E^{perp_G} lies in every channel kernel and
is J-stable. On real dimension 26, the common kernel has complex dimension
at least 5. This is an exact conditional linear-algebra bound. The support
rank 4 is certified with exact rational arithmetic below. Saturation and
all period-dependent comparisons remain double-precision diagnostics.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import sympy as sp


def fro(A):
    return float(np.linalg.norm(A))


def realify(A):
    return np.block([[A.real, -A.imag], [A.imag, A.real]])


def exact_support_certificate(W, edge_indices):
    A = W[:, edge_indices]
    _, column_pivots = A.rref()
    H = A[:, column_pivots]
    _, row_pivots = H.T.rref()
    minor = H[list(row_pivots), :]
    coordinates = minor.inv() * A[list(row_pivots), :]
    reconstruction_exact = H * coordinates == A
    assert len(column_pivots) == 4 and minor.det() != 0
    assert reconstruction_exact
    return H, {
        "epistemic": "EXACT rational arithmetic on packaged harmonic matrix",
        "edge_indices": list(edge_indices),
        "independent_edge_indices": [edge_indices[i] for i in column_pivots],
        "independent_row_indices": list(row_pivots),
        "rank": len(column_pivots),
        "nonzero_minor_determinant": str(minor.det()),
        "reconstruction_exact": reconstruction_exact,
        "coordinates_4x8": [[str(x) for x in row] for row in coordinates.tolist()],
    }


def svd_kernel(matrices, antilinear=False, tol=1e-9):
    _, singular_values, vh = np.linalg.svd(np.vstack(matrices), full_matrices=False)
    nullity = int(np.count_nonzero(singular_values < tol * singular_values[0]))
    V = vh.conj().T[:, -nullity:]
    # Kernel of B K consists of conjugates of kernel vectors of B.
    if antilinear:
        V = V.conj()
    P = V @ V.conj().T
    return P, {
        "nullity_complex": nullity,
        "singular_values": singular_values.tolist(),
        "relative_threshold": tol,
        "largest_discarded": float(singular_values[-nullity]),
        "smallest_retained": float(singular_values[-nullity - 1]),
        "gap": float(singular_values[-nullity - 1] / singular_values[-nullity]),
        "rank_by_relative_threshold": {
            str(t): int(np.count_nonzero(singular_values > t * singular_values[0]))
            for t in (1e-7, 1e-9, 1e-11, 1e-13)
        },
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=Path(__file__).parent)
    parser.add_argument("--dps", type=int, default=50)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    import mtft
    from mtft import hecke as H, liealg as L
    from mtft.periods import data_path, physics as PH

    # Legend discovery occurred before authoring; these module-local methods
    # are absent from the registry in 0.26.2, so source definitions were read.
    harmonic_path = data_path("X0_143_m7_harmonic_basis.json")
    record = json.loads(harmonic_path.read_text())
    Wq = sp.Matrix([[sp.Rational(a, b) for a, b in row]
                    for row in record["basis_26x84"]])
    W = np.array(Wq, dtype=float)
    model = H.model()
    edges = [(model["tri_of"][model["erep"][k]],
              model["tri_of"][model["sS"][model["erep"][k]]])
             for k in range(model["E"])]
    triangles = (1, 11, 12)
    support = sorted(k for k, (a, b) in enumerate(edges)
                     if a in triangles or b in triangles)
    Hq, certificate = exact_support_certificate(Wq, support)
    certificate["per_triangle"] = {
        str(t): {
            "edges": [k for k, (a, b) in enumerate(edges) if a == t or b == t],
            "rank": Wq[:, [k for k, (a, b) in enumerate(edges)
                             if a == t or b == t]].rank(),
        } for t in triangles
    }

    R, Ri, B = L.x0143_fixed_channels(dps=args.dps, seed=143)
    C = [B[i] @ B[j].conj() - B[j] @ B[i].conj()
         for i in range(3) for j in range(i + 1, 3)]
    P, svd_report = svd_kernel(C)
    PB, B_report = svd_kernel(B, antilinear=True)
    energy = sum(A.conj().T @ A for A in C)
    eigenvalues, eigenvectors = np.linalg.eigh((energy + energy.conj().T) / 2)
    gram_nullity = int(np.count_nonzero(np.abs(eigenvalues) < 1e-10 * eigenvalues[-1]))
    assert gram_nullity == 5
    Vgram = eigenvectors[:, :gram_nullity]
    Pgram = Vgram @ Vgram.conj().T

    # Independent localization route, using a 26 x 16 real matrix only.
    G = PH.hodge_metric_hecke(args.dps)
    J = PH.hodge_structure_hecke(args.dps)
    Ja = R @ J @ Ri
    Hfloat = np.array(Hq, dtype=float)
    Y = R @ np.linalg.solve(W @ W.T, Hfloat)
    Z = R @ np.linalg.solve(G, Hfloat)
    K = np.column_stack([Y, Z, Ja @ Y, Ja @ Z])
    U, bound_singular_values, _ = np.linalg.svd(K, full_matrices=False)
    assert bound_singular_values[-1] > 1e-3
    P_localization = np.eye(26) - U @ U.T
    P_real = realify(P)
    localization = {
        "epistemic": "DIAGNOSTIC; exact rank upper bound, numerical saturation",
        "rank_growth_real": [int(np.linalg.matrix_rank(M, tol=1e-10))
                             for M in (Y, np.column_stack([Y, Z]), K)],
        "span_singular_values": bound_singular_values.tolist(),
        "span_condition_number": float(bound_singular_values[0] / bound_singular_values[-1]),
        "projector_distance_to_Cij_svd_real": fro(P_localization - P_real),
        "fixed_dimension_real": 26 - K.shape[1],
        "fixed_dimension_complex": (26 - K.shape[1]) // 2,
        "projector_J_commutator": fro(P_localization @ Ja - Ja @ P_localization),
    }

    # Comparison with the existing Lie closure is downstream, not used to
    # obtain the SVD, Gram, or localization projectors.
    closure = L.close_lie(C)
    rep = L.rep_summary(closure["basis"])
    STARa = R @ L.x0143_symmetry_ops()["STAR"] @ Ri
    Ustar = STARa[:13, :13] + 1j * STARa[13:, :13]
    Ustar_i = np.linalg.inv(Ustar)
    star_report = {
        "antilinear_structure_residual": fro(STARa @ Ja + Ja @ STARa) / fro(STARa),
        "unitarity_residual": fro(Ustar @ Ustar.conj().T - np.eye(13)),
        "involution_residual": fro(Ustar @ Ustar.conj() - np.eye(13)),
        "fixed_projector_preservation_residual": fro(Ustar @ P.conj() @ Ustar.conj().T - P),
        "localized_real_projector_commutator": fro(STARa @ P_localization - P_localization @ STARa),
        "Cij_identity_action_relative_residuals": [
            fro(Ustar @ A.conj() @ Ustar_i - A) / fro(A) for A in C],
        "interpretation": "Antiunitary real structure preserving both subspaces; not a spacetime signature.",
    }

    source_root = Path(L.__file__).resolve().parents[2]
    source_files = [Path(L.__file__), Path(PH.__file__), harmonic_path]
    results = {
        "epistemic": "DIAGNOSTIC except the exact rational support certificate and stated conditional rank bound",
        "version": mtft.__version__,
        "period_input_dps": args.dps,
        "operator_arithmetic": "numpy complex128 / float64, not interval-certified",
        "frame_seed": 143,
        "exact_support": certificate,
        "Cij_svd": svd_report,
        "B_antilinear_svd": B_report,
        "Cij_energy_gram": {
            "nullity": gram_nullity,
            "eigenvalues": eigenvalues.tolist(),
            "projector_distance_to_Cij_svd_complex": fro(Pgram - P),
            "note": "Forming C* C squares conditioning; small signed eigenvalues are floating-point roundoff.",
        },
        "localization_reconstruction": localization,
        "cross_checks": {
            "B_C_projector_distance_complex": fro(PB - P),
            "Cij_kernel_relative_residuals": [fro(A @ P) / fro(A) for A in C],
            "B_antilinear_kernel_relative_residuals": [fro(A @ P.conj()) / fro(A) for A in B],
            "closure_dimension": len(closure["basis"]),
            "closure_projector_distance_complex": fro(rep["P_fixed"] - P),
            "Cij_antihermitian_relative_residuals": [fro(A + A.conj().T) / fro(A) for A in C],
            "B_symmetric_relative_residuals": [fro(A - A.T) / fro(A) for A in B],
        },
        "STAR": star_report,
        "source_sha256": {
            str(p.relative_to(source_root)) if p.is_relative_to(source_root) else str(p):
            hashlib.sha256(p.read_bytes()).hexdigest() for p in source_files
        },
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    assert svd_report["nullity_complex"] == B_report["nullity_complex"] == 5
    assert localization["projector_distance_to_Cij_svd_real"] < 1e-9
    assert star_report["fixed_projector_preservation_residual"] < 1e-9
    (args.output / "kernel_check_results.json").write_text(json.dumps(results, indent=2) + "\n")
    np.savez_compressed(args.output / "kernel_check_projectors.npz",
                        Cij_fixed_complex=P, localization_fixed_real=P_localization,
                        localization_active_real=np.eye(26) - P_localization,
                        localization_span=K, R=R, Ri=Ri)
    print(json.dumps({"support_rank": certificate["rank"],
                      "fixed_dimension_complex": svd_report["nullity_complex"],
                      "localization_rank_growth": localization["rank_growth_real"],
                      "localization_projector_distance": localization["projector_distance_to_Cij_svd_real"],
                      "STAR_preservation": star_report["fixed_projector_preservation_residual"]}, indent=2))


if __name__ == "__main__":
    main()
