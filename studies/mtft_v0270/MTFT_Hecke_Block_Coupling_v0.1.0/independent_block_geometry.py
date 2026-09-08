#!/usr/bin/env python3
"""Independent Hecke-block geometry from exact bases and localization.

This route never uses Chinese-remainder projectors or the primary script's
helpers. The block bases are exact package input; all Hodge geometry here
is float64 DIAGNOSTIC, using the preceding study's shared period/frame data.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np


BLOCKS = ("ell", "old", "q4", "q6")
PRIMES = (2, 3, 5, 7, 17, 19)
ANGLE_TOL = 1e-8
ANGLE_AMBIGUOUS_LIMIT = 1e-6
SINGULAR_NULL_LIMIT = 1e-9
SINGULAR_NONZERO_LIMIT = 1e-7


def norm(a):
    return float(np.linalg.norm(a, "fro"))


def run(previous: Path, source: Path):
    sys.path.insert(0, str(source / "src"))
    import mtft
    from mtft import hecke as H
    from mtft.periods.involutions import al_matrix

    projector_path = previous / "kernel_check_projectors.npz"
    data = np.load(projector_path)
    r, ri = data["R"], data["Ri"]
    ploc = data["localization_active_real"]
    peig, pu = np.linalg.eigh((ploc + ploc.T) / 2)
    active_basis = pu[:, peig > 0.5]
    fixed_basis = pu[:, peig <= 0.5]
    p = active_basis @ active_basis.T
    n = len(p)
    q = fixed_basis @ fixed_basis.T
    # R is the frozen Hodge-orthonormal, complex-adapted real frame.
    # This J uses its documented realification convention, independent
    # of the primary result matrices.
    j = np.block([[np.zeros((13, 13)), -np.eye(13)],
                  [np.eye(13), np.zeros((13, 13))]])
    exact_bases = H.blocks()
    qs, blocks = {}, {}
    for name in BLOCKS:
        transported = r @ np.asarray(exact_bases[name], dtype=float).T
        u, triangular = np.linalg.qr(transported, mode="reduced")
        qs[name] = u
        s = np.linalg.svd(transported, compute_uv=False)
        compressed = u.T @ p @ u
        ev = np.linalg.eigvalsh((compressed + compressed.T) / 2)
        overlap = float(np.trace(compressed))
        blocks[name] = {
            "real_dimension": u.shape[1],
            "complex_dimension": u.shape[1] // 2,
            "transported_basis_singular_values": s.tolist(),
            "transported_basis_condition_number": float(s[0] / s[-1]),
            "basis_orthonormality_residual": norm(u.T @ u - np.eye(u.shape[1])),
            "active_trace_overlap_real": overlap,
            "active_trace_overlap_complex": overlap / 2,
            "fixed_trace_overlap_real": u.shape[1] - overlap,
            "active_fraction_of_block": overlap / u.shape[1],
            "squared_principal_cosines_real": ev.tolist(),
            "squared_principal_cosine_pairing_max_difference":
                float(np.max(np.abs(ev[::2] - ev[1::2]))),
            "active_intersection_dimension_real_at_frozen_tolerance":
                int(np.count_nonzero(np.abs(ev - 1) <= ANGLE_TOL)),
            "fixed_intersection_dimension_real_at_frozen_tolerance":
                int(np.count_nonzero(np.abs(ev) <= ANGLE_TOL)),
            "ambiguous_endpoint_eigenvalue_count": int(np.count_nonzero(
                (np.minimum(np.abs(ev), np.abs(ev - 1)) > ANGLE_TOL)
                & (np.minimum(np.abs(ev), np.abs(ev - 1)) <= ANGLE_AMBIGUOUS_LIMIT))),
            "mixed_eigenvalue_count": int(np.count_nonzero(
                np.minimum(np.abs(ev), np.abs(ev - 1)) > ANGLE_AMBIGUOUS_LIMIT)),
            "block_projector_J_commutator": norm((u @ u.T) @ j - j @ (u @ u.T)),
            "compressed_projector_idempotence_residual": norm(compressed @ compressed - compressed),
        }

    all_u = np.column_stack([qs[name] for name in BLOCKS])
    pdiag = sum((u @ (u.T @ p @ u) @ u.T for u in qs.values()), np.zeros_like(p))
    cross_p_squared = [[norm(qs[b].T @ p @ qs[c]) ** 2 for c in BLOCKS] for b in BLOCKS]

    exact_ops = {f"T{k}": H.cuspidal_hecke(k) for k in PRIMES}
    exact_ops.update({f"W{k}": al_matrix(k) for k in (11, 13, 143)})
    exact_ops["STAR"] = H.star_involution()
    operators = {}
    for name, exact in exact_ops.items():
        a = r @ np.asarray(exact, dtype=float) @ ri
        coupling = fixed_basis.T @ a @ active_basis
        coupling_singular_values = np.linalg.svd(coupling, compute_uv=False)
        coupling_singular_normalized = coupling_singular_values / norm(a)
        # Independently form each component in the 10 x 16 rectangular
        # fixed-to-active coordinate matrix, rather than using the
        # primary script's full 26 x 26 K_b matrices.
        coupling_components = [fixed_basis.T @ u @ (u.T @ a @ u) @ u.T @ active_basis
                               for u in qs.values()]
        component_rows = np.vstack([k.ravel() / norm(a) for k in coupling_components])
        signed_gram = component_rows @ component_rows.T
        leakage_energy = norm(coupling) ** 2 / norm(a) ** 2
        comm = p @ a - a @ p
        matrix_pairs, formula_pairs = [], []
        offblock_a_squared = 0.0
        max_formula_difference = 0.0
        for b in BLOCKS:
            ub = qs[b]
            abb = ub.T @ a @ ub
            norm_row, formula_row = [], []
            for c in BLOCKS:
                uc = qs[c]
                if b != c:
                    offblock_a_squared += norm(ub.T @ a @ uc) ** 2
                block = ub.T @ comm @ uc
                pbc = ub.T @ p @ uc
                # Valid when A preserves each arithmetic block. The
                # discrepancy is retained, never assumed to vanish.
                formula = pbc @ (uc.T @ a @ uc) - abb @ pbc
                norm_row.append(norm(block) ** 2 / norm(a) ** 2)
                formula_row.append(norm(formula) ** 2 / norm(a) ** 2)
                max_formula_difference = max(max_formula_difference, norm(block - formula) / norm(a))
            matrix_pairs.append(norm_row)
            formula_pairs.append(formula_row)
        total = norm(comm) ** 2 / norm(a) ** 2
        operators[name] = {
            "coupling_singular_values": coupling_singular_values.tolist(),
            "singular_values_over_operator_norm": coupling_singular_normalized.tolist(),
            "coupling_rank_real_at_frozen_bands": int(np.sum(coupling_singular_normalized >= SINGULAR_NONZERO_LIMIT)),
            "coupling_null_singular_count_at_frozen_bands": int(np.sum(coupling_singular_normalized <= SINGULAR_NULL_LIMIT)),
            "coupling_rank_ambiguous": bool(np.any((coupling_singular_normalized > SINGULAR_NULL_LIMIT)
                                                  & (coupling_singular_normalized < SINGULAR_NONZERO_LIMIT))),
            "signed_coupling_Gram_over_operator_squared": signed_gram.tolist(),
            "signed_Gram_sum": float(np.sum(signed_gram)),
            "signed_Gram_diagonal_sum": float(np.trace(signed_gram)),
            "signed_Gram_offdiagonal_sum": float(np.sum(signed_gram) - np.trace(signed_gram)),
            "coupling_norm_squared_over_operator_squared": leakage_energy,
            "signed_Gram_sum_absolute_error": abs(float(np.sum(signed_gram)) - leakage_energy),
            "signed_Gram_minimum_eigenvalue": float(np.linalg.eigvalsh(signed_gram)[0]),
            "coupling_component_sum_relative_to_operator_error": norm(sum(coupling_components) - coupling) / norm(a),
            "signed_Gram_caveat": "Negative off-diagonal entries encode cancellation. The full Gram matrix is positive semidefinite mathematically; its sum is nonnegative mathematically but may round slightly below zero. Absolute errors are reported, including for STAR.",
            "commutator_squared_over_operator_squared": total,
            "block_pair_commutator_squared_over_operator_squared": matrix_pairs,
            "block_pair_formula_squared_over_operator_squared": formula_pairs,
            "maximum_pair_formula_difference_over_operator": max_formula_difference,
            "sum_pair_energy_minus_full_energy": float(np.sum(matrix_pairs) - total),
            "off_arithmetic_block_operator_fraction_squared": offblock_a_squared / norm(a) ** 2,
            "note": "Every ordered pair, including diagonal blocks, is retained. Off-diagonal arithmetic-block fractions concern A itself; commutator fractions concern [P,A].",
        }

    # Read the primary outputs only after our independent geometry,
    # operators, singular spectra, and component Grams are computed.
    primary_path = Path(__file__).with_name("block_coupling_results.json")
    primary = json.loads(primary_path.read_text())
    comparisons = []
    for reference in primary["runs"]:
        geometry_ref = reference["geometry"]["blocks"]
        per_block = {}
        for name in BLOCKS:
            own, ref = blocks[name], geometry_ref[name]
            per_block[name] = {
                "active_trace_overlap_real_absolute_difference": abs(own["active_trace_overlap_real"] - ref["active_overlap_real"]),
                "principal_cosine_spectrum_max_absolute_difference": float(np.max(np.abs(
                    np.asarray(own["squared_principal_cosines_real"]) - ref["principal_cosines_squared"]))),
                "active_intersection_counts_agree": own["active_intersection_dimension_real_at_frozen_tolerance"] == ref["active_intersection_real_diagnostic"],
                "fixed_intersection_counts_agree": own["fixed_intersection_dimension_real_at_frozen_tolerance"] == ref["fixed_intersection_real_diagnostic"],
            }
        per_operator = {}
        for name, own in operators.items():
            ref = reference["operators"][name]
            per_operator[name] = {
                "commutator_pair_attributions_max_absolute_difference": float(np.max(np.abs(
                    np.asarray(own["block_pair_commutator_squared_over_operator_squared"])
                    - ref["commutator_block_squared_fractions"]))),
                "singular_values_normalized_max_absolute_difference": float(np.max(np.abs(
                    np.asarray(own["singular_values_over_operator_norm"]) - ref["singular_values_over_operator_norm"]))),
                "signed_Gram_max_absolute_difference": float(np.max(np.abs(
                    np.asarray(own["signed_coupling_Gram_over_operator_squared"]) - ref["signed_coupling_Gram_over_operator_squared"]))),
                "signed_Gram_sum_absolute_difference": abs(own["signed_Gram_sum"] - ref["coupling_Gram_sum"]),
                "coupling_rank_agrees": own["coupling_rank_real_at_frozen_bands"] == ref["coupling_rank_real"],
                "rank_ambiguity_agrees": own["coupling_rank_ambiguous"] == ref["rank_ambiguous"],
            }
        differences = [value for group in (per_block, per_operator)
                       for row in group.values() for key, value in row.items()
                       if key.endswith("difference")]
        agreements = [value for group in (per_block, per_operator)
                      for row in group.values() for key, value in row.items()
                      if key.endswith("agree") or key.endswith("agrees")]
        comparisons.append({"frame_seed": reference["frame_seed"],
                            "period_input_dps": reference["period_input_dps"],
                            "blocks": per_block, "operators": per_operator,
                            "maximum_absolute_difference": max(differences),
                            "all_discrete_counts_agree": all(agreements)})

    checks = {
        "active_dimension_is_16_real": active_basis.shape[1] == 16,
        "basis_dimensions_are_2_4_8_12": [qs[k].shape[1] for k in BLOCKS] == [2, 4, 8, 12],
        "combined_block_basis_orthonormal_within_1e_9": norm(all_u.T @ all_u - np.eye(n)) < 1e-9,
        "active_trace_sum_is_16_within_1e_9": abs(sum(b["active_trace_overlap_real"] for b in blocks.values()) - 16) < 1e-9,
        "fixed_trace_sum_is_10_within_1e_9": abs(sum(b["fixed_trace_overlap_real"] for b in blocks.values()) - 10) < 1e-9,
        "principal_cosines_within_unit_interval_at_1e_9": all(
            min(b["squared_principal_cosines_real"]) >= -1e-9
            and max(b["squared_principal_cosines_real"]) <= 1+1e-9 for b in blocks.values()),
        "squared_principal_cosines_paired_within_1e_9": all(
            b["squared_principal_cosine_pairing_max_difference"] < 1e-9 for b in blocks.values()),
        "all_pair_commutator_formulas_within_1e_9": all(
            a["maximum_pair_formula_difference_over_operator"] < 1e-9 for a in operators.values()),
        "all_pair_energy_sums_within_1e_9": all(
            abs(a["sum_pair_energy_minus_full_energy"]) < 1e-9 for a in operators.values()),
        "all_signed_Gram_sums_within_absolute_1e_9": all(
            a["signed_Gram_sum_absolute_error"] < 1e-9 for a in operators.values()),
        "all_coupling_component_reconstructions_within_1e_9": all(
            a["coupling_component_sum_relative_to_operator_error"] < 1e-9 for a in operators.values()),
        "all_signed_Grams_numerically_positive_semidefinite_at_1e_12": all(
            a["signed_Gram_minimum_eigenvalue"] >= -1e-12 for a in operators.values()),
        "all_primary_comparisons_within_absolute_1e_9": all(
            row["maximum_absolute_difference"] < 1e-9 for row in comparisons),
        "all_primary_discrete_counts_agree": all(row["all_discrete_counts_agree"] for row in comparisons),
    }
    result = {
        "epistemic": "DIAGNOSTIC numerical Hodge geometry; exact input bases do not make numerical intersections exact",
        "version": mtft.__version__,
        "route": "Exact H.blocks basis columns transported by frozen R, independent QR per block; prior localization projector independently diagonalized; no CRT projector construction",
        "shared_inputs": "mtft arithmetic source, frozen period data and the previous study's R, Ri, localization projector",
        "block_order": list(BLOCKS),
        "frozen_squared_principal_cosine_distance_tolerance": ANGLE_TOL,
        "frozen_ambiguous_endpoint_distance_upper_bound": ANGLE_AMBIGUOUS_LIMIT,
        "frozen_coupling_singular_bands_normalized_by_operator_Frobenius_norm": {
            "null_at_or_below": SINGULAR_NULL_LIMIT, "nonzero_at_or_above": SINGULAR_NONZERO_LIMIT,
            "between": "ambiguous"},
        "angle_caveat": "Tolerance applies to squared cosines near 0 or 1, not to angles in radians; counts are numerical diagnostics.",
        "input_sha256": {"kernel_check_projectors.npz": hashlib.sha256(projector_path.read_bytes()).hexdigest()},
        "source_sha256": {"src/mtft/hecke.py": hashlib.sha256(Path(H.__file__).read_bytes()).hexdigest()},
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "protocol_sha256": hashlib.sha256(Path(__file__).with_name("PROTOCOL.md").read_bytes()).hexdigest(),
        "projector_reconstruction_residual": norm(p - ploc),
        "active_projector_J_commutator": norm(p @ j - j @ p),
        "combined_block_basis_orthonormality_residual": norm(all_u.T @ all_u - np.eye(n)),
        "blocks": blocks,
        "block_pair_active_projector_norm_squared": cross_p_squared,
        "block_diagonal_truncation": {
            "formula": "P_diag = sum_b E_b P E_b, E_b = U_b U_b^T",
            "trace": float(np.trace(pdiag)),
            "idempotence_residual": norm(pdiag @ pdiag - pdiag),
            "discarded_off_block_norm": norm(p - pdiag),
            "discarded_norm_squared": norm(p - pdiag) ** 2,
            "trace_defect_trace_Pdiag_minus_Pdiag_squared": float(np.trace(pdiag - pdiag @ pdiag)),
            "discarded_energy_trace_identity_difference": norm(p - pdiag) ** 2 - float(np.trace(pdiag - pdiag @ pdiag)),
            "interpretation": "Block-diagonal truncation is generally a positive contraction rather than a projector; its failure of idempotence measures lost cross-block coherence, not failure of the original projector.",
        },
        "operators": operators,
        "primary_comparison": {"input_sha256": hashlib.sha256(primary_path.read_bytes()).hexdigest(),
                               "runs": comparisons,
                               "maximum_absolute_difference_all_runs": max(row["maximum_absolute_difference"] for row in comparisons)},
        "checks": checks,
        "all_checks_passed": all(checks.values()),
    }
    assert result["all_checks_passed"], checks
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--previous", type=Path, required=True)
    args = parser.parse_args()
    result = run(args.previous.resolve(), args.source.resolve())
    target = Path(__file__).with_name("independent_block_geometry_results.json")
    target.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"all_checks_passed": result["all_checks_passed"],
                      "blocks": result["blocks"],
                      "block_diagonal_truncation": result["block_diagonal_truncation"]}, indent=2))
