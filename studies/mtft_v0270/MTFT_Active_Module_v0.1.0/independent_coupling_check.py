#!/usr/bin/env python3
"""Independent replay from geometric localization, without coupling helpers.

Consumes the localization reconstruction from kernel_check.py, not the
main study's primitive-commutator kernel. The shared period/frame inputs
mean this is a distinct computation route, not independent source data.
"""
from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from pathlib import Path

import numpy as np
from mtft import hecke
from mtft.periods.involutions import al_matrix


def fnorm(a):
    return float(np.linalg.norm(a, "fro"))


def exact_commutator(a, b):
    n = len(a)
    return [[sum((Fraction(a[i][k])*Fraction(b[k][j])
                  -Fraction(b[i][k])*Fraction(a[k][j]) for k in range(n)),
                 Fraction(0)) for j in range(n)] for i in range(n)]


def run(root):
    path = root/"kernel_check_projectors.npz"
    source = np.load(path)
    # E comes only from the localization projector; its range is
    # orthonormalized independently by a symmetric eigendecomposition.
    ploc = source["localization_active_real"]
    ev, u = np.linalg.eigh((ploc+ploc.T)/2)
    e = u[:, ev > .5]
    p = e @ e.T
    q = np.eye(26)-p
    r, ri = source["R"], source["Ri"]
    exact = {f"T{k}": hecke.cuspidal_hecke(k) for k in (2,3,5,7,17,19)}
    exact.update({f"W{k}": al_matrix(k) for k in (11,13,143)})
    exact["STAR"] = hecke.star_involution()
    ops = {k: r @ np.asarray(v, float) @ ri for k,v in exact.items()}
    # The main result file is read only for numerical comparison, after
    # the independent projector and operators have been constructed.
    main = json.loads((root/"coupling_results.json").read_text())
    primary = next(x for x in main["runs"] if x["frame_seed"] == 143
                   and x["dps_period_inputs"] == 50)
    rows = {}
    for name, o in ops.items():
        af = fnorm(q @ o @ p)/fnorm(o)
        fa = fnorm(p @ o @ q)/fnorm(o)
        comm = fnorm(p @ o-o @ p)/fnorm(o)
        ref = primary["operators"][name]
        rows[name] = {
            "active_to_fixed_over_operator": af,
            "fixed_to_active_over_operator": fa,
            "commutator_over_operator": comm,
            "difference_from_main_active_to_fixed": abs(af-ref["active_to_fixed_over_operator"]),
            "difference_from_main_fixed_to_active": abs(fa-ref["fixed_to_active_over_operator"]),
            "difference_from_main_commutator": abs(comm-ref["commutator_over_operator"]),
        }
    # One augmentation, no iterative rank decisions and no Lie closure.
    scale = fnorm(ops["T2"])
    enlarged = np.column_stack((e, ops["T2"] @ e/scale))
    sv = np.linalg.svd(enlarged, compute_uv=False)
    exact_rows = {}
    exploratory = {}
    for other in ("T3", "W11"):
        cexact = exact_commutator(exact["T2"], exact[other])
        exact_rows[f"T2_{other}"] = {
            "all_entries_exactly_zero": all(x == 0 for row in cexact for x in row),
            "nonzero_entries": sum(x != 0 for row in cexact for x in row),
        }
        a,b = ops["T2"], ops[other]
        ac,bc = e.T @ a @ e, e.T @ b @ e
        c = ac @ bc-bc @ ac
        through_complement = e.T @ (b @ q @ a-a @ q @ b) @ e
        exploratory[f"T2_{other}"] = {
            "compressed_commutator_norm": fnorm(c),
            "compressed_commutator_over_product_norms": fnorm(c)/(fnorm(ac)*fnorm(bc)),
            "off_sector_identity_residual_relative_to_commutator": fnorm(c-through_complement)/fnorm(c),
        }
    differences = [v[k] for v in rows.values() for k in v if k.startswith("difference_")]
    assert e.shape == (26,16)
    assert sv[-1]/sv[0] > 1e-3
    assert max(differences) < 1e-10
    assert all(v["all_entries_exactly_zero"] for v in exact_rows.values())
    return {
        "epistemic": "DIAGNOSTIC except exact Fraction arithmetic commutator identities",
        "route": "independent localization projector and fresh eigendecomposition; no coupling_study helpers or main projector used",
        "shared_inputs": "frozen periods, frame matrices, package arithmetic operators",
        "localization_input_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        "active_projector_rank": e.shape[1],
        "projector_spectrum": ev.tolist(),
        "projector_reconstruction_error": fnorm(ploc-p),
        "operators": rows,
        "maximum_absolute_difference_from_main": max(differences),
        "T2_one_augmentation": {
            "construction": "column_stack(E, T2_adapted E / ||T2_adapted||_F)",
            "T2_scale": scale,
            "singular_values": sv.tolist(),
            "smallest_normalized_singular_value": float(sv[-1]/sv[0]),
            "rank_at_relative_1e_9": int(np.count_nonzero(sv/sv[0] > 1e-9)),
            "real_dimension": 26,
        },
        "exact_full_commutators": exact_rows,
        "exploratory_compression": {
            "status": "EXPLORATORY follow-up; not a preregistered primary endpoint",
            "identity": "For [A,B]=0, [E^T A E,E^T B E]=E^T(B Q A-A Q B)E, Q=I-EE^T",
            "interpretation": "Compression to a noninvariant space can make commuting operators fail to commute. This is a leakage identity, not a new exact representation of the Hecke algebra.",
            "pairs": exploratory,
        },
        "closure_scope": {
            "note": "Main study reports ambiguous package Lie-closure gates at frame seeds 2026 and 31415. This independent localization calculation does not use Lie closure and cannot repair or upgrade that closure claim.",
            "reported_runs": [{"seed": x["frame_seed"], "dps": x["dps_period_inputs"], "closure_status": x["checks"]["closure_status"]} for x in main["runs"]],
        },
        "checks_passed": True,
    }


if __name__ == "__main__":
    root = Path(__file__).resolve().parent
    result = run(root)
    (root/"independent_coupling_check_results.json").write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps(result, indent=2))
