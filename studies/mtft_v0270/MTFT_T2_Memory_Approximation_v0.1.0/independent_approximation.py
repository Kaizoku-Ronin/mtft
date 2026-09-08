#!/usr/bin/env python3
"""Independent real-frame Galerkin study; no primary matrices or scores are read."""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from pathlib import Path
import time

import numpy as np
from scipy.linalg import expm, null_space


TOL = 1e-9
FAMILIES = ("svd", "poles", "snapshots")
BUDGETS = (0.01, 0.05, 0.10)


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save_json(path, data):
    Path(path).write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")


def complex_branch(m):
    n = m.shape[0] // 2
    return (m[:n, :n] + m[n:, n:] + 1j * (m[n:, :n] - m[:n, n:])) / 2


def realify(m):
    return np.block([[m.real, -m.imag], [m.imag, m.real]])


def hermitian(m):
    return (m + m.conj().T) / 2


def norm(m):
    return float(np.linalg.norm(m))


def opnorm(m):
    return float(np.linalg.norm(m, 2)) if min(m.shape) else 0.0


def complex_json(m):
    return {"real": m.real.tolist(), "imag": m.imag.tolist()}


def descending_gaps(values):
    v = np.asarray(values)
    scale = max(float(np.max(np.abs(v))), np.finfo(float).tiny)
    return [float(abs(v[j] - v[j+1]) / scale) for j in range(len(v)-1)]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parent)
    args = parser.parse_args()
    root = args.root.resolve()
    started = time.perf_counter()
    input_path = root / "inputs" / "T2_independent_geometry.npz"
    protocol_path = root / "PROTOCOL.md"
    raw = np.load(input_path)
    hr, pr = raw["H_real"], raw["P_real"]
    h_raw, p_raw = complex_branch(hr), complex_branch(pr)
    h, p = hermitian(h_raw), hermitian(p_raw)
    pe, pv = np.linalg.eigh(p)
    u, f = pv[:, pe > .5], pv[:, pe <= .5]
    assert u.shape == (13, 8) and f.shape == (13, 5)
    a, b, d = u.conj().T @ h @ u, u.conj().T @ h @ f, f.conj().T @ h @ f
    s = np.column_stack((u, f))
    full_blocks = np.block([[a, b], [b.conj().T, d]])
    defects = {
        "H_real_hermitian": norm(hr - hr.T),
        "P_real_hermitian": norm(pr - pr.T),
        "H_complex_linearity": norm(hr - realify(h_raw)),
        "P_complex_linearity": norm(pr - realify(p_raw)),
        "H_hermitian_cleanup": norm(h - h_raw),
        "P_hermitian_cleanup": norm(p - p_raw),
        "projector_idempotency": norm(p @ p - p),
        "sector_orthonormality": norm(s.conj().T @ s - np.eye(13)),
        "projector_reconstruction": norm(u @ u.conj().T - p),
        "block_reconstruction": norm(s @ full_blocks @ s.conj().T - h),
    }
    assert max(defects.values()) < TOL, defects

    training_z = np.array([x + 1j*eta for eta in (.05, .2, .7)
                           for x in np.linspace(-3, 3, 61)])
    i8, i13 = np.eye(8), np.eye(13)

    def full_response(z):
        # Reference solves the ambient 13-dimensional problem, with fresh U.
        return u.conj().T @ np.linalg.solve(z * i13 - h, u)

    training_full = np.array([full_response(z) for z in training_z])
    training_norms = np.linalg.norm(training_full, axis=(1, 2))

    def model(q):
        br = b @ q
        dr = q.conj().T @ d @ q
        return {
            "q": q, "br": br, "dr": dr,
            "hr": np.block([[a, br], [br.conj().T, dr]]),
        }

    def schur_response(m, z):
        r = m["q"].shape[1]
        if r:
            sigma = m["br"] @ np.linalg.solve(z*np.eye(r) - m["dr"], m["br"].conj().T)
            return np.linalg.solve(z*i8 - a - sigma, i8)
        return np.linalg.solve(z*i8 - a, i8)

    def train(m):
        return np.array([norm(schur_response(m, z) - g)/gn
                         for z, g, gn in zip(training_z, training_full, training_norms)])

    _, singular, vh = np.linalg.svd(b, full_matrices=True)
    svd_vectors = vh.conj().T
    de, dv = np.linalg.eigh(d)
    gram = np.zeros((5, 5), dtype=complex)
    for z in training_z:
        x = np.linalg.solve(z*np.eye(5) - d, b.conj().T)
        gram += x @ x.conj().T / len(training_z)
    gram = hermitian(gram)
    ge, gv = np.linalg.eigh(gram)
    ge, gv = ge[::-1], gv[:, ::-1]

    models, training_errors, subset_rows = {}, {}, []
    selected_subsets = {}
    for rank in range(6):
        candidates = []
        for indices in itertools.combinations(range(5), rank):
            q = dv[:, list(indices)]
            m = model(q)
            errors = train(m)
            score = float(np.max(errors))
            subset_rows.append({"rank": rank, "indices": list(indices), "max_relative_error": score})
            candidates.append((indices, m, errors, score))
        best_score = min(row[3] for row in candidates)
        tied = [row for row in candidates if row[3] <= best_score + 1e-12]
        indices, m, errors, score = min(tied, key=lambda row: row[0])
        mid = f"poles_r{rank}"
        models[mid], training_errors[mid] = m, errors
        selected_subsets[str(rank)] = list(indices)
        for family, vectors in (("svd", svd_vectors), ("snapshots", gv)):
            mid = f"{family}_r{rank}"
            models[mid] = model(vectors[:, :rank])
            training_errors[mid] = train(models[mid])

    choices = {}
    for family in FAMILIES:
        choices[family] = {}
        for budget in BUDGETS:
            qualifying = [r for r in range(6)
                          if np.max(training_errors[f"{family}_r{r}"]) <= budget]
            choices[family][str(budget)] = min(qualifying) if qualifying else None

    gaps = {
        "B_singular_values": descending_gaps(singular),
        "D_eigenvalues_increasing": descending_gaps(de),
        "snapshot_eigenvalues": descending_gaps(ge),
    }
    ambiguity = {}
    for mid in models:
        family, rs = mid.rsplit("_r", 1)
        rank = int(rs)
        if rank in (0, 5):
            ambiguity[mid] = False
        elif family == "svd":
            ambiguity[mid] = bool(gaps["B_singular_values"][rank-1] <= TOL)
        elif family == "snapshots":
            ambiguity[mid] = bool(gaps["snapshot_eigenvalues"][rank-1] <= TOL)
        else:
            selected = set(selected_subsets[str(rank)])
            scale = max(max(abs(de)), np.finfo(float).tiny)
            ambiguity[mid] = any(abs(de[i]-de[j])/scale <= TOL
                                 for i in selected for j in set(range(5))-selected)

    # This file is fully written and hashed BEFORE any held-out arrays or calls.
    selection = {
        "epistemic": "DIAGNOSTIC; independent algorithm shares upstream geometry",
        "input_sha256": sha256(input_path), "protocol_sha256": sha256(protocol_path),
        "script_sha256": sha256(__file__),
        "route": "real26 complex branch, independent localization bases, Schur solves",
        "training_point_count": len(training_z),
        "training_grid": {"x": "linspace(-3,3,61)", "eta": [.05, .2, .7]},
        "raw_geometry_defects": defects,
        "singular_values_B": singular.tolist(), "eigenvalues_D": de.tolist(),
        "snapshot_gram_eigenvalues_descending": ge.tolist(),
        "relative_spectral_gaps": gaps,
        "truncated_subspace_ambiguous": ambiguity,
        "all_pole_subset_training_scores": subset_rows,
        "selected_pole_indices": selected_subsets,
        "training_maximum_by_model": {k: float(np.max(v)) for k, v in training_errors.items()},
        "training_selected_rank_by_budget": choices,
        "selected_basis_coordinates_in_fresh_complement": {k: complex_json(v["q"]) for k,v in models.items()},
        "freeze_order": "Written before generating or evaluating held-out response and dynamics points",
    }
    selection_path = root / "independent_selection.json"
    save_json(selection_path, selection)
    frozen_sha = sha256(selection_path)
    print(json.dumps({"stage": "training_selection_frozen", "sha256": frozen_sha, "choices": choices}), flush=True)

    heldout_z = np.array([x + 1j*eta for eta in (.05, .2, .7)
                          for x in -2.995 + .01*np.arange(600)])
    times = np.array(sorted(set(((np.arange(256)+.5)*2*np.pi/256).tolist()
                               + [0, .25, 1, float(np.pi), float(2*np.pi)])))
    assert len(heldout_z) == 1800 and len(times) == 261
    heldout_full = np.array([full_response(z) for z in heldout_z])
    heldout_norms = np.linalg.norm(heldout_full, axis=(1,2))
    time_full = np.array([u.conj().T @ expm(-1j*t*h) @ u for t in times])
    time_norms = np.linalg.norm(time_full, axis=(1,2))
    full_transfer = 1 - time_norms**2/8
    arrays = {
        "training_z": training_z, "heldout_z": heldout_z, "times": times,
        "H_independent_complex": h, "P_independent_complex": p,
        "U_independent": u, "F_independent": f,
        "full_average_transfer": full_transfer,
    }
    metrics = {}
    for family in FAMILIES:
        for rank in range(6):
            mid = f"{family}_r{rank}"
            m, q = models[mid], models[mid]["q"]
            gr = np.array([schur_response(m, z) for z in heldout_z])
            delta = gr-heldout_full
            relative = np.linalg.norm(delta, axis=(1,2))/heldout_norms
            absolute = np.array([opnorm(v) for v in delta])
            ur = np.array([expm(-1j*t*m["hr"])[:8, :8] for t in times])
            time_delta = ur-time_full
            time_relative = np.linalg.norm(time_delta, axis=(1,2))/time_norms
            time_absolute = np.array([opnorm(v) for v in time_delta])
            transfer = 1 - np.linalg.norm(ur, axis=(1,2))**2/8
            transfer_error = abs(transfer-full_transfer)
            qd = np.eye(5, dtype=complex) if rank == 0 else (np.empty((5,0), dtype=complex) if rank == 5 else null_space(q.conj().T))
            discarded_direct = b @ qd
            mixing = q.conj().T @ d @ qd
            c = np.vstack((discarded_direct, mixing))
            cn = opnorm(c)
            spectral_bound = cn/heldout_z.imag**2
            time_bound = np.minimum(2, abs(times)*cn)
            # r=5 has exactly zero embedding C; permit floating solve noise in gate.
            spectral_excess = float(np.max(absolute-spectral_bound))
            time_excess = float(np.max(time_absolute-time_bound))
            metrics[mid] = {
                "rank": rank, "training_maximum_relative": float(np.max(training_errors[mid])),
                "response_maximum_relative": float(np.max(relative)),
                "response_maximum_absolute_operator": float(np.max(absolute)),
                "time_maximum_relative": float(np.max(time_relative)),
                "time_maximum_absolute_operator": float(np.max(time_absolute)),
                "transfer_maximum_absolute_difference": float(np.max(transfer_error)),
                "per_eta_response_maximum_relative": {str(eta): float(np.max(relative[np.isclose(heldout_z.imag,eta)])) for eta in (.05,.2,.7)},
                "per_eta_response_maximum_absolute_operator": {str(eta): float(np.max(absolute[np.isclose(heldout_z.imag,eta)])) for eta in (.05,.2,.7)},
                "time_prefix_maximum_relative": {str(horizon): float(np.max(time_relative[times<=horizon])) for horizon in (.25,1,float(2*np.pi))},
                "time_prefix_transfer_maximum_absolute_difference": {str(horizon): float(np.max(transfer_error[times<=horizon])) for horizon in (.25,1,float(2*np.pi))},
                "discarded_direct_coupling_operator_norm": opnorm(discarded_direct),
                "retained_discarded_mixing_operator_norm": opnorm(mixing),
                "embedding_C_operator_norm": cn,
                "spectral_bound_maximum_excess": spectral_excess,
                "time_bound_maximum_excess": time_excess,
                "bounds_pass_at_tolerance": spectral_excess <= TOL and time_excess <= TOL,
                "retained_basis_orthonormality": norm(q.conj().T@q - np.eye(rank)),
                "reduced_H_hermitian_defect": norm(m["hr"]-m["hr"].conj().T),
                "raw_average_transfer_range": [float(np.min(transfer)),float(np.max(transfer))],
            }
            arrays.update({
                f"projector_{mid}": f @ q @ q.conj().T @ f.conj().T,
                f"training_relative_{mid}": training_errors[mid],
                f"response_relative_{mid}": relative,
                f"response_absolute_{mid}": absolute,
                f"time_relative_{mid}": time_relative,
                f"time_absolute_{mid}": time_absolute,
                f"transfer_error_{mid}": transfer_error,
                f"average_transfer_{mid}": transfer,
                f"spectral_bound_{mid}": spectral_bound,
                f"time_bound_{mid}": time_bound,
            })

    verdicts = {}
    for family in FAMILIES:
        verdicts[family] = {}
        for budget in BUDGETS:
            rank = choices[family][str(budget)]
            row = metrics[f"{family}_r{rank}"] if rank is not None else None
            spectral_pass = bool(row and row["response_maximum_relative"] <= budget)
            time_pass = bool(row and row["time_maximum_relative"] <= budget)
            verdicts[family][str(budget)] = {
                "training_selected_rank": rank, "heldout_response_pass": spectral_pass,
                "heldout_time_pass": time_pass, "joint_validated": spectral_pass and time_pass,
                "response_maximum_relative": row["response_maximum_relative"] if row else None,
                "time_maximum_relative": row["time_maximum_relative"] if row else None,
            }
    full_rank_checks = {family: {
        "response": metrics[f"{family}_r5"]["response_maximum_relative"] < TOL,
        "dynamics": metrics[f"{family}_r5"]["time_maximum_relative"] < TOL,
    } for family in FAMILIES}
    assert all(all(row.values()) for row in full_rank_checks.values()), full_rank_checks
    assert all(row["bounds_pass_at_tolerance"] for row in metrics.values())
    assert sha256(selection_path) == frozen_sha
    result = {
        "epistemic": "DIAGNOSTIC; floating point maxima on finite deterministic grids",
        "route": selection["route"],
        "input_sha256": sha256(input_path), "protocol_sha256": sha256(protocol_path),
        "selection_sha256": frozen_sha, "selection_unchanged_after_holdout": True,
        "script_sha256": sha256(__file__),
        "training_point_count": len(training_z), "heldout_response_point_count": len(heldout_z), "heldout_time_point_count": len(times),
        "raw_geometry_defects": defects, "selected_pole_indices": selected_subsets,
        "training_selected_rank_by_budget": choices,
        "models": metrics, "heldout_verdicts": verdicts,
        "full_rank_controls": full_rank_checks,
        "all_identity_and_bound_checks_passed": True,
        "elapsed_seconds": time.perf_counter()-started,
        "limitations": [
            "Independent implementation shares the frozen upstream numerical geometry.",
            "Coordinate phases and within-sector bases may differ from the primary route.",
            "A failed reduced-rank budget remains a failure; r=5 is an exactness control.",
            "No continuous-domain, all-time, physical, interval-certified, or probabilistic error claim.",
        ],
    }
    save_json(root/"independent_approximation_results.json",result)
    np.savez_compressed(root/"independent_approximation_matrices.npz",**arrays)
    print(json.dumps({"stage":"completed", "elapsed_seconds":result["elapsed_seconds"], "heldout_verdicts":verdicts}),flush=True)


if __name__ == "__main__":
    main()
