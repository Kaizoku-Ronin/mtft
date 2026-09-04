"""mtft.surface.frozen — frozen certified data for X0(143), re-verified at call time without PARI/GP.

The GP layers (transport, periods) are the only parts of mtft.surface that need
PARI.  Their N = 143 outputs are frozen here with SHA-1 provenance so that a
runtime without gp can still verify every downstream gate exactly: integrality,
AL involutions and commutation, Hecke self-adjointness for the intersection
form, and the Riemann bilinear relations + Hecke/AL commutators of J_true.
The deterministic cycle basis is frozen too and compared with the live one.
"""
from __future__ import annotations

import hashlib
import json
from importlib import resources

import numpy as np

_FILE = "x0143_certified.json"


def _load() -> dict:
    with resources.files("mtft.surface._data").joinpath(_FILE).open() as fh:
        return json.load(fh)


def x0143(verify: bool = True) -> dict:
    raw = _load()
    d = {k: np.array(v) for k, v in raw["arrays"].items()}
    d["provenance"] = raw["provenance"]
    if verify:
        d["gates"] = verify_gates(d)
        bad = [k for k, v in d["gates"].items() if not v]
        if bad:
            raise AssertionError(f"frozen X0(143) data failed gates: {bad}")
    return d


def verify_gates(d: dict, tol: float = 1e-11) -> dict:
    from .manin import cell_complex
    from .cycles import tree_cotree
    I = np.eye(26, dtype=np.int64)
    Jint = d["intersection_cycles"].astype(np.int64)
    ops = {k: d[k].astype(np.int64) for k in ("T2", "T3", "W11", "W13", "W143")}
    Q = d["period_Q"]
    g = 13
    S = np.block([[np.zeros((g, g)), -np.eye(g)], [np.eye(g), np.zeros((g, g))]])
    J = np.linalg.solve(Q, S @ Q)
    G = Jint @ J
    if np.linalg.eigvalsh((G + G.T) / 2)[0] < 0:
        J, G = -J, -G
    rel = lambda a, b: np.linalg.norm(a - b) / max(np.linalg.norm(b), 1e-300)
    gates = {
        "intersection_skew_unimodular": bool(np.array_equal(Jint.T, -Jint) and abs(round(np.linalg.det(Jint))) == 1),
        "AL_involutions": all(np.array_equal(ops[w] @ ops[w], I) for w in ("W11", "W13", "W143")),
        "W11_W13_equals_W143": bool(np.array_equal(ops["W11"] @ ops["W13"], ops["W143"])),
        "Hecke_commutes_with_AL": all(np.array_equal(ops[t] @ ops[w], ops[w] @ ops[t]) for t in ("T2", "T3") for w in ("W11", "W13")),
        "AL_preserve_intersection": all(np.array_equal(ops[w].T @ Jint @ ops[w], Jint) for w in ("W11", "W13", "W143")),
        "Hecke_selfadjoint": all(np.array_equal(ops[t].T @ Jint, Jint @ ops[t]) for t in ("T2", "T3")),
        "AL_traces_2_m2_m18": (int(np.trace(ops["W11"])), int(np.trace(ops["W13"])), int(np.trace(ops["W143"]))) == (2, -2, -18),
        "riemann_I_symplectic": rel(J.T @ Jint @ J, Jint) < tol,
        "riemann_II_positive": bool(np.linalg.eigvalsh((G + G.T) / 2)[0] > 0.05),
        "J_commutes_with_Hecke_AL": all(rel(J @ A, A @ J) < tol for A in ops.values()),
        "frozen_J_matches_recomputed": rel(d["J_true"], J) < 1e-9,
    }
    cx = cell_complex(143)
    cb = tree_cotree(cx)
    gates["live_cycle_basis_equals_frozen"] = bool(np.array_equal(cb.basis_matrix, d["basis_matrix"]))
    if "U13" in d:                                   # v0.26.0 extension
        U11, U13 = d["U11"].astype(np.int64), d["U13"].astype(np.int64)
        gates["U_commute_with_good_Hecke"] = all(np.array_equal(U @ T, T @ U) for U in (U11, U13) for T in (ops["T2"], ops["T3"]))
        gates["U_commute_with_J"] = all(rel(J @ U, U @ J) < tol for U in (U11, U13))
        Gs = (G + G.T) / 2
        gates["AL_adjoint_identity_U13"] = rel(np.linalg.solve(Gs, U13.T @ Gs), ops["W13"] @ U13 @ ops["W13"]) < 1e-10
        gates["AL_adjoint_identity_U11"] = rel(np.linalg.solve(Gs, U11.T @ Gs), ops["W11"] @ U11 @ ops["W11"]) < 1e-10
        gates["U11_squared_identity"] = bool(np.array_equal(U11 @ U11, I))
        x = ops["T2"]
        polys = {"ell": x, "ghost": x + 2 * I, "q4": x @ x @ x @ x - 3 * x @ x @ x - x @ x + 5 * x + I,
                 "q6": np.linalg.matrix_power(x, 6) - 10 * np.linalg.matrix_power(x, 4) + 2 * x @ x @ x + 24 * x @ x - 7 * x - 12 * I}
        ok_blocks, ok_cov = True, True
        for k, f in polys.items():
            B = d[f"block_{k}"].astype(np.int64)
            ok_blocks = ok_blocks and not np.any(f @ B) and np.linalg.matrix_rank(B.astype(float)) == B.shape[1]
            P = B.T @ Jint @ B
            cov = float(np.sqrt(abs(np.linalg.det(B.T @ Gs @ B))))
            ok_cov = ok_cov and abs(cov * cov - abs(round(np.linalg.det(P.astype(float))))) < 1e-6 * cov * cov
        gates["hecke_block_lattices_annihilated_and_full_rank"] = ok_blocks
        gates["hodge_covolume_squared_equals_polarization_det"] = ok_cov
    return gates


def block_invariants() -> dict:
    """EXACT polarization types of the saturated Hecke-block lattices (elementary divisors of the
    intersection form restricted to H_1(Z) ∩ V_block) and the Hodge covolumes, which equal the
    square roots of the polarization determinants (J-stable sublattice: det J|_B = 1)."""
    return _load()["block_invariants"]


def sha1_of(obj) -> str:
    return hashlib.sha1(json.dumps(obj, sort_keys=True).encode()).hexdigest()
