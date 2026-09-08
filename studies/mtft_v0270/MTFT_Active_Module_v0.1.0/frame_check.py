#!/usr/bin/env python3
"""Independent coordinate/adjoint audit for the active-module experiment.

Numbers are float64 diagnostics using the release's frozen period stage;
`dps` controls upstream period evaluation, not the numpy linear algebra.
No source files are changed. Run with the intended release on PYTHONPATH.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from mtft import hecke, liealg
from mtft.periods import physics


def norm(a):
    return float(np.linalg.norm(a, "fro"))


def linear_realification(a):
    return np.block([[a.real, -a.imag], [a.imag, a.real]])


def antilinear_realification(b):
    return np.block([[b.real, b.imag], [b.imag, -b.real]])


def complex_parts(o):
    """For O(x,y), return L,B such that z maps to L z + B conjugate(z)."""
    n = o.shape[0] // 2
    a, b, c, d = o[:n, :n], o[:n, n:], o[n:, :n], o[n:, n:]
    return (a+d + 1j*(c-b))/2, (a-d + 1j*(c+b))/2


def run(dps=50, seed=143):
    r, ri, fixed = liealg.x0143_fixed_channels(dps=dps, seed=seed)
    g = physics.hodge_metric_hecke(dps)
    j = physics.hodge_structure_hecke(dps)
    n = r.shape[0]//2
    ident = np.eye(2*n)
    js = np.block([[np.zeros((n,n)), -np.eye(n)],
                   [np.eye(n), np.zeros((n,n))]])
    sym = liealg.x0143_symmetry_ops()
    ops = {f"T{p}": np.array(hecke.cuspidal_hecke(p), float)
           for p in (2,3,5,7)}
    ops.update(sym)
    rows = {}
    for name, o in ops.items():
        oa = r @ o @ ri
        scale = norm(oa)
        lin, anti = complex_parts(oa)
        adj = np.linalg.solve(g, o.T @ g)
        rows[name] = {
            "complex_linear_commutator_rel": norm(oa @ js-js @ oa)/scale,
            "complex_antilinear_anticommutator_rel": norm(oa @ js+js @ oa)/scale,
            "linear_component_fraction": np.sqrt(2)*norm(lin)/scale,
            "antilinear_component_fraction": np.sqrt(2)*norm(anti)/scale,
            "complex_reconstruction_rel": norm(oa-linear_realification(lin)-antilinear_realification(anti))/scale,
            "adapted_selfadjoint_defect_rel": norm(oa.T-oa)/scale,
            "transported_hodge_adjoint_residual_rel": norm(r @ adj @ ri-oa.T)/scale,
            "adapted_isometry_defect_rel": norm(oa.T @ oa-ident)/norm(ident),
        }
    star_direct = np.array(hecke.star_involution(), float)
    # Frame conjugacy should preserve all operator residuals once measured
    # in the adapted Hodge-orthonormal frame. The transition is complex linear.
    r2, ri2, _ = liealg.x0143_fixed_channels(dps=dps, seed=seed+1)
    transition = r2 @ ri
    out = {
        "epistemic": "DIAGNOSTIC: float64, frozen certified period inputs; no interval bounds",
        "dps": dps, "frame_seed": seed,
        "definitions": {
            "frame": "adapted real coordinates a=R h, h=Ri a; z=a[:13]+i a[13:]",
            "operator": "O_adapted=R O_Hecke Ri",
            "adjoint": "O_dagger_G=G^-1 O^T G; transported adjoint is O_adapted^T",
            "complex_parts": "O(z)=L z+B conjugate(z); L=(A+D+i(C-B))/2; B=(A-D+i(C+B))/2",
            "primitive_channels": "fixed[i] is B in z->B conjugate(z), not a complex-linear action",
            "antilinear_kernel": "kernel(B conjugate(z))=conjugate(kernel(B)); conjugate its right nullspace/projector",
        },
        "frame": {
            "R_Ri_residual": norm(r @ ri-ident),
            "Ri_R_residual": norm(ri @ r-ident),
            "R_transpose_R_vs_G_rel": norm(r.T @ r-g)/norm(g),
            "Ri_transpose_G_Ri_vs_I_rel": norm(ri.T @ g @ ri-ident)/norm(ident),
            "adapted_J_vs_standard_rel": norm(r @ j @ ri-js)/norm(js),
            "J_squared_residual_rel": norm(j @ j+ident)/norm(ident),
            "G_condition_number": float(np.linalg.cond(g)),
            "STAR_bridge_vs_direct_rel": norm(sym["STAR"]-star_direct)/norm(star_direct),
            "next_seed_transition_isometry_rel": norm(transition.T @ transition-ident)/norm(ident),
            "next_seed_transition_complex_linear_rel": norm(transition @ js-js @ transition)/norm(transition),
            "next_seed_operator_conjugacy_max_rel": max(norm(r2 @ o @ ri2-transition @ (r @ o @ ri) @ transition.T)/norm(r2 @ o @ ri2) for o in ops.values()),
        },
        "primitive_complex_symmetry_residuals": [norm(b-b.T)/norm(b) for b in fixed],
        "operators": rows,
    }
    return out


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=Path(__file__).with_name("frame_check_results.json"))
    parser.add_argument("--dps", type=int, default=50)
    args = parser.parse_args()
    result = run(args.dps)
    args.output.write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps(result, indent=2))
