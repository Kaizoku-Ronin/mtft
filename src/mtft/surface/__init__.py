"""mtft.surface — the Modular Surface Laboratory as an mtft subpackage (v0.26.0).

Layers, in the order the computation runs, each with its claim class:

  manin            EXACT      Manin complex of X0(N): CRT P^1(Z/N), S/R/T, cells, gates
  cycles           EXACT      tree/cotree integral basis of H_1, unimodular gate
  hodge            CERTIFIED  intersection form from period-dual Whitney forms;
                   DIAGNOSTIC Whitney masses, refinement, polar star, branch split
  transport        EXACT      Hecke/AL on the cycle lattice (route A) vs PARI (route B)  [needs gp]
  hodge_structure  CERTIFIED  J_true from periods; Riemann + Hecke gates; Siegel
                              distance of the Whitney family; elliptic-block j E2  [needs gp]
  gauge            EXACT      gauge theory ON the surface: flux, YM, cusp holonomy,
                   /OVERLAY   AST line-operator census (theorem-gated), Riemann-Roch
  ising            EXACT      Ising model on the dual Manin graph: Fisher decoration +
                              Cimasoni-Reshetikhin spin-structure Pfaffian sum; brute-force
                              two-route gate through genus 5; 4^13 job checkpointed at 143
  frozen           EXACT      N=143 transport + periods frozen with SHA-1 provenance,
                              all gates re-verified at call time without PARI/GP
  dynamics         CERTIFIED  Hamiltonian flows on H_1(R), Lie closure with mandatory
                   /DIAGNOSTIC genericity controls (351 generic / 4 abelian / 127 block)
  bimodule         EXACT      doubled-space (real spectral triple) census: order-zero,
                   /CERTIFIED first-order, one-forms with absolute scales, sector support;
                              AF-09 result at 143 (AL twist fixes order-zero, calculus zero)

Provenance: Modular Surface Laboratory v0.3.0 (Sol) audited and extended to
v0.4.0 (Claude, 2026-09-02); CC-MSL-01 (nu3 at p=2) and CC-MSL-02 (nu2 guard)
carried over.  ``report(N)`` runs the exact and certified layers; the GP
layers run only when PARI/GP is found.
"""
from __future__ import annotations

from typing import Dict

import numpy as np

from . import cycles, gauge, hodge, manin
from .cycles import CycleBasis, tree_cotree
from .manin import Invariants, ManinComplex, cell_complex, invariants, complex_gates, assert_gates

__all__ = ["invariants", "cell_complex", "tree_cotree", "report", "Invariants", "ManinComplex",
           "CycleBasis", "manin", "cycles", "hodge", "gauge", "transport", "hodge_structure",
           "ising", "frozen", "dynamics", "bimodule"]


def report(N: int, max_refinement: int = 2, gp_layers: bool = True, primes=(2, 3)) -> Dict:
    """Exact + certified layers; GP layers (transport, periods) when PARI/GP is present."""
    cx = cell_complex(N)
    g_exact = complex_gates(cx)
    cb = tree_cotree(cx)
    uh = hodge.unweighted_hodge(cx, cb)
    Jint = uh["intersection_cycles"]
    fam = hodge.refinement_family(cx, cb, Jint, max_refinement)
    out = {
        "N": N, "invariants": cx.inv, "cells": cx.counts,
        "gates_exact": g_exact + cb.gates, "gates_certified": uh["gates"],
        "intersection_cycles": Jint, "unweighted": uh, "whitney_family": fam,
        "branch_bracket_ratio": [lv["branch"]["bracket_ratio"] for lv in fam],
    }
    if gp_layers:
        from ..gprun import find_gp
        if find_gp() is not None:
            from . import hodge_structure as HS, transport as TR
            tr = TR.run(cx, cb, Jint, primes)
            ops = {f"T{p}": M for p, M in tr.hecke.items()}
            ops.update({f"W{Q}": M for Q, M in tr.atkin_lehner.items()})
            hs = HS.from_periods(N, tr.period_Q, Jint, ops)
            out.update({"transport": tr, "gates_transport": tr.gates, "hodge_structure": hs,
                        "gates_hodge_structure": hs.gates,
                        "family_distances": HS.family_distances(hs, fam)})
    return out


def all_pass(rep: Dict) -> bool:
    keys = ["gates_exact", "gates_certified", "gates_transport"]
    ok = all(g["status"] == "PASS" for k in keys if k in rep for g in rep[k])
    if "hodge_structure" in rep:
        from .hodge_structure import gates_pass
        ok = ok and gates_pass(rep["hodge_structure"])
    return ok
