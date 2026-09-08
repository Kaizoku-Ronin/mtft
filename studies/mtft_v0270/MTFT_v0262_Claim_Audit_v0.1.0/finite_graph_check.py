#!/usr/bin/env python3
"""Finite Manin-dual Ising and diffusion audit; never infer a critical exponent.

Run against mtft 0.26.2. Exact integer density of states uses package route C;
independent small-level brute force and graph isomorphism verify the negative
control. Laplacian = B B^T for the oriented dual incidence matrix. Self-loops
give zero columns (diffusion) but a constant ferromagnetic energy (Ising).
Numeric diffusion and heat capacities are DIAGNOSTIC, not continuum claims.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path

import mtft
import numpy as np
import sympy as sp
from scipy.optimize import minimize_scalar
from mtft.surface.ising import density_of_states, even_subgraph_polynomial
from mtft.surface.manin import cell_complex


def make_graph(level):
    cx = cell_complex(level)
    n = len(cx.faces)
    edges = [(cx.face_of[d], cx.face_of[cx.S(d)]) for d in cx.edges]
    incidence = np.zeros((n, len(edges)), dtype=np.int64)
    adjacency = np.zeros((n, n), dtype=np.int64)
    for e, (u, v) in enumerate(edges):
        incidence[u, e] += 1
        incidence[v, e] -= 1
        adjacency[u, v] += 1
        adjacency[v, u] += 1
    return cx, edges, incidence @ incidence.T, adjacency


def brute_dos(n, edges):
    answer = [0] * (len(edges) + 1)
    for bits in itertools.product((0, 1), repeat=n):
        answer[sum(bits[u] != bits[v] for u, v in edges)] += 1
    return answer


def thermo(D, edges, vertices, beta):
    k = np.flatnonzero(np.array(D, dtype=float))
    log_weights = np.log(np.array(D, dtype=float)[k]) - 2 * beta * k
    weights = np.exp(log_weights - max(log_weights))
    weights /= sum(weights)
    energies = 2 * k - edges
    mean = float(weights @ energies)
    variance = float(weights @ ((energies - mean) ** 2))
    return {"beta": float(beta), "energy_per_spin": mean / vertices,
            "heat_capacity_per_spin": beta * beta * variance / vertices}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, default=Path(__file__).with_name('finite_graph_results.json'))
    args = parser.parse_args()
    times = np.logspace(-3, 3, 241)
    result = {
        'package_version': mtft.__version__,
        'source_sha256': {},
        'conventions': {
            'Ising': 'J=1, zero field, H=-sum_edges sigma_u sigma_v; beta dimensionless',
            'Laplacian': 'L=B B^T unnormalized symmetric dual incidence Laplacian; self-loops cancel',
            'heat_trace': 'K(t)=Tr exp(-t L), P(t)=K(t)/n',
            'running_dimension': 'd_eff(t)=2t sum lambda exp(-t lambda)/sum exp(-t lambda)',
            'zero_mode_removed': 'Same diagnostic after omitting the single exact zero eigenvalue',
            'times': '241 log-spaced points from 0.001 through 1000, no power-law fit',
            'claim_class': 'EXACT integer DOS/graph certificates; DIAGNOSTIC floating spectral/thermal values',
        },
        'levels': {},
        'sources': [
            'https://arxiv.org/html/1205.2723v1',
            'https://davidtong.org/pdfs/teaching/statistical-physics/statphys5.pdf',
        ],
    }
    source = Path(mtft.__file__).parent
    for rel in ['surface/ising.py', 'surface/manin.py']:
        result['source_sha256'][rel] = hashlib.sha256((source / rel).read_bytes()).hexdigest()
    adjacencies = {}
    for N in [6, 11, 143]:
        cx, edges, L, adjacency = make_graph(N)
        n, E = len(cx.faces), len(edges)
        D = density_of_states(N)
        width = density_of_states.last_width
        A = even_subgraph_polynomial(D, n, E)
        assert sum(D) == 1 << n
        assert D[0] == 2
        assert A[0] == 1 and sum(A) == 1 << (E - n + 1)
        assert np.array_equal(L, cx.boundary_2.T @ cx.boundary_2)
        evals = np.linalg.eigvalsh(L.astype(float))
        assert abs(evals[0]) < 1e-12 and evals[1] > 1e-12
        evals[0] = 0.0
        x = sp.Symbol('x')
        charpoly = sp.Matrix(L.tolist()).charpoly(x).as_expr()
        factors = sp.factor_list(charpoly)[1]
        named_multiplicities = {}
        for lam in [0, 1, 2, 4, 5]:
            count = 0
            for polynomial, multiplicity in factors:
                if polynomial.subs(x, lam) == 0:
                    count += multiplicity
            named_multiplicities[str(lam)] = count
        heat = []
        for t in times:
            w = np.exp(-t * evals)
            # Stable even if exp(-t lambda_min) underflows in the primed trace.
            wp = np.exp(-t * (evals[1:] - evals[1]))
            heat.append({
                't': float(t), 'K': float(w.sum()), 'P': float(w.mean()),
                'd_eff': float(2 * t * (evals @ w) / w.sum()),
                'd_eff_without_zero': float(2 * t * (evals[1:] @ wp) / wp.sum()),
            })
        peak = minimize_scalar(lambda beta: -thermo(D, E, n, beta)['heat_capacity_per_spin'],
                               bounds=(0.001, 4.0), method='bounded', options={'xatol': 1e-12})
        rec = {
            'genus': cx.inv.genus, 'vertices': n, 'edges': E,
            'self_loops': sum(u == v for u, v in edges),
            'edge_list': edges, 'laplacian': L.tolist(),
            'DOS_exact': D, 'even_subgraph_polynomial_exact': A,
            'elimination_width': width,
            'sum_D': sum(D), 'sum_A': sum(A),
            'A_at_honeycomb_t_diagnostic': float(sum(a / 3 ** (j / 2) for j, a in enumerate(A))),
            'laplacian_charpoly_factorization_exact': [[str(p), int(m)] for p, m in factors],
            'integer_eigenvalue_multiplicities_exact': named_multiplicities,
            'eigenvalues_diagnostic': evals.tolist(),
            'gap_diagnostic': float(evals[1]),
            'heat_trace': heat,
            'heat_capacity_peak_diagnostic': thermo(D, E, n, float(peak.x)),
            'finite_sum_has_no_real_beta_singularity': True,
            'partition_formula_exact': f'Z(beta)=exp({E} beta) sum_k D_k exp(-2 beta k)',
        }
        if n <= 4:
            rec['DOS_matches_independent_bruteforce'] = D == brute_dos(n, edges)
            assert rec['DOS_matches_independent_bruteforce']
            adjacencies[N] = adjacency
        result['levels'][str(N)] = rec
    A6, A11 = adjacencies[6], adjacencies[11]
    permutations = [list(p) for p in itertools.permutations(range(4))
                    if np.array_equal(A6, A11[np.ix_(p, p)])]
    assert permutations
    assert result['levels']['6']['DOS_exact'] == result['levels']['11']['DOS_exact']
    result['negative_control'] = {
        'levels': [6, 11], 'genera': [0, 1],
        'adjacency_isomorphism_permutation': permutations[0],
        'identical_exact_DOS': result['levels']['6']['DOS_exact'],
        'conclusion': 'Distinct surface genera, isomorphic dual multigraphs: identical scalar zero-field Ising thermodynamics at every beta.',
    }
    result['finite_limit_conclusions'] = {
        'including_zero': 'For every fixed connected finite graph, d_eff(t) tends to zero as t tends to zero and as t tends to infinity.',
        'excluding_zero': 'd_eff_primed(t) tends to zero as t tends to zero and grows as 2 lambda_gap t at large t; this is exponential decay, not a dimension.',
        'spectral_zeta': 'The finite sum over positive eigenvalues is entire in s and has no pole from which to infer a continuum spectral dimension.',
        'Ising': 'Z is a finite sum of positive exponentials for real beta; log Z is real analytic. No thermodynamic critical exponent is determined by this one finite graph.',
        'next_gate': 'Register a growing/refining graph family, metric and coupling scaling, and finite-size/continuum window before estimating universal critical behavior.',
    }
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({'output': str(args.output),
                      'negative_control': result['negative_control'],
                      'N143_gap': result['levels']['143']['gap_diagnostic'],
                      'N143_named_eigenvalues': result['levels']['143']['integer_eigenvalue_multiplicities_exact'],
                      'N143_heat_capacity_peak': result['levels']['143']['heat_capacity_peak_diagnostic'],
                      'N143_heat_samples': [result['levels']['143']['heat_trace'][j] for j in [0, 40, 80, 120, 160, 200, 240]]}, indent=2))


if __name__ == '__main__':
    main()
