"""Exact zero-field Ising energy counts by packed-polynomial variable elimination.

For a multigraph, D[k] counts spin assignments with k disagreeing edges.
Z(beta) = exp(beta * E) * sum_k D[k] * exp(-2 * beta * k).
Self loops are retained in E and always contribute zero disagreements.

All counting arithmetic uses Python integers. Base 2**(n+1) prevents carries
between polynomial coefficients: no partial count can exceed 2**n.
The min-fill order is a heuristic; its width is an upper bound, not a proof
of optimal treewidth. No Pfaffian or spin-structure code is used here.

Run with mtft==0.26.0 and numpy installed:
    python ising_exact.py --levels 6 11 15 35 55 77 105 143 --out results.json
"""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import platform
import time

import numpy as np


def min_fill_order(n, edges, reverse_ties=False):
    adj = {v: set() for v in range(n)}
    for u, v in edges:
        if u != v:
            adj[u].add(v)
            adj[v].add(u)
    order, widths = [], []
    while adj:
        def score(v):
            nb = sorted(adj[v])
            fill = sum(b not in adj[a] for i, a in enumerate(nb) for b in nb[i+1:])
            return fill, len(nb), -v if reverse_ties else v
        v = min(adj, key=score)
        nb = adj[v]
        order.append(v)
        widths.append(len(nb))
        for a in nb:
            adj[a].update(nb - {a})
            adj[a].remove(v)
        del adj[v]
    return order, widths


def packed_dos(n, edges, reverse_ties=False, max_width=18):
    """Return the exact integer density of states and an elimination certificate."""
    edges = [tuple(map(int, e)) for e in edges]
    if any(min(u, v) < 0 or max(u, v) >= n for u, v in edges):
        raise ValueError("edge endpoint outside 0..n-1")
    start = time.perf_counter()
    order, widths = min_fill_order(n, edges, reverse_ties)
    if max(widths, default=0) > max_width:
        raise RuntimeError(f"elimination width exceeds budget {max_width}")
    digit_bits = n + 1
    base = 1 << digit_bits
    # Every variable has a neutral unary factor, including isolated spins.
    factors = [((v,), np.array([1, 1], dtype=object)) for v in range(n)]
    for u, v in edges:
        if u != v:
            factors.append((tuple(sorted((u, v))), np.array([[1, base], [base, 1]], dtype=object)))
    peak_entries = 1
    for v in order:
        selected = [(s, a) for s, a in factors if v in s]
        factors = [(s, a) for s, a in factors if v not in s]
        scope = tuple(sorted(set().union(*(s for s, _ in selected))))
        peak_entries = max(peak_entries, 1 << len(scope))
        product = np.ones((2,) * len(scope), dtype=object)
        for s, a in selected:
            product *= a.reshape(tuple(2 if w in s else 1 for w in scope))
        axis = scope.index(v)
        factors.append((tuple(w for w in scope if w != v), product.sum(axis=axis)))
    packed = 1
    for scope, a in factors:
        assert not scope
        packed *= int(a)
    mask = base - 1
    counts = [(packed >> (digit_bits * k)) & mask for k in range(len(edges) + 1)]
    assert packed >> (digit_bits * (len(edges) + 1)) == 0
    assert sum(counts) == 1 << n
    return counts, {
        "order": order, "neighbor_widths": widths,
        "width_upper_bound": max(widths, default=0),
        "peak_spin_table_entries": peak_entries,
        "packing_base_power_of_two": digit_bits,
        "seconds": time.perf_counter() - start,
    }


def even_subgraph_counts(counts, n, E):
    """Exact MacWilliams/Fourier transform: spin cut counts -> even-subgraph counts.

    A(t) = 2**(-n) * sum_k D[k]*(1+t)**(E-k)*(1-t)**k.
    The coefficient of t**j in the product is a Krawtchouk polynomial.
    """
    acc = [0] * (E + 1)
    for k, dk in enumerate(counts):
        if not dk:
            continue
        prev, cur = 0, 1
        for j in range(E + 1):
            acc[j] += dk * cur
            if j < E:
                numerator = (E - 2*k) * cur - (E - j + 1) * prev
                assert numerator % (j + 1) == 0
                prev, cur = cur, numerator // (j + 1)
    den = 1 << n
    assert all(a % den == 0 for a in acc)
    out = [a // den for a in acc]
    assert all(a >= 0 for a in out)
    return out


def thermodynamics(counts, n, E, beta):
    """Dimensionless J=k_B=1 observables; numerical evaluation of exact counts."""
    k = np.array([k for k, a in enumerate(counts) if a], dtype=float)
    d = np.array([a for a in counts if a], dtype=float)
    energy = 2*k - E
    logw = np.log(d) - beta*energy
    m = float(logw.max())
    w = np.exp(logw - m)
    p = w / w.sum()
    mean_energy = float(p @ energy)
    variance = float(p @ ((energy - mean_energy)**2))
    return {"beta": float(beta), "log_Z": m + math.log(float(w.sum())),
            "energy_per_spin": mean_energy/n,
            "heat_capacity_per_spin": beta*beta*variance/n,
            "mean_disagreeing_edges": float(p @ k),
            "energy_variance": variance}


def brute_dos(n, edges, chunk=1 << 18):
    """Independent direct enumeration, used only for n<=24 verification."""
    if n > 24:
        raise ValueError("enumeration limited to 24 spins")
    out = np.zeros(len(edges) + 1, dtype=np.int64)
    for start in range(0, 1 << n, chunk):
        cfg = np.arange(start, min(start + chunk, 1 << n), dtype=np.int64)
        cut = np.zeros(len(cfg), dtype=np.int64)
        for u, v in edges:
            cut += ((cfg >> u) ^ (cfg >> v)) & 1
        out += np.bincount(cut, minlength=len(out))
    return out.tolist()


def float_partition(n, edges, beta):
    """Independent floating evaluation, without packed counts or their transform."""
    order, _ = min_fill_order(n, edges, reverse_ties=True)
    factors = [((v,), np.ones(2)) for v in range(n)]
    loop_factor = 1.0
    for u, v in edges:
        if u == v:
            loop_factor *= math.exp(beta)
        else:
            a, b = math.exp(beta), math.exp(-beta)
            factors.append((tuple(sorted((u, v))), np.array([[a, b], [b, a]])))
    for v in order:
        selected = [(s, a) for s, a in factors if v in s]
        factors = [(s, a) for s, a in factors if v not in s]
        scope = tuple(sorted(set().union(*(s for s, _ in selected))))
        product = np.ones((2,) * len(scope))
        for s, a in selected:
            product *= a.reshape(tuple(2 if w in s else 1 for w in scope))
        factors.append((tuple(w for w in scope if w != v), product.sum(axis=scope.index(v))))
    return loop_factor * math.prod(float(a) for _, a in factors)


def inspect_level(N):
    from mtft.surface.manin import cell_complex
    cx = cell_complex(N)
    n, E = len(cx.faces), len(cx.edges)
    edges = [(cx.face_of[d], cx.face_of[cx.S(d)]) for d in cx.edges]
    counts, cert = packed_dos(n, edges)
    second, cert2 = packed_dos(n, edges, reverse_ties=True)
    even = even_subgraph_counts(counts, n, E)
    loops = sum(a == b for a, b in edges)
    mult = Counter(tuple(sorted(e)) for e in edges if e[0] != e[1])
    total = 1 << n
    mu = Fraction(sum(k*d for k, d in enumerate(counts)), total)
    var = Fraction(sum(k*k*d for k, d in enumerate(counts)), total) - mu**2
    gates = {"distinct_elimination_orders": cert["order"] != cert2["order"],
             "two_orders_identical_integers": counts == second,
             "total_spin_states": sum(counts) == total,
             "global_flip_even_counts": all(d % 2 == 0 for d in counts),
             "connected_ground_degeneracy_two": counts[0] == 2,
             "infinite_temperature_mean_cut": mu == Fraction(E-loops, 2),
             "infinite_temperature_cut_variance": var == Fraction(sum(m*m for m in mult.values()), 4),
             "cycle_space_total": sum(even) == 1 << (E-n+1),
             "empty_even_subgraph": even[0] == 1,
             "single_edge_even_subgraphs_are_loops": even[1] == loops}
    if n <= 24:
        gates["exact_brute_force_histogram"] = counts == brute_dos(n, edges)
    numerical = []
    for beta in [-0.7, 0.0, 0.2, math.atanh(1/math.sqrt(3)), 1.1]:
        a = thermodynamics(counts, n, E, beta)
        Z2 = float_partition(n, edges, beta)
        a["floating_contraction_relative_difference"] = abs(math.expm1(math.log(Z2) - a["log_Z"]))
        numerical.append(a)
    gates["independent_float_contraction"] = max(r["floating_contraction_relative_difference"] for r in numerical) < 1e-11
    assert all(gates.values()), (N, gates)
    nonzero = [k for k, d in enumerate(counts) if d]
    return {"N": N, "genus": cx.inv.genus, "spins": n, "edges_count": E,
            "edges": edges, "self_loops": loops, "density_of_states": counts,
            "even_subgraph_counts": even, "elimination": cert, "alternate_elimination": cert2,
            "gates": gates, "infinite_temperature_mean_cut_exact": str(mu),
            "infinite_temperature_cut_variance_exact": str(var),
            "min_cut": nonzero[0], "max_cut": nonzero[-1],
            "max_cut_degeneracy": counts[nonzero[-1]],
            "numerical_checks": numerical}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--levels", nargs="+", type=int, default=[6, 11, 15, 35, 55, 77, 105, 143])
    parser.add_argument("--out", type=Path, default=Path("results.json"))
    args = parser.parse_args()
    import mtft
    report = {"mtft_version": mtft.__version__, "python": platform.python_version(),
              "numpy": np.__version__, "count_class": "Cert / EXACT integer finite-graph enumeration",
              "thermal_curve_class": "DIAGNOSTIC float evaluation of exact finite sums",
              "levels": []}
    for N in args.levels:
        row = inspect_level(N)
        report["levels"].append(row)
        args.out.write_text(json.dumps(report, indent=2) + "\n")
        print(f"N={N}: {row['spins']} spins, width<={row['elimination']['width_upper_bound']}, "
              f"exact DOS in {row['elimination']['seconds']:.3f}s, "
              f"max cut={row['max_cut']}, degeneracy={row['max_cut_degeneracy']}, "
              f"{sum(row['gates'].values())} gates PASS", flush=True)


if __name__ == "__main__":
    main()
