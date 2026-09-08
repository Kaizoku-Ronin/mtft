"""MTFT Ising companion v0.2.0: fields, signs, two layers, and exact sampling.

All polynomial coefficients and rational-weight sampling probabilities use
unbounded integers. Numerical thermodynamic evaluation is a separate operation.
Spins use bits: 0 -> -1, 1 -> +1. A preferred edge bit b has energy cost
(spin_u XOR spin_v) XOR b; b=0 is ferro and b=1 is antiferro.
"""
from __future__ import annotations
from fractions import Fraction
import math
import random
import time
import numpy as np

from ising_exact import min_fill_order


def graph(N):
    from mtft.surface.manin import cell_complex
    cx = cell_complex(N)
    edges = [(cx.face_of[d], cx.face_of[cx.S(d)]) for d in cx.edges]
    return cx, len(cx.faces), edges


def eliminate_sum(n, edges, factors, q=2, reverse_ties=False, keep_conditionals=False,
                  max_entries=1 << 20):
    """Exact sum-product elimination; optional tables support reverse sampling."""
    order, widths = min_fill_order(n, edges, reverse_ties)
    peak = q ** (max(widths, default=0) + 1)
    if peak > max_entries:
        raise RuntimeError(f"{peak} entries exceeds the declared contraction budget")
    factors = list(factors)
    records = []
    started = time.perf_counter()
    for v in order:
        selected = [(s, a) for s, a in factors if v in s]
        factors = [(s, a) for s, a in factors if v not in s]
        scope = tuple(sorted(set().union(*(s for s, _ in selected))))
        table = np.ones((q,) * len(scope), dtype=object)
        for s, a in selected:
            table *= np.asarray(a, dtype=object).reshape(tuple(q if w in s else 1 for w in scope))
        rest = tuple(w for w in scope if w != v)
        axis = scope.index(v)
        if keep_conditionals:
            records.append((v, rest, np.moveaxis(table, axis, -1)))
        factors.append((rest, table.sum(axis=axis)))
    total = 1
    for scope, a in factors:
        assert not scope
        total *= int(a)
    return total, {"order": order, "width_upper_bound": max(widths, default=0),
                   "peak_spin_table_entries": peak, "seconds": time.perf_counter() - started}, records


def signed_counts(n, edges, preferred=None, joint=False, reverse_ties=False):
    """D[k] or D[k,m], m=#up spins. Self-loops count as frustrated iff b=1."""
    E = len(edges)
    preferred = [0] * E if preferred is None else list(map(int, preferred))
    if len(preferred) != E or any(b not in (0, 1) for b in preferred):
        raise ValueError("one preferred bit per edge is required")
    bits = n + 1
    B = 1 << bits
    stride = E + 1
    marker = 1 << (bits * stride) if joint else 1
    factors = [((v,), np.array([1, marker], dtype=object)) for v in range(n)]
    loop_cost = 0
    for (u, v), b in zip(edges, preferred):
        if u == v:
            loop_cost += b
        else:
            a = np.array([[B**b, B**(1-b)], [B**(1-b), B**b]], dtype=object)
            factors.append((tuple(sorted((u, v))), a))
    packed, cert, _ = eliminate_sum(n, edges, factors, reverse_ties=reverse_ties)
    packed <<= bits * loop_cost
    mask = B - 1
    count_columns = n + 1 if joint else 1
    out = [[(packed >> (bits*(k+stride*m))) & mask for m in range(count_columns)] for k in range(E+1)]
    assert packed >> (bits*stride*count_columns) == 0
    assert sum(map(sum, out)) == 1 << n
    cert.update({"coefficient_bits": bits, "magnetization_stride": stride if joint else None,
                 "count_class": "Cert / EXACT integers"})
    return (out if joint else [row[0] for row in out]), cert


def brute_joint(n, edges, preferred=None):
    if n > 20:
        raise ValueError("full joint brute-force gate limited to 20 spins")
    preferred = [0]*len(edges) if preferred is None else preferred
    cfg = np.arange(1 << n, dtype=np.int64)
    m = np.zeros(len(cfg), dtype=np.int64)
    k = np.zeros(len(cfg), dtype=np.int64)
    for v in range(n):
        m += (cfg >> v) & 1
    for (u, v), b in zip(edges, preferred):
        k += (((cfg >> u) ^ (cfg >> v)) & 1) ^ b
    return np.bincount(k*(n+1)+m, minlength=(len(edges)+1)*(n+1)).reshape(len(edges)+1, n+1).tolist()


class IntegerSampler:
    """Independent exact equilibrium samples at rational local weight ratios.

    Each satisfied edge gets weight a, each frustrated edge weight d. For
    uniform coupling beta=log(a/d)/2. Field weights (w0,w1) mean
    eta=log(w1/w0)/2, where p(s) is proportional to exp(beta*S+eta*M).
    Random integer draws implement conditional probabilities without rounding.
    """
    def __init__(self, n, edges, a=4, d=1, field_weights=(1, 1), preferred=None):
        if any(not isinstance(x, int) or x <= 0 for x in (a, d, *field_weights)):
            raise ValueError("strictly positive integer weights required")
        self.n, self.edges = n, list(edges)
        self.a, self.d, self.field_weights = a, d, tuple(field_weights)
        self.preferred = [0]*len(edges) if preferred is None else list(preferred)
        factors = [((v,), np.array(field_weights, dtype=object)) for v in range(n)]
        loop_factor = 1
        for (u, v), b in zip(edges, self.preferred):
            equal, unequal = (a, d) if b == 0 else (d, a)
            if u == v:
                loop_factor *= equal
            else:
                factors.append((tuple(sorted((u, v))), np.array([[equal, unequal], [unequal, equal]], dtype=object)))
        total, self.certificate, self.records = eliminate_sum(n, edges, factors, keep_conditionals=True)
        self.Z_integer = total * loop_factor
        self.beta = math.log(a/d)/2
        self.eta = math.log(field_weights[1]/field_weights[0])/2

    def sample(self, count, seed=143):
        rng = random.Random(seed)
        out = np.zeros((count, self.n), dtype=np.uint8)
        for row in out:
            for v, rest, table in reversed(self.records):
                weights = table[tuple(int(row[w]) for w in rest)]
                w0, w1 = map(int, weights)
                row[v] = rng.randrange(w0+w1) >= w0
        return out

    def probability(self, assignment):
        p = Fraction(1)
        for v, rest, table in reversed(self.records):
            w = table[tuple(assignment[i] for i in rest)]
            p *= Fraction(int(w[assignment[v]]), int(w[0]+w[1]))
        return p

    def weight(self, assignment):
        k = sum((assignment[u] ^ assignment[v]) ^ b for (u,v), b in zip(self.edges, self.preferred))
        m = sum(assignment)
        return self.a**(len(self.edges)-k)*self.d**k*self.field_weights[0]**(self.n-m)*self.field_weights[1]**m


def signed_minimum(n, edges, preferred, reference=None):
    """Exact lexicographic minimum (frustrated edges, Hamming distance to reference)."""
    scale = n + 1
    ref = [0]*n if reference is None else reference
    factors = [((v,), np.array([int(ref[v] != 0), int(ref[v] != 1)], dtype=np.int64)) for v in range(n)]
    loop_cost = 0
    for (u,v), b in zip(edges, preferred):
        if u == v:
            loop_cost += b
        else:
            factors.append((tuple(sorted((u,v))), scale*np.array([[b,1-b],[1-b,b]], dtype=np.int64)))
    order, widths = min_fill_order(n, edges)
    records = []
    for v in order:
        selected = [(s,a) for s,a in factors if v in s]
        factors = [(s,a) for s,a in factors if v not in s]
        scope = tuple(sorted(set().union(*(s for s,_ in selected))))
        tab = np.zeros((2,)*len(scope), dtype=np.int64)
        for s,a in selected:
            tab += a.reshape(tuple(2 if w in s else 1 for w in scope))
        rest = tuple(w for w in scope if w != v)
        axis = scope.index(v)
        records.append((v, rest, tab.argmin(axis=axis)))
        factors.append((rest, tab.min(axis=axis)))
    cost = sum(int(a) for _,a in factors) + scale*loop_cost
    spins = [0]*n
    for v,rest,choice in reversed(records):
        spins[v] = int(choice[tuple(spins[w] for w in rest)])
    k = sum((spins[u]^spins[v])^b for (u,v),b in zip(edges,preferred))
    distance = sum(a != b for a,b in zip(spins,ref))
    assert cost == scale*k+distance
    return {"min_frustrated": k, "minimum_hamming_distance": distance, "spin_bits": spins}


def bilayer_counts(n, edges):
    """Exact B[k,j]: total intra-layer disagreements k, inter-layer disagreements j."""
    E = len(edges)
    bits, stride = 2*n+1, 2*E+1
    B = 1 << bits
    states = [(v & 1, (v >> 1) & 1) for v in range(4)]
    factors = [((v,), np.array([1 << (bits*stride*(s^t)) for s,t in states], dtype=object)) for v in range(n)]
    pair = np.array([[B**((s^u)+(t^v)) for u,v in states] for s,t in states], dtype=object)
    for u,v in edges:
        if u != v:
            factors.append((tuple(sorted((u,v))), pair))
    packed, cert, _ = eliminate_sum(n, edges, factors, q=4)
    mask = B-1
    out = [[(packed >> (bits*(k+stride*j))) & mask for j in range(n+1)] for k in range(2*E+1)]
    assert packed >> (bits*stride*(n+1)) == 0
    assert sum(map(sum,out)) == 1 << (2*n)
    return out, cert


def joint_observables(counts, n, E, beta, eta=0., higher=False):
    """Numerical cumulants of exact counts, in natural coordinates (beta, eta)."""
    D = np.asarray(counts, dtype=float)
    k,m = np.nonzero(D)
    S, M = E-2*k, 2*m-n
    logw = np.log(D[k,m]) + beta*S + eta*M
    shift = float(logw.max())
    weights = np.exp(logw-shift)
    p = weights/weights.sum()
    stats = np.column_stack((S,M)).astype(float)
    mean = p @ stats
    centered = stats - mean
    g = np.einsum('i,ij,ik->jk',p,centered,centered)
    determinant = float(np.linalg.det(g))
    out = {"beta": float(beta), "eta": float(eta), "log_Z": shift+math.log(float(weights.sum())),
           "magnetization_per_spin": float(mean[1]/n),
           "absolute_magnetization_per_spin": float(p @ np.abs(M)/n),
           "interaction_energy_per_spin": float(-mean[0]/n),
           "susceptibility_per_spin": float(beta*g[1,1]/n),
           "heat_capacity_zero_field_convention": float(beta*beta*g[0,0]/n),
           "fisher_metric": g.tolist(), "fisher_determinant": determinant}
    if higher:
        c3 = np.einsum('i,ij,ik,il->jkl',p,centered,centered,centered)
        # Gaussian curvature convention: categorical simplex has K=+1/4.
        mat = np.array([[g[0,0],g[0,1],g[1,1]],
                        [c3[0,0,0],c3[0,0,1],c3[0,1,1]],
                        [c3[0,0,1],c3[0,1,1],c3[1,1,1]]])
        out["gaussian_curvature"] = float(-np.linalg.det(mat)/(4*determinant**2)) if determinant > 1e-14 else None
        out["third_cumulants"] = c3.tolist()
    return out
