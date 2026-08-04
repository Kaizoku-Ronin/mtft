#!/usr/bin/env python3
"""
x0143_particle_box.py — MTFT particle box, stage 1 (Python prototype)
======================================================================
Roger Tano / MTFT Research Program — built with Claude, August 2026

Pipeline (independent route; mtft package used ONLY as final oracle):

  Stage A  Farey ideal triangulation of X0(N), N squarefree, torsion-free case.
           Combinatorics from P^1(Z/N): triangles = tau-orbits (168/3 = 56),
           edges = sigma-orbits (168/2 = 84), vertices = cusps (4).
           Certificates: orbit sizes, cusp widths, Euler characteristic,
           genus by a route that never touches the dimension formula.

  Stage B  Dual trivalent graph; unitary wavepacket evolution by the graph
           Laplacian (CURVE) with a Cayley / Crank-Nicolson step.
           Certificate: norm conservation to machine precision.

  Stage C  Modular symbols via Manin symbols on the SAME coset space.
           Relations x + x.sigma = 0, x + x.tau + x.tau^2 = 0.
           Boundary map -> cuspidal subspace (dim 2g, third genus route).
           Hecke T_p by the p+1 matrices + Manin's continued-fraction trick.
           Certificate: char poly of T_2 on the cuspidal space factors as
           (old)^2 x [f1]^2 x [f2]^2 x [f3]^2 and matches the mtft/LMFDB
           oracle exactly (E2: two independent routes).

  Stage D  Orbit blocks (generalized eigenspaces of T_2) = the generations.
           Spectral signature of any drawn cycle = weight vector over
           [old, f1(dim1), f2(dim4), f3(dim6)].  DIAGNOSTIC normalization:
           Euclidean in Manin coordinates (Petersson upgrade is staged).

Epistemic classes on every reported number: EXACT / Cert / DIAGNOSTIC.
"""

from __future__ import annotations
from fractions import Fraction
from math import gcd
import json
import sys

import numpy as np
import os as _os
_HERE = _os.path.dirname(_os.path.abspath(__file__))


# ----------------------------------------------------------------------
# Stage A — P^1(Z/N) and the Farey tessellation of X0(N)
# ----------------------------------------------------------------------

SIGMA = ((0, -1), (1, 0))    # order 2 in PSL: z -> -1/z, reverses {0,oo}
TAU   = ((0, -1), (1, -1))   # order 3 in PSL: cycles 0 -> 1 -> oo -> 0


class P1:
    """The projective line P^1(Z/N) with the right SL2(Z) action on cosets
    Gamma_0(N)\\SL2(Z):  (c:d) . g = (c*g00 + d*g10 : c*g01 + d*g11)."""

    def __init__(self, N: int):
        self.N = N
        units = [u for u in range(1, N) if gcd(u, N) == 1]
        canon = {}
        reps = []
        for c in range(N):
            for d in range(N):
                if gcd(gcd(c, d), N) != 1:
                    continue
                key = (c, d)
                if key in canon:
                    continue
                orbit = sorted(((u * c) % N, (u * d) % N) for u in units)
                rep = orbit[0]
                for p in orbit:
                    canon[p] = rep
                if rep not in reps:
                    reps.append(rep)
        self.canon = canon
        self.reps = sorted(set(reps))
        self.index = {x: i for i, x in enumerate(self.reps)}

    def act(self, x, g):
        (c, d), ((a, b), (cc, dd)) = x, g
        return self.canon[((c * a + d * cc) % self.N, (c * b + d * dd) % self.N)]

    def lift(self, x):
        """Lift (c:d) to a matrix [[a,b],[c0,d0]] in SL2(Z)."""
        c, d = x
        if c % self.N == 0:
            c0, d0 = 0, 1
        else:
            c0 = c
            d0 = d
            k = 0
            while gcd(c0, d0) != 1:
                k += 1
                d0 = d + k * self.N
                if k > 10 * self.N:
                    raise RuntimeError(f"lift failed for {x}")
        # a*d0 - b*c0 = 1
        g, a, b = extended_gcd(d0, -c0)
        assert g == 1, (x, c0, d0)
        return ((a, b), (c0, d0))


def extended_gcd(p, q):
    """returns (g, x, y) with x*p + y*q = g"""
    if q == 0:
        return (abs(p), 1 if p >= 0 else -1, 0)
    g, x1, y1 = extended_gcd(q, p % q)
    return (g, y1, x1 - (p // q) * y1)


def orbits(p1: P1, g, order):
    seen, out = set(), []
    for x in p1.reps:
        if x in seen:
            continue
        orb, y = [x], x
        seen.add(x)
        for _ in range(order - 1):
            y = p1.act(y, g)
            if y not in seen:
                orb.append(y)
                seen.add(y)
        out.append(tuple(orb))
    return out


def cusp_divisor(N, num, den):
    """Cusp class of num/den (lowest terms; den=0 means oo) for squarefree N:
    classes correspond to gcd(den, N) | N."""
    return gcd(den % N if den != 0 else 0, N) if den != 0 else N


def tessellation(N):
    """Build the descended Farey ideal triangulation of X0(N) (N squarefree,
    Gamma_0(N) torsion-free) and certify its combinatorics."""
    p1 = P1(N)
    n_cosets = len(p1.reps)

    tris = orbits(p1, TAU, 3)
    edges = orbits(p1, SIGMA, 2)
    T_orb = orbits(p1, ((1, 1), (0, 1)), n_cosets)   # parabolic T, big order cap

    cert = {}
    cert["index (EXACT)"] = n_cosets
    cert["triangles (EXACT)"] = len(tris)
    cert["tau orbit sizes all 3 (EXACT)"] = all(len(t) == 3 for t in tris)
    cert["edges (EXACT)"] = len(edges)
    cert["sigma orbit sizes all 2 (EXACT)"] = all(len(e) == 2 for e in edges)
    cert["cusps (EXACT)"] = len(T_orb)
    widths = sorted(len(o) for o in T_orb)
    cert["cusp widths (EXACT)"] = widths
    chi = len(T_orb) - len(edges) + len(tris)
    cert["Euler characteristic V-E+F (EXACT)"] = chi
    g = (2 - chi) // 2
    cert["genus from tessellation (EXACT)"] = g
    return p1, tris, edges, cert


# ----------------------------------------------------------------------
# Stage B — dual trivalent graph and unitary wavepacket dynamics
# ----------------------------------------------------------------------

def dual_graph(p1, tris, edges):
    """Vertices = triangles; each unoriented primal edge {x, x.sigma}
    joins triangle(x) and triangle(x.sigma)."""
    tri_of = {}
    for i, t in enumerate(tris):
        for x in t:
            tri_of[x] = i
    n = len(tris)
    A = np.zeros((n, n))
    edge_list = []           # (tri_i, tri_j, manin_rep_x) per unoriented edge
    for e in edges:
        x = e[0]
        xs = p1.act(x, SIGMA)
        i, j = tri_of[x], tri_of[xs]
        A[i, j] += 1
        A[j, i] += 1
        edge_list.append((i, j, x))
    return A, edge_list, tri_of


def cayley_evolution(A, steps=2000, dt=0.05, start=0):
    """Crank-Nicolson (Cayley) unitary step for i dpsi/dt = L psi, L = D - A."""
    n = A.shape[0]
    L = np.diag(A.sum(axis=1)) - A
    I = np.eye(n)
    M = np.linalg.solve(I + 0.5j * dt * L, I - 0.5j * dt * L)
    psi = np.zeros(n, dtype=complex)
    psi[start] = 1.0
    norms, snaps = [], {}
    snap_at = {0, steps // 3, 2 * steps // 3, steps - 1}
    for s in range(steps):
        if s in snap_at:
            snaps[s] = np.abs(psi) ** 2
        psi = M @ psi
        norms.append(np.linalg.norm(psi))
    drift = max(abs(x - 1.0) for x in norms)
    return snaps, drift, L


# ----------------------------------------------------------------------
# Stage C — modular symbols, boundary, Hecke
# ----------------------------------------------------------------------

class Cusp:
    """A point of P^1(Q) as a reduced integer pair (num, den), den >= 0;
    (1,0) is infinity."""
    __slots__ = ("n", "d")

    def __init__(self, n, d):
        if d < 0:
            n, d = -n, -d
        g = gcd(abs(n), d) if d != 0 else abs(n)
        g = g if g != 0 else 1
        self.n, self.d = n // g, d // g

    def apply(self, m):
        (a, b), (c, d) = m
        return Cusp(a * self.n + b * self.d, c * self.n + d * self.d)

    def __repr__(self):
        return "oo" if self.d == 0 else f"{self.n}/{self.d}"


def convergent_chain(c: Cusp):
    """Manin's trick: the chain {oo -> c} as a list of Manin generators.
    Each segment {p_{k-1}/q_{k-1} -> p_k/q_k} is g_k{0,oo} with
    g_k = [[p_k, s*p_{k-1}], [q_k, s*q_{k-1}]], s = (-1)^{k-1}, det = +1.
    Returns list of (bottom_row) integer pairs (before reduction mod N)."""
    if c.d == 0:
        return []
    # continued fraction of n/d
    n, d = c.n, c.d
    quots = []
    while d != 0:
        a = n // d           # floor division (works for negatives)
        quots.append(a)
        n, d = d, n - a * d
    # convergents
    p_prev, q_prev = 1, 0    # p_{-1}/q_{-1} = oo
    p_cur, q_cur = quots[0], 1
    segs = [(q_cur, (1) * q_prev)]  # k = 0: s = (-1)^{-1} = -1? fix below
    # Recompute carefully with explicit signs and a det assertion:
    segs = []
    p_m1, q_m1 = 1, 0
    p_0, q_0 = quots[0], 1
    convs = [(p_m1, q_m1), (p_0, q_0)]
    for a in quots[1:]:
        p_m1, q_m1, p_0, q_0 = p_0, q_0, a * p_0 + p_m1, a * q_0 + q_m1
        convs.append((p_0, q_0))
    for k in range(1, len(convs)):
        (pk, qk), (pk1, qk1) = convs[k], convs[k - 1]
        s = 1 if (k % 2 == 0) else -1     # sign making det = +1
        det = pk * (s * qk1) - (s * pk1) * qk
        if det != 1:
            s = -s
            det = pk * (s * qk1) - (s * pk1) * qk
        assert det == 1, ("det", det, convs)
        segs.append((qk, s * qk1))        # bottom row (c, d) of g_k
    return segs


def path_chain(alpha: Cusp, beta: Cusp):
    """{alpha -> beta} = {oo -> beta} - {oo -> alpha} as Manin bottom rows."""
    out = []
    for row in convergent_chain(beta):
        out.append((row, +1))
    for row in convergent_chain(alpha):
        out.append((row, -1))
    return out


def hecke_matrices(p):
    """The p+1 matrices for T_p, p prime not dividing N."""
    ms = [((p, 0), (0, 1))]
    for j in range(p):
        ms.append(((1, j), (0, p)))
    return ms


class ModularSymbols:
    def __init__(self, p1: P1):
        self.p1 = p1
        self.n = len(p1.reps)
        self._build_quotient()
        self._build_boundary()

    # -- quotient by Manin relations, exact over Q -----------------------
    def _build_quotient(self):
        import sympy as sp
        n, p1 = self.n, self.p1
        rows = []
        seen = set()
        for x in p1.reps:
            i = p1.index[x]
            j = p1.index[p1.act(x, SIGMA)]
            key = tuple(sorted((i, j)))
            if ("s", key) not in seen:
                seen.add(("s", key))
                r = [0] * n
                r[i] += 1
                r[j] += 1
                rows.append(r)
            a = p1.index[p1.act(x, TAU)]
            b = p1.index[p1.act(p1.act(x, TAU), TAU)]
            key = tuple(sorted((i, a, b)))
            if ("t", key) not in seen:
                seen.add(("t", key))
                r = [0] * n
                r[i] += 1
                r[a] += 1
                r[b] += 1
                rows.append(r)
        R = sp.Matrix(rows)
        Rr, piv = R.rref()
        piv = list(piv)
        free = [j for j in range(n) if j not in piv]
        self.free = free
        self.dim = len(free)
        # projection F -> Q^dim : e_free[j] -> unit; e_piv -> -sum coeffs
        proj = {}
        for j_idx, j in enumerate(free):
            proj[j] = {j_idx: Fraction(1)}
        for r_idx, p_col in enumerate(piv):
            comb = {}
            for f_idx, f_col in enumerate(free):
                v = Rr[r_idx, f_col]
                if v != 0:
                    comb[f_idx] = Fraction(-sp.Rational(v).p, sp.Rational(v).q)
            proj[p_col] = comb
        self.proj = proj

    def project(self, chain):
        """chain: dict coset_index -> Fraction. Returns vector in Q^dim."""
        v = [Fraction(0)] * self.dim
        for j, c in chain.items():
            for k, w in self.proj[j].items():
                v[k] += c * w
        return v

    def lift_to_F(self, vec):
        """Canonical lift of a quotient vector back to the 168 Manin coords
        (free coords as given, pivot coords determined by the relations).
        Used only for drawing edge currents (DIAGNOSTIC)."""
        F = [Fraction(0)] * self.n
        for j_idx, j in enumerate(self.free):
            F[j] = Fraction(vec[j_idx])
        # pivots: e_p = sum over free of proj coeffs -> value
        for col, comb in self.proj.items():
            if col in self.free:
                continue
            F[col] = sum((w * Fraction(vec[k]) for k, w in comb.items()),
                         Fraction(0))
        return F

    # -- Manin symbol of a bottom row -----------------------------------
    def manin_index(self, row):
        c, d = row
        N = self.p1.N
        return self.p1.index[self.p1.canon[(c % N, d % N)]]

    def chain_of_path(self, alpha: Cusp, beta: Cusp):
        ch = {}
        for row, sgn in path_chain(alpha, beta):
            j = self.manin_index(row)
            ch[j] = ch.get(j, Fraction(0)) + sgn
        return ch

    # -- boundary --------------------------------------------------------
    def _build_boundary(self):
        import sympy as sp
        p1, N = self.p1, self.p1.N
        divs = sorted(d for d in range(1, N + 1) if N % d == 0)
        self.divisors = divs
        didx = {d: i for i, d in enumerate(divs)}
        cols = []
        for j_idx, j in enumerate(self.free):
            x = p1.reps[j]
            (a, b), (c0, d0) = p1.lift(x)
            # symbol = path from b/d0 (= gamma.0) to a/c0 (= gamma.oo)
            col = [Fraction(0)] * len(divs)
            col[didx[cusp_divisor(N, a, c0)]] += 1
            col[didx[cusp_divisor(N, b, d0)]] -= 1
            cols.append(col)
        B = sp.Matrix([[sp.Rational(cols[j][i]) for j in range(self.dim)]
                       for i in range(len(divs))])
        self.boundary = B
        ns = B.nullspace()
        self.cuspidal_basis = ns          # list of dim-vectors
        self.cuspidal_dim = len(ns)

    # -- Hecke -----------------------------------------------------------
    def hecke_on_quotient(self, p):
        """Exact matrix of T_p on the quotient (dim x dim)."""
        import sympy as sp
        cols = []
        for j_idx, j in enumerate(self.free):
            x = self.p1.reps[j]
            (a, b), (c0, d0) = self.p1.lift(x)
            alpha, beta = Cusp(b, d0), Cusp(a, c0)
            total = {}
            for m in hecke_matrices(p):
                ch = self.chain_of_path(alpha.apply(m), beta.apply(m))
                for k, c in ch.items():
                    total[k] = total.get(k, Fraction(0)) + c
            cols.append(self.project(total))
        M = sp.Matrix([[sp.Rational(cols[j][i].numerator,
                                    cols[j][i].denominator)
                        for j in range(self.dim)] for i in range(self.dim)])
        return M

    def restrict_to_cuspidal(self, M):
        import sympy as sp
        B = sp.Matrix.hstack(*self.cuspidal_basis)     # dim x 2g
        BtB = B.T * B
        A = BtB.LUsolve(B.T * (M * B))
        # certificate: invariance
        assert sp.simplify(B * A - M * B) == sp.zeros(*B.shape)
        return A, B


# ----------------------------------------------------------------------
# Stage D — orbit blocks and spectral signatures
# ----------------------------------------------------------------------

def orbit_blocks(A, factors):
    """factors: list of (sympy poly in x, name, expected_block_dim).
    Returns list of (name, block basis matrix)."""
    import sympy as sp
    x = sp.symbols("x")
    blocks = []
    for poly, name, expdim in factors:
        P = sp.Poly(poly, x)
        M = sp.zeros(A.rows, A.cols)
        Apow = sp.eye(A.rows)
        coeffs = P.all_coeffs()[::-1]
        for c in coeffs:
            M += sp.Rational(c) * Apow
            Apow = Apow * A
        ker = sp.Matrix.hstack(*M.nullspace()) if M.nullspace() else None
        blocks.append((name, ker, expdim))
    return blocks


def signature(vec_np, block_mats_np):
    """DIAGNOSTIC Euclidean weights of a cuspidal vector across blocks.
    Solve v = sum block components in the least-squares (here exact) sense."""
    Ball = np.hstack(block_mats_np)
    coef, *_ = np.linalg.lstsq(Ball, vec_np, rcond=None)
    out, start = [], 0
    for Bm in block_mats_np:
        k = Bm.shape[1]
        comp = Bm @ coef[start:start + k]
        out.append(float(np.dot(comp, comp)))
        start += k
    total = sum(out)
    return [w / total for w in out] if total > 0 else out


# ----------------------------------------------------------------------
# main
# ----------------------------------------------------------------------

def run_level(N, verbose=True):
    p1, tris, edges, cert = tessellation(N)
    ms = ModularSymbols(p1)
    cert["modular symbols dim = 2g + cusps - 1 (EXACT)"] = ms.dim
    cert["cuspidal dim = 2g (EXACT)"] = ms.cuspidal_dim
    if verbose:
        for k, v in cert.items():
            print(f"  {k}: {v}")
    return p1, tris, edges, ms, cert


def main():
    import sympy as sp
    x = sp.symbols("x")
    ledger = {}

    print("=" * 70)
    print("STAGE 0 — engine self-test at N = 11 (oracle: a2 = -2, a3 = -1)")
    print("=" * 70)
    p1_11, _, _, ms11, c11 = run_level(11)
    T2_11, _ = ms11.restrict_to_cuspidal(ms11.hecke_on_quotient(2))
    cp2_11 = sp.factor(T2_11.charpoly(x).as_expr())
    T3_11, _ = ms11.restrict_to_cuspidal(ms11.hecke_on_quotient(3))
    cp3_11 = sp.factor(T3_11.charpoly(x).as_expr())
    print(f"  charpoly T2 on cuspidal: {cp2_11}   [expect (x+2)^2]")
    print(f"  charpoly T3 on cuspidal: {cp3_11}   [expect (x+1)^2]")
    ok11 = (sp.expand(cp2_11 - (x + 2) ** 2) == 0 and
            sp.expand(cp3_11 - (x + 1) ** 2) == 0)
    ledger["N=11 engine self-test (Cert)"] = bool(ok11)
    print(f"  N=11 self-test: {'PASS' if ok11 else 'FAIL'}")
    if not ok11:
        sys.exit("Engine failed the level-11 oracle. Stop.")

    print()
    print("=" * 70)
    print("STAGE A/C — X0(143): tessellation + modular symbols")
    print("=" * 70)
    p1, tris, edges, ms, certA = run_level(143)
    ledger["tessellation"] = {k: (v if not isinstance(v, list) else v)
                              for k, v in certA.items()}

    print()
    print("STAGE C — Hecke T2 on the 26-dim cuspidal space (exact)")
    T2q = ms.hecke_on_quotient(2)
    A2, Bc = ms.restrict_to_cuspidal(T2q)
    cp2 = A2.charpoly(x).as_expr()
    cp2f = sp.factor(cp2)
    print(f"  charpoly(T2) factors:\n    {cp2f}")

    # ---- oracle comparison (E2: first contact with the mtft package) ----
    import mtft.x0_143 as ox
    q2 = sum(int(round(float(c))) * x ** (4 - i)
             for i, c in enumerate(ox.hecke_polynomial_f2_T2()))
    q3 = sum(int(round(float(c))) * x ** (6 - i)
             for i, c in enumerate(ox.hecke_polynomial_f3_T2()))
    a2_f1 = ox.CURVE_143A1.hecke_eigenvalues[2]          # 0
    a2_old = -2                                          # level-11 curve 11a
    expected = sp.expand(((x - a2_old) ** 4) * ((x - a2_f1) ** 2)
                         * q2 ** 2 * q3 ** 2)
    match = sp.expand(cp2 - expected) == 0
    ledger["charpoly(T2) matches mtft/LMFDB oracle (Cert)"] = bool(match)
    print(f"  ORACLE MATCH (old^2 * f1^2 * f2^2 * f3^2): "
          f"{'PASS' if match else 'FAIL'}")
    if not match:
        sys.exit("Char poly disagrees with the oracle. Stop.")

    # ---- orbit blocks ---------------------------------------------------
    print()
    print("STAGE D — orbit blocks (generalized eigenspaces of T2)")
    blocks = orbit_blocks(A2, [
        (x - a2_old, "old (level 11, flavor-universal)", 4),
        (x - a2_f1, "f1  electron orbit (dim 1)", 2),
        (q2, "f2  muon orbit (dim 4)", 8),
        (q3, "f3  tau orbit (dim 6)", 12),
    ])
    dims_ok = True
    block_np = []
    for name, ker, expdim in blocks:
        d = 0 if ker is None else ker.cols
        print(f"  {name}: block dim {d} (expected {expdim})")
        dims_ok &= (d == expdim)
        block_np.append(np.array(ker.tolist(), dtype=float))
    ledger["block dims [4,2,8,12] (Cert)"] = bool(dims_ok)
    print(f"  block dimension certificate: {'PASS' if dims_ok else 'FAIL'}")

    # ---- signatures -----------------------------------------------------
    print()
    print("STAGE D — spectral signatures (DIAGNOSTIC Euclidean weights)")
    names = ["old", "f1 e", "f2 mu", "f3 tau"]
    demo = {}
    for i, (name, _, _) in enumerate(blocks):
        v = block_np[i][:, 0]
        sig = signature(v, block_np)
        demo[f"eigencycle of {names[i]}"] = [round(s, 6) for s in sig]
    rng = np.random.default_rng(143)
    vrand = sum(rng.integers(-3, 4) * block_np[i][:, k]
                for i in range(4) for k in range(min(2, block_np[i].shape[1])))
    demo["random drawn cycle"] = [round(s, 6)
                                  for s in signature(vrand, block_np)]
    for k, v in demo.items():
        print(f"  {k}:  [old, e, mu, tau] = {v}")
    ledger["signatures (DIAGNOSTIC)"] = demo

    # ---- Stage B: dynamics ---------------------------------------------
    print()
    print("STAGE B — unitary wavepacket on the dual trivalent graph")
    A, edge_list, tri_of = dual_graph(p1, tris, edges)
    deg_ok = np.allclose(A.sum(axis=1), 3)
    print(f"  dual graph 3-regular: {'PASS' if deg_ok else 'FAIL'} "
          f"({A.shape[0]} nodes, {int(A.sum() / 2)} edges)")
    snaps, drift, L = cayley_evolution(A)
    print(f"  norm drift over 2000 Cayley steps: {drift:.3e} (EXACT-class "
          f"unitarity to machine precision)")
    evals = np.linalg.eigvalsh(L)
    print(f"  dual-graph Laplacian spectral gap lambda_1 = {evals[1]:.6f} "
          f"(DIAGNOSTIC)")
    ledger["dual graph"] = {"nodes": int(A.shape[0]),
                            "edges": int(A.sum() / 2),
                            "three_regular (EXACT)": bool(deg_ok),
                            "norm_drift (EXACT-class)": float(drift),
                            "lambda_1 (DIAGNOSTIC)": float(evals[1])}

    # ---- figures --------------------------------------------------------
    print()
    print("FIGURES")
    make_figures(p1, tris, edges, ms, blocks, block_np, A, edge_list,
                 tri_of, snaps, L, demo)

    with open(_os.path.join(_HERE, "certificates.json"), "w") as f:
        json.dump(ledger, f, indent=2, default=str)
    print("\nledger written: certificates.json")
    print("ALL CERTIFICATES PASS" if (ok11 and match and dims_ok and deg_ok)
          else "SOME CERTIFICATES FAILED")


# ----------------------------------------------------------------------
# figures
# ----------------------------------------------------------------------

def layout_positions(L):
    """Spectral layout from Laplacian eigenvectors (CURVE-derived), with a
    short deterministic spring refinement for readability."""
    w, V = np.linalg.eigh(L)
    pos = V[:, 1:3].copy()
    pos /= (np.abs(pos).max() + 1e-12)
    A = np.diag(np.diag(L)) - L
    n = L.shape[0]
    for _ in range(300):
        disp = np.zeros_like(pos)
        diff = pos[:, None, :] - pos[None, :, :]
        dist2 = (diff ** 2).sum(-1) + 1e-6
        rep = diff / dist2[..., None] * 0.002
        disp += rep.sum(axis=1)
        attr = -diff * A[..., None] * 0.02
        disp += attr.sum(axis=1)
        pos += np.clip(disp, -0.03, 0.03)
    pos /= np.abs(pos).max()
    return pos


def make_figures(p1, tris, edges, ms, blocks, block_np, A, edge_list,
                 tri_of, snaps, L, demo):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    pos = layout_positions(L)
    Bc = np.array(np.hstack([b for b in
                             [np.array(k[1].tolist(), dtype=float)
                              for k in blocks]]))
    # map: block vector (26) -> quotient (29) -> lift to F (168) -> currents
    Bcusp = np.array(
        np.hstack([np.array(v.tolist(), dtype=float)
                   for v in ms.cuspidal_basis]))

    def edge_currents(block_vec26):
        q = Bcusp @ block_vec26                    # 29-dim quotient coords
        F = ms.lift_to_F([Fraction(t).limit_denominator(10 ** 9)
                          for t in q])
        cur = []
        for (i, j, xrep) in edge_list:
            a = float(F[p1.index[xrep]])
            b = float(F[p1.index[p1.act(xrep, SIGMA)]])
            cur.append(abs(a - b))
        cur = np.array(cur)
        return cur / (cur.max() + 1e-12)

    # -- figure 1: the three generations drawn on the map -----------------
    fig, axes = plt.subplots(1, 3, figsize=(16.5, 5.6))
    titles = [("f1 — electron orbit (dim 1)", 1, "#e4572e"),
              ("f2 — muon orbit (dim 4)", 2, "#17bebb"),
              ("f3 — tau orbit (dim 6)", 3, "#9b5de5")]
    for ax, (title, bi, color) in zip(axes, titles):
        cur = edge_currents(block_np[bi][:, 0])
        for (i, j, _), w in zip(edge_list, cur):
            ax.plot([pos[i, 0], pos[j, 0]], [pos[i, 1], pos[j, 1]],
                    color=color, alpha=0.12 + 0.88 * w,
                    linewidth=0.4 + 4.6 * w, zorder=1)
        ax.scatter(pos[:, 0], pos[:, 1], s=14, c="#222", zorder=2)
        ax.set_title(title, fontsize=11)
        ax.set_aspect("equal")
        ax.axis("off")
    fig.suptitle("Particles drawn from the map — Hecke eigencycles on the "
                 "56-triangle Farey tessellation of $X_0(143)$ "
                 "(dual graph; edge weight = homology current, DIAGNOSTIC)",
                 fontsize=12)
    fig.tight_layout()
    fig.savefig("fig_particles_on_map.png", dpi=170)
    plt.close(fig)
    print("  wrote fig_particles_on_map.png")

    # -- figure 2: wavepacket --------------------------------------------
    fig, axes = plt.subplots(1, len(snaps), figsize=(16.5, 4.4))
    for ax, (step, prob) in zip(axes, sorted(snaps.items())):
        for (i, j, _) in edge_list:
            ax.plot([pos[i, 0], pos[j, 0]], [pos[i, 1], pos[j, 1]],
                    color="#bbb", linewidth=0.5, zorder=1)
        ax.scatter(pos[:, 0], pos[:, 1], s=3000 * prob + 4,
                   c=prob, cmap="inferno", zorder=2)
        ax.set_title(f"step {step}", fontsize=10)
        ax.set_aspect("equal")
        ax.axis("off")
    fig.suptitle("Unitary wavepacket under the graph Laplacian (CURVE) — "
                 "Cayley step, norm conserved to machine precision; spread = "
                 "purely-a.c.-style propagation (no measurement, no freezing)",
                 fontsize=11)
    fig.tight_layout()
    fig.savefig("fig_wavepacket.png", dpi=170)
    plt.close(fig)
    print("  wrote fig_wavepacket.png")

    # -- figure 3: signatures ---------------------------------------------
    fig, ax = plt.subplots(figsize=(9, 4.6))
    labels = list(demo.keys())
    cats = ["old", "e", "mu", "tau"]
    colors = ["#777", "#e4572e", "#17bebb", "#9b5de5"]
    bottoms = np.zeros(len(labels))
    for ci, cat in enumerate(cats):
        vals = np.array([demo[k][ci] for k in labels])
        ax.bar(range(len(labels)), vals, bottom=bottoms,
               color=colors[ci], label=cat)
        bottoms += vals
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels([l.replace("eigencycle of ", "") for l in labels],
                       rotation=15, ha="right", fontsize=9)
    ax.set_ylabel("spectral weight (DIAGNOSTIC)")
    ax.set_title("Spectral signature of drawn cycles — weight over "
                 "[old, f1, f2, f3] orbit blocks")
    ax.legend(loc="upper right", fontsize=9)
    fig.tight_layout()
    fig.savefig("fig_signatures.png", dpi=170)
    plt.close(fig)
    print("  wrote fig_signatures.png")


if __name__ == "__main__":
    main()
