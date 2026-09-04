"""mtft.surface.ising — the Ising model on the dual Manin graph: a first QFT on X0(N).

Spins live on the Manin triangles (vertices of the trivalent dual graph G*),
couplings on the 84 edges; the four cusps are the faces of G*.  With the
Fisher triangle decoration each vertex of G* becomes a triangle, and the
high-temperature expansion Z = 2^F cosh(beta)^E sum_P t^{|P|} over even
subgraphs P becomes a dimer sum on the decorated graph G^F (dimer
configurations <-> complements of even subgraphs).  On a closed surface of
genus g the dimer sum is the Cimasoni-Reshetikhin formula

    Z_dimer = 2^-g  sum_{eps in (Z/2)^{2g}} (-1)^{Arf(q_eps)} Pf_norm(K_eps),

one Pfaffian per spin structure, with K_eps the 4^g Kasteleyn orientations
(a base orientation solving the clockwise-odd face conditions mod 2, twisted
by representatives of H^1(Sigma; Z/2)) and q_eps the quadratic form whose
values on a lifted dual cycle C~ are n_K(C~) + 1 (Kasteleyn's lemma), Arf
computed on a mod-2 symplectic basis.  Pf_norm divides by the sign of the
reference matching (all external edges) so that q_eps(0) = 0.

Gates (EXACT, two routes): at levels with few faces the Pfaffian sum is
compared with brute-force enumeration of all 2^F spin configurations
(N = 6, 11, 15, 35: genus 0, 1, 1, 3).  Structural gates at every level:
Kasteleyn system solvable, H^1 representative count 2g, cocycle/cycle
pairing nondegenerate (so eps -> q_eps is a bijection onto the 4^g forms),
Fourier orthogonality of the Arf signs.

At N = 143 the full sum has 4^13 = 67,108,864 Pfaffians of size 168 (about a
CPU-day); ``partition_function`` refuses it unless ``allow_long=True`` and
``full_sum_job`` checkpoints.  ``sample`` returns per-spin-structure
Pfaffians for a random sample, split by parity; in the continuum limit the
odd-parity Pfaffians vanish (theta-nulls of odd characteristics are zero),
which is the first quantitative fermionic prediction to test under
refinement (DIAGNOSTIC at r = 0).

Critical coupling for the honeycomb-like refinement family: tanh(beta_c) =
1/sqrt(3).  Claim class of numbers: EXACT (sums of exactly computed
Pfaffians at double precision; brute-force cross-check at small genus).
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Dict, List, Sequence, Tuple

import numpy as np

from .cycles import CycleBasis, tree_cotree
from .hodge import unweighted_hodge
from .manin import ManinComplex, cell_complex

T_CRITICAL_HONEYCOMB = 1 / math.sqrt(3)
BETA_CRITICAL_HONEYCOMB = math.atanh(T_CRITICAL_HONEYCOMB)


# ----------------------------------------------------------- GF(2) algebra
def gf2_rref(M: np.ndarray) -> Tuple[np.ndarray, List[int]]:
    M = (M % 2).astype(np.uint8).copy()
    rows, cols = M.shape
    r, piv = 0, []
    for c in range(cols):
        if r == rows:
            break
        nz = np.nonzero(M[r:, c])[0]
        if len(nz) == 0:
            continue
        p = r + nz[0]
        if p != r:
            M[[r, p]] = M[[p, r]]
        others = np.nonzero(M[:, c])[0]
        others = others[others != r]
        M[others] ^= M[r]
        piv.append(c)
        r += 1
    return M[:r], piv


def gf2_solve(A: np.ndarray, b: np.ndarray) -> np.ndarray | None:
    aug = np.hstack([A % 2, (b % 2)[:, None]]).astype(np.uint8)
    R, piv = gf2_rref(aug)
    n = A.shape[1]
    if n in piv:
        return None
    x = np.zeros(n, dtype=np.uint8)
    for i, c in enumerate(piv):
        x[c] = R[i, -1]
    return x


def gf2_nullspace(A: np.ndarray) -> np.ndarray:
    R, piv = gf2_rref(A)
    n = A.shape[1]
    free = [c for c in range(n) if c not in piv]
    out = []
    for f in free:
        x = np.zeros(n, dtype=np.uint8)
        x[f] = 1
        for i, c in enumerate(piv):
            x[c] = R[i, f]
        out.append(x)
    return np.array(out, dtype=np.uint8).reshape(len(out), n)


def gf2_complement_reps(vectors: np.ndarray, subspace: np.ndarray) -> np.ndarray:
    """Independent representatives of span(vectors) modulo row space of subspace."""
    S, piv = gf2_rref(subspace)
    red = []
    for v in vectors:
        v = v.copy()
        for i, c in enumerate(piv):
            if v[c]:
                v ^= S[i]
        red.append(v)
    R, _ = gf2_rref(np.array(red, dtype=np.uint8))
    return R


# ---------------------------------------------------------------- Pfaffian
def pfaffian(A: np.ndarray) -> float:
    """Pfaffian of a real antisymmetric matrix (Parlett-Reid elimination, O(n^3))."""
    A = np.array(A, dtype=float, copy=True)
    n = A.shape[0]
    if n % 2:
        return 0.0
    pf = 1.0
    for k in range(0, n - 1, 2):
        kp = k + 1 + int(np.argmax(np.abs(A[k, k + 1:])))
        if kp != k + 1:
            A[[k + 1, kp], :] = A[[kp, k + 1], :]
            A[:, [k + 1, kp]] = A[:, [kp, k + 1]]
            pf = -pf
        if A[k, k + 1] == 0.0:
            return 0.0
        pf *= A[k, k + 1]
        if k + 2 < n:
            tau = A[k, k + 2:] / A[k, k + 1]
            col = A[k + 2:, k + 1].copy()
            A[k + 2:, k + 2:] += np.outer(tau, col) - np.outer(col, tau)
    return float(pf)


# --------------------------------------------------------- decorated graph
@dataclass
class FisherGraph:
    N: int
    cx: ManinComplex
    n_vertices: int
    edges: List[Tuple[int, int]]          # (u, v) reference direction u -> v
    kind: List[str]                       # 'ext' or 'gad'
    faces: List[List[Tuple[int, int]]]    # ccw boundaries: (edge, +1 if traversed u->v)
    ext_of_manin_edge: List[int]          # Manin edge -> decorated edge index
    vid: Dict[Tuple[int, Tuple[int, int]], int]


def fisher_graph(cx: ManinComplex) -> FisherGraph:
    vid = {(f, d): 3 * f + i for f, o in enumerate(cx.faces) for i, d in enumerate(o)}
    edges, kind = [], []
    ext_of = []
    for e, d in enumerate(cx.edges):
        u, v = vid[(cx.face_of[d], d)], vid[(cx.face_of[cx.S(d)], cx.S(d))]
        ext_of.append(len(edges))
        edges.append((u, v))
        kind.append("ext")
    gad_index: Dict[Tuple[int, int], int] = {}
    for f, o in enumerate(cx.faces):
        for i in range(3):
            gad_index[(f, i)] = len(edges)
            edges.append((vid[(f, o[i])], vid[(f, o[(i + 1) % 3])]))
            kind.append("gad")
    faces: List[List[Tuple[int, int]]] = []
    for f in range(len(cx.faces)):                       # gadget triangles, R order = ccw
        faces.append([(gad_index[(f, i)], 1) for i in range(3)])
    local = {(f, d): i for f, o in enumerate(cx.faces) for i, d in enumerate(o)}
    def gadget_between(f: int, d_from, d_to) -> Tuple[int, int]:
        u, v = vid[(f, d_from)], vid[(f, d_to)]
        for i in range(3):
            a, b = edges[gad_index[(f, i)]]
            if (a, b) == (u, v):
                return gad_index[(f, i)], 1
            if (b, a) == (u, v):
                return gad_index[(f, i)], -1
        raise ValueError("darts not in the same gadget")

    cusp_faces = []
    for orbit in cx.vertices:                            # cusp faces: edge(y_k) separates f_{k-1} and f_k
        bd = []
        w = len(orbit)
        for k, y in enumerate(orbit):
            f = cx.face_of[y]
            y_next = orbit[(k + 1) % w]
            depart = cx.S(y_next)                        # dart of edge(y_next) lying in f
            if cx.face_of[depart] != f:
                raise ValueError("cusp walk not closed at a gadget")
            bd.append(gadget_between(f, y, depart))
            e = cx.edge_of[y_next][0]
            bd.append((ext_of[e], 1 if cx.edges[e] == depart else -1))
        cusp_faces.append(bd)
    # orientation sense: gadget edges must be traversed oppositely by triangle and cusp face
    same = sum(1 for bd in cusp_faces for e, sgn in bd if kind[e] == "gad" and sgn == 1)
    if same:                                             # cusp walks run clockwise: reverse them
        cusp_faces = [[(e, -sgn) for e, sgn in reversed(bd)] for bd in cusp_faces]
    faces.extend(cusp_faces)
    return FisherGraph(cx.N, cx, 3 * len(cx.faces), edges, kind, faces, ext_of, vid)


# ------------------------------------------------------- Kasteleyn setup
@dataclass
class KasteleynData:
    graph: FisherGraph
    base_bits: np.ndarray                 # 0: reference direction, 1: reversed
    twists: np.ndarray                    # 2g x E_F mod-2 cocycle representatives
    dual_cycles: List[List[Tuple[int, int]]]   # lifted cycles: (decorated edge, traversal dir)
    qvec0: np.ndarray                     # q_0 on dual cycles
    pairing: np.ndarray                   # 2g x 2g: |twist_i ∩ cycle_j| mod 2
    bilinear: np.ndarray                  # 2g x 2g intersection form mod 2 on dual cycles
    symp_a: np.ndarray                    # g x 2g coefficient vectors
    symp_b: np.ndarray
    ref_sign0: float
    ref_twist_parity: np.ndarray          # 2g: |twist_i ∩ external edges| mod 2
    gates: List[dict]


def _face_edge_system(G: FisherGraph):
    F, E = len(G.faces), len(G.edges)
    A = np.zeros((F, E), dtype=np.uint8)
    b = np.zeros(F, dtype=np.uint8)
    for f, bd in enumerate(G.faces):
        against = 0
        for e, s in bd:
            A[f, e] ^= 1
            if s == -1:
                against += 1
        b[f] = (1 + against) % 2          # sum bits = 1 + #(dir=-1)  <=>  #against ccw odd
    return A, b


def _coboundaries(G: FisherGraph) -> np.ndarray:
    V, E = G.n_vertices, len(G.edges)
    D = np.zeros((V, E), dtype=np.uint8)
    for e, (u, v) in enumerate(G.edges):
        D[u, e] ^= 1
        D[v, e] ^= 1
    return D


def _dual_fundamental_cycles(cx: ManinComplex, cb: CycleBasis) -> List[List[Tuple[int, int, int]]]:
    """Each cycle as a list of (face, entering dart, leaving dart); loops handled."""
    dual_adj: Dict[int, List[Tuple[int, int]]] = {f: [] for f in range(len(cx.faces))}
    for e in cb.dual_tree:
        d = cx.edges[e]
        a, b = cx.face_of[d], cx.face_of[cx.S(d)]
        dual_adj[a].append((b, e))
        dual_adj[b].append((a, e))
    from collections import deque

    def tree_path(src: int, dst: int) -> List[Tuple[int, int]]:
        par = {src: None}
        q = deque([src])
        while q:
            x = q.popleft()
            if x == dst:
                break
            for nb, e in dual_adj[x]:
                if nb not in par:
                    par[nb] = (x, e)
                    q.append(nb)
        path = []
        x = dst
        while x != src:
            px, e = par[x]
            path.append((px, e, x))
            x = px
        return path[::-1]                      # (from_face, edge, to_face)

    out = []
    for lo in cb.leftover_edges:
        d = cx.edges[lo]
        a, b = cx.face_of[d], cx.face_of[cx.S(d)]     # leftover goes a -> b
        seq = [(a, lo, b)] + tree_path(b, a)
        # convert to (face, entering dart, leaving dart)
        cyc = []
        m = len(seq)
        for k in range(m):
            f_prev, e_in, f = seq[(k - 1) % m]
            _, e_out, _ = seq[k]
            rin, rout = cx.edges[e_in], cx.edges[e_out]
            din = rin if cx.face_of[rin] == f else cx.S(rin)
            dout = rout if cx.face_of[rout] == f else cx.S(rout)
            if f_prev == f and e_in == e_out:            # single loop-edge cycle
                din, dout = cx.S(rout), rout
            cyc.append((f, din, dout))
        out.append(cyc)
    return out


def kasteleyn_setup(cx: ManinComplex, cb: CycleBasis, Jint: np.ndarray) -> KasteleynData:
    G = fisher_graph(cx)
    g2 = cb.rank
    A, b = _face_edge_system(G)
    cover = np.zeros(len(G.edges), dtype=int)
    dsum = np.zeros(len(G.edges), dtype=int)
    for bd in G.faces:
        for e, sgn in bd:
            cover[e] += 1
            dsum[e] += sgn
    gates = [{"name": "faces_cover_each_edge_twice_oppositely",
              "status": "PASS" if (np.all(cover == 2) and np.all(dsum == 0)) else "FAIL"}]
    base = gf2_solve(A, b)
    gates.append({"name": "kasteleyn_system_solvable", "status": "PASS" if base is not None else "FAIL"})
    if base is None:
        raise ArithmeticError("no Kasteleyn orientation: face system inconsistent")
    Z1 = gf2_nullspace(A)
    B1 = _coboundaries(G)
    twists = gf2_complement_reps(Z1, B1)
    gates.append({"name": "H1_mod2_rank_is_2g", "status": "PASS" if twists.shape[0] == g2 else "FAIL",
                  "evidence": f"{twists.shape[0]} vs {g2}"})
    # lifted dual cycles
    cycles = _dual_fundamental_cycles(cx, cb)
    local = {(f, d): i for f, o in enumerate(cx.faces) for i, d in enumerate(o)}
    lifted = []
    for cyc in cycles:
        L = []
        for f, din, dout in cyc:
            i, j = local[(f, din)], local[(f, dout)]
            # gadget edge between local i and j: edge (f,i) if j == i+1, else (f,j) traversed backwards
            gad_ij = None
            for k in range(3):
                u = G.vid[(f, cx.faces[f][k])]
                v = G.vid[(f, cx.faces[f][(k + 1) % 3])]
                idx = len(cx.edges) + 3 * f + k
                if (u, v) == (G.vid[(f, din)], G.vid[(f, dout)]):
                    gad_ij = (idx, 1)
                elif (v, u) == (G.vid[(f, din)], G.vid[(f, dout)]):
                    gad_ij = (idx, -1)
            L.append(gad_ij)
            e = cx.edge_of[dout][0]
            L.append((G.ext_of_manin_edge[e], 1 if cx.edges[e] == dout else -1))
        lifted.append(L)

    def n_against(bits, L):
        return sum(1 for e, s in L if (bits[e] == 0) != (s == 1)) % 2

    qvec0 = np.array([(n_against(base, L) + 1) % 2 for L in lifted], dtype=np.uint8)
    pairing = np.array([[sum(t[e] for e, _ in L) % 2 for L in lifted] for t in twists], dtype=np.uint8).reshape(len(twists), len(lifted)).T
    piv = gf2_rref(pairing)[1] if g2 else []
    gates.append({"name": "cocycle_cycle_pairing_nondegenerate", "status": "PASS" if len(piv) == g2 else "FAIL"})
    # classes of dual cycles in the primal basis and the mod-2 intersection form
    Bm = np.abs(cb.basis_matrix) % 2
    J2 = (np.abs(Jint) % 2).astype(np.uint8)
    X = []
    for cyc in cycles:
        edges_in = [cx.edge_of[dout][0] for _, _, dout in cyc]
        counts = np.array([sum(Bm[e, k] for e in edges_in) % 2 for k in range(g2)], dtype=np.uint8)
        x = gf2_solve(J2, counts)
        X.append(x)
    X = np.array(X, dtype=np.uint8).reshape(len(X), g2)
    bil = (X.astype(int) @ J2.astype(int) @ X.T.astype(int)) % 2
    gates.append({"name": "dual_cycle_intersection_form_nondegenerate", "status": "PASS" if (g2 == 0 or len(gf2_rref(bil)[1]) == g2) else "FAIL"})
    sa, sb = _symplectic_basis(bil.astype(np.uint8)) if g2 else (np.zeros((0, 0), np.uint8), np.zeros((0, 0), np.uint8))
    # reference matching sign (all external edges), base orientation
    ext = [e for e, k in enumerate(G.kind) if k == "ext"]
    ref_sign0 = _matching_sign(G, base, ext)
    ref_parity = np.array([sum(t[e] for e in ext) % 2 for t in twists], dtype=np.uint8)
    return KasteleynData(G, base, twists, lifted, qvec0, pairing, bil.astype(np.uint8), sa, sb, ref_sign0, ref_parity, gates)


def _symplectic_basis(B: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    n = B.shape[0]
    vecs = [np.eye(n, dtype=np.uint8)[i] for i in range(n)]
    A_, B_ = [], []
    form = lambda x, y: int(x.astype(int) @ B.astype(int) @ y.astype(int)) % 2
    while vecs:
        a = vecs.pop(0)
        partner = next((i for i, v in enumerate(vecs) if form(a, v)), None)
        if partner is None:
            continue
        bvec = vecs.pop(partner)
        A_.append(a)
        B_.append(bvec)
        new = []
        for v in vecs:
            v = (v + form(v, bvec) * a + form(v, a) * bvec) % 2
            if v.any():
                new.append(v.astype(np.uint8))
        vecs = new
    return np.array(A_, dtype=np.uint8), np.array(B_, dtype=np.uint8)


def _matching_sign(G: FisherGraph, bits: np.ndarray, matching_edges: Sequence[int]) -> float:
    pairs = []
    prod = 1.0
    for e in matching_edges:
        u, v = G.edges[e]
        if bits[e]:
            u, v = v, u                    # oriented u -> v
        if u < v:
            pairs.append((u, v))
        else:
            pairs.append((v, u))
            prod = -prod
    pairs.sort()
    perm = [x for p in pairs for x in p]
    # parity of permutation
    seen = [False] * len(perm)
    pos = {val: i for i, val in enumerate(perm)}
    sign = 1
    order = sorted(perm)
    for i in range(len(perm)):
        if seen[i]:
            continue
        j, length = i, 0
        while not seen[j]:
            seen[j] = True
            j = pos[order[j]]
            length += 1
        if length % 2 == 0:
            sign = -sign
    return sign * prod


# --------------------------------------------------------- evaluation
class IsingSurface:
    def __init__(self, N: int):
        self.N = N
        self.cx = cell_complex(N)
        self.cb = tree_cotree(self.cx)
        self.Jint = (unweighted_hodge(self.cx, self.cb)["intersection_cycles"] if self.cb.rank
                     else np.zeros((0, 0), dtype=np.int64))
        self.K = kasteleyn_setup(self.cx, self.cb, self.Jint)
        self.genus = self.cx.inv.genus
        self.F, self.E = len(self.cx.faces), len(self.cx.edges)

    # quadratic form and Arf for a twist vector eps
    def q_values(self, eps: np.ndarray) -> np.ndarray:
        return (self.K.qvec0 + (self.K.pairing.astype(int) @ eps.astype(int)) % 2) % 2

    def _q_of_vector(self, qvec: np.ndarray, v: np.ndarray) -> int:
        lin = int(v.astype(int) @ qvec.astype(int)) % 2
        B = self.K.bilinear.astype(int)
        quad = 0
        idx = np.nonzero(v)[0]
        for a in range(len(idx)):
            for b in range(a + 1, len(idx)):
                quad += B[idx[a], idx[b]]
        return (lin + quad) % 2

    def arf(self, eps: np.ndarray) -> int:
        qv = self.q_values(eps)
        return sum(self._q_of_vector(qv, a) * self._q_of_vector(qv, b) for a, b in zip(self.K.symp_a, self.K.symp_b)) % 2

    def matrix(self, eps: np.ndarray, t: float) -> np.ndarray:
        G = self.K.graph
        bits = (self.K.base_bits.astype(int) + self.K.twists.astype(int).T @ eps.astype(int)) % 2
        n = G.n_vertices
        A = np.zeros((n, n))
        for e, (u, v) in enumerate(G.edges):
            w = 1.0 if G.kind[e] == "ext" else t
            if bits[e]:
                u, v = v, u
            A[u, v] += w
            A[v, u] -= w
        return A

    def pfaffian_normalized(self, eps: np.ndarray, t: float) -> float:
        pf = pfaffian(self.matrix(eps, t))
        ref = self.K.ref_sign0 * (-1) ** int(self.K.ref_twist_parity.astype(int) @ eps.astype(int) % 2)
        return pf / ref

    def all_eps(self):
        g2 = 2 * self.genus
        for k in range(1 << g2):
            yield np.array([(k >> i) & 1 for i in range(g2)], dtype=np.uint8)

    def dimer_sum(self, t: float, allow_long: bool = False) -> float:
        if self.genus > 7 and not allow_long:
            raise RuntimeError(f"4^{self.genus} Pfaffians; pass allow_long=True or use full_sum_job")
        total = 0.0
        for eps in self.all_eps():
            total += (-1) ** self.arf(eps) * self.pfaffian_normalized(eps, t)
        return total / 2 ** self.genus

    def partition_function(self, beta: float, allow_long: bool = False) -> float:
        t = math.tanh(beta)
        return 2 ** self.F * math.cosh(beta) ** self.E * self.dimer_sum(t, allow_long)

    def brute_force(self, beta: float, chunk: int = 1 << 20) -> float:
        if self.F > 24:
            raise RuntimeError("brute force limited to F <= 24 faces")
        us = np.array([self.cx.face_of[d] for d in self.cx.edges])
        vs = np.array([self.cx.face_of[self.cx.S(d)] for d in self.cx.edges])
        total = 0.0
        n = 1 << self.F
        for start in range(0, n, chunk):
            cfg = np.arange(start, min(n, start + chunk), dtype=np.int64)
            agree = np.zeros(len(cfg))
            for u, v in zip(us, vs):
                agree += 1 - 2 * (((cfg >> u) ^ (cfg >> v)) & 1)
            total += float(np.sum(np.exp(beta * agree)))
        return total

    def fourier_gate(self, trials: int = 3, seed: int = 0) -> dict:
        """sum_eps (-1)^{Arf(q_eps) + q_eps(h)} = 2^g for every h (small genus only): the identity
        that makes the Cimasoni-Reshetikhin sum count each homology class exactly once."""
        rng = np.random.default_rng(seed)
        g2 = 2 * self.genus
        eps_list = list(self.all_eps())
        arfs = np.array([self.arf(e) for e in eps_list])
        results = {}
        for trial in range(trials + 1):
            h = np.zeros(g2, dtype=np.uint8) if trial == 0 else rng.integers(0, 2, g2).astype(np.uint8)
            if trial and not h.any():
                continue
            s = sum((-1) ** (int(arfs[i]) + self._q_of_vector(self.q_values(e), h)) for i, e in enumerate(eps_list))
            results[str(h.tolist())] = int(s)
        ok = all(v == 2 ** self.genus for v in results.values())
        return {"status": "PASS" if ok else "FAIL", "values": results}

    def sample(self, t: float, n_samples: int = 32, seed: int = 0) -> dict:
        rng = np.random.default_rng(seed)
        even, odd = [], []
        for _ in range(n_samples):
            eps = rng.integers(0, 2, 2 * self.genus).astype(np.uint8)
            pf = self.pfaffian_normalized(eps, t)
            (odd if self.arf(eps) else even).append(pf)
        return {"class": "DIAGNOSTIC", "t": t, "even": even, "odd": odd,
                "median_abs_even": float(np.median(np.abs(even))) if even else None,
                "median_abs_odd": float(np.median(np.abs(odd))) if odd else None}


def gate_report(N: int, beta: float = BETA_CRITICAL_HONEYCOMB) -> dict:
    """Two-route gate: spin-structure Pfaffian sum vs brute force (small genus)."""
    S = IsingSurface(N)
    Zpf = S.partition_function(beta)
    Zbf = S.brute_force(beta)
    rel = abs(Zpf - Zbf) / Zbf
    fg = S.fourier_gate()
    return {"N": N, "genus": S.genus, "faces": S.F, "beta": beta, "Z_pfaffian_sum": Zpf, "Z_brute_force": Zbf,
            "relative_difference": rel, "fourier_gate": fg["status"],
            "structural_gates": S.K.gates,
            "status": "PASS" if rel < 1e-10 and fg["status"] == "PASS" and all(g["status"] == "PASS" for g in S.K.gates) else "FAIL"}


def full_sum_job(N: int, t: float, checkpoint: str, chunk: int = 1 << 16) -> float:
    """Resumable 4^g spin-structure sum (CPU-day at N=143). Appends partial sums to ``checkpoint``."""
    import json, os
    S = IsingSurface(N)
    g2 = 2 * S.genus
    total, done = 0.0, 0
    if os.path.exists(checkpoint):
        with open(checkpoint) as fh:
            for line in fh:
                rec = json.loads(line)
                total, done = rec["total"], rec["done"]
    while done < (1 << g2):
        part = 0.0
        for k in range(done, min(done + chunk, 1 << g2)):
            eps = np.array([(k >> i) & 1 for i in range(g2)], dtype=np.uint8)
            part += (-1) ** S.arf(eps) * S.pfaffian_normalized(eps, t)
        total += part
        done = min(done + chunk, 1 << g2)
        with open(checkpoint, "a") as fh:
            fh.write(json.dumps({"total": total, "done": done}) + "\n")
    return total / 2 ** S.genus
