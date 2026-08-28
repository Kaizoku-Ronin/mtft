"""Theta characteristics / spin structures mod 2 for X0(143).

Frame-agnostic machinery for the action of symplectic (mod 2) operators on
theta characteristics, together with the certified census for the group
G2 = <W11, W13, STAR> acting on the 2^26 characteristics of X0(143)
(Cert EXACT, E2, 2026-08-27: mtft.periods route vs an independent PARI/GP
mslattice route agree on every integer).

Dictionary.  A characteristic is t in F_2^{2g}, split t = (t', t'') in the
standard mod-2 symplectic frame J2 = [[0, I], [I, 0]].  It labels the
quadratic form q_t = q0 + <t, .> with q0(x) = x'.x'', and equivalently the
theta characteristic with half-integer shifts a = t'/2, b = t''/2.  The
parity of theta[a; b] equals Arf(q_t) = t'.t'' + (correction determined by
the frame); in the standard frame used here parity(t) = t'.t'' mod 2, and
the census reproduces the textbook counts 2^{g-1}(2^g +- 1).

A symplectic mod-2 operator M acts on characteristics affinely,
t -> A t + d with A = (M^{-1})^T and d_j = q0(M^{-1} e_j); this is the
unique affine action for which q_t -> q_t o M^{-1}, and it preserves Arf.
The self-check :func:`verify_action` re-derives this on random samples.

Main entry points:
  * ``SpinAction(ops)``        - action machinery for named mod-2 symplectic
                                 operators in the standard J2 frame;
  * ``x0143_periods_frame()``  - the G2 action in the frame of the frozen
                                 Riemann matrix tau0 (mtft.periods);
  * ``x0143_gp_frame()``       - the same action in the independent PARI/GP
                                 integral frame (mtft.homology);
  * ``census()``               - the packaged certified census;
  * ``SpinAction.full_census`` - recompute everything from scratch;
  * ``SpinAction.invariant_characteristics`` - the joint fixed locus
                                 (for X0(143): affine F_2^7, 96 even / 32
                                 odd) as explicit vectors in the acting
                                 frame.

Counts are frame-invariant; the explicit characteristic vectors are not.
Theta-function work at tau0 must use the periods frame.
"""
from __future__ import annotations

import itertools
import json
from importlib import resources

import numpy as np

__all__ = [
    "gf2_rank", "gf2_inv", "gf2_solve_affine", "standard_J2",
    "SpinAction", "census", "x0143_periods_frame", "x0143_gp_frame",
    "verify_action",
]


# ---------------------------------------------------------------- GF(2) core
def gf2_rank(A):
    A = A.copy() % 2
    m, nc = A.shape
    r = 0
    for c in range(nc):
        p = next((i for i in range(r, m) if A[i, c]), None)
        if p is None:
            continue
        A[[r, p]] = A[[p, r]]
        for i in range(m):
            if i != r and A[i, c]:
                A[i] ^= A[r]
        r += 1
    return r


def gf2_inv(A):
    m = A.shape[0]
    M = np.concatenate([A.copy() % 2, np.eye(m, dtype=np.uint8)], axis=1)
    r = 0
    for c in range(m):
        p = next((i for i in range(r, m) if M[i, c]), None)
        if p is None:
            raise ValueError("singular over F_2")
        M[[r, p]] = M[[p, r]]
        for i in range(m):
            if i != r and M[i, c]:
                M[i] ^= M[r]
        r += 1
    return M[:, m:]


def gf2_solve_affine(A, b):
    """Solve A t = b over F_2.  Returns (t0, B): solution set t0 + col span
    of B, or (None, None) if inconsistent."""
    A = A.copy() % 2
    b = b.copy() % 2
    m, nc = A.shape
    Mx = np.concatenate([A, b.reshape(-1, 1)], axis=1)
    piv = []
    r = 0
    for c in range(nc):
        p = next((i for i in range(r, m) if Mx[i, c]), None)
        if p is None:
            continue
        Mx[[r, p]] = Mx[[p, r]]
        for i in range(m):
            if i != r and Mx[i, c]:
                Mx[i] ^= Mx[r]
        piv.append(c)
        r += 1
    for i in range(r, m):
        if not Mx[i, :nc].any() and Mx[i, nc]:
            return None, None
    free = [c for c in range(nc) if c not in piv]
    t0 = np.zeros(nc, dtype=np.uint8)
    for i, c in enumerate(piv):
        t0[c] = Mx[i, nc]
    B = np.zeros((nc, len(free)), dtype=np.uint8)
    for j, f in enumerate(free):
        B[f, j] = 1
        for i, c in enumerate(piv):
            if Mx[i, f]:
                B[c, j] = 1
    return t0, B


def standard_J2(g):
    J = np.zeros((2 * g, 2 * g), dtype=np.uint8)
    J[:g, g:] = np.eye(g, dtype=np.uint8)
    J[g:, :g] = np.eye(g, dtype=np.uint8)
    return J


# ------------------------------------------------------------ action object
class SpinAction:
    """Affine action of named mod-2 symplectic operators on characteristics.

    ``ops``: dict name -> (2g x 2g) uint8 matrix, symplectic for the
    standard J2 (asserted).  All counting is exact.
    """

    def __init__(self, ops: dict, g: int | None = None):
        names = list(ops)
        M0 = np.asarray(ops[names[0]]) % 2
        n = M0.shape[0]
        self.n = n
        self.g = g if g is not None else n // 2
        J2 = standard_J2(self.g)
        self.ops = {}
        for k, M in ops.items():
            M = np.asarray(M, dtype=np.uint8) % 2
            assert np.array_equal((M.T @ J2 @ M) % 2, J2), \
                f"{k}: not symplectic mod 2 in the standard frame"
            self.ops[k] = M

    # -- parity ------------------------------------------------------------
    def parity(self, T):
        """Arf parity of characteristic rows T (0 = even, 1 = odd)."""
        T = np.atleast_2d(np.asarray(T, dtype=np.uint8))
        g = self.g
        return ((T[:, :g] * T[:, g:]).sum(axis=1) % 2).astype(np.uint8)

    # -- affine action -----------------------------------------------------
    def affine(self, M):
        """(A, d) with the action t -> A t + d, A = (M^{-1})^T,
        d_j = q0(M^{-1} e_j)."""
        Mi = gf2_inv(np.asarray(M, dtype=np.uint8))
        A = Mi.T % 2
        g = self.g
        d = np.array([(Mi[:g, j] @ Mi[g:, j]) % 2 for j in range(self.n)],
                     dtype=np.uint8)
        return A, d

    def _group(self):
        """All words in the (commuting, involutive) generators; returns
        dict word-name -> matrix, including the identity ''."""
        names = list(self.ops)
        out = {"": np.eye(self.n, dtype=np.uint8)}
        for r in range(1, len(names) + 1):
            for combo in itertools.combinations(names, r):
                M = np.eye(self.n, dtype=np.uint8)
                for c in combo:
                    M = (M @ self.ops[c]) % 2
                out["*".join(combo)] = M
        return out

    # -- fixed loci --------------------------------------------------------
    def fixed_locus(self, M):
        """(t0, B) for the affine fixed locus of a single operator, or
        (None, None) if empty."""
        A, d = self.affine(M)
        I = np.eye(self.n, dtype=np.uint8)
        return gf2_solve_affine((I + A) % 2, d)

    def joint_fixed_locus(self, subset=None):
        subset = list(self.ops) if subset is None else list(subset)
        AA, dd = [], []
        I = np.eye(self.n, dtype=np.uint8)
        for k in subset:
            A, d = self.affine(self.ops[k])
            AA.append((I + A) % 2)
            dd.append(d)
        return gf2_solve_affine(np.vstack(AA), np.concatenate(dd))

    # -- exact parity counting on affine subspaces -------------------------
    def parity_count(self, t0, B, chunk=20):
        d = B.shape[1]
        g = self.g
        if d == self.n:
            return 2 ** (g - 1) * (2 ** g + 1), 2 ** (g - 1) * (2 ** g - 1)
        even = odd = 0
        lo_bits = min(d, chunk)
        k = np.arange(2 ** lo_bits, dtype=np.uint32)
        coeff = ((k[:, None] >> np.arange(lo_bits)[None, :]) & 1
                 ).astype(np.uint8)
        for hi in range(2 ** max(0, d - chunk)):
            if d > chunk:
                hb = np.array([(hi >> j) & 1 for j in range(d - chunk)],
                              dtype=np.uint8)
                base = (t0 + (B[:, chunk:] @ hb)) % 2
            else:
                base = t0
            pts = (base[None, :] + coeff @ B[:, :lo_bits].T) % 2
            a = self.parity(pts)
            odd += int(a.sum())
            even += int(len(a) - a.sum())
        return even, odd

    def invariant_characteristics(self):
        """The joint fixed locus as explicit vectors (rows), with parities.

        For the X0(143) G2 action: 128 rows, 96 even / 32 odd."""
        t0, B = self.joint_fixed_locus()
        if t0 is None:
            return np.zeros((0, self.n), np.uint8), np.zeros(0, np.uint8)
        d = B.shape[1]
        k = np.arange(2 ** d, dtype=np.uint32)
        coeff = ((k[:, None] >> np.arange(d)[None, :]) & 1).astype(np.uint8)
        pts = (t0[None, :] + coeff @ B.T) % 2
        return pts, self.parity(pts)

    def parity_polarization(self, t0, B):
        """Rank and radical dimension of the parity function restricted to
        the affine locus t0 + span B, as a quadratic form on F_2^d.
        The polarization bilinear form is
        beta(u, v) = P(u + v) - P(u) - P(v) + P(0) with P(x) =
        parity(t0 + B x)."""
        d = B.shape[1]
        P0 = int(self.parity(t0)[0])
        Pe = [int(self.parity((t0 + B[:, i]) % 2)[0]) for i in range(d)]
        beta = np.zeros((d, d), dtype=np.uint8)
        for i in range(d):
            for j in range(d):
                x = (t0 + B[:, i] + B[:, j]) % 2
                beta[i, j] = (int(self.parity(x)[0]) + Pe[i] + Pe[j] + P0) % 2
        rank = gf2_rank(beta)
        return rank, d - rank

    # -- orbit census ------------------------------------------------------
    def full_census(self, chunk=20):
        """Recompute the complete census (per-element fixed counts, joint
        locus, parity polarization, Burnside orbit counts, orbit-size
        census).  Exact; a few minutes for 2g = 26."""
        grp = self._group()
        out = {"element_fixed": {}, "joint": {}, "burnside": {},
               "orbit_sizes": {}}
        fixed = {}
        for name, M in grp.items():
            key = name if name else "I"
            t0, B = self.fixed_locus(M) if name else (
                np.zeros(self.n, np.uint8), np.eye(self.n, dtype=np.uint8))
            if t0 is None:
                out["element_fixed"][key] = {"affine_dim": None,
                                             "even": 0, "odd": 0}
                fixed[name] = (0, 0)
                continue
            ev, od = self.parity_count(t0, B, chunk)
            out["element_fixed"][key] = {"affine_dim": int(B.shape[1]),
                                         "even": ev, "odd": od}
            fixed[name] = (ev, od)

        t0, B = self.joint_fixed_locus()
        ev, od = self.parity_count(t0, B, chunk)
        rank, rad = self.parity_polarization(t0, B)
        out["joint"] = {"affine_dim": int(B.shape[1]), "even": ev,
                        "odd": od, "parity_quadratic_rank": rank,
                        "parity_radical_dim": rad}

        # Burnside
        G = len(grp)
        se = sum(v[0] for v in fixed.values())
        so = sum(v[1] for v in fixed.values())
        assert se % G == 0 and so % G == 0
        out["burnside"] = {"even_orbits": se // G, "odd_orbits": so // G}

        # orbit sizes: exact-stabilizer Moebius inversion over the FULL
        # subgroup lattice of the (commuting involution) group.  Elements
        # are bitmasks over the generators; every subgroup arises as a span
        # of <= k nonzero masks.  Divisibility asserts certify exactness.
        names = list(self.ops)
        k = len(names)
        Gsize = 2 ** k
        name_of = {}
        for mask in range(1, Gsize):
            combo = [names[i] for i in range(k) if (mask >> i) & 1]
            name_of[mask] = combo
        subgroups = set()
        for r in range(k + 1):
            for gsel in itertools.combinations(range(1, Gsize), r):
                H = {0}
                for x in gsel:
                    H |= {h ^ x for h in list(H)}
                subgroups.add(frozenset(H))
        subgroups = sorted(subgroups, key=lambda h: (len(h), sorted(h)))

        Icnt = (2 ** (self.g - 1) * (2 ** self.g + 1),
                2 ** (self.g - 1) * (2 ** self.g - 1))

        def fixed_of_subgroup(H):
            if len(H) == 1:
                return Icnt
            AA, dd = [], []
            I = np.eye(self.n, dtype=np.uint8)
            for mask in sorted(H):
                if mask == 0:
                    continue
                M = np.eye(self.n, dtype=np.uint8)
                for nm in name_of[mask]:
                    M = (M @ self.ops[nm]) % 2
                A, d = self.affine(M)
                AA.append((I + A) % 2)
                dd.append(d)
            tj, Bj = gf2_solve_affine(np.vstack(AA), np.concatenate(dd))
            if tj is None:
                return (0, 0)
            return self.parity_count(tj, Bj, chunk)

        fixH = {H: fixed_of_subgroup(H) for H in subgroups}
        exact = {}
        for H in sorted(subgroups, key=lambda h: -len(h)):
            ev, od = fixH[H]
            for K, (kev, kod) in exact.items():
                if H < K:
                    ev -= kev
                    od -= kod
            exact[H] = (ev, od)
        bysize = {}
        for H in subgroups:
            ev, od = exact[H]
            osz = Gsize // len(H)
            assert ev % osz == 0 and od % osz == 0, \
                ("orbit divisibility", sorted(H), ev, od)
            se_, so_ = bysize.get(osz, (0, 0))
            bysize[osz] = (se_ + ev // osz, so_ + od // osz)
        for size in sorted(bysize):
            ev, od = bysize[size]
            out["orbit_sizes"][str(size)] = {"even_orbits": ev,
                                             "odd_orbits": od}
        assert sum(v[0] for v in bysize.values()) == \
            out["burnside"]["even_orbits"]
        assert sum(v[1] for v in bysize.values()) == \
            out["burnside"]["odd_orbits"]
        return out


# ------------------------------------------------------------- constructors
def x0143_periods_frame():
    """The G2 = <W11, W13, STAR> action in the mtft.periods symplectic
    frame -- the frame of the frozen Riemann matrix tau0.  This is the
    frame in which theta-function computations at tau0 must be phrased."""
    from mtft import homology
    ops = homology.periods_frame_ops()
    J = homology.mod2(ops["J"])
    act_ops = {k: homology.mod2(ops[k]) for k in ("W11", "W13", "STAR")}
    # bring to standard J2 if the frame's J differs
    g = J.shape[0] // 2
    J2 = standard_J2(g)
    if not np.array_equal(J % 2, J2):
        raise ValueError("periods symplectic form is not standard mod 2; "
                         "conjugate before constructing the action")
    return SpinAction(act_ops, g=g)


def x0143_gp_frame():
    """The same action in the independent PARI/GP integral frame
    (mtft.homology.symplectic_frame).  Counts must agree with the periods
    frame; explicit vectors will differ."""
    from mtft import homology
    _, ops = homology.symplectic_frame()
    return SpinAction({k: homology.mod2(v) for k, v in ops.items()},
                      g=homology.standard_J().shape[0] // 2)


def census():
    """The packaged certified census (Cert EXACT, E2)."""
    with resources.files(__package__).joinpath(
            "_data/x0143_theta_census.json").open() as fh:
        return json.load(fh)


def verify_action(action: SpinAction, samples: int = 200, seed: int = 143):
    """Self-check: the affine action intertwines q_t -> q_t o M^{-1} and
    preserves Arf, on random samples.  Raises on failure."""
    rng = np.random.default_rng(seed)
    n = action.n
    for name, M in action.ops.items():
        A, d = action.affine(M)
        Mi = gf2_inv(M)
        for _ in range(samples):
            t = rng.integers(0, 2, n).astype(np.uint8)
            x = rng.integers(0, 2, n).astype(np.uint8)
            t2 = (A @ t + d) % 2
            g = action.g
            q_t_Mx = (int((Mi @ x)[:g] @ (Mi @ x)[g:] % 2)
                      + int(t @ ((Mi @ x) % 2) % 2)) % 2
            q_t2_x = (int(x[:g] @ x[g:] % 2) + int(t2 @ x % 2)) % 2
            assert q_t_Mx == q_t2_x, f"{name}: affine action law fails"
            assert action.parity(t)[0] == action.parity(t2)[0], \
                f"{name}: Arf not preserved"
    return True
