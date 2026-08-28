"""Numerical Lie-algebra fingerprinting for subalgebras of u(n).

Library form of the route-B machinery that E2-certified the D4 / so(8)
fingerprint of the STAR-fixed triangle algebra of X0(143)
(CERT(tol, E2), 2026-08-27).  Generic over n; the X0(143) constructors
reproduce the certified pipeline in one call (:func:`d4_report`).

Two hard-won gate designs are baked in (ledgered as method notes):

  * :func:`close_lie` uses an ABSOLUTE-residual two-tier gate with an
    ambiguity assert.  Relative gates (residual / bracket norm) run away
    on small-norm brackets: in the original session a relative gate at
    1e-7 leaked past dim 28 and closed on all of u(13).  Accept above
    ``tol_hi``, reject below ``tol_lo``, and REFUSE to decide silently in
    between -- an ambiguous residual raises :class:`LieGateAmbiguous`.
  * :func:`spectral_gap_kernel` reads kernel dimensions off the largest
    log-gap in the singular spectrum, never a fixed tolerance.  A fixed
    cut undercounted the X0(143) normalizer at an impossible 33 (it must
    contain g + Z(g) of dimension 54); the true gap, once located, is
    1.9e7 wide.

All state is explicit; nothing is cached except the X0(143) fixed
channels.
"""
from __future__ import annotations

import json

import numpy as np

__all__ = [
    "LieGateAmbiguous", "inner", "antiherm", "vec_u", "unvec_u", "u_basis",
    "close_lie", "structure_constants", "invariants", "rep_summary",
    "spectral_gap_kernel", "normalizer_in_u", "cartan_and_roots",
    "cosine_bucket_dev", "symmetry_screen",
    "x0143_fixed_channels", "x0143_symmetry_ops", "d4_report",
]


class LieGateAmbiguous(RuntimeError):
    """Raised when a closure residual falls between the reject and accept
    thresholds: the algebra dimension cannot be certified at this noise
    level.  Widen the data precision or the gate consciously -- never
    silently."""


# ------------------------------------------------------------- basic algebra
def inner(X, Y):
    """Real Hilbert-Schmidt inner product on gl(n, C) viewed over R."""
    return float(np.real(np.sum(np.conj(X) * Y)))


def antiherm(X):
    return (X - X.conj().T) / 2


def vec_u(X):
    """Isometric real coordinates of an anti-hermitian n x n matrix
    (dimension n^2)."""
    n = X.shape[0]
    d = X.diagonal().imag.copy()
    iu = np.triu_indices(n, 1)
    return np.concatenate([d, np.sqrt(2) * X[iu].real,
                           np.sqrt(2) * X[iu].imag])


def unvec_u(v, n):
    X = np.zeros((n, n), complex)
    X[np.diag_indices(n)] = 1j * v[:n]
    iu = np.triu_indices(n, 1)
    m = len(iu[0])
    a = v[n:n + m] / np.sqrt(2)
    b = v[n + m:] / np.sqrt(2)
    X[iu] = a + 1j * b
    X[iu[1], iu[0]] = -a + 1j * b
    return X


def u_basis(n):
    """Orthonormal (Hilbert-Schmidt) basis of u(n), n^2 elements."""
    out = []
    for i in range(n):
        E = np.zeros((n, n), complex)
        E[i, i] = 1j
        out.append(E)
    for i in range(n):
        for j in range(i + 1, n):
            E = np.zeros((n, n), complex)
            E[i, j] = 1
            E[j, i] = -1
            out.append(E / np.sqrt(2))
            E = np.zeros((n, n), complex)
            E[i, j] = 1j
            E[j, i] = 1j
            out.append(E / np.sqrt(2))
    return out


# ------------------------------------------------------------------- closure
def close_lie(seeds, tol_hi=1e-5, tol_lo=1e-7, max_rounds=10):
    """Lie closure of anti-hermitized seeds by greedy orthogonalization in
    isometric vec coordinates, with the absolute two-tier gate.

    Returns dict with keys ``basis`` (list of unit anti-hermitian
    matrices), ``growth``, ``min_accepted``, ``max_rejected``,
    ``separation``.  Raises :class:`LieGateAmbiguous` if any residual
    lands in (tol_lo, tol_hi]."""
    n = seeds[0].shape[0]
    Vrows = np.zeros((0, n * n))
    growth = []
    state = {"acc": np.inf, "rej": 0.0}
    ambiguous: list[float] = []

    def try_add(C):
        w = vec_u(antiherm(C))
        if np.linalg.norm(w) < 1e-12:
            return False
        for _ in range(3):
            if len(Vrows):
                w = w - Vrows.T @ (Vrows @ w)
        r = float(np.linalg.norm(w))            # ABSOLUTE residual
        if r > tol_hi:
            nonlocal_rows.append(w / r)
            state["acc"] = min(state["acc"], r)
            return True
        if r > tol_lo:
            ambiguous.append(r)
        state["rej"] = max(state["rej"], r)
        return False

    nonlocal_rows = []

    def flushrows():
        nonlocal Vrows, nonlocal_rows
        if nonlocal_rows:
            Vrows = np.vstack([Vrows] + nonlocal_rows)
            nonlocal_rows = []

    for C in seeds:
        try_add(C)
        flushrows()
    growth.append(len(Vrows))
    for _ in range(max_rounds):
        cur = [unvec_u(v, n) for v in Vrows]
        changed = False
        for i in range(len(cur)):
            for j in range(i + 1, len(cur)):
                if try_add(cur[i] @ cur[j] - cur[j] @ cur[i]):
                    changed = True
                flushrows()
        growth.append(len(Vrows))
        assert len(Vrows) <= n * n, "closure exceeded dim u(n)"
        if not changed:
            break
    if ambiguous:
        raise LieGateAmbiguous(
            f"{len(ambiguous)} residuals in ({tol_lo:g}, {tol_hi:g}]; "
            f"worst {sorted(ambiguous)[:5]}")
    basis = [unvec_u(v, n) for v in Vrows]
    return {"basis": basis, "growth": growth,
            "min_accepted": float(state["acc"]),
            "max_rejected": float(state["rej"]),
            "separation": float(state["acc"] / max(state["rej"], 1e-300))}


def structure_constants(basis):
    """(ad, maxres): ad[i][:, j] holds the coefficients of [b_i, b_j] in
    the basis; maxres is the worst relative reconstruction residual (the
    closure certificate)."""
    m = len(basis)
    ad = np.zeros((m, m, m))
    maxres = 0.0
    for i in range(m):
        for j in range(m):
            C = basis[i] @ basis[j] - basis[j] @ basis[i]
            c = np.array([inner(bk, C) for bk in basis])
            rec = sum(ck * bk for ck, bk in zip(c, basis))
            nC = np.sqrt(max(inner(C, C), 1e-300))
            maxres = max(maxres,
                         np.sqrt(inner(C - rec, C - rec)) / nC)
            ad[i][:, j] = c
    return ad, float(maxres)


def invariants(basis, ad, seed=2026):
    """center_dim, derived_dim, rank, Cartan basis, cartan_comm,
    killing_signature (neg, zero, pos)."""
    m = len(basis)
    adm = [ad[i] for i in range(m)]
    Cen = np.vstack([np.array([[ad[j][k, i] for j in range(m)]
                               for k in range(m)]) for i in range(m)])
    scen = np.linalg.svd(Cen, compute_uv=False)
    center_dim = int((scen < 1e-8 * scen[0]).sum())
    Der = np.array([ad[i][:, j] for i in range(m) for j in range(m)])
    sder = np.linalg.svd(Der, compute_uv=False)
    derived_dim = int((sder > 1e-8 * sder[0]).sum())

    rng = np.random.default_rng(seed)
    coef = rng.standard_normal(m)
    adX = sum(c * A for c, A in zip(coef, adm))
    _, s, vt = np.linalg.svd(adX)
    rank = int((s < 1e-7 * s[0]).sum())
    Hcoef = vt[-rank:]
    Hs = [sum(Hcoef[a2, k] * basis[k] for k in range(m))
          for a2 in range(rank)]
    cart_comm = max(
        np.sqrt(inner(Hs[i] @ Hs[j] - Hs[j] @ Hs[i],
                      Hs[i] @ Hs[j] - Hs[j] @ Hs[i]))
        for i in range(rank) for j in range(i + 1, rank)) if rank > 1 else 0.0

    K = np.array([[np.trace(adm[a2] @ adm[b2]) for b2 in range(m)]
                  for a2 in range(m)])
    kev = np.linalg.eigvalsh((K + K.T) / 2)
    thr = 1e-8 * np.max(np.abs(kev))
    ksig = (int((kev < -thr).sum()), int((np.abs(kev) <= thr).sum()),
            int((kev > thr).sum()))
    return {"center_dim": center_dim, "derived_dim": derived_dim,
            "rank": rank, "cartan": Hs, "Hcoef": Hcoef,
            "cartan_comm": float(cart_comm),
            "killing_signature": list(ksig)}


def rep_summary(basis):
    """Decomposition data of the defining representation on C^n: common
    fixed dimension, active dimension, invariance leak, commutant of the
    algebra in u(n), commutant restricted to the active block."""
    n = basis[0].shape[0]
    Mstack = np.vstack([b for b in basis])
    sv = np.linalg.svd(Mstack, compute_uv=False)
    fix_dim = int((sv < 1e-9 * sv[0]).sum())
    _, _, vh = np.linalg.svd(Mstack)
    Vfix = vh.conj().T[:, n - fix_dim:]
    Pfix = Vfix @ Vfix.conj().T
    Pact = np.eye(n) - Pfix
    leak = max(np.linalg.norm(Pfix @ b @ Pact) for b in basis)

    def flat(X):
        return np.concatenate([X.real.ravel(), X.imag.ravel()])

    Un = u_basis(n)
    Comm = np.array([np.concatenate([flat(e @ b - b @ e) for b in basis])
                     for e in Un]).T
    scomm = np.linalg.svd(Comm, compute_uv=False)
    comm_dim = int((scomm < 1e-8 * scomm[0]).sum())

    W = np.linalg.svd(Pact)[0][:, :n - fix_dim]
    bact = [W.conj().T @ b @ W for b in basis]
    Ua = u_basis(n - fix_dim)
    rows = [np.concatenate([flat(e @ b - b @ e) for b in bact])
            for e in Ua]
    sa = np.linalg.svd(np.array(rows).T, compute_uv=False)
    act_comm = int((sa < 1e-8 * sa[0]).sum())
    return {"common_fixed_dim": fix_dim, "active_dim": n - fix_dim,
            "invariance_leak": float(leak), "u_n_commutant_dim": comm_dim,
            "active_commutant_dim": act_comm, "P_active": Pact,
            "P_fixed": Pfix}


def spectral_gap_kernel(s):
    """Kernel dimension of a singular spectrum by the largest log-gap.
    Returns (kernel_dim, gap, sv_below/s_max, sv_above/s_max)."""
    s = np.asarray(s, float)
    sv = np.sort(s)
    ratios = sv[1:] / np.maximum(sv[:-1], 1e-300)
    k = int(np.argmax(ratios)) + 1
    return (k, float(ratios[k - 1]), float(sv[k - 1] / sv[-1]),
            float(sv[k] / sv[-1]))


def normalizer_in_u(basis):
    """Dimension data of N_{u(n)}(g) via the spectral-gap cut."""
    n = basis[0].shape[0]

    def flat(X):
        return np.concatenate([X.real.ravel(), X.imag.ravel()])

    Bflat = np.array([flat(b) for b in basis])
    Q, _ = np.linalg.qr(Bflat.T)
    Pg = Q @ Q.T
    Un = u_basis(n)
    rows = []
    I = np.eye(Pg.shape[0])
    for e in Un:
        rows.append(np.concatenate([(I - Pg) @ flat(e @ b - b @ e)
                                    for b in basis]))
    snorm = np.linalg.svd(np.array(rows).T, compute_uv=False)
    k, gap, below, above = spectral_gap_kernel(snorm)
    return {"kernel_dim": k, "gap": gap, "sv_below": below,
            "sv_above": above, "projector_g": Pg}


def cartan_and_roots(basis, ad, inv=None):
    """Roots of the algebra with respect to the generic Cartan of
    :func:`invariants`.  Returns dict with the root array and
    diagnostics (count, length_ratio, plus_minus_paired, zero_modes)."""
    if inv is None:
        inv = invariants(basis, ad)
    m = len(basis)
    rank = inv["rank"]
    Hcoef = inv["Hcoef"]
    adm = [ad[i] for i in range(m)]
    adH = [sum(Hcoef[a2, k] * adm[k] for k in range(m))
           for a2 in range(rank)]
    cgen = np.array([1.0, np.e / 2, np.pi / 3, np.sqrt(7) / 2,
                     np.sqrt(2), np.log(3), 5 / 7, np.sqrt(11) / 3,
                     np.euler_gamma, 8 / 13, np.sqrt(5) / 4, 3 / np.e,
                     np.pi / 7][:rank])
    A0 = sum(c * A for c, A in zip(cgen, adH))
    ev, V = np.linalg.eigh(1j * A0)
    roots = []
    for k2, lam in enumerate(ev):
        if abs(lam) > 1e-6 * np.max(np.abs(ev)):
            v = V[:, k2]
            roots.append(np.array(
                [float(np.real(v.conj() @ (1j * A) @ v)) for A in adH]))
    roots = np.array(roots)
    lens = np.linalg.norm(roots, axis=1)
    pair_ok = True
    for r in roots:
        d = np.min(np.linalg.norm(roots + r, axis=1))
        pair_ok &= bool(d < 1e-7 * lens.max())
    return {"roots": roots,
            "count": int(len(roots)),
            "length_ratio": float(lens.max() / lens.min()),
            "plus_minus_paired": pair_ok,
            "zero_modes": int(m - len(roots))}


def cosine_bucket_dev(roots, buckets=(0.0, 0.5, -0.5, 1.0, -1.0)):
    """Max deviation of pairwise root cosines from the nearest bucket."""
    lens = np.linalg.norm(roots, axis=1)
    G = (roots @ roots.T) / np.outer(lens, lens)
    dev = 0.0
    for i in range(len(roots)):
        for j in range(i + 1, len(roots)):
            dev = max(dev, min(abs(G[i, j] - b) for b in buckets))
    return float(dev)


def symmetry_screen(ops, R, Ri, basis, Pg=None):
    """Test named real 2n x 2n operators for normalization of g and for
    identity action.  ``R, Ri`` map the real frame to the adapted frame;
    linear vs antilinear is detected by (anti)commutation with the
    standard complex structure."""
    n = basis[0].shape[0]

    def flat(X):
        return np.concatenate([X.real.ravel(), X.imag.ravel()])

    if Pg is None:
        Bflat = np.array([flat(b) for b in basis])
        Q, _ = np.linalg.qr(Bflat.T)
        Pg = Q @ Q.T
    Jf = np.zeros((2 * n, 2 * n))
    Jf[:n, n:] = -np.eye(n)
    Jf[n:, :n] = np.eye(n)
    screen = {}
    for name, O in ops.items():
        Oa = R @ O @ Ri
        comm = np.linalg.norm(Oa @ Jf - Jf @ Oa) / np.linalg.norm(Oa)
        anti = np.linalg.norm(Oa @ Jf + Jf @ Oa) / np.linalg.norm(Oa)
        Uc = Oa[:n, :n] + 1j * Oa[n:, :n]
        Uci = np.linalg.inv(Uc)
        if comm < anti:
            phi = lambda b, Uc=Uc, Uci=Uci: Uc @ b @ Uci
            kind, struct = "linear", comm
        else:
            phi = lambda b, Uc=Uc, Uci=Uci: Uc @ np.conj(b) @ Uci
            kind, struct = "antilinear", anti
        res = iden = 0.0
        for b in basis:
            im = phi(b)
            f = flat(im)
            res = max(res, np.linalg.norm(f - Pg @ f) / np.linalg.norm(f))
            iden = max(iden, np.sqrt(inner(im - b, im - b)))
        screen[name] = {"kind": kind, "structure_residual": float(struct),
                        "normalization_residual": float(res),
                        "identity_deviation_on_g": float(iden)}
    return screen


# --------------------------------------------------------- X0(143) pipeline
_X143: dict = {}


def x0143_fixed_channels(dps=50, seed=143):
    """The three primitive STAR-fixed channels of the X0(143) stage as
    complex 13 x 13 matrices, plus the adapted-frame maps (R, Ri).

    Faithful port of the certified construction: Hodge metric and complex
    structure at ``dps`` digits, a J-adapted orthonormal frame seeded by
    RNG(``seed``), the m7 harmonic basis, per-triangle averaging masks at
    the three primitive triangles (1, 11, 12), hermitian and hamiltonian
    splits, and the complexification to B-coordinates."""
    key = (dps, seed)
    if key in _X143:
        return _X143[key]
    from mtft import hecke as H
    from mtft.periods import physics as PH, data_path
    from mtft.periods.hamiltonian import hermitian_split, hamiltonian_split

    N = 26
    rng = np.random.default_rng(seed)
    G = PH.hodge_metric_hecke(dps)
    J = PH.hodge_structure_hecke(dps)
    w, Q = np.linalg.eigh((G + G.T) / 2)
    Gh = Q @ np.diag(np.sqrt(w)) @ Q.T
    Ghi = Q @ np.diag(1 / np.sqrt(w)) @ Q.T
    Jt = Gh @ J @ Ghi
    Jt = (Jt - Jt.T) / 2
    Us, Ws = [], []
    for _ in range(13):
        v = rng.standard_normal(N)
        for u2, z in zip(Us, Ws):
            v -= (u2 @ v) * u2 + (z @ v) * z
        v /= np.linalg.norm(v)
        z = Jt @ v
        for u2, wv in zip(Us, Ws):
            z -= (u2 @ z) * u2 + (wv @ z) * wv
        z /= np.linalg.norm(z)
        Us.append(v)
        Ws.append(z)
    Qb = np.array(Us + Ws).T
    R, Ri = Qb.T @ Gh, Ghi @ Qb

    def toB(A):
        Aa = R @ A @ Ri
        Bm = np.empty((13, 13), complex)
        for k in range(13):
            x = np.zeros(N)
            x[k] = 1
            y = Aa @ x
            Bm[:, k] = y[:13] + 1j * y[13:]
        return Bm

    m = H.model()
    rec = json.loads(data_path("X0_143_m7_harmonic_basis.json").read_text())
    Wh = np.array([[a / b for a, b in row] for row in rec["basis_26x84"]],
                  float)
    Ggi = np.linalg.inv(Wh @ Wh.T)
    edges = [(m["tri_of"][m["erep"][k]], m["tri_of"][m["sS"][m["erep"][k]]])
             for k in range(m["E"])]
    fixed = []
    for t in (1, 11, 12):
        gav = np.array([.5 * ((a2 == t) + (b2 == t)) for a2, b2 in edges])
        V = Ggi @ (Wh @ np.diag(gav) @ Wh.T)
        A, _ = hermitian_split(V, dps)
        _, Am = hamiltonian_split(A, dps)
        fixed.append(toB(Am))
    _X143[key] = (R, Ri, fixed)
    return _X143[key]


def x0143_symmetry_ops():
    """The known arithmetic operators in the (real, Hecke-frame) 26 x 26
    representation used by the screen: W11, W13, W143 and STAR."""
    from mtft import periods as P
    from mtft.periods.involutions import al_matrix
    Csym = np.array(P.hecke_to_symplectic_change(), float)
    Cinv = np.linalg.inv(Csym)
    ops = {"W11": np.array(al_matrix(11), float),
           "W13": np.array(al_matrix(13), float)}
    try:
        ops["W143"] = np.array(al_matrix(143), float)
    except Exception:
        ops["W143"] = ops["W11"] @ ops["W13"]
    ops["STAR"] = Cinv @ np.array(P.star_symplectic(), float) @ Csym
    return ops


def d4_report(screen=True, dps=50):
    """One-call reproduction of the certified D4 fingerprint of the
    STAR-fixed triangle algebra: closure, structure, representation,
    normalizer, roots, and (optionally) the arithmetic symmetry screen.
    Expected values (CERT(tol, E2)): dim 28, growth [3, 6, 17, 28, 28],
    center 0, derived 28, rank 4, Killing (28, 0, 0), rep 8 + 1^5 with
    active commutant 1, normalizer 54, 24 equal-length roots on the D4
    cosine buckets, STAR = identity on g, W11/W13/W143 non-normalizing."""
    R, Ri, fixed = x0143_fixed_channels(dps=dps)
    seeds = [fixed[i] @ np.conj(fixed[j]) - fixed[j] @ np.conj(fixed[i])
             for i in range(3) for j in range(i + 1, 3)]
    cl = close_lie(seeds)
    basis = cl["basis"]
    ad, maxres = structure_constants(basis)
    inv = invariants(basis, ad)
    rep = rep_summary(basis)
    nz = normalizer_in_u(basis)
    rt = cartan_and_roots(basis, ad, inv)
    out = {
        "closure": {k: cl[k] for k in
                    ("growth", "min_accepted", "max_rejected", "separation")},
        "dim": len(basis),
        "structure_closure_maxres": maxres,
        "structure": {k: inv[k] for k in
                      ("center_dim", "derived_dim", "rank", "cartan_comm",
                       "killing_signature")},
        "representation": {k: rep[k] for k in
                           ("common_fixed_dim", "active_dim",
                            "invariance_leak", "u_n_commutant_dim",
                            "active_commutant_dim")},
        "normalizer": {k: nz[k] for k in
                       ("kernel_dim", "gap", "sv_below", "sv_above")},
        "roots": {k: rt[k] for k in ("count", "length_ratio",
                                     "plus_minus_paired", "zero_modes")},
        "roots_cosine_dev": cosine_bucket_dev(rt["roots"]),
    }
    if screen:
        out["symmetry_screen"] = symmetry_screen(
            x0143_symmetry_ops(), R, Ri, basis, Pg=nz["projector_g"])
    return out
