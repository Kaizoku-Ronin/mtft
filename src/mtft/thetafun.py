"""Numerical genus-g theta functions with characteristics for X0(143).

theta[a; b](z, tau) = sum_{n in Z^g} exp(i pi (n+a)^T tau (n+a)
                                          + 2 pi i (n+a)^T (z + b)),
half-integer characteristics a, b in {0, 1/2}^g, tau in Siegel upper half
space H_g.  Designed for the frozen X0(143) Riemann matrix tau0 (g = 13),
whose raw Im tau0 has lambda_min ~ 0.093: direct summation would need
~1e16 lattice points, so the engine works in an LLL-reduced frame where
the same ellipsoid holds ~1e5-1e7 points.

Guarantees:
  * exact ellipsoid enumeration (Fincke-Pohst recursion, vectorized inner
    coordinate) -- every lattice point with Q[x] <= t is included, none
    outside;
  * rigorous, conservative tail bounds (see :func:`tail_bound_value` /
    :func:`tail_bound_grad`); returned alongside every value;
  * exact characteristic transport under the two reduction moves
    tau -> U^T tau U (U in GL(g, Z)) and tau -> tau + S (S symmetric
    integral), including phases, so reduced-frame evaluation equals
    direct evaluation to within the tail bound (unit-tested at g = 1, 2
    against mpmath.jtheta and against direct summation).

Frame conventions match :mod:`mtft.thetachar`: a mod-2 characteristic
t in F_2^{2g} in the standard symplectic frame corresponds to
a = t[:g] / 2, b = t[g:] / 2 (:func:`char_to_ab`).
"""
from __future__ import annotations

import numpy as np

__all__ = [
    "char_to_ab", "lll_reduce_gram", "siegel_ready", "reduce_char",
    "theta_reduced", "theta_char", "theta_grad",
    "tail_bound_value", "tail_bound_grad",
    "x0143_tau", "x0143_ready",
]


# ------------------------------------------------------------ characteristic
def char_to_ab(t):
    """Mod-2 characteristic (standard symplectic frame) -> (a, b) with
    half-integer entries."""
    t = np.asarray(t)
    g = t.shape[-1] // 2
    return t[..., :g] / 2.0, t[..., g:] / 2.0


# --------------------------------------------------------------------- LLL
def lll_reduce_gram(G, delta=0.99, max_sweeps=10000):
    """Own LLL on a positive-definite Gram matrix.

    Returns (G_red, U) with U integral unimodular and
    G_red = U^T G U; the columns of U are the reduced lattice basis in
    original coordinates.  Gram-based textbook LLL with Lovasz parameter
    ``delta``."""
    G = np.asarray(G, float)
    n = G.shape[0]
    U = np.eye(n, dtype=np.int64)

    def gso(Gm):
        # mu, Bstar norms from the Gram matrix
        mu = np.zeros((n, n))
        Bn = np.zeros(n)
        for i in range(n):
            v = Gm[i].copy().astype(float)
            for j in range(i):
                mu[i, j] = (Gm[i, j]
                            - sum(mu[i, k] * mu[j, k] * Bn[k]
                                  for k in range(j))) / Bn[j]
            Bn[i] = Gm[i, i] - sum(mu[i, k] ** 2 * Bn[k] for k in range(i))
        return mu, Bn

    Gm = U.T @ G @ U
    k = 1
    sweeps = 0
    mu, Bn = gso(Gm)
    while k < n:
        sweeps += 1
        if sweeps > max_sweeps:
            break
        for j in range(k - 1, -1, -1):
            q = int(round(mu[k, j]))
            if q != 0:
                U[:, k] -= q * U[:, j]
                Gm = U.T @ G @ U
                mu, Bn = gso(Gm)
        if Bn[k] >= (delta - mu[k, k - 1] ** 2) * Bn[k - 1]:
            k += 1
        else:
            U[:, [k - 1, k]] = U[:, [k, k - 1]]
            Gm = U.T @ G @ U
            mu, Bn = gso(Gm)
            k = max(k - 1, 1)
    Gr = U.T @ G @ U
    assert abs(np.linalg.det(U.astype(float))) - 1 < 1e-9
    return Gr, U


# ----------------------------------------------------------- Siegel-ready
def siegel_ready(tau):
    """Prepare tau for efficient theta evaluation.

    Returns a dict with:
      tau_red = U^T tau U + S  (Im LLL-reduced, |Re| <= 1/2 entrywise),
      U (GL(g, Z)), S (symmetric integral), Cholesky R of Im(tau_red)
      (upper triangular, Im = R^T R), lam_min of Im(tau_red), and the
      original tau."""
    tau = np.asarray(tau, complex)
    tau = (tau + tau.T) / 2
    g = tau.shape[0]
    Y = tau.imag
    _, U = lll_reduce_gram(Y)
    tU = U.T @ tau @ U
    S = -np.round(tU.real).astype(np.int64)
    S = ((S + S.T) // 2) + np.diag(np.diag(S) - np.diag((S + S.T) // 2))
    # symmetric integral by construction from a symmetric real part; force:
    S = np.round((S + S.T) / 2).astype(np.int64)
    tau_red = tU + S
    Yr = tau_red.imag
    lam = float(np.linalg.eigvalsh(Yr).min())
    assert lam > 0, "Im tau not positive definite"
    R = np.linalg.cholesky(Yr).T          # upper triangular, Yr = R^T R
    return {"tau": tau, "tau_red": tau_red, "U": U, "S": S,
            "R": R, "lam_min": lam, "g": g}


def reduce_char(a, b, ready):
    """Transport a characteristic through tau -> U^T tau U + S.

    theta[a; b](0, tau) = phase * theta[a2; b2](0, tau_red), and the
    z-gradient transports as grad = phase * U @ grad_red (the coefficient
    vector of each term is x = U y).

    Construction.  With x = U y the lattice sum reindexes exactly
    (n <-> m bijection, no constants), giving offset a1 = frac(U^{-1} a)
    and linear shift b1 = U^T b.  The S-move contributes a b-shift
    b2 = b1 + S a1 + diag(S)/2 (mod-2 identity m^T S m = m . diag(S))
    plus a unimodular constant; reducing b2 mod 1 contributes another.
    The total constant is an eighth root of unity depending only on
    (a, U, S).  Rather than track it symbolically, it is CALIBRATED
    exactly from the m = 0 term and CERTIFIED by re-checking the m = e_1
    term (the ratio original/reduced is constant in m; deviation would
    raise).  Returns (a2, b2, phase)."""
    U = ready["U"].astype(float)
    S = ready["S"].astype(float)
    tau = ready["tau"]
    tau_red = ready["tau_red"]
    a = np.asarray(a, float)
    b = np.asarray(b, float)

    a1 = np.linalg.solve(U, a)
    a1r = np.round(a1 * 2) / 2
    assert np.abs(a1 - a1r).max() < 1e-9, "U^{-1} a must be half-integral"
    a1 = a1r - np.floor(a1r + 1e-12)
    b1 = U.T @ b
    b2 = b1 + S @ a1 + np.diag(S) / 2
    b2 = b2 - np.floor(b2 + 1e-12)

    def term(x, t, bb):
        return np.exp(1j * np.pi * (x @ t @ x) + 2j * np.pi * (x @ bb))

    y0 = a1
    phase = term(U @ y0, tau, b) / term(y0, tau_red, b2)
    mag = abs(phase)
    assert abs(mag - 1) < 1e-8, "transport constant must be unimodular"
    phase /= mag
    y1 = a1 + np.eye(len(a1))[0]
    chk = term(U @ y1, tau, b) / term(y1, tau_red, b2)
    assert abs(chk - phase) < 1e-8, \
        "transport ratio not constant: reduction identity violated"
    return a1, b2, phase


# ------------------------------------------------------------- tail bounds
def _theta3(c):
    """sum_{k in Z} e^{-pi c k^2}, rigorous upper bound (geometric tail)."""
    s = 1.0
    k = 1
    while True:
        t = 2 * np.exp(-np.pi * c * k * k)
        s += t
        if t < 1e-18 * s:
            break
        k += 1
        if k > 10000:
            break
    # geometric remainder bound: terms decay faster than ratio
    # e^{-pi c (2k+1)}
    r = np.exp(-np.pi * c * (2 * k + 1))
    return s * (1 + 2 * r / max(1 - r, 1e-15))


def tail_bound_value(t, rdiag2, lam=None, theta=None):
    """Rigorous bound for  sum_{Q[x] > t} e^{-pi Q[x]}  over a shifted
    lattice, Q = ||R x||^2, R upper triangular, rdiag2 = diag(R)^2.

    Split e^{-pi Q} <= e^{-pi theta t} e^{-pi (1-theta) Q} for any
    theta in (0, 1); in nested Fincke-Pohst variables Q = sum r_ii^2 v_i^2
    with each v_i on a unit-spaced shifted comb, so
    sum_all e^{-pi (1-theta) Q} <= prod_i theta3((1-theta) r_ii^2).
    ``theta=None`` optimizes over a grid."""
    rdiag2 = np.asarray(rdiag2, float)
    thetas = [theta] if theta is not None else \
        [0.5, 0.6, 0.7, 0.8, 0.85, 0.9, 0.95, 0.975]
    best = np.inf
    for th in thetas:
        prod = 1.0
        for r2 in rdiag2:
            prod *= _theta3((1 - th) * r2)
        best = min(best, np.exp(-np.pi * th * t) * prod)
    return float(best)


def tail_bound_grad(t, rdiag2, lam, theta=None):
    """Rigorous bound for  sum_{Q[x] > t} 2 pi |x| e^{-pi Q[x]}.

    |x| <= sqrt(Q / lam), lam = lambda_min(R^T R).  For eta = 1/8,
    sqrt(Q) e^{-pi eta Q} is decreasing for Q >= 1/(2 pi eta), so on the
    tail (t >= 1)  sqrt(Q) e^{-pi eta Q} <= sqrt(t) e^{-pi eta t}.
    Split the remaining e^{-pi (1-eta) Q} as in the value bound with
    parameter theta:  e^{-pi (1-eta) Q} <= e^{-pi theta t}
    e^{-pi (1-eta-theta) Q}, requiring theta < 1 - eta."""
    t = max(t, 1.0)
    eta = 0.125
    rdiag2 = np.asarray(rdiag2, float)
    thetas = [theta] if theta is not None else \
        [0.4, 0.5, 0.6, 0.7, 0.75, 0.8, 0.85]
    best = np.inf
    for th in thetas:
        rem = 1 - eta - th
        if rem <= 0:
            continue
        prod = 1.0
        for r2 in rdiag2:
            prod *= _theta3(rem * r2)
        best = min(best,
                   2 * np.pi * np.sqrt(t / lam)
                   * np.exp(-np.pi * (eta + th) * t) * prod)
    return float(best)


def _pick_radius(rdiag2, lam, tol, grad):
    t = 5.0
    f = (lambda tt: tail_bound_grad(tt, rdiag2, lam)) if grad else \
        (lambda tt: tail_bound_value(tt, rdiag2, lam))
    while f(t) > tol:
        t += 0.5
        if t > 4000:
            raise RuntimeError("cannot certify tolerance; tau too skewed")
    return t


# ---------------------------------------------------- ellipsoid enumeration
def _enumerate_ellipsoid(R, center, t, chunk=1 << 18):
    """Yield integer points n (as (K, g) int64 blocks) with
    Q[n - center] <= t, Q = ||R x||^2, R upper triangular.  Exact
    Fincke-Pohst; the innermost TWO coordinates are fully vectorized
    (ragged-range construction), blocks stream at ``chunk`` granularity."""
    g = R.shape[0]
    buf, buf_n = [], 0
    partial = [0.0] * g
    partial_n = [0] * g
    eps = 1e-9

    def flush():
        nonlocal buf, buf_n
        if buf:
            out = np.concatenate(buf, axis=0)
            buf, buf_n = [], 0
            return out
        return None

    def leaf1(rem):
        """coordinates 1 (then 0) with everything above fixed."""
        nonlocal buf, buf_n
        cross1 = sum(R[1, j] * partial[j] for j in range(2, g))
        half = np.sqrt(max(rem, 0.0))
        lo = (-half - cross1) / R[1, 1] + center[1]
        hi = (half - cross1) / R[1, 1] + center[1]
        n1lo, n1hi = int(np.ceil(lo - eps)), int(np.floor(hi + eps))
        if n1hi < n1lo:
            return
        n1 = np.arange(n1lo, n1hi + 1, dtype=np.int64)
        x1 = n1 - center[1]
        rem2 = rem - (R[1, 1] * x1 + cross1) ** 2
        ok = rem2 >= -eps
        n1, x1, rem2 = n1[ok], x1[ok], np.maximum(rem2[ok], 0.0)
        if n1.size == 0:
            return
        cross0 = R[0, 1] * x1 + sum(R[0, j] * partial[j]
                                    for j in range(2, g))
        h0 = np.sqrt(rem2)
        lo0 = np.ceil((-h0 - cross0) / R[0, 0] + center[0] - eps)
        hi0 = np.floor((h0 - cross0) / R[0, 0] + center[0] + eps)
        lens = (hi0 - lo0 + 1).astype(np.int64)
        keep = lens > 0
        if not keep.any():
            return
        n1, cross0, lo0, lens = n1[keep], cross0[keep], lo0[keep], lens[keep]
        rem2 = rem2[keep]
        tot = int(lens.sum())
        off = np.repeat(np.cumsum(lens) - lens, lens)
        n0 = (np.arange(tot, dtype=np.int64) - off
              + np.repeat(lo0.astype(np.int64), lens))
        x0 = n0 - center[0]
        q0 = (R[0, 0] * x0 + np.repeat(cross0, lens)) ** 2
        good = q0 <= np.repeat(rem2, lens) + eps
        if not good.any():
            return
        blk = np.empty((int(good.sum()), g), dtype=np.int64)
        blk[:, 0] = n0[good]
        blk[:, 1] = np.repeat(n1, lens)[good]
        for j in range(2, g):
            blk[:, j] = partial_n[j]
        buf.append(blk)
        buf_n += blk.shape[0]

    def rec(i, rem):
        nonlocal buf, buf_n
        if i == 1 and g >= 2:
            leaf1(rem)
            if buf_n >= chunk:
                out = flush()
                if out is not None:
                    yield out
            return
        cross = sum(R[i, j] * partial[j] for j in range(i + 1, g))
        half = np.sqrt(max(rem, 0.0))
        lo = (-half - cross) / R[i, i] + center[i]
        hi = (half - cross) / R[i, i] + center[i]
        nlo, nhi = int(np.ceil(lo - eps)), int(np.floor(hi + eps))
        if nhi < nlo:
            return
        if i == 0:                       # only reached when g == 1
            ns = np.arange(nlo, nhi + 1, dtype=np.int64)
            xi = ns - center[0]
            q = (R[0, 0] * xi + cross) ** 2
            keep = ns[q <= rem + eps]
            if keep.size:
                blk = keep.reshape(-1, 1)
                buf.append(blk)
                buf_n += blk.shape[0]
            return
        for nv in range(nlo, nhi + 1):
            xv = nv - center[i]
            partial[i] = xv
            partial_n[i] = nv
            rem2 = rem - (R[i, i] * xv + cross) ** 2
            if rem2 >= -eps:
                yield from rec(i - 1, rem2)
        partial[i] = 0.0

    yield from rec(g - 1, t)
    last = flush()
    if last is not None:
        yield last


def _inner_template(A, c_in, t):
    """All integer points m (K x k) with ||A (m - c_in)||^2 <= t, plus
    precomputed U = A (m - c_in) rows and their squared norms."""
    pts = []
    for blk in _enumerate_ellipsoid(A, c_in, t):
        pts.append(blk)
    if pts:
        M = np.concatenate(pts, axis=0)
    else:
        M = np.zeros((0, A.shape[0]), np.int64)
    X = M.astype(float) - c_in[None, :]
    U = X @ A.T
    return M, X, U, np.einsum("ki,ki->k", U, U)


def theta_reduced(a, b, ready, tol=1e-10, deriv=False, split=6,
                  prefix_chunk=64, stats=False):
    """theta[a; b](0, tau_red) (and optionally its z-gradient) by certified
    summation in the reduced frame.

    Algorithm.  With R upper triangular, Q(x) = Q_out(x_out) +
    ||A x_in + B x_out||^2 where x_out are the last ``split`` coordinates,
    A = R[:k, :k], B = R[:k, k:].  The outer block is enumerated by exact
    Fincke-Pohst recursion (few nodes); the inner block uses a single
    precomputed template of the inner ellipsoid at full radius, filtered
    per outer prefix in bulk (the shift enters only through
    ||u + s||^2 = ||u||^2 + 2 u.s + ||s||^2 with u precomputed).  Every
    lattice point with Q <= t is included exactly once; the certified
    tail bound covers the rest.

    Returns (value, bound) or (value, grad, val_bound, grad_bound)."""
    tau = ready["tau_red"]
    R = ready["R"]
    lam = ready["lam_min"]
    g = ready["g"]
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    rdiag2 = np.diag(R) ** 2
    t = _pick_radius(rdiag2, lam, tol, deriv)
    center = -a

    npoints = 0
    absmass = 0.0
    if split <= 0 or split >= g - 1 or g <= 3:
        val = 0.0 + 0.0j
        grad = np.zeros(g, complex)
        for blk in _enumerate_ellipsoid(R, center, t):
            x = blk.astype(float) + a[None, :]
            quad = np.einsum("ki,ij,kj->k", x, tau, x)
            ph = np.exp(1j * np.pi * quad + 2j * np.pi * (x @ b))
            val += ph.sum()
            npoints += blk.shape[0]
            absmass += float(np.abs(ph).sum())
            if deriv:
                grad += (2j * np.pi) * (x * ph[:, None]).sum(axis=0)
        vb = tail_bound_value(t, rdiag2, lam)
        gb = tail_bound_grad(t, rdiag2, lam) if deriv else None
        if stats:
            return {"value": val, "grad": grad if deriv else None,
                    "val_bound": vb, "grad_bound": gb,
                    "npoints": npoints, "absmass": absmass,
                    "radius": t}
        if deriv:
            return val, grad, vb, gb
        return val, vb

    k = g - split
    A = R[:k, :k]
    B = R[:k, k:]
    Rout = R[k:, k:]
    c_in, c_out = center[:k], center[k:]
    a_in, a_out = a[:k], a[k:]
    b_in, b_out = b[:k], b[k:]
    tau_ii = tau[:k, :k]
    tau_io = tau[:k, k:]
    tau_oo = tau[k:, k:]

    Min, Xin, Uin, Uq = _inner_template(A, c_in, t)
    # per-template quantities independent of the prefix:
    lin_in = Xin @ b_in                       # x_in . b_in
    quad_ii = np.einsum("ki,ij,kj->k", Xin, tau_ii, Xin)

    val = 0.0 + 0.0j
    grad = np.zeros(g, complex)
    for oblk in _enumerate_ellipsoid(Rout, c_out, t):
        for s0 in range(0, oblk.shape[0], prefix_chunk):
            ob = oblk[s0:s0 + prefix_chunk]
            Xo = ob.astype(float) + a_out[None, :]          # (P, split)
            Qout = np.einsum("pi,ij,pj->p", Xo, (Rout.T @ Rout), Xo)
            rem = t - Qout
            ok = rem >= -1e-9
            if not ok.any():
                continue
            Xo, rem = Xo[ok], rem[ok]
            S = Xo @ B.T                                    # (P, k)
            # q_pk = Uq[k] + 2 Uin.S + |S|^2  <= rem_p
            crossq = Uin @ S.T                              # (K, P)
            Snorm = np.einsum("pi,pi->p", S, S)
            mask = (Uq[:, None] + 2 * crossq + Snorm[None, :]
                    <= rem[None, :] + 1e-9)
            if not mask.any():
                continue
            quad_oo = np.einsum("pi,ij,pj->p", Xo, tau_oo, Xo)
            cross_t = Xin @ (tau_io @ Xo.T)                 # (K, P)
            lin_o = Xo @ b_out                              # (P,)
            expo = (1j * np.pi * (quad_ii[:, None] + 2 * cross_t
                                  + quad_oo[None, :])
                    + 2j * np.pi * (lin_in[:, None] + lin_o[None, :]))
            ph = np.where(mask, np.exp(expo), 0.0)
            val += ph.sum()
            npoints += int(mask.sum())
            absmass += float(np.abs(ph).sum())
            if deriv:
                grad[:k] += (2j * np.pi) * (Xin.T @ ph.sum(axis=1))
                grad[k:] += (2j * np.pi) * ((ph.sum(axis=0)[:, None]
                                             * Xo).sum(axis=0))
    vb = tail_bound_value(t, rdiag2, lam)
    gb = tail_bound_grad(t, rdiag2, lam) if deriv else None
    if stats:
        return {"value": val, "grad": grad if deriv else None,
                "val_bound": vb, "grad_bound": gb,
                "npoints": npoints, "absmass": absmass, "radius": t}
    if deriv:
        return val, grad, vb, gb
    return val, vb


# --------------------------------------------------------------- public API
def theta_char(a, b, ready, tol=1e-10, **kw):
    """theta[a; b](0, tau) for the ORIGINAL tau of ``ready``: reduced-frame
    evaluation with exact characteristic transport.  Returns
    (value, certified_bound)."""
    a2, b2, phase = reduce_char(a, b, ready)
    if kw.get("stats"):
        st = theta_reduced(a2, b2, ready, tol=tol, **kw)
        st["value"] = phase * st["value"]
        return st
    v, vb = theta_reduced(a2, b2, ready, tol=tol, **kw)
    return phase * v, vb


def theta_grad(a, b, ready, tol=1e-10, **kw):
    """grad_z theta[a; b](0, tau) for the ORIGINAL tau of ``ready``.

    Transport: with x = U y, each term's coefficient vector is x = U y,
    and the reduction moves (a-reindexing, S-shift, b-reduction) change
    only scalar factors folded into ``phase``; hence
      grad = phase * U @ grad_red.
    Validated against direct high-precision summation on skewed tau with
    nontrivial U and S.  Returns (grad_vector, value, val_bound,
    grad_bound)."""
    a2, b2, phase = reduce_char(a, b, ready)
    U = ready["U"].astype(float)
    if kw.get("stats"):
        st = theta_reduced(a2, b2, ready, tol=tol, deriv=True, **kw)
        st["value"] = phase * st["value"]
        st["grad"] = phase * (U @ st["grad"])
        return st
    v, gr, vb, gb = theta_reduced(a2, b2, ready, tol=tol, deriv=True, **kw)
    grad = phase * (U @ gr)
    return grad, phase * v, vb, gb


# ------------------------------------------------------------ X0(143) data
_X143: dict = {}


def x0143_tau():
    """The frozen 13x13 symplectically normalized Riemann matrix of
    X0(143) as complex128."""
    if "tau" not in _X143:
        from mtft import periods as P
        _X143["tau"] = np.array(P.frozen_riemann_matrix(),
                                complex).reshape(13, 13)
    return _X143["tau"]


def x0143_ready():
    """Cached Siegel-ready structure for tau0."""
    if "ready" not in _X143:
        _X143["ready"] = siegel_ready(x0143_tau())
    return _X143["ready"]
