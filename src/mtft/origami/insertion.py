"""mtft.origami.insertion — the observable insertion calculus.

One formula underlies the Tano moments, the Dirichlet log-n metric, dimer edge
occupations, the two-spin Ising model, the prism, and differentiation under
the integral sign:

    Z(lambda) = sum_x w(x) exp(sum_i lambda_i A_i(x))
    d^k/dlambda^k  log Z  =  kappa_k(A)          (cumulants, not moments)

The hierarchy this module exposes:

    psi   = log Z                     scalar potential
    grad psi                          conservative expectation field
    Hess psi = g                      Fisher metric
    third derivatives = T             cubic Amari tensor
    Brioschi(g) = K                   Gaussian curvature

Everything is done in MULTIPLICATIVE coordinates X_i = exp(theta_i) with the
operator D_i = X_i d/dX_i, which keeps expressions rational and avoids the
sympy stalls that plague exp/log forms.
"""
from __future__ import annotations

import sympy as sp

__all__ = [
    "D_log", "cumulants", "fisher_metric", "cubic_tensor",
    "brioschi_curvature", "cumulant_curvature", "fisher_pack",
    "path_independence",
]


def D_log(f, X):
    """The Euler operator X d/dX, i.e. d/d(log X)."""
    return sp.cancel(X * sp.diff(f, X))


def cumulants(Z, coords, order=4):
    """Joint cumulants of the sufficient statistics, as rational functions.

    Returns dict: tuple of coordinate indices (sorted, with repetition) ->
    the corresponding derivative of log Z.  Order 1 = means, 2 = covariance,
    3 = cubic Amari tensor, and so on.
    """
    psi = sp.log(Z)
    out = {}
    cur = {(): psi}
    for _ in range(order):
        nxt = {}
        for key, val in cur.items():
            for i, X in enumerate(coords):
                nk = tuple(sorted(key + (i,)))
                if nk in nxt or nk in out:
                    continue
                nxt[nk] = D_log(val, X)
        out.update(nxt)
        cur = nxt
    return out


def fisher_metric(Z, coords):
    """g_ij = Cov(A_i, A_j) = D_i D_j log Z, as a sympy Matrix."""
    psi = sp.log(Z)
    n = len(coords)
    first = [D_log(psi, X) for X in coords]
    return sp.Matrix(n, n, lambda i, j: D_log(first[i], coords[j]))


def cubic_tensor(Z, coords):
    """T_ijk = D_i D_j D_k log Z (third cumulants), as a nested dict."""
    g = fisher_metric(Z, coords)
    n = len(coords)
    return {(i, j, k): D_log(g[i, j], coords[k])
            for i in range(n) for j in range(n) for k in range(n)}


def brioschi_curvature(Z, coords):
    """Gaussian curvature of the 2-D Fisher metric via the Brioschi formula."""
    if len(coords) != 2:
        raise ValueError("Brioschi curvature needs exactly 2 coordinates")
    a, b = coords
    psi = sp.log(Z)
    E = D_log(D_log(psi, a), a)
    F = D_log(D_log(psi, a), b)
    G = D_log(D_log(psi, b), b)
    Ea, Eb = D_log(E, a), D_log(E, b)
    Fa, Fb = D_log(F, a), D_log(F, b)
    Ga, Gb = D_log(G, a), D_log(G, b)
    Ebb, Fab, Gaa = D_log(Eb, b), D_log(Fa, b), D_log(Ga, a)
    M1 = sp.Matrix([[-Ebb / 2 + Fab - Gaa / 2, Ea / 2, Fa - Eb / 2],
                    [Fb - Ga / 2, E, F],
                    [Gb / 2, F, G]])
    M2 = sp.Matrix([[0, Eb / 2, Ga / 2],
                    [Eb / 2, E, F],
                    [Ga / 2, F, G]])
    return sp.cancel((M1.det() - M2.det()) / (E * G - F ** 2) ** 2)


def cumulant_curvature(Z, coords):
    """Independent (E2) curvature route: Amari 0-curvature from third cumulants.

    SIGN CONVENTION (documented, calibrated on the trinomial where Brioschi is
    symbolically proven to give +1/4): the index arrangement
    R = (1/4) g^{mn} (T_00m T_11n - T_01m T_01n) returns -K under the
    orientation Brioschi fixes, hence the leading minus below.
    """
    if len(coords) != 2:
        raise ValueError("cumulant curvature needs exactly 2 coordinates")
    g = fisher_metric(Z, coords)
    T = cubic_tensor(Z, coords)
    gi = g.inv()
    R = sp.Rational(1, 4) * sum(
        gi[m, n] * (T[(0, 0, m)] * T[(1, 1, n)] - T[(0, 1, m)] * T[(0, 1, n)])
        for m in range(2) for n in range(2))
    return sp.cancel(-R / g.det())


def fisher_pack(Z, coords):
    """Convenience: (E, F, G, K) for a two-parameter family."""
    a, b = coords
    psi = sp.log(Z)
    E = D_log(D_log(psi, a), a)
    F = D_log(D_log(psi, a), b)
    G = D_log(D_log(psi, b), b)
    return E, F, G, brioschi_curvature(Z, coords)


def path_independence(Z, coords, theta_A, theta_B, paths, dps=30):
    """POTENTIAL PATH INDEPENDENCE gate.

    d psi = sum_e <N_e> d theta_e is an EXACT differential, so the line
    integral of the expectation field depends only on the endpoints, and every
    closed loop integrates to zero.  This is the statistical-mechanical form of
    the conservative electric field.

    Parameters
    ----------
    paths : list of callables t in [0,1] -> list of theta values.  A path whose
            endpoints coincide is treated as a loop (target 0).

    Returns list of (value, target, residual).  Raises on failure.
    """
    import mpmath as mp
    mp.mp.dps = dps
    psi_f = sp.lambdify(coords, sp.log(Z), "mpmath")
    means = [sp.lambdify(coords, D_log(sp.log(Z), X), "mpmath") for X in coords]
    nvar = len(coords)

    def integrate(path):
        def integrand(t):
            th = path(t)
            dth = [mp.diff(lambda s, i=i: path(s)[i], t) for i in range(nvar)]
            w = [mp.e ** c for c in th]
            return sum(means[i](*w) * dth[i] for i in range(nvar))
        return mp.quad(integrand, [0, 1])

    dpsi = (psi_f(*[mp.e ** c for c in theta_B])
            - psi_f(*[mp.e ** c for c in theta_A]))
    out = []
    tol = mp.mpf(10) ** (-(dps - 8))
    for path in paths:
        start, end = path(0), path(1)
        loop = all(abs(mp.mpf(start[i]) - mp.mpf(end[i])) < tol
                   for i in range(nvar))
        target = mp.mpf(0) if loop else dpsi
        val = integrate(path)
        res = abs(val - target)
        assert res < tol, f"path independence FAILED: {val} vs {target}"
        out.append((val, target, res))
    return out
