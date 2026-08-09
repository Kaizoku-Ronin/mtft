#!/usr/bin/env python3
"""
Curvature of the Tano statistical manifold
==========================================

MIT License — Copyright (c) 2026 Roger Tano
See LICENSE file for full terms.

The two-parameter exponential family

    p_{beta,lambda}(n) = n^-beta e^{lambda w_n} / Z(beta, lambda),
    w_n = sum_{d|n} (log d)/d,

is a statistical manifold for beta > 1 and all real lambda (w_n grows
like (log log n)^2 on champions, so Z converges).  Its Fisher metric
and Amari-Chentsov tensor come from mtft.moments; this module supplies
the Levi-Civita curvature.

THE CANCELLATION THEOREM (EXACT).  In (beta, lambda) coordinates the
metric is a Hessian, g_ij = d_i d_j psi with psi = log Z, so
Gamma_{ij,k} = (1/2) d_i d_j d_k psi and the second-derivative block
of the Brioschi form is

    -1/2 d_vv E + d_uv F - 1/2 d_uu G = -1/2 k4 + k4 - 1/2 k4 = 0,

identically, because all three terms are the SAME totally symmetric
fourth derivative of psi (the fourth cumulant).  THE FOURTH CUMULANT
CONTRIBUTES NOTHING TO THE CURVATURE OF ANY EXPONENTIAL FAMILY.
Curvature therefore needs only the certified second and third
cumulants.  Verified analytically, on the Gaussian family, and by
direct finite-difference measurement (block ~ 1e-11 against O(1)).

Contents and epistemic classes
------------------------------
1.  brioschi(...) — the reduced Brioschi form (EXACT given inputs);
    convention locked by gaussian_family_curvature() = -1/2.

2.  gaussian_curvature(beta) — K along lambda = 0 (Pr).  The profile
    is a POSITIVE DOME on (1, beta_0) and a hyperbolic cold tail
    beyond:
        leaves the Hagedorn wall flat, K -> 0+ with slope
            A = (zeta''(2) kappa3_cold - kappa_wwl_cold chi_cold)
                / (2 chi_cold^2)  =  0.4236574637970935 ,
        summits at beta* = 4.593591164956, K* = 1.1956959819919385,
        crosses zero at beta_0 = 8.8565170425,
        and dives as K ~ -c (6/5)^beta with c = 0.2701264653054.

3.  finite_atom_curvature(beta, atoms) — the same geometry on a
    finite support (EXACT).  Rigidity locks: {1,2,3} and {1,2,3,4}
    both give K = 1/4 IDENTICALLY in beta (atom 4 is inert, since
    (l_4, w_4) = 2 (l_2, w_2) makes {1,2,4} collinear and equally
    spaced in statistic space).  Atom 5 flips the sign; with atoms 5
    and 6 present the cold asymptote 6/5 already appears — the deep-
    cold geometry of the Tano ensemble is the geometry of the first
    six integers.

Provenance: studies/curvature_tano_manifold.py (11 gates).  Note that
the sieve route is numerically honest only to beta ~ 18: the Brioschi
numerator cancels about 0.653*beta digits, so use the exact route (or
raise precision) in the cold tail.
"""

from __future__ import annotations

from mpmath import mp

from . import moments as _mom

__all__ = [
    "brioschi", "gaussian_curvature", "finite_atom_curvature",
    "gaussian_family_curvature", "hagedorn_slope", "CURVATURE",
]


def brioschi(E, F, G, Eu, Ev, Fu, Fv, Gu, Gv, block=0):
    """Reduced Brioschi form of a 2D metric.

    E, F, G are the metric components and the rest their first
    derivatives (u = first coordinate, v = second).  `block` is the
    second-derivative combination -1/2 E_vv + F_uv - 1/2 G_uu, which
    vanishes identically for a Hessian metric — leave it at 0 unless
    deliberately testing the cancellation theorem.
    """
    h = mp.mpf("0.5")
    M1 = mp.matrix([[block, h * Eu, Fu - h * Ev],
                    [Fv - h * Gu, E, F],
                    [h * Gv, F, G]])
    M2 = mp.matrix([[0, h * Ev, h * Gu],
                    [h * Ev, E, F],
                    [h * Gu, F, G]])
    return (mp.det(M1) - mp.det(M2)) / (E * G - F * F) ** 2


def gaussian_family_curvature(a=mp.mpf("0.3"), b=mp.mpf("-0.7")):
    """Convention lock: the Gaussian family has K = -1/2 everywhere."""
    psi = lambda x, y: -x * x / (4 * y) + mp.mpf("0.5") * mp.log(-mp.pi / y)
    d = lambda i, j: mp.diff(psi, (mp.mpf(a), mp.mpf(b)), (i, j))
    return brioschi(d(2, 0), d(1, 1), d(0, 2), d(3, 0), d(2, 1),
                    d(2, 1), d(1, 2), d(1, 2), d(0, 3))


def metric_components(beta, **kw):
    """The nine Brioschi inputs at (beta, lambda = 0).

    E = Var(log n), F = -Cov(log n, w), G = chi_w, with first
    derivatives supplied by the third cumulants.  Uses the sign rule
    d_beta^a d_lambda^b psi = (-1)^a kappa(l^a, w^b), so the mixed
    partials Ev == Fu and Fv == Gu are the cumulant identities.
    """
    b = mp.mpf(beta)
    logz = lambda x: mp.log(mp.zeta(x))
    E = mp.diff(logz, b, 2)
    Eu = mp.diff(logz, b, 3)
    F = -mp.diff(mp.zeta, b + 1, 2)
    Fu = -mp.diff(mp.zeta, b + 1, 3)
    Ev = Fu
    chi = lambda x: _mom.weight_susceptibility(x, **kw)
    G = chi(b)
    Gu = mp.diff(chi, b, 1)
    Fv = Gu
    zp = mp.diff(mp.zeta, b + 1, 1)
    T = _mom.weight_second_moment(b, **kw)
    U = _mom.weight_third_moment(b, **kw)
    Gv = U + 3 * T * zp - 2 * zp ** 3
    return [E, F, G, Eu, Ev, Fu, Fv, Gu, Gv]


def gaussian_curvature(beta, **kw):
    """K(beta) along lambda = 0.  Pr.

    Positive on (1, beta_0), zero at beta_0 = 8.8565170425, negative
    beyond with asymptotic law -c (6/5)^beta.  In the cold tail raise
    mp.dps: the Brioschi numerator cancels ~0.653*beta digits.
    """
    return brioschi(*metric_components(beta, **kw))


def hagedorn_slope(**kw):
    """dK/dbeta at the Hagedorn edge, in closed form.

    A = (zeta''(2) kappa3_cold - kappa_wwl_cold chi_cold)
        / (2 chi_cold^2) = 0.423657463797093...

    A pure combination of four constants already certified in
    mtft.moments.COLD; K -> 0+ linearly from positive curvature, so
    the Hagedorn wall is approached FLAT.
    """
    one = mp.mpf(1)
    chi1 = _mom.weight_susceptibility(one, **kw)
    chip1 = mp.diff(lambda x: _mom.weight_susceptibility(x, **kw), one, 1)
    z2 = mp.diff(mp.zeta, mp.mpf(2), 2)
    zp2 = mp.diff(mp.zeta, mp.mpf(2), 1)
    T1 = _mom.weight_second_moment(one, **kw)
    U1 = _mom.weight_third_moment(one, **kw)
    k3 = U1 + 3 * T1 * zp2 - 2 * zp2 ** 3
    return (chip1 * chi1 + z2 * k3) / (2 * chi1 ** 2)


def _weight(n):
    return sum(mp.log(d) / d for d in range(2, n + 1) if n % d == 0)


def finite_atom_curvature(beta, atoms=(1, 2, 3), dps=None):
    """K(beta) for the ensemble restricted to a finite atom set.

    EXACT (arbitrary precision).  Rigidity locks: atoms (1,2,3) and
    (1,2,3,4) both give exactly 1/4 for every beta; atom 5 flips the
    sign; (1,...,6) already carries the 6/5 cold asymptote.
    """
    old = mp.dps
    mp.dps = dps or max(50, int(0.7 * float(beta)) + 45)
    try:
        b = mp.mpf(beta)
        ls = [mp.log(n) for n in atoms]
        ws = [_weight(n) for n in atoms]
        M = len(atoms)
        p = [mp.e ** (-b * ls[i]) for i in range(M)]
        Z = sum(p)
        p = [x / Z for x in p]
        mL = sum(p[i] * ls[i] for i in range(M))
        mW = sum(p[i] * ws[i] for i in range(M))
        cL = [ls[i] - mL for i in range(M)]
        cW = [ws[i] - mW for i in range(M)]
        m = lambda f: sum(p[i] * f(i) for i in range(M))
        out = brioschi(
            m(lambda i: cL[i] ** 2), -m(lambda i: cL[i] * cW[i]),
            m(lambda i: cW[i] ** 2), -m(lambda i: cL[i] ** 3),
            m(lambda i: cL[i] ** 2 * cW[i]),
            m(lambda i: cL[i] ** 2 * cW[i]),
            -m(lambda i: cL[i] * cW[i] ** 2),
            -m(lambda i: cL[i] * cW[i] ** 2),
            m(lambda i: cW[i] ** 3))
    finally:
        mp.dps = old
    return out


CURVATURE = {
    "K_2p5": "0.559191364467921087",
    "K_3p5": "0.9978851502096039",
    "hagedorn_slope": "0.423657463797093480081718158187",
    "beta_star": "4.593591164956",
    "K_star": "1.19569598199193852905",
    "beta_zero": "8.8565170425",
    "cold_rate": "6/5",
    "cold_amplitude": "0.270126465305424759517602",
    "K_three_atom": "1/4 identically",
    "K_four_atom": "1/4 identically (atom 4 inert)",
    "K_12346_limit": "1.3549368866",
}
