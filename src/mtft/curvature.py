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
    finite support (EXACT).  The K = 1/4 locks are now a THEOREM
    rather than an observation (studies/arithmetic_area_geometry.py):

      RIGIDITY THEOREM.  A two-dimensional discrete exponential family
      whose support consists of any number (>= 2) of statistic points
      on a single affine line, together with EXACTLY ONE point off
      that line, with arbitrary positive base weights, has K = 1/4
      identically on the whole parameter manifold.  Proof: normalize
      the line to y = 0 and the lone point to (0, 1); then
      Z = A(u) + e^v, the metric diagonalizes in (u, q) with
      q = v - log A(u), and the warped form is
      ds^2 = dy^2 + cos^2(y/2) d(rho)^2, whence K = -f''/f = 1/4.

    The hypothesis is sharp: zero off-line points makes the two
    statistics affinely dependent and det g vanishes; two off-line
    points destroys the lock.  Use rigidity_class() to test a support.
    {1,2,3,4} is in the class because w_4 = (1/2) log 4 gives
    X_4 = 2 X_2 EXACTLY, so {1,2,4} are collinear and 3 is the lone
    off-line point; the theorem was confirmed OUT OF SAMPLE on
    {1,2,4,8} and {1,2,4,16}.

    COLD LAW, CORRECTED.  The cold core is {1, 2, 3, 5}, NOT the first
    six integers: 6/5 is 36/30, the squared metric triangle (1,2,3)
    against the leading curvature triangle (2,3,5).  The intermediate
    candidate (2,3,4) with product 24 would give a faster (3/2)^beta
    mode, but its factor cross(X_2, X_4) vanishes by the same
    collinearity — a dyadic shield.  Perturbing w_4 by delta restores
    that mode with amplitude 3 delta / log 3.

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
    "gaussian_family_curvature", "hagedorn_slope", "cold_amplitude",
    "rigidity_class", "CURVATURE",
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


def rigidity_class(points, tol=1e-30):
    """Is a planar support in the K = 1/4 rigidity class?

    Returns the number of points lying off the affine line through the
    support's dominant direction, or None when fewer than three points
    are supplied.  The class is EXACTLY: at least two points on one
    affine line, plus exactly ONE point off it (see the theorem in the
    module docstring).  Zero off-line points is degenerate — the two
    statistics are affinely dependent and det g vanishes identically.
    """
    pts = [(mp.mpf(x), mp.mpf(y)) for x, y in points]
    if len(pts) < 3:
        return None
    best = None
    for i in range(len(pts)):
        for j in range(i + 1, len(pts)):
            dx = pts[j][0] - pts[i][0]
            dy = pts[j][1] - pts[i][1]
            n = 0
            for k in range(len(pts)):
                cr = (dx * (pts[k][1] - pts[i][1])
                      - dy * (pts[k][0] - pts[i][0]))
                if abs(cr) > mp.mpf(tol):
                    n += 1
            best = n if best is None else min(best, n)
    return best


def _curvature_at_dps(beta, atoms, dps):
    """Single-precision evaluation; see finite_atom_curvature."""

    old = mp.dps
    mp.dps = dps
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


def finite_atom_curvature(beta, atoms=(1, 2, 3), dps=None,
                          rtol="1e-20", max_doublings=6):
    """K(beta) for the ensemble restricted to a finite atom set.

    EXACT (arbitrary precision), with ADAPTIVE precision: the value is
    recomputed at doubled precision until two successive evaluations
    agree to `rtol`, so the caller cannot silently receive a
    cancellation artifact.

    CORRECTION CC-03 (append-only).  Earlier versions took a fixed
    floor of 0.7*beta + 45 digits and honoured a caller's `dps`
    override even when it was lower.  That is unsafe, and not merely
    because of the override: the required precision is
    SUPPORT-dependent, not a function of beta alone.  When the support
    contains a COLLINEAR TRIPLE the leading metric triangle vanishes,
    det g drops to the next triangle product, and the cancellation
    deepens far past any beta-based estimate — for {1,2,4,8} at
    beta = 80 the old floor of 101 digits returned ~-2e+19 instead of
    1/4, and about 200 digits are needed.  `dps` is now a STARTING
    precision that may only be raised.

    Rigidity locks: (1,2,3) and (1,2,3,4) give exactly 1/4 for every
    beta, by the theorem in the module docstring; atom 5 flips the
    sign; the cold core is {1,2,3,5}.
    """
    d = max(dps or 0, 50, int(0.7 * float(beta)) + 45)
    prev = None
    for _ in range(max_doublings):
        cur = _curvature_at_dps(beta, atoms, d)
        if prev is not None:
            if cur == 0 and prev == 0:
                return cur
            scale = max(abs(cur), abs(prev))
            if scale > 0 and abs(cur - prev) <= scale * mp.mpf(rtol):
                return cur
        prev = cur
        d *= 2
    raise ValueError(
        f"finite_atom_curvature did not stabilise for atoms={atoms} "
        f"at beta={beta} within {max_doublings} doublings; the support "
        f"may contain a near-degenerate triangle")


def cold_amplitude():
    """The exact cold amplitude c in K(beta) ~ -c (6/5)^beta.

    Closed form in log 2, log 3, log 5, from the leading curvature
    triangle (2, 3, 5) against the squared metric triangle (1, 2, 3):

        c = (9 L5^2)/(25 L2 L3) [1 - (9/5) L5/L3 + (4/5) L5/L2]
          = 0.27012646530542495706433719670365...

    CORRECTION CC-04 (append-only).  v0.14.0 ledgered
    0.270126465305424759517602, which is wrong in the 16th digit: it
    was extracted from the SIX-atom model at beta = 200, which still
    carries an atom-6 contamination of relative size (5/6)^beta —
    1.46e-16 at that beta, against an observed error of 1.98e-16.  The
    retracted value is preserved in CURVATURE for the record.
    """
    L2, L3, L5 = mp.log(2), mp.log(3), mp.log(5)
    return (9 * L5 ** 2) / (25 * L2 * L3) * (
        1 - mp.mpf(9) / 5 * L5 / L3 + mp.mpf(4) / 5 * L5 / L2)


CURVATURE = {
    "K_2p5": "0.559191364467921087",
    "K_3p5": "0.9978851502096039",
    "hagedorn_slope": "0.423657463797093480081718158187",
    "beta_star": "4.593591164956",
    "K_star": "1.19569598199193852905",
    "beta_zero": "8.8565170425",
    "cold_rate": "6/5 = 36/30 (squared metric triangle 1*2*3 over "
                 "leading curvature triangle 2*3*5)",
    "cold_core": "{1, 2, 3, 5}",
    "cold_amplitude": "0.27012646530542495706433719670365",
    "cold_amplitude_RETRACTED_CC04": "0.270126465305424759517602",
    "K_three_atom": "1/4 identically",
    "K_four_atom": "1/4 identically (rigidity theorem: {1,2,4} "
                    "collinear since X_4 = 2 X_2, and 3 is the lone "
                    "off-line point)",
    "K_12346_limit": "1.3549368866",
}
