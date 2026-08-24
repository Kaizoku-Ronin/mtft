"""mtft.exception_spectrum — the Exception-Spacing Curvature Law (v0.19.0).

Fisher curvature of the marked primon ensemble Z(beta, lambda) =
sum n^-beta e^{lambda F(n)} for a binary arithmetic property F, via the
exact cumulant/Brioschi reduction (third cumulants only).  Certified in
the 2026-08-24 arc (curvature gate + certificates v3-v6):

  * rigidity: line-plus-one-marked-point => K = 1/4 exactly.
  * law: F(1) = F(2) = 0, first exceptions e1 < e2, rho = 2 e1/e2:
        K = 1/4 - C rho^beta + o(rho^beta),
        C = log^2(e2/e1) / (4 log^2 2).
    Phases: rho < 1 spherical (K -> 1/4); rho = 1 critical (K -> 0,
    exact product-family flatness); rho > 1 hyperbolic (K -> -inf).
  * line-atom lemma (rigid base) and flat-base lemma, with the
    lambda-dressing V -> V e^lambda; complete fibers are curvature-
    invisible at finite weight.
  * signed defect spectrum with EXACT cancellation when e_j = e1 * m;
    the inverse problem hears the support modulo multiplicative fibers.

Truncation rule (hard-won): drop marked atom m only when the RELATIVE
scale (2 e1/m)^beta is below target; absolute weight is the wrong
criterion and silently deletes the phenomenon.
"""
from __future__ import annotations

from mpmath import exp, ln, matrix, mpf, det

__all__ = [
    "K_from_AB", "K_atoms", "K_marked_set", "two_exception_C", "rho",
    "classify_phase", "line_atom_variation", "flat_base_variation",
    "defect_spectrum",
]


def K_from_AB(A, B, lam=mpf(0)):
    """Gaussian curvature of the (beta, lambda) family from raw moments.

    A[a] = sum over ALL atoms of w * x^a; B[a] = same over marked atoms;
    a = 0..3.  Exact third-cumulant Brioschi reduction.
    """
    el = exp(lam)
    Z = A[0] + (el - 1) * B[0]
    m = {}
    for a in range(4):
        m[(a, 0)] = (A[a] + (el - 1) * B[a]) / Z
        m[(a, 1)] = el * B[a] / Z
    m10, m01 = m[(1, 0)], m[(0, 1)]
    k20 = m[(2, 0)] - m10 ** 2
    k11 = m[(1, 1)] - m10 * m01
    k02 = m[(0, 1)] - m01 ** 2
    k30 = m[(3, 0)] - 3 * m[(2, 0)] * m10 + 2 * m10 ** 3
    k21 = m[(2, 1)] - m[(2, 0)] * m01 - 2 * m[(1, 1)] * m10 \
        + 2 * m10 ** 2 * m01
    k12 = m[(1, 1)] - m[(0, 1)] * m10 - 2 * m[(1, 1)] * m01 \
        + 2 * m01 ** 2 * m10
    k03 = m[(0, 1)] - 3 * m[(0, 1)] * m01 + 2 * m01 ** 3
    M1 = matrix([[0, k30 / 2, k21 / 2], [k12 / 2, k20, k11],
                 [k03 / 2, k11, k02]])
    M2 = matrix([[0, k21 / 2, k12 / 2], [k21 / 2, k20, k11],
                 [k12 / 2, k11, k02]])
    return (det(M1) - det(M2)) / (k20 * k02 - k11 ** 2) ** 2


def K_atoms(xs, marks, wts, lam=mpf(0)):
    A = [sum(w * x ** a for x, w in zip(xs, wts)) for a in range(4)]
    B = [sum(w * x ** a for x, mk, w in zip(xs, marks, wts) if mk)
         for a in range(4)]
    return K_from_AB(A, B, lam)


def K_marked_set(beta, marked, e1, rel_target=mpf(10) ** -40, lam=mpf(0)):
    """Curvature of the full-integer ensemble with the given marked
    iterable (increasing), honoring the RELATIVE truncation rule."""
    from mpmath import zeta
    b = mpf(beta)
    A = [(-1) ** a * zeta(b, derivative=a) for a in range(4)]
    B = [mpf(0)] * 4
    kept = 0
    for n in marked:
        if kept > 1 and (mpf(2 * e1) / n) ** b < rel_target:
            break
        t = mpf(n) ** (-b)
        lg = ln(n)
        for a in range(4):
            B[a] += t * lg ** a
        kept += 1
    return K_from_AB(A, B, lam)


def two_exception_C(e1, e2):
    return ln(mpf(e2) / e1) ** 2 / (4 * ln(2) ** 2)


def rho(e1, e2):
    return mpf(2 * e1) / e2


def classify_phase(e1, e2):
    r = rho(e1, e2)
    if r < 1:
        return "spherical (K -> 1/4)"
    if r == 1:
        return "critical (K -> 0, product-family flat)"
    return "hyperbolic (K -> -infinity)"


def line_atom_variation(line_xw, b, v, c):
    """dK/dw at w = 0 for a second marked atom at c over the RIGID base:
    baseline cloud line_xw = [(x_i, q_i)], first marked (b, v).
    Certificates v4-v5 (machine-exact)."""
    L0 = sum(w for _, w in line_xw)
    L1 = sum(w * x for x, w in line_xw)
    L2 = sum(w * x * x for x, w in line_xw)
    return -(L0 + v) ** 2 * (c - b) ** 2 / (4 * v * (L0 * L2 - L1 ** 2))


def flat_base_variation(cloud, V, r, lam=mpf(0), marked=False):
    """dK/dw at w = 0 over the FLAT base (cloud doubled with weights V q_i
    at translation t; t drops out).  lambda-dressed: V -> V e^lambda in
    the prefactor; the marked/unmarked raw-weight ratio is -1/V at any
    lambda.  Certificate v6."""
    Q = sum(q for _, q in cloud)
    mu = sum(q * x for x, q in cloud) / Q
    s2 = sum(q * (x - mu) ** 2 for x, q in cloud) / Q
    Ve = V * exp(lam)
    base = (1 - Ve) / (1 + Ve) * ((r - mu) ** 2 - s2) / (4 * Q * s2)
    return -base / V if marked else base


def defect_spectrum(marked, e1, rmax):
    """Surviving channels of the multiplicative defect measure
    mu_F = sum_{unmarked} delta_n - sum_{marked} delta_{n/e1}, as
    (r, sign, amplitude) with amplitude log^2 r / (4 log^2 2); channels
    with e_j = e1 * m cancel exactly (certificates v5-v6)."""
    marked = set(marked)
    chans = {}
    n = 3
    while n <= rmax:
        if n not in marked:
            chans[mpf(n)] = chans.get(mpf(n), 0) + 1
        n += 1
    for e in sorted(marked):
        r = mpf(e) / e1
        if r > rmax:
            break
        if r > 2:
            chans[r] = chans.get(r, 0) - 1
    out = []
    for r, s in sorted(chans.items()):
        if s:
            out.append((r, s, s * ln(r) ** 2 / (4 * ln(2) ** 2)))
    return out
