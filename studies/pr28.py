# --- stage 5 (Integration Plan v0.1): re-pointed to the mtft package. ---
# internal() delegates to mtft.chain.internal (certified spectrally
# identical: g to 9e-15, full H(u) spectrum to 2.5e-14; B differs only
# by gauge inside near-degenerate blocks).  Where a local gsq() remains
# it is retained VERBATIM: it evaluates complex u in the f64 backend,
# which mtft.ep deliberately does not offer (PR-20 floors) -- S5-1.
# ------------------------------------------------------------------------
#!/usr/bin/env python3
"""
pr28.py — PR-28(a): the eigenvector coefficient, by monomial extraction
=======================================================================
MIT License — (c) 2026 Roger Tano — MTFT Research Program

CONVENTION PINNED BEFORE ARITHMETIC (Addendum AY's flag).  I count
MONOMIALS in the coupling symbols, not bond walks.  With

    T_45 = sqrt(rho_4 rho_5) R,  T_56 = sqrt(rho_5 rho_6) S,
    T_46 = sqrt(rho_4 rho_6) Q,  M_45 = R, M_56 = S, M_46 = Q,

R = r^kappa, S = r'^kappa, Q = q^kappa and r r' = q, the monomials
reaching q^{2 kappa} are exactly those with R^a S^b Q^c, a + c = 2,
b = a — i.e. THREE:  Q^2,  R S Q,  R^2 S^2.  Predicted count: 3.
(The walk count is larger — walks of the same monomial differ by
denominator structure — which is the ambiguity AY flagged; this is the
convention that makes the count well-defined.)

METHOD.  The three-level block {4,5,6} carries every contributing walk
(AY's exhaustive enumeration), so the coefficient is the sum of three
multivariate Taylor coefficients of the EXACT 3x3 problem:

    coeff = [Q^2] + [R S Q] + [R^2 S^2]  of  (V^T M V)_44 - 1,

each taken at R = S = Q = 0.  Nothing is truncated: the 3x3 eigenproblem
is solved exactly at each node and the derivatives are exact to the
working precision.

Same extraction applied to delta-g_4 cross-checks PR-26's closed form.
"""
from __future__ import annotations
import mpmath as mp

mp.mp.dps = 60
OUT = open("pr28_run.txt", "w", buffering=1)


def say(s):
    print(s, flush=True)
    OUT.write(s + "\n")


def rho(n):
    return mp.log(n) / n ** 3


r4, r5, r6 = rho(6), rho(7), rho(8)          # levels 4,5,6 = sites 6,7,8


def blocks(R, S, Q):
    t45 = mp.sqrt(r4 * r5) * R
    t56 = mp.sqrt(r5 * r6) * S
    t46 = mp.sqrt(r4 * r6) * Q
    T = mp.matrix([[r4, t45, t46], [t45, r5, t56], [t46, t56, r6]])
    M = mp.matrix([[1, R, Q], [R, 1, S], [Q, S, 1]])
    return T, M


def pick(T):
    """eigenpair continuing from level 4 (largest rho)."""
    lam, V = mp.eigsy(T)
    k = max(range(3), key=lambda j: abs(V[0, j]))
    return lam[k], mp.matrix([V[0, k], V[1, k], V[2, k]])


def fV(R, S, Q):
    T, M = blocks(R, S, Q)
    _, v = pick(T)
    return (v.T * M * v)[0] / (v.T * v)[0] - 1


def fg(R, S, Q):
    T, _ = blocks(R, S, Q)
    lam, _ = pick(T)
    return -(lam - r4) / r4                    # dg_4 = -dlambda/lambda


def coeffs(f):
    c_Q2 = mp.diff(lambda q: f(0, 0, q), 0, 2) / 2
    c_RSQ = mp.diff(lambda a, b, c: f(a, b, c), (0, 0, 0), (1, 1, 1))
    c_R2S2 = mp.diff(lambda a, b: f(a, b, 0), (0, 0), (2, 2)) / 4
    return c_Q2, c_RSQ, c_R2S2


if __name__ == "__main__":
    say("=" * 92)
    say("  PR-28(a) — EIGENVECTOR COEFFICIENT BY MONOMIAL EXTRACTION")
    say("=" * 92)
    say("\n  convention pinned in advance: monomials, not walks.  "
        "predicted count = 3")
    say("  monomials at q^{2k}:  Q^2,  R*S*Q,  R^2*S^2\n")

    say("  [cross-check] delta-g_4 by the same extraction")
    a, b, c = coeffs(fg)
    say(f"    [Q^2]     {mp.nstr(a, 12)}")
    say(f"    [R S Q]   {mp.nstr(b, 12)}")
    say(f"    [R^2 S^2] {mp.nstr(c, 12)}")
    say(f"    sum       {mp.nstr(a + b + c, 12)}")
    say(f"    PR-26 closed form -9.601605567 ; ratio "
        f"{mp.nstr((a + b + c) / mp.mpf('-9.601605567'), 10)}\n")

    say("  [PR-28(a)] delta-(V^T M V)_44")
    a2, b2, c2 = coeffs(fV)
    say(f"    [Q^2]     {mp.nstr(a2, 12)}")
    say(f"    [R S Q]   {mp.nstr(b2, 12)}")
    say(f"    [R^2 S^2] {mp.nstr(c2, 12)}")
    tot = a2 + b2 + c2
    say(f"    sum       {mp.nstr(tot, 12)}")
    say(f"    Kimi target +47.25156 ; ratio "
        f"{mp.nstr(tot / mp.mpf('47.25156'), 10)}")
    say(f"    factor vs naive 2.74191: {mp.nstr(tot / mp.mpf('2.74191'), 8)}"
        f"  (target 17.2330)")
    say("\n" + "=" * 92)
    OUT.close()
