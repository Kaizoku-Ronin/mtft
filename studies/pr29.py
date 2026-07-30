# --- stage 5 (Integration Plan v0.1): re-pointed to the mtft package. ---
# internal() delegates to mtft.chain.internal (certified spectrally
# identical: g to 9e-15, full H(u) spectrum to 2.5e-14; B differs only
# by gauge inside near-degenerate blocks).  Where a local gsq() remains
# it is retained VERBATIM: it evaluates complex u in the f64 backend,
# which mtft.ep deliberately does not offer (PR-20 floors) -- S5-1.
# ------------------------------------------------------------------------
#!/usr/bin/env python3
"""
pr29.py — PR-29(a): the repulsion channel, and K_4 direct
==========================================================
MIT License — (c) 2026 Roger Tano — MTFT Research Program

AZ-F1 accepted: the enhancement factor is 17.2330613, not 17.233083
(seventh-digit transcription slip; both operands verified).

METHOD.  Rather than deriving X_rep separately in H-space — where the
denominators are u-dependent and the object is a crossing condition, not
a spectrum — extract K_4 DIRECTLY by the same monomial method, which
captures all three channels at once.

Isolate the q^{2 kappa} sector exactly: keep only the (4,5), (5,6),
(4,6) couplings as symbols R, S, Q and set every other coupling to zero.
Level 3 (the crossing's lower member) is then exactly decoupled, so

    eps_3(u) = g_3 - u B_33   with B_33 = e^{-g_3} exactly,
    eps_4(u) = lowest branch of the coupled {4,5,6} block,

and the crossing eps_3 = eps_4 gives |u_cross|.  Every ingredient — the
g's from the block's eigenvalues, B from e^{-h/2}(V^T M V)e^{-h/2}, and
the off-diagonal repulsion inside the block — is present, so the three
monomial coefficients of |u_cross| - L sum to K_4 with no channel
decomposition required and no closure equation assumed.

Predicted count: 3 (Q^2, RSQ, R^2S^2), same convention as PR-28.
Target: K_4 = -54.7822 (Kimi, clean extraction); X_rep = -729.31.
"""
from __future__ import annotations
import mpmath as mp

mp.mp.dps = 60
OUT = open("pr29_run.txt", "w", buffering=1)


def say(s):
    print(s, flush=True)
    OUT.write(s + "\n")


def rho(n):
    return mp.log(n) / n ** 3


r2 = rho(2)
r3, r4, r5, r6 = rho(5), rho(6), rho(7), rho(8)   # levels 3,4,5,6


def block(R, S, Q):
    """{4,5,6} sub-block: g's, B's from the symbolised T."""
    t45 = mp.sqrt(r4 * r5) * R
    t56 = mp.sqrt(r5 * r6) * S
    t46 = mp.sqrt(r4 * r6) * Q
    T = mp.matrix([[r4, t45, t46], [t45, r5, t56], [t46, t56, r6]])
    M = mp.matrix([[1, R, Q], [R, 1, S], [Q, S, 1]])
    lam, V = mp.eigsy(T)
    order = sorted(range(3), key=lambda j: -lam[j])
    lam = [lam[j] for j in order]
    Vd = mp.matrix(3, 3)
    for c, j in enumerate(order):
        for row in range(3):
            Vd[row, c] = V[row, j]
    g = [mp.log(r2 / lam[j]) for j in range(3)]
    Mt = Vd.T * M * Vd
    B = mp.matrix(3, 3)
    for a in range(3):
        for b in range(3):
            B[a, b] = mp.e ** (-(g[a] + g[b]) / 2) * Mt[a, b]
    return g, B


g3 = mp.log(r2 / r3)
B33 = mp.e ** (-g3)


def eps4(u, g, B):
    H = mp.matrix(3, 3)
    for a in range(3):
        for b in range(3):
            H[a, b] = (g[a] if a == b else 0) - u * B[a, b]
    return min(mp.eigsy(H, eigvals_only=True))


def ucross(R, S, Q):
    g, B = block(R, S, Q)
    f = lambda u: (g3 - u * B33) - eps4(u, g, B)
    lo, hi = mp.mpf(-40), mp.mpf(-1)
    flo = f(lo)
    for _ in range(300):
        mid = (lo + hi) / 2
        if f(mid) * flo > 0:
            lo = mid
        else:
            hi = mid
        if hi - lo < mp.mpf(10) ** (-55):
            break
    return -(lo + hi) / 2


L = ucross(mp.mpf(0), mp.mpf(0), mp.mpf(0))


def dev(R, S, Q):
    return ucross(R, S, Q) - L


if __name__ == "__main__":
    say("=" * 92)
    say("  PR-29(a) — K_4 DIRECT, ALL THREE CHANNELS AT ONCE")
    say("=" * 92)
    say(f"\n  [AZ-F1] factor corrected: 17.2330613 (was printed 17.233083)")
    say(f"  unperturbed crossing L = {mp.nstr(L, 14)}   "
        f"(u_34(inf) = 8.3166465032)")
    say(f"  convention: 3 monomials, as PR-28\n")
    # explicit finite differences: mp.diff's internal steps fall below
    # the root-finder's output resolution, which is why the first pass
    # returned exact zeros.
    h = mp.mpf(10) ** (-6)
    a = (dev(0, 0, h) + dev(0, 0, -h)) / (2 * h ** 2)
    b = sum(sx * sy * sz * dev(sx * h, sy * h, sz * h)
            for sx in (1, -1) for sy in (1, -1) for sz in (1, -1)) \
        / (8 * h ** 3)
    hh = mp.mpf(10) ** (-4)
    d2 = lambda y: (dev(hh, y, 0) - 2 * dev(0, y, 0) + dev(-hh, y, 0)) \
        / hh ** 2
    c = (d2(hh) - 2 * d2(0) + d2(-hh)) / (4 * hh ** 2)   # BA: /4, not /2
    say(f"    [Q^2]     {mp.nstr(a, 12)}")
    say(f"    [R S Q]   {mp.nstr(b, 12)}")
    say(f"    [R^2 S^2] {mp.nstr(c, 12)}")
    tot = a + b + c
    say(f"    K_4 = sum {mp.nstr(tot, 12)}")
    say(f"    target -54.7822 ; ratio "
        f"{mp.nstr(tot / mp.mpf('-54.7822'), 10)}")
    say("\n" + "=" * 92)
    OUT.close()
