# --- stage 5 (Integration Plan v0.1): re-pointed to the mtft package. ---
# internal() delegates to mtft.chain.internal (certified spectrally
# identical: g to 9e-15, full H(u) spectrum to 2.5e-14; B differs only
# by gauge inside near-degenerate blocks).  Where a local gsq() remains
# it is retained VERBATIM: it evaluates complex u in the f64 backend,
# which mtft.ep deliberately does not offer (PR-20 floors) -- S5-1.
# ------------------------------------------------------------------------
#!/usr/bin/env python3
"""
pr23.py — PR-23(b): the amplitude derived, term by term
========================================================
MIT License — (c) 2026 Roger Tano — MTFT Research Program
Prints to stdout AND pr23_run.txt.

Addendum AT put the full partition on record (pair (0,1), kappa = 56):
    delta-g        = -7.2119335e-15
    delta-(VtMV)   = +3.9999846e-14
    B_12 repulsion = -1.4385005e-14
    B_01           = -9.0e-18
    sum            =  1.8402817e-14
and AT.2 showed all three channels are one coupling's bookkeeping, the
site ratio r = (i+3)/(i+4).  The bar (AS.3): reproduce the CONSTANTS AND
THE SPLIT, not the sum alone.

THE DERIVATION.  Write a = g_i, b = g_{i+1}, c = g_{i+2} (limit gaps),
p_j = B_jj -> e^{-g_j}, D = e^{-a} - e^{-b}, m = b - a, u0 = m/D, and
R2 = r^{2 kappa}.  With rho_j the weight at site j+2:

 (1) T-EIGENVALUE.  Second order in T: d(lambda_{i+1}) =
     rho_{i+1} rho_{i+2} R2/(rho_{i+1} - rho_{i+2}), and since
     g = log(lambda_0/lambda), dg_{i+1} = -rho_{i+2} R2/(rho_{i+1} -
     rho_{i+2}).  Enters u through d|u|/db = (D - m e^{-b})/D^2.

 (2) T-EIGENVECTOR.  V = I + W with W_{i+2,i+1} = T_{i+1,i+2}/
     (lambda_{i+1} - lambda_{i+2}), so (V^T M V)_{i+1,i+1} - 1 =
     2 sqrt(rho_{i+1} rho_{i+2}) R2/(rho_{i+1} - rho_{i+2}).  This
     multiplies e^{-b} in B_{i+1,i+1}, and enters u through
     d|u|/dp_{i+1} = m/D^2.  NO delta-lambda computation contains it.

 (3) B_12 REPULSION.  Level i+2 sits above the crossing, pushing the
     upper member down by u0^2 B_{i+1,i+2}^2/(eps_{i+2} - eps_{i+1}),
     which advances the crossing: divide by D.  NEGATIVE, as PR-22 §3
     established and AT.2 confirmed.

Every coefficient is closed form; nothing is fitted.
"""
from __future__ import annotations
import mpmath as mp

mp.mp.dps = 40
OUT = open("pr23_run.txt", "w", buffering=1)


def say(s=""):
    print(s, flush=True)
    OUT.write(s + "\n")


def rho(site):
    return mp.log(site) / site ** 3


def gaps(nlev=5):
    return [mp.log(rho(2) / rho(2 + j)) for j in range(nlev)]


def terms(i, kappa):
    g = gaps(i + 3)
    a, b, c = g[i], g[i + 1], g[i + 2]
    ea, eb, ec = mp.e ** (-a), mp.e ** (-b), mp.e ** (-c)
    D, m = ea - eb, b - a
    u0 = m / D
    r = mp.mpf(i + 3) / (i + 4)
    R2 = r ** (2 * kappa)
    r1, r2 = rho(i + 3), rho(i + 4)          # sites of levels i+1, i+2

    dudb = (D - m * eb) / D ** 2
    dg = -r2 * R2 / (r1 - r2)
    T1 = dudb * dg

    dMt = 2 * mp.sqrt(r1 * r2) * R2 / (r1 - r2)
    T2 = (m / D ** 2) * (eb * dMt)

    Bc = mp.e ** (-(b + c) / 2) * r ** kappa
    eps_gap = (c + u0 * ec) - (b + u0 * eb)
    T3 = (-u0 ** 2 * Bc ** 2 / eps_gap) / D

    return T1, T2, T3, T1 + T2 + T3


if __name__ == "__main__":
    say("=" * 96)
    say("  PR-23(b) — THE AMPLITUDE DERIVED, TERM BY TERM")
    say("=" * 96)
    sealed = {
        (0, 56): dict(dg=mp.mpf('-7.2119335e-15'),
                      dV=mp.mpf('3.9999846e-14'),
                      B12=mp.mpf('-1.4385005e-14'),
                      tot=mp.mpf('1.8402817e-14')),
    }
    say("\n  AT.2 sealed partition, pair (0,1), kappa = 56 — "
        "derived vs measured")
    T1, T2, T3, tot = terms(0, mp.mpf(56))
    s = sealed[(0, 56)]
    for lbl, der, meas in (("(1) T-eigenvalue  dg", T1, s['dg']),
                           ("(2) T-eigenvector dV", T2, s['dV']),
                           ("(3) B_12 repulsion  ", T3, s['B12']),
                           ("    TOTAL           ", tot, s['tot'])):
        say(f"    {lbl}  derived {mp.nstr(der, 8):>15}   measured "
            f"{mp.nstr(meas, 8):>15}   ratio {mp.nstr(der/meas, 6)}")

    say("\n  the same derivation at other (pair, kappa) — measured "
        "totals from the record")
    known = {(0, 72): mp.mpf('1.8486e-18'), (1, 88): mp.mpf('3.913976586e-17'),
             (1, 112): mp.mpf('8.72938e-22'), (2, 96): mp.mpf('4.83822e-15'),
             (2, 112): mp.mpf('1.42133e-17')}
    say(f"    {'pair':>6} {'kappa':>6} {'derived total':>16} "
        f"{'measured total':>16} {'ratio':>9}")
    for (i, kap), meas in sorted(known.items()):
        _, _, _, tt = terms(i, mp.mpf(kap))
        say(f"    {f'({i},{i+1})':>6} {kap:6d} {mp.nstr(tt, 8):>16} "
            f"{mp.nstr(meas, 8):>16} {mp.nstr(tt/meas, 6):>9}")

    say("\n  amplitudes A_i (dev = A_i r^{2 kappa}), derived vs AS.2 "
        "sealed values")
    seal_A = {0: mp.mpf('1.811457162'), 1: mp.mpf('4.454780532'),
              2: mp.mpf('7.749863218')}
    for i in (0, 1, 2):
        kap = mp.mpf(80)
        _, _, _, tt = terms(i, kap)
        r = mp.mpf(i + 3) / (i + 4)
        A = tt / r ** (2 * kap)
        say(f"    A_{i}: derived {mp.nstr(A, 10):>14}   sealed "
            f"{mp.nstr(seal_A[i], 10):>14}   ratio "
            f"{mp.nstr(A/seal_A[i], 8)}")
    say("\n" + "=" * 96)
    OUT.close()
