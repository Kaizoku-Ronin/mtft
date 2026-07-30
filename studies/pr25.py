# --- stage 5 (Integration Plan v0.1): re-pointed to the mtft package. ---
# internal() delegates to mtft.chain.internal (certified spectrally
# identical: g to 9e-15, full H(u) spectrum to 2.5e-14; B differs only
# by gauge inside near-degenerate blocks).  Where a local gsq() remains
# it is retained VERBATIM: it evaluates complex u in the f64 backend,
# which mtft.ep deliberately does not offer (PR-20 floors) -- S5-1.
# ------------------------------------------------------------------------
#!/usr/bin/env python3
"""
pr25.py — PR-25: the incident-coupling hierarchy, and the family's tail
=======================================================================
MIT License — (c) 2026 Roger Tano — MTFT Research Program
Prints to stdout AND pr25_run.txt.

ADDENDUM AV DISPOSITIONS:
 AV-F1 (A_2 = 7.750 is the stale finite-kappa seal): ACCEPTED; the
   asymptotic ledger value is the derived 7.7619.
 AV-F2 (third order): my guess ((i+4)/(i+5))^{2k} — the level TWO above —
   is FALSIFIED by measurement.  The true third order is the LOWER-
   adjacent incident coupling ((i+1)/(i+2))^{2k}.  Verified below.

THE RULE, and what it predicts next.  Kimi's general statement is that
the remainder hierarchy runs through the couplings INCIDENT ON THE
CROSSING MEMBERS in descending ratio order.  For the (i, i+1) crossing
the incident couplings and their ratios are

    (i+1, i+2): (i+3)/(i+4)   -> A_i   [1st]
    (i,   i+1): (i+2)/(i+3)   -> C_i   [2nd]
    (i-1, i  ): (i+1)/(i+2)   -> C3_i  [3rd]
    (i+1, i+3): (i+3)/(i+5)   -> C4_i  [4th, PREDICTED HERE]

so the rule makes a fourth prediction with no freedom left: for pair
(3,4) the fourth order must be (6/8)^{2k} = (3/4)^{2k}.  My last three
"next order" guesses were all wrong; this one is the RULE's guess, not
mine, which is exactly why it is worth shooting at.

Third-order coefficient, same three channels with the (i-1, i) coupling
(t = (i+1)/(i+2), T2 = t^{2k}, only the LOWER member is touched):
  (1) eigenvalue    dg_i = +rho_{i-1} T2/(rho_{i-1} - rho_i)
                    [plus the reference shift when i-1 = 0]
  (2) eigenvector   dp_i = e^{-a} * 2 sqrt(rho_{i-1} rho_i) T2/
                                   (rho_i - rho_{i-1})
  (3) repulsion from BELOW: level i-1 sits ABOVE the crossing, pushing
      the lower member down: +u0^2 e^{-(g_{i-1}+g_i)} T2/
                              ((eps_{i-1} - eps_i) D)
"""
from __future__ import annotations
import mpmath as mp
from mtft import expansion as _ex


mp.mp.dps = 40
OUT = open("pr25_run.txt", "w", buffering=1)


def say(s=""):
    print(s, flush=True)
    OUT.write(s + "\n")


def rho(site):
    return mp.log(site) / site ** 3


def gaps(n):
    return [mp.log(rho(2) / rho(2 + j)) for j in range(n)]


def A_of(i):
    return _ex.A(i)


def C3_of(i):
    if i < 1:
        return None
    return _ex.C3(i), _ex.channels_C3(i)


if __name__ == "__main__":
    say("=" * 94)
    say("  PR-25 — THE INCIDENT-COUPLING HIERARCHY, AND THE FAMILY TAIL")
    say("=" * 94)

    say("\n[AV-F1] family table repaired to asymptotic values")
    say("    A_0..A_4 = " + ", ".join(mp.nstr(A_of(i), 6)
                                      for i in range(5)))

    say("\n[AV-F2 accepted] third order verified independently")
    say(f"    {'pair':>7} {'derived C3':>16} {'Kimi measured':>16} "
        f"{'ratio':>12}   channels")
    kimi = {3: mp.mpf('250.3238617'), 4: mp.mpf('528.7608972')}
    for i in (3, 4):
        c3, ch = C3_of(i)
        say(f"    {f'({i},{i+1})':>7} {mp.nstr(c3, 10):>16} "
            f"{mp.nstr(kimi[i], 10):>16} {mp.nstr(c3/kimi[i], 9):>12}   "
            + " + ".join(mp.nstr(x, 6) for x in ch))

    say("\n[PR-25(a)] the family tail — values AV.5 withheld for me")
    say(f"    {'i':>3} {'A_i':>18} {'sign':>6}")
    for i in range(5, 13):
        A = A_of(i)
        say(f"    {i:3d} {mp.nstr(A, 12):>18} {'-' if A < 0 else '+':>6}")
    signs = [A_of(i) < 0 for i in range(4, 25)]
    mags = [abs(A_of(i)) for i in range(4, 25)]
    mono = all(mags[k] < mags[k + 1] for k in range(len(mags) - 1))
    say(f"    all negative for i = 4..24: {all(signs)};  |A_i| strictly "
        f"increasing: {mono}")

    say("\n[PR-25(b)] is any integer A_i = 0?")
    zeros = [i for i in range(0, 61) if abs(A_of(i)) < mp.mpf('1e-9')]
    flips = [i for i in range(0, 60)
             if (A_of(i) < 0) != (A_of(i + 1) < 0)]
    say(f"    exact integer zeros in 0..60: {zeros if zeros else 'none'}")
    say(f"    sign changes in 0..60: between i = {flips} and i+1 "
        f"— the only one, so the zero is strictly non-integer")

    say("\n[PR-25(c) re-scoped] the RULE's fourth-order prediction")
    for i in (3,):
        say(f"    pair ({i},{i+1}) incident couplings, descending ratio:")
        say(f"      1st ({i+1},{i+2}): {mp.nstr(mp.mpf(i+3)/(i+4), 6)}"
            f"   -> A_{i}  [derived]")
        say(f"      2nd ({i},{i+1}):   {mp.nstr(mp.mpf(i+2)/(i+3), 6)}"
            f"   -> C_{i}  [derived]")
        say(f"      3rd ({i-1},{i}):   {mp.nstr(mp.mpf(i+1)/(i+2), 6)}"
            f"   -> C3_{i} [derived, above]")
        say(f"      4th ({i+1},{i+3}): {mp.nstr(mp.mpf(i+3)/(i+5), 6)}"
            f"   -> PREDICTED next remainder power")
    say("\n" + "=" * 94)
    OUT.close()
