#!/usr/bin/env python3
"""li_density_surrogate.py — Table 4.5 of "Three Ensembles on One Arithmetic".

The density-first demonstration (Df 4.4 / Measured 4.5): exact Li
coefficients lambda_n (certified series method, mtft.critical_ensemble)
versus the surrogate

    lam_hat_n(K) = sum_{k<=K} pair(gamma_k, n)
                   + integral_Gamma^inf pair(t, n) * (1/2pi) log(t/2pi) dt,

pair(t, n) = 4 sin^2(n arctan(1/(2t)))   — the on-line zero-pair term,

K = 100 exact pairs, Gamma = gamma_100. The residual r_n = lambda_n -
lam_hat_n isolates what the smooth Riemann-von Mangoldt density cannot
see: |r_n|/lambda_n < 1.5e-4 for n <= 12, constant to 2% across n
(-1.47 -> -1.50 e-4), i.e. an O(1) boundary count offset, DeltaN_eff ~ 0.19
(Heur 4.6 — interpretation only; the measured facts are the table).

First run computes gamma_1..gamma_100 with mpmath.zetazero (minutes) and
caches to zeros100.json alongside this script. Requires: mtft, mpmath.
Two-engine verified July 2026 (Claude Fable 5 / Kimi K3).
"""
import json
import os

import mpmath as mp

mp.mp.dps = 25

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, "zeros100.json")
K = 100
NS = (1, 2, 3, 4, 6, 8, 10, 12)


def load_zeros(k_max):
    have = json.load(open(CACHE)) if os.path.exists(CACHE) else []
    for k in range(len(have) + 1, k_max + 1):
        have.append(mp.nstr(mp.zetazero(k).imag, 20))
        if k % 25 == 0:
            json.dump(have, open(CACHE, "w"))
            print(f"  ...cached {k} ordinates", flush=True)
    json.dump(have, open(CACHE, "w"))
    return [mp.mpf(s) for s in have]


def pair(t, n):
    return 4 * mp.sin(n * mp.atan(1 / (2 * t))) ** 2


def main():
    import mtft
    gam = load_zeros(K)
    G = gam[-1]
    print(f"K = {K} exact pairs, Gamma = gamma_100 = {mp.nstr(G, 12)}\n")
    density = lambda t: mp.log(t / (2 * mp.pi)) / (2 * mp.pi)
    hdr = f"{'n':>3} {'lambda_n':>14} {'head':>14} {'tail':>14} {'lam_hat':>14} {'r_n':>11} {'r_n/lam':>11}"
    print(hdr)
    rows = []
    for n in NS:
        head = mp.fsum([pair(t, n) for t in gam])
        tail = mp.quad(lambda t: pair(t, n) * density(t), [G, mp.inf])
        lam = mp.mpf(mtft.li_lambda(n))
        hat = head + tail
        r = lam - hat
        rows.append((n, lam, r, r / lam))
        print(f"{n:>3} {mp.nstr(lam, 11):>14} {mp.nstr(head, 11):>14} "
              f"{mp.nstr(tail, 11):>14} {mp.nstr(hat, 11):>14} "
              f"{mp.nstr(r, 3):>11} {mp.nstr(r / lam, 3):>11}")
    dN1 = -rows[0][2] / pair(G, 1)
    dN12 = -rows[-1][2] / pair(G, 12)
    print(f"\nboundary offset (Heur 4.6): DeltaN_eff = {mp.nstr(dN1, 4)} "
          f"(n=1), {mp.nstr(dN12, 4)} (n=12) — kernel-proportional residual")
    ratios = [row[3] for row in rows]
    spread = (max(ratios) - min(ratios)) / abs(sum(ratios) / len(ratios))
    print(f"r_n/lambda_n constancy: {mp.nstr(min(ratios), 3)} .. "
          f"{mp.nstr(max(ratios), 3)}, spread {mp.nstr(spread, 2)} of mean")


if __name__ == "__main__":
    main()
