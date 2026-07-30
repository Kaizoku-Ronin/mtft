# --- stage 5 (Integration Plan v0.1): re-pointed to the mtft package. ---
# internal() delegates to mtft.chain.internal (certified spectrally
# identical: g to 9e-15, full H(u) spectrum to 2.5e-14; B differs only
# by gauge inside near-degenerate blocks).  Where a local gsq() remains
# it is retained VERBATIM: it evaluates complex u in the f64 backend,
# which mtft.ep deliberately does not offer (PR-20 floors) -- S5-1.
# ------------------------------------------------------------------------
#!/usr/bin/env python3
"""
pr18.py — PR-18: the owed item, and the mechanism's free second test
=====================================================================
MIT License — (c) 2026 Roger Tano — MTFT Research Program
Prints to stdout AND pr18_run.txt.

ADDENDUM AO DISPOSITIONS:
 AO-F1 (overshoot is structure): ACCEPTED on the auditor's evidence.
   Their engine resolves R(kappa) to 10 digits with nb-stability 1e-14;
   mine floors at ~2e-4.  I cannot independently resolve a 1.5e-4
   overshoot and will not pretend to.  Their table is the ledger.
 AO-F2 (lost pre-registered item): DELIVERED here — the four-parameter
   fit extended to kappa != 5 with the k-window analysis AN-F1 required.
 AO-F3 (gate classes): applied in the note.

PR-18(b) — THE FREE SECOND TEST.  The diagonal-limit construction is not
special to levels 0,1.  Levels i,j cross at
    u_ij = -(g_j - g_i)/(e^{-g_j} - e^{-g_i}),
and AO.3 decodes the limiting gaps as g_i(inf) = ln(rho_2/rho_{i+2}) —
pure arithmetic weight ratios, with g_2(inf) = ln 4 EXACTLY (because
rho_2/rho_4 = (ln2/8)/(ln4/64) = 4).  Hence a parameter-free prediction
for the levels-1/2 exceptional point PR-14 measured at 2.481614:

    u_12(inf) = (ln4 - m_inf)/(e^{-m_inf} - 1/4) = 2.8707...

The same mechanism, a different level pair, tested against data already
in the record.  Falsifiable: the measured levels-1/2 EP must climb to
2.8707 and stop.
"""
from __future__ import annotations
import math
import numpy as np
import mpmath as mp
from mtft.chain import internal as _chain_internal


BETA = 2.0
OUT = open("pr18_run.txt", "w", buffering=1)


def say(s=""):
    print(s, flush=True)
    OUT.write(s + "\n")


def internal(N=1600, kappa=5.0, nb=30, gcap=200.0):
    ic = _chain_internal(kappa, nb=nb, backend="f64", N=N, gcap=gcap)
    return np.asarray(ic.g), np.asarray(ic.B)


def pair_winding(r, g, B, i0=0, i1=1, npts=900):
    th = 2 * math.pi * np.arange(npts + 1) / npts
    ev = np.linalg.eigvalsh(np.diag(g) - r * B)
    a, b = complex(ev[i0]), complex(ev[i1])
    ap, bp = a, b
    tot, qprev = 0.0, b - a
    for t in th[1:]:
        u = r * np.exp(1j * t)
        cur = np.linalg.eigvals(np.diag(g) - u * B)
        pa, pb = 2 * a - ap, 2 * b - bp
        ia = int(np.argmin(np.abs(cur - pa)))
        rest = np.delete(cur, ia)
        ib = int(np.argmin(np.abs(rest - pb)))
        ap, bp = a, b
        a, b = cur[ia], rest[ib]
        q = b - a
        tot += np.angle(q / qprev)
        qprev = q
    return 2.0 * tot / (2 * math.pi)


def R_pair(kap, i0, i1, lo, hi, nb=30, npts=900, iters=40):
    g, B = internal(kappa=kap, nb=nb)
    if pair_winding(lo, g, B, i0, i1, npts) > 1.0:
        return float('nan')
    if pair_winding(hi, g, B, i0, i1, npts) < 1.0:
        return float('nan')
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        if pair_winding(mid, g, B, i0, i1, npts) < 1.0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def cauchy(kap, rho, npts=512, kmax=40, nb=30, lev=1):
    g, B = internal(kappa=kap, nb=nb)
    th = 2 * math.pi * np.arange(npts) / npts
    vals = np.empty(npts, dtype=complex)
    cur = complex(np.linalg.eigvalsh(np.diag(g) - rho * B)[lev])
    for j, t in enumerate(th):
        u = rho * np.exp(1j * t)
        ev = np.linalg.eigvals(np.diag(g) - u * B)
        cur = ev[int(np.argmin(np.abs(ev - cur)))]
        vals[j] = cur
    return np.array([np.sum(vals * np.exp(-1j * k * th)) / npts / rho ** k
                     for k in range(kmax + 1)]).real


def fit_R(cs, kmin, kmax=40, Rlo=1.0, Rhi=1.8):
    ks = np.arange(kmin, kmax + 1)
    y = cs[ks]
    best = None
    for R in np.arange(Rlo, Rhi, 0.0004):
        base = ks ** -1.5 * R ** (-ks.astype(float))
        for thh in np.arange(2.30, 3.10, 0.004):
            A = np.column_stack([base * np.cos(ks * thh),
                                 -base * np.sin(ks * thh)])
            sol, *_ = np.linalg.lstsq(A, y, rcond=None)
            r2 = float(np.sum((A @ sol - y) ** 2))
            if best is None or r2 < best[0]:
                best = (r2, R, thh)
    return best[1], best[2]


if __name__ == "__main__":
    say("=" * 94)
    say("  PR-18 — THE OWED ITEM (fit at kappa != 5), AND THE FREE "
        "SECOND TEST")
    say("=" * 94)

    with mp.workdps(30):
        m_inf = float(mp.log(27 * mp.log(2) / (8 * mp.log(3))))
    e1 = math.exp(-m_inf)
    u12 = (math.log(4) - m_inf) / (e1 - 0.25)
    say(f"\n[b] parameter-free prediction for the LEVELS-1/2 crossing")
    say(f"    g_1(inf) = m_inf = ln(rho2/rho3) = {m_inf:.7f}")
    say(f"    g_2(inf) = ln(rho2/rho4) = ln 4  = {math.log(4):.7f}  "
        f"(AO.3: rho2/rho4 = 4 exactly)")
    say(f"    u_12(inf) = (ln4 - m)/(e^-m - 1/4) = {u12:.7f}")
    say(f"    PR-14 measured the levels-1/2 EP at kappa=5: 2.481614")
    say(f"\n{'kappa':>7} {'|u| levels 1-2':>16} {'u12(inf) - |u|':>16}")
    for kap in (5.0, 8.0, 16.0, 24.0):
        v = R_pair(kap, 1, 2, 1.60, 4.20)
        say(f"{kap:7.1f} {v:16.6f} {u12 - v:16.6f}")

    say(f"\n[AO-F2, owed] four-parameter fit at kappa != 5, with the "
        f"k-window analysis")
    say(f"{'kappa':>7} {'window':>10} {'R fit':>11} {'R ledger':>11} "
        f"{'dev':>10}")
    ledger = {4.0: 1.2568164063, 8.0: 1.4244407112}
    for kap, rho in ((4.0, 0.95), (8.0, 1.05)):
        cs = cauchy(kap, rho)
        for kmin in (15, 25, 35):
            R, th = fit_R(cs, kmin)
            say(f"{kap:7.1f} {('k>=' + str(kmin)):>10} {R:11.4f} "
                f"{ledger[kap]:11.4f} {abs(R-ledger[kap]):10.4f}")
    say("\n" + "=" * 94)
    OUT.close()
