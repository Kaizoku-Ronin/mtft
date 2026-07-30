# --- stage 5 (Integration Plan v0.1): re-pointed to the mtft package. ---
# internal() delegates to mtft.chain.internal (certified spectrally
# identical: g to 9e-15, full H(u) spectrum to 2.5e-14; B differs only
# by gauge inside near-degenerate blocks).  Where a local gsq() remains
# it is retained VERBATIM: it evaluates complex u in the f64 backend,
# which mtft.ep deliberately does not offer (PR-20 floors) -- S5-1.
# ------------------------------------------------------------------------
#!/usr/bin/env python3
"""
rung5_v1_postmortem_ARCHIVE.py — the FAILED v1 design, archived
================================================================
MIT License — (c) 2026 Roger Tano — MTFT Research Program

PROVENANCE NOTE (Addendum Z, F4).  The original v1 script was deleted
during the rung-5 build before the failure was written up.  This file is
a FAITHFUL REPRODUCTION of the v1 design — raw Mellin hopping, no vacuum
dressing — written to make the postmortem numbers executable, not the
original artifact.  Labeled as such per corpus provenance discipline.

THE FAILED MODEL:   H(k) = h(kappa) - 2 tau cos(k) M(kappa)
with M(n,m) = (min/max)^kappa the RAW Mellin kernel.

WHY IT DIES.  The level set x_n = log n is log-dense: a vector flat over
a cluster of levels has Rayleigh quotient growing with cluster size, so
||M_N|| ~ N/kappa (measured: ||M_N||*kappa/N = 1.122, constant across
N = 400/800/1600).  M is therefore unbounded, H is UNBOUNDED BELOW in
infinite volume, and the band structure is meaningless: the "band
bottom" runs off to -infinity with N.

THE CURE (shipped as v2):  B = e^{-h/2} M e^{-h/2}, Kato class,
||B|| = 1.07543049 identical to 8 digits at N = 400/800/1600.

Run:  py rung5_v1_postmortem_ARCHIVE.py
"""
import math
import numpy as np
import mpmath as mp

BETA, KAPPA, TAU = 2.0, 5.0, 0.05


def raw_v1(N, nb=300, gcap=200.0):
    n = np.arange(2, N + 1, dtype=np.float64)
    with mp.workdps(30):
        Z2 = float(-mp.zeta(BETA + 1, derivative=1))
    rho = np.log(n) * n ** (-(BETA + 1.0)) / Z2
    x = np.log(n)
    M = np.exp(-KAPPA * np.abs(np.subtract.outer(x, x)))
    D = np.sqrt(rho)
    T = (D[:, None] * M) * D[None, :]
    lam, V = np.linalg.eigh(T)
    lam = lam[::-1]; V = V[:, ::-1]
    g = np.log(lam[0] / np.maximum(lam, lam[0] * math.exp(-gcap)))
    nb = min(nb, N - 1)
    Mtil = V[:, :nb].T @ M @ V[:, :nb]          # RAW — no dressing
    return g[:nb], Mtil, M


if __name__ == "__main__":
    print("=" * 84)
    print("  RUNG-5 v1 POSTMORTEM (ARCHIVED FAILURE) — raw Mellin hopping")
    print("=" * 84)
    print(f"{'N':>6} {'||M_N||':>10} {'||M||k/N':>10} {'band bottom':>14} "
          f"{'edge drift':>12}")
    prev = None
    for N in (400, 800, 1600, 2400):
        g, Mtil, M = raw_v1(N)
        nrm = float(np.linalg.eigvalsh(M)[-1])
        w = np.linalg.eigvalsh(np.diag(g) - 2 * TAU * 1.0 * Mtil)[:10]
        drift = "-" if prev is None else f"{float(np.max(np.abs(w - prev))):.2e}"
        prev = w
        print(f"{N:>6} {nrm:>10.2f} {nrm*KAPPA/N:>10.3f} {w[0]:>14.3f} "
              f"{drift:>12}")
    print("-" * 84)
    print("  DIAGNOSIS: ||M||*kappa/N constant => ||M_N|| ~ N/kappa "
          "(unbounded);")
    print("  band bottom runs to -infinity and edges drift O(10) with N.")
    print("  The gates (BG2/BG3) killed this design on its first run. "
          "Cure: vacuum dressing (v2).")
    print("=" * 84)
