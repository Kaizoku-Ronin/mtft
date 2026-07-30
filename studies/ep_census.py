# --- stage 5 (Integration Plan v0.1): re-pointed to the mtft package. ---
# internal() delegates to mtft.chain.internal (certified spectrally
# identical: g to 9e-15, full H(u) spectrum to 2.5e-14; B differs only
# by gauge inside near-degenerate blocks).  Where a local gsq() remains
# it is retained VERBATIM: it evaluates complex u in the f64 backend,
# which mtft.ep deliberately does not offer (PR-20 floors) -- S5-1.
# ------------------------------------------------------------------------
#!/usr/bin/env python3
"""
ep_census.py — certified EP census by the argument principle (PR-15a)
======================================================================
MIT License — (c) 2026 Roger Tano — MTFT Research Program

Addendum AJ-F2: five separate Newton nets (two engines) each missed
roots the other found.  EP *existence* is certifiable by Newton; EP
*census* is not.  Kimi proposed the discriminant polynomial; its
numerically stable form is the ARGUMENT PRINCIPLE on the discriminant's
logarithmic derivative, which never forms the (astronomically large)
discriminant itself:

    D(u) = prod_{i<j} (eps_i - eps_j)^2   (the discriminant, a
                                           polynomial in u)
    D'/D = sum_{i<j} 2 (eps_i' - eps_j')/(eps_i - eps_j),
    eps_i' = -(v_i^T B v_i)/(v_i^T v_i)   (complex symmetric: left
                                           eigenvector = v^T)
    N(|u| < r) = (1/2 pi i) contour-integral of D'/D.

A simple EP makes (eps_a-eps_b)^2 ~ (u-u_c), so D has a SIMPLE zero
there: the integral counts exceptional points exactly, with multiplicity,
and with no seed, basin, or net anywhere.  Trapezoid on a circle is
spectrally accurate for the analytic integrand.
"""
import math
import numpy as np
import mpmath as mp
from mtft.chain import internal as _chain_internal


BETA, KSTAR = 2.0, 5.0


def internal(N=1600, kappa=KSTAR, nb=40, gcap=200.0):
    ic = _chain_internal(kappa, nb=nb, backend="f64", N=N, gcap=gcap)
    return np.asarray(ic.g), np.asarray(ic.B)


def dlogD(u, g, B):
    ev, V = np.linalg.eig(np.diag(g) - u * B)
    nrm = np.einsum('ij,ij->j', V, V)              # v^T v (no conjugate)
    dep = -np.einsum('ij,jk,ki->i', V.T, B, V) / nrm
    d = ev[:, None] - ev[None, :]
    dd = dep[:, None] - dep[None, :]
    np.fill_diagonal(d, 1.0); np.fill_diagonal(dd, 0.0)
    return float(np.sum(dd / d).real) + 1j * float(np.sum(dd / d).imag)


def count(r, g, B, npts=800):
    th = 2 * math.pi * np.arange(npts) / npts
    u = r * np.exp(1j * th)
    f = np.array([dlogD(uu, g, B) for uu in u])
    integ = np.sum(f * 1j * u) * (2 * math.pi / npts)
    return integ / (2j * math.pi)


if __name__ == "__main__":
    g, B = internal()
    print("== argument-principle EP census, kappa = 5, nb = 40 ==",
          flush=True)
    print(f"{'radius':>8} {'count (should be ~integer)':>32}")
    for r in (1.5, 2.0, 2.6, 3.0, 3.5, 4.6, 4.8, 5.5, 9.0):
        c = count(r, g, B)
        print(f"{r:8.2f}   {c.real:12.4f} {c.imag:+.2e}i", flush=True)
    print("\n== nb-convergence of the |u|<5.5 count ==", flush=True)
    for nb in (20, 30, 40, 50):
        gg, BB = internal(nb=nb)
        print(f"  nb={nb}: {count(5.5, gg, BB).real:.4f}", flush=True)
