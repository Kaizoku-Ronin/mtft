"""
mtft.chain — the internal (rung-4) model, de-duplicated
========================================================
MIT License — (c) 2026 Roger Tano — MTFT Research Program

Absorbs the SIXTEEN hand-rolled copies of `internal()` written across the
rung-4 → PR-36 arc.  One definition, two precision backends, and the
crossover documented rather than discovered again.

THE MODEL.  Sites n = 2, 3, ... carry arithmetic weights
    rho(n) = (log n) n^{-3}
(the normalisation -zeta'(3) cancels in every gap, so it is not carried).
The Mellin kernel and transfer operator are
    K_ij = (min/max)^kappa,   T = D K D,  D = diag(sqrt(rho)),
the internal gaps and vacuum-dressed hopping are
    g_i = log(lambda_0/lambda_i),   B = e^{-h/2} (V^T K V) e^{-h/2},
and the strong-coupling limits are pure weight ratios
    g_i(inf) = log(rho_2/rho_{i+2}).

BACKENDS — the choice is in the API, not in convention (PR-20 sec.2).
  'f64' : N-site float64 diagonalisation.  Fast; floors at ~2e-12 on EP
          locations and ~2e-4 on winding bisections.
  'mp'  : closed-form T at working precision, few sites.  The ONLY valid
          backend below ~1e-12, because g and B INHERITED from a float64
          diagonalisation cap any measurement at ~1e-15 regardless of the
          working precision downstream.  Truncation is legitimate at large
          kappa precisely because the couplings die: (4/6)^80 ~ 8e-15.

STANDING RULE ENFORCED HERE (rung-5): spatial couplings must be
vacuum-dressed.  `Internal.B` is the dressed hopping; the raw kernel is
exposed only as `Internal.K_raw` and is documented as UV-sick on the
log-dense level set (||K_N|| ~ N/kappa).
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any, List, Sequence

import numpy as np

try:
    import mpmath as mp
except ImportError:                                       # pragma: no cover
    mp = None

__all__ = ["rho", "gap", "Internal", "internal", "limit_gaps", "mass_gap",
           "crossing_limit", "selftest"]

F64_EP_FLOOR = 2e-12
F64_WINDING_FLOOR = 2e-4


def rho(n, backend: str = "f64"):
    """Arithmetic weight at site n (unnormalised: -zeta'(3) cancels in g)."""
    if backend == "mp":
        n = mp.mpf(n)
        return mp.log(n) / n ** 3
    return math.log(n) / n ** 3


@dataclass
class Internal:
    """Internal model at one coupling.  `B` is the DRESSED hopping."""
    kappa: float
    g: Any                      # gaps g_i = log(lambda_0/lambda_i)
    B: Any                      # vacuum-dressed hopping
    backend: str
    nb: int
    sites: Sequence
    K_raw: Any = field(default=None, repr=False)
    lam: Any = field(default=None, repr=False)
    V: Any = field(default=None, repr=False)

    @property
    def m(self):
        """Mass gap m(kappa) = g_1."""
        return self.g[1]

    def diag(self, i):
        """B_ii — the level's dressed diagonal (= e^{-g_i} in the limit)."""
        return self.B[i, i]


def _internal_f64(kappa, nb, N, gcap):
    n = np.arange(2, N + 1, dtype=np.float64)
    r = np.log(n) * n ** -3.0
    x = np.log(n)
    K = np.exp(-kappa * np.abs(np.subtract.outer(x, x)))
    D = np.sqrt(r)
    T = (D[:, None] * K) * D[None, :]
    lam, V = np.linalg.eigh(T)
    lam = lam[::-1]
    V = V[:, ::-1]
    g = np.log(lam[0] / np.maximum(lam, lam[0] * math.exp(-gcap)))
    nb = min(nb, N - 1)
    VB = V[:, :nb]
    Kt = VB.T @ K @ VB
    eh = np.exp(-0.5 * g[:nb])
    B = (eh[:, None] * Kt) * eh[None, :]
    return g[:nb], B, K, lam, V, list(range(2, N + 1))


def _internal_mp(kappa, nb, nsite, gcap):
    if mp is None:                                        # pragma: no cover
        raise RuntimeError("backend='mp' requires mpmath")
    kappa = mp.mpf(kappa)
    sites = [mp.mpf(s) for s in range(2, 2 + nsite)]
    r = [mp.log(s) / s ** 3 for s in sites]
    K = mp.matrix(nsite, nsite)
    T = mp.matrix(nsite, nsite)
    for i in range(nsite):
        for j in range(nsite):
            q = min(sites[i], sites[j]) / max(sites[i], sites[j])
            K[i, j] = q ** kappa
            T[i, j] = mp.sqrt(r[i] * r[j]) * K[i, j]
    lam, V = mp.eigsy(T)
    order = sorted(range(nsite), key=lambda j: -lam[j])
    lam = [lam[j] for j in order]
    Vd = mp.matrix(nsite, nsite)
    for c, j in enumerate(order):
        for row in range(nsite):
            Vd[row, c] = V[row, j]
    g = [mp.log(lam[0] / lam[j]) for j in range(nsite)]
    nb = min(nb, nsite)
    Kt = Vd.T * K * Vd
    B = mp.matrix(nb, nb)
    for i in range(nb):
        for j in range(nb):
            B[i, j] = mp.e ** (-(g[i] + g[j]) / 2) * Kt[i, j]
    return g[:nb], B, K, lam, Vd, [int(s) for s in sites]


def internal(kappa, nb: int = 60, backend: str = "f64", N: int = 1600,
             nsite: int = 14, gcap: float = 200.0) -> Internal:
    """The internal model at one coupling.

    backend='f64': N-site numerical diagonalisation (fast, floors ~1e-15).
    backend='mp' : nsite-site closed-form T at working precision (exact).
    """
    if backend not in ("f64", "mp"):
        raise ValueError("backend must be 'f64' or 'mp'")
    if backend == "f64":
        g, B, K, lam, V, sites = _internal_f64(kappa, nb, N, gcap)
    else:
        g, B, K, lam, V, sites = _internal_mp(kappa, nb, nsite, gcap)
    return Internal(kappa=kappa, g=g, B=B, backend=backend, nb=nb,
                    sites=sites, K_raw=K, lam=lam, V=V)


def gap(i, backend: str = "mp"):
    """g_i(infinity) = log(rho_2/rho_{i+2}) at REAL i.

    The continuous-i form is the primitive; `limit_gaps` is the integer
    view of it.  Keeping them separate matters: an integer-indexed helper
    silently breaks at real i, which is how a whole-family zero scan once
    reported "no zeros" for a family whose zero was known (Add. BE)."""
    r2 = rho(2, backend)
    if backend == "mp":
        return mp.log(r2 / (mp.log(mp.mpf(i) + 2) / (mp.mpf(i) + 2) ** 3))
    return math.log(r2 / (math.log(i + 2) / (i + 2) ** 3))


def limit_gaps(nlev: int = 6, backend: str = "mp"):
    """g_i(infinity) = log(rho_2/rho_{i+2}) — pure weight ratios.

    g_1(inf) = m_inf = log(27 ln2/(8 ln3)); g_2(inf) = log 4 exactly."""
    return [gap(j, backend) for j in range(nlev)]


def mass_gap(kappa, **kw):
    """m(kappa) = g_1.  Non-monotone: interior minimum near kappa* = 5."""
    return internal(kappa, nb=kw.pop("nb", 30), **kw).m


def crossing_limit(i: int, backend: str = "mp"):
    """|u_{i,i+1}(inf)| = (g_{i+1}-g_i)/(e^{-g_i} - e^{-g_{i+1}}).

    The diagonal limit's real level crossing (rung-4 Pr G): as kappa -> inf
    the hopping becomes diagonal, levels are linear in u, and members i,i+1
    cross on the negative real axis.  i=0 gives R(inf) = m_inf/(1-e^{-m_inf}).
    """
    g = limit_gaps(i + 2, backend)
    a, b = g[i], g[i + 1]
    if backend == "mp":
        return abs((b - a) / (mp.e ** (-a) - mp.e ** (-b)))
    return abs((b - a) / (math.exp(-a) - math.exp(-b)))



def _L(name):
    """Ledger scalar as float.  Real-ledger API: LEDGER[name] is an Entry
    whose .value is a decimal STRING (exact representation stored)."""
    from . import ledger as _lg
    if name not in _lg.LEDGER:
        raise KeyError(
            f"{name!r} is not a ledger entry — a selftest may not certify "
            f"against an unregistered number (BI-F2).  Register it first.")
    return float(_lg.LEDGER[name].value)


def _LF(family, i):
    """Family member as float (A_FAMILY / C_FAMILY hold Entry lists)."""
    from . import ledger as _lg
    return float(getattr(_lg, family).entries[i].value)

def selftest(verbose: bool = True):
    """Reproduce ledger constants from the model definition.

    Targets are LEDGER LOOKUPS (BI-F2).  The previous version made this
    claim in prose while hardcoding every literal -- the claim is now
    true as implemented.  mu_0 and mu_1 were certifying against numbers
    the ledger did not carry; they are registered ledger entries now."""
    if mp is not None:
        mp.mp.dps = 40
    checks: List = []

    def chk(name, got, want, tol):
        ok = abs(float(got) - want) <= tol
        checks.append((name, float(got), want, tol, ok))

    g = limit_gaps(6)
    chk("m_inf = g_1(inf)", g[1], _L("m_inf"), 1e-11)
    chk("g_2(inf) = log 4", g[2], math.log(4), 1e-11)
    chk("u_ep_01 = R_inf", crossing_limit(0), _L("u_ep_01"), 2e-9)
    chk("u_ep_12", crossing_limit(1), _L("u_ep_12"), 2e-9)
    chk("u_ep_23", crossing_limit(2), _L("u_ep_23"), 2e-9)
    chk("u_ep_34", crossing_limit(3), _L("u_ep_34"), 2e-9)
    chk("u_ep_45", crossing_limit(4), _L("u_ep_45"), 2e-9)

    R = float(crossing_limit(0))
    chk("R_inf/m_inf", R / float(g[1]), _L("R_inf_over_m_inf"), 5e-7)

    ms = mass_gap(_L("kappa_star"), backend="f64", N=1600)
    chk("m_star = m(kappa*)", ms, _L("m_star"), 1e-5)
    chk("m(kappa=200) -> m_inf", mass_gap(200.0, backend="mp", nsite=12),
        _L("m_inf"), 1e-9)

    ic = internal(5.0, nb=30, backend="f64")
    chk("mu_0 = B_00", ic.diag(0), _L("mu_0"), 1e-5)
    chk("mu_1 = B_11", ic.diag(1), _L("mu_1"), 1e-5)
    ev = np.linalg.eigvalsh(ic.B)
    chk("B strictly PD (min eig > 0)", 1.0 if ev[0] > 0 else 0.0, 1.0, 0.0)

    Kn = [float(np.linalg.eigvalsh(
        internal(5.0, nb=8, backend="f64", N=n).K_raw[:n - 1, :n - 1])[-1])
        * 5.0 / n for n in (400, 800)]
    chk("raw K is UV-sick: ||K||k/N const", Kn[1] / Kn[0], 1.0, 0.01)

    if verbose:
        for name, got, want, tol, ok in checks:
            print(f"[{'PASS' if ok else 'FAIL'}] {name:<32} "
                  f"{got:>18.10g}  vs {want:<16.10g} tol {tol:g}")
        n_ok = sum(c[4] for c in checks)
        print(f"  {n_ok}/{len(checks)} chain.py self-checks green")
    return all(c[4] for c in checks)


if __name__ == "__main__":
    selftest()
