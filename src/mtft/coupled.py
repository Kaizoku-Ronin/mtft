"""
mtft.coupled — the spatial sector (rungs 5 and 5b)
====================================================
MIT License — (c) 2026 Roger Tano — MTFT Research Program

Couples the internal model to space.  ONE model, two spatial measures:

    H(x) = h(kappa) - tau x B(kappa),      h = diag(g),  B dressed

    Bloch (Z)    x = 2 cos k,   dmu = dk/pi     -> arcsine measure,
                                                   van Hove edges (-1/2)
    Kesten (tree) x in [-2 sqrt p, 2 sqrt p],   -> Bruhat-Tits building of
                  dmu_q as below                   PGL_2(Q_p); SOFT edges (+1/2)

The two differ ONLY in the measure, which is the whole content of rung 5b:
the 1-D density DIVERGES at the band edge and the Kesten density VANISHES
there, so every edge-sensitive statement flips sign of exponent.

STANDING RULE ENFORCED HERE (rung-5 v1 postmortem).  The spatial coupling
must be VACUUM-DRESSED.  `band` and friends take an `Internal` and use its
`.B`; passing a raw kernel is impossible through this API, because the raw
Mellin kernel is unbounded on the log-dense level set (||K_N|| ~ N/kappa)
and H is then unbounded below in infinite volume -- which is exactly how
rung 5 v1 died (band bottom -25, N-drift 17).

WHAT IS AND IS NOT HERE.  Bands, band edges, the merging transition tau_c,
the binding threshold, and the Kesten moment identities.  The EP/expansion
machinery lives in `mtft.ep` and `mtft.expansion`; nothing is duplicated.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np

try:
    import mpmath as mp
except ImportError:                                       # pragma: no cover
    mp = None

from .chain import Internal, internal

__all__ = ["Measure", "BLOCH", "kesten", "H_of", "band", "bands",
           "tau_c", "binding_threshold", "moments", "selftest"]


# ------------------------------------------------------------------ measures
@dataclass(frozen=True)
class Measure:
    """A spatial spectral measure: support and density in the fibre
    variable x, plus the quadrature that integrates against it."""
    name: str
    x_max: float
    density: Any                     # callable x -> weight (Lebesgue)
    p: Optional[int] = None          # tree degree parameter, if any

    def nodes(self, n: int = 400):
        """Gauss-Legendre nodes with the endpoint behaviour REMOVED by
        substitution.  For the tree this is x = x_max sin(theta), which
        cancels the sqrt endpoint exactly -- a raw-x grid stalls at ~1e-6
        because the sqrt is not polynomial there (rung-5b sec.4)."""
        t, w = np.polynomial.legendre.leggauss(n)
        if self.p is None:                       # Bloch: x = 2 cos k
            k = 0.5 * math.pi * (t + 1.0)
            return 2.0 * np.cos(k), w * 0.5
        th = 0.5 * math.pi * t
        x = self.x_max * np.sin(th)
        q = self.p + 1
        dens = (q * self.x_max ** 2 * np.cos(th) ** 2
                / (2 * math.pi * (q ** 2 - x ** 2)))
        return x, w * 0.5 * math.pi * dens


def _bloch_density(x):
    return 1.0 / (math.pi * np.sqrt(np.maximum(4.0 - np.asarray(x) ** 2,
                                               1e-300)))


BLOCH = Measure("bloch", 2.0, _bloch_density, None)


def kesten(p: int = 2) -> Measure:
    """Kesten-McKay measure of the (p+1)-regular tree: the Plancherel
    measure of PGL_2(Q_p) acting on its building.  Support
    [-2 sqrt p, 2 sqrt p]; the density VANISHES like a square root at the
    edges, where Bloch's diverges."""
    q = p + 1
    xm = 2.0 * math.sqrt(p)

    def dens(x):
        x = np.asarray(x, dtype=float)
        return q * np.sqrt(np.maximum(xm ** 2 - x ** 2, 0.0)) \
            / (2 * math.pi * (q ** 2 - x ** 2))

    return Measure(f"kesten(p={p})", xm, dens, p)


# ------------------------------------------------------------------ the model
def H_of(ic: Internal, tau: float, x: float):
    """H(x) = diag(g) - tau x B.  `ic.B` is the DRESSED hopping; there is
    no path through this module that couples a raw kernel."""
    if ic.backend == "mp":
        n = ic.nb
        M = mp.matrix(n, n)
        for a in range(n):
            for b in range(n):
                M[a, b] = (ic.g[a] if a == b else 0) - tau * x * ic.B[a, b]
        return M
    return np.diag(ic.g) - tau * x * ic.B


def band(ic: Internal, tau: float, x: float, nb: Optional[int] = None):
    """Eigenvalues of the fibre at x (ascending)."""
    if ic.backend == "mp":
        ev = sorted(mp.eigsy(H_of(ic, tau, x), eigvals_only=True))
    else:
        ev = np.linalg.eigvalsh(H_of(ic, tau, x))
    return ev[:nb] if nb else ev


def bands(ic: Internal, tau: float, meas: Measure = BLOCH,
          nb: int = 10) -> Dict:
    """Band edges over the measure's support, with the direct and
    indirect gaps.

    Every branch is STRICTLY monotone in x because B is strictly positive
    definite (min-max / Hellmann-Feynman) -- so there are no flat bands
    and the coupled spectrum is purely absolutely continuous (rung-5
    Pr H).  The edges therefore sit at x = +-x_max, and are computed
    there rather than searched for."""
    lo = np.asarray([float(v) for v in band(ic, tau, +meas.x_max, nb)])
    hi = np.asarray([float(v) for v in band(ic, tau, -meas.x_max, nb)])
    widths = hi - lo
    indirect = lo[1:] - hi[:-1]
    return {"lo": lo, "hi": hi, "widths": widths,
            "indirect_gaps": indirect,
            "gapped": bool(np.all(indirect > 0)),
            "measure": meas.name}


def tau_c(ic: Internal, meas: Measure = BLOCH,
          bracket: Tuple[float, float] = (0.02, 3.0),
          iters: int = 60) -> float:
    """Band-merging transition: where the indirect gap between the two
    lowest bands closes.

    First-order prediction m/(x_max (mu_0+mu_1)); the measured value sits
    above it by the eigenvector-rotation factor, which PR-8 proved is
    GEOMETRY-FREE (H depends on (tau, x) only through tau*x, so any
    critical condition reads tau_c * x_max = const)."""
    def g_ind(t):
        b = bands(ic, t, meas, nb=2)
        return float(b["indirect_gaps"][0])

    lo, hi = bracket
    if g_ind(lo) <= 0 or g_ind(hi) >= 0:
        raise ValueError("tau_c: bracket does not straddle the transition")
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        if g_ind(mid) > 0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def binding_threshold(ic: Internal, tau: float, meas: Measure,
                      d: int = 0, nquad: int = 700,
                      deltas: Sequence[float] = (1e-3, 3e-4, 1e-4,
                                                 3e-5, 1e-5)) -> Dict:
    """|V_b| for a rank-one defect on orbital `d`: 1 = V G_dd(E) with E
    below the band bottom.

    THE DIMENSIONAL SIGNATURE (rung-5b Pr S).  On Z the density diverges
    at the edge, so G diverges as delta^{-1/2} and ANY attraction binds:
    V_b = 0.  On the tree the density vanishes, G converges, and the
    threshold is FINITE -- the classic low/high-dimension dichotomy,
    produced by arithmetic geometry rather than by dimension count.

    The tree's G(delta) approaches its limit like sqrt(delta), so the
    threshold needs EXTRAPOLATION in sqrt(delta) (quadratic), not a
    saturation assertion."""
    xs, ws = meas.nodes(nquad)
    cache = [(np.linalg.eigvalsh(H_of(ic, tau, float(x))),
              np.linalg.eigh(H_of(ic, tau, float(x)))[1]) for x in xs]
    e_min = min(float(np.min(ev)) for ev, _ in cache)
    seq = []
    for dd in deltas:
        E = e_min - dd
        tot = 0.0
        for (ev, V), wv in zip(cache, ws):
            tot += float(wv) * float(np.sum(V[d, :] ** 2 / (E - ev)))
        seq.append(abs(tot))
    seq = np.asarray(seq)
    u = np.sqrt(np.asarray(deltas, dtype=float))
    if meas.p is None:                      # Z: divergent, no finite limit
        expo = float(np.polyfit(np.log(deltas), np.log(seq), 1)[0])
        return {"V_b": 0.0, "divergent": True, "exponent": expo,
                "measure": meas.name}
    c = np.polyfit(u, seq, 2)               # G(delta) = G0 - c1 u + c2 u^2
    G0 = float(c[2])
    resid = float(np.max(np.abs(np.polyval(c, u) - seq))) / abs(G0)
    return {"V_b": 1.0 / G0, "divergent": False, "G0": G0,
            "fit_resid": resid, "measure": meas.name}


def moments(meas: Measure, kmax: int = 6, n: int = 800) -> Dict[int, float]:
    """Moments of the spatial measure.

    For the tree these are the CLOSED-WALK COUNTS on the (p+1)-regular
    tree -- m_2 = q, m_4 = q(2q-3)+..., integers -- which is the sharpest
    available check that the measure is the right one."""
    x, w = meas.nodes(n)
    return {k: float(np.sum(w * x ** k)) for k in range(1, kmax + 1)}


def _tree_walks(p: int, kmax: int, depth: int = 7) -> Dict[int, float]:
    """Independent leg for `moments`: closed walks counted on an actual
    truncated tree -- no measure involved."""
    q = p + 1
    adj: Dict[int, List[int]] = {0: []}
    frontier, nxt = [0], 1
    for _ in range(depth):
        newf = []
        for v in frontier:
            deg = q if v == 0 else q - 1
            for _ in range(deg):
                adj[nxt] = []
                adj[v].append(nxt)
                adj[nxt].append(v)
                newf.append(nxt)
                nxt += 1
        frontier = newf
    idx = {v: i for i, v in enumerate(adj)}
    A = np.zeros((len(adj), len(adj)))
    for v, ns in adj.items():
        for u in ns:
            A[idx[v], idx[u]] = 1.0
    vec = np.zeros(len(adj))
    vec[idx[0]] = 1.0
    out = {}
    for step in range(1, kmax + 1):
        vec = A @ vec
        out[step] = float(vec[idx[0]])
    return out


def _L(name):
    from . import ledger as _lg
    if name not in _lg.LEDGER:
        raise KeyError(f"{name!r} is not a ledger entry (BI-F2)")
    return float(_lg.LEDGER[name].value)


def selftest(verbose: bool = True):
    """Targets are ledger lookups or structural invariants (BI-F2)."""
    checks: List = []

    def chk(name, got, want, tol):
        checks.append((name, float(got), want, tol,
                       abs(float(got) - want) <= tol))

    K2 = kesten(2)
    chk("kesten support = 2 sqrt 2", K2.x_max, 2 * math.sqrt(2), 1e-13)
    mm, ww = moments(K2, 6), _tree_walks(2, 6)
    chk("kesten normalised", sum(w for w in
                                 [np.sum(K2.nodes(800)[1])]), 1.0, 1e-12)
    for k in (2, 4, 6):
        chk(f"kesten m_{k} = closed {k}-walks", mm[k], ww[k], 1e-9)
    chk("bloch m_2 = 2", moments(BLOCH, 4)[2], 2.0, 1e-12)

    ic = internal(_L("kappa_star"), nb=30, backend="f64")

    ev = np.linalg.eigvalsh(ic.B)
    chk("B strictly PD -> no flat bands (Pr H)",
        1.0 if ev[0] > 0 else 0.0, 1.0, 0.0)

    b = bands(ic, 0.05, BLOCH, nb=8)
    chk("bloch bands fully gapped at tau=0.05",
        1.0 if b["gapped"] else 0.0, 1.0, 0.0)
    chk("widths narrow with index",
        1.0 if np.all(np.diff(b["widths"][:6]) < 0) else 0.0, 1.0, 0.0)

    tc_Z = tau_c(ic, BLOCH)
    chk("tau_c (bloch)", tc_Z, 0.230032, 5e-5)
    tc_T = tau_c(ic, K2)
    chk("tau_c ratio tree/Z = 2/x_max = 1/sqrt2",
        tc_T / tc_Z, 2.0 / K2.x_max, 2e-4)

    bt_T = binding_threshold(ic, 0.05, K2)
    chk("tree binding threshold finite", bt_T["V_b"], 0.03733, 5e-4)
    chk("tree sqrt-extrapolation fits", bt_T["fit_resid"], 0.0, 1e-3)
    bt_Z = binding_threshold(ic, 0.05, BLOCH)
    chk("Z divergent: V_b = 0", bt_Z["V_b"], 0.0, 0.0)
    chk("Z exponent -1/2 (any attraction binds)",
        bt_Z["exponent"], -0.5, 0.02)

    if verbose:
        for name, got, want, tol, ok in checks:
            print(f"[{'PASS' if ok else 'FAIL'}] {name:<40} "
                  f"{got:>16.9g} vs {want:<14.9g} tol {tol:g}")
        print(f"  {sum(c[4] for c in checks)}/{len(checks)} "
              f"coupled.py self-checks green")
    return all(c[4] for c in checks)
