"""
mtft.expansion — the remainder hierarchy, and safe extraction
==============================================================
MIT License — (c) 2026 Roger Tano — MTFT Research Program

The deviation of an exceptional point from its strong-coupling limit is
known to four derived orders, each governed by a coupling INCIDENT ON THE
CROSSING MEMBERS, taken in descending ratio order (PR-19 → PR-34):

    dev(i;kappa) = A_i r^{2k} - C_i s^{2k} + C3_i t^{2k} + K_i q^{2k} + ...
      r = (i+3)/(i+4)   (i+1,i+2)  upward from the UPPER member
      s = (i+2)/(i+3)   (i,i+1)    within the pair
      t = (i+1)/(i+2)   (i-1,i)    downward from the LOWER member
      q = (i+3)/(i+5)   (i+1,i+3)  upward from the upper member, skipping
      (5th) (i+2)/(i+4) (i,i+2)    upward from the LOWER member

Rung 4's bridge law (5/6)^{2k} is the i = 2 member of the first family.

WHY THIS MODULE IS DEFENSIVE.  Both engines' recent errors concentrated
here: a stencil divided by 2 instead of 4 (BA), an mp.diff that returned
exact zeros on a root-defined observable (PR-29), a scan that returned
"no zeros" because an exception was swallowed (BE).  Three rules are
therefore enforced in the API rather than documented:

  1. `extract()` always validates against a control whose coefficients are
     known -- a well-designed default if none is supplied -- runs it on
     EVERY call, and validates EACH STENCIL SEPARATELY.  A validation
     that covers one of three components is not a validation of the sum
     (PR-30 sec.3).
  2. `scan_zeros()` REQUIRES a positive control that must be found, or it
     raises.  A vacuous scan and a negative scan are indistinguishable
     without one (BE.0).
  3. Nothing here differentiates an implicitly defined observable with
     `mp.diff`; explicit finite differences with Richardson only.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from fractions import Fraction
from typing import Callable, Dict, List, Optional, Sequence, Tuple

try:
    import mpmath as mp
except ImportError:                                       # pragma: no cover
    mp = None

from .chain import gap, limit_gaps, rho                   # noqa: F401

__all__ = ["incident_couplings", "channels_A", "channels_C", "channels_C3",
           "A", "C", "C3", "cancellation_ratio", "scan_zeros",
           "extract", "richardson", "completeness", "power_fit",
           "separation_hazard", "MONOMIAL_COUNT", "selftest"]

MONOMIAL_COUNT = 3   # X^2, XYZ, X^2Y^2 — see `monomials.__doc__`
HAZARD_RATIO = 0.95  # power separations slower than this contaminate any kappa


# ----------------------------------------------------------------- hierarchy
def incident_couplings(i: int) -> List[Tuple[str, Fraction, str]]:
    """The couplings incident on the crossing members {i, i+1}, in
    DESCENDING ratio order — which is the order of the remainder."""
    raw = [(f"({i+1},{i+2})", Fraction(i + 3, i + 4), "A  up-from-upper"),
           (f"({i},{i+1})", Fraction(i + 2, i + 3), "C  within-pair"),
           (f"({i-1},{i})", Fraction(i + 1, i + 2), "C3 down-from-lower"),
           (f"({i+1},{i+3})", Fraction(i + 3, i + 5), "K  up-from-upper"),
           (f"({i},{i+2})", Fraction(i + 2, i + 4), "K5 up-from-lower")]
    return sorted(raw, key=lambda z: -z[1])


def _mpf(x):
    return mp.mpf(x) if mp is not None else float(x)


def _rho_at(x):
    """rho at a REAL site index — the continuous-i primitive."""
    x = mp.mpf(x)
    return mp.log(x) / x ** 3


def channels_A(i):
    """(eigenvalue, eigenvector, repulsion) constituents of A_i."""
    a, b, c = gap(i), gap(i + 1), gap(i + 2)
    ea, eb, ec = mp.e ** (-a), mp.e ** (-b), mp.e ** (-c)
    D, m = ea - eb, b - a
    u0 = m / D
    r1, r2 = _rho_at(i + 3), _rho_at(i + 4)
    return (((D - m * eb) / D ** 2) * (-r2 / (r1 - r2)),
            (m / D ** 2) * eb * (2 * mp.sqrt(r1 * r2) / (r1 - r2)),
            (-u0 ** 2 * mp.e ** (-(b + c))
             / ((c + u0 * ec) - (b + u0 * eb))) / D)


def channels_C(i):
    """Constituents of C_i.  For i = 0 the LOWER member is the reference
    level lambda_0, so the within-pair coupling shifts lambda_0 itself and
    g_1 picks up d(lambda_0)/lambda_0 as well (PR-24 sec.2)."""
    a, b = gap(i), gap(i + 1)
    ea, eb = mp.e ** (-a), mp.e ** (-b)
    D, m = ea - eb, b - a
    r0, r1 = _rho_at(i + 2), _rho_at(i + 3)
    if abs(i) < 1e-12:
        dg1, dg0 = (r0 + r1) / (r0 - r1), _mpf(0)
    else:
        dg1, dg0 = r0 / (r0 - r1), -r1 / (r0 - r1)
    E1 = ((D - m * eb) / D ** 2) * dg1 + ((m * ea - D) / D ** 2) * dg0
    X = 2 * mp.sqrt(r0 * r1) / (r0 - r1)
    return (-E1, m * X * (ea + eb) / D ** 2,
            2 * m * mp.e ** (-(a + b)) / D ** 3)


def channels_C3(i):
    """Constituents of C3_i.  Level i-1 sits ABOVE the crossing, so the
    repulsion acts from below."""
    if i < 1:
        raise ValueError("C3 undefined for i < 1 (no lower-adjacent level)")
    am, a, b = gap(i - 1), gap(i), gap(i + 1)
    eam, ea, eb = mp.e ** (-am), mp.e ** (-a), mp.e ** (-b)
    D, m = ea - eb, b - a
    u0 = m / D
    rm, r0 = _rho_at(i + 1), _rho_at(i + 2)
    dga = ((rm + r0) / (rm - r0) if abs(i - 1) < 1e-12
           else rm / (rm - r0))
    X = 2 * mp.sqrt(rm * r0) / (r0 - rm)
    return (((m * ea - D) / D ** 2) * dga,
            -m * (ea * X) / D ** 2,
            (u0 ** 2 * mp.e ** (-(am + a))
             / ((am + u0 * eam) - (a + u0 * ea))) / D)


def A(i):
    return sum(channels_A(i))


def C(i):
    """Note the channel convention: `channels_C` returns the constituents
    already carrying the sign they contribute to dev = ... - C s^{2k},
    so C is their plain sum.  Flipping twice was this module's first bug,
    and the cancellation-ratio checks did NOT catch it — ratios are
    sign-invariant.  Per-quantity checks, not per-derived-quantity."""
    return sum(channels_C(i))


def C3(i):
    return sum(channels_C3(i))


# --------------------------------------------------------- zero-capability
def cancellation_ratio(channels: Sequence) -> float:
    """|sum| / max|constituent|.  A coefficient is ZERO-CAPABLE exactly
    when this dips below 1 — opposition that actually bites.  Mixed signs
    are necessary, not sufficient (BF.3: C is mixed-sign with no zero)."""
    s = abs(float(sum(channels)))
    return s / max(abs(float(x)) for x in channels)


def scan_zeros(f: Callable, lo: float, hi: float, step: float = 0.25,
               control: Optional[Tuple[Callable, float, float]] = None):
    """Bracket sign changes of f on [lo, hi].

    `control` = (g, clo, chi) with a KNOWN zero of g in (clo, chi); it must
    be found or this raises.  Without it a vacuous scan (an exception
    swallowed, a wrong signature) is indistinguishable from a negative
    result — which is exactly how BE's first scan reported 'no zeros'
    for every family including one whose zero was known.

    BLIND SPOT (inherent, BI-F5): sign-change bracketing cannot see an
    even-multiplicity (touching) zero that falls between grid points.
    Zeros landing exactly on the grid are caught by the `prev == 0` clause.
    """
    if control is None:
        raise ValueError("scan_zeros requires a positive control "
                         "(standing rule: a vacuous scan looks negative)")
    g, clo, chi = control
    if not _brackets(g, clo, chi, step):
        raise RuntimeError("positive control FAILED — scan is not trustworthy")
    return _brackets(f, lo, hi, step)


def _brackets(f, lo, hi, step):
    out, x = [], lo
    prev = f(x)
    while x < hi:
        x2 = min(x + step, hi)
        cur = f(x2)
        if prev == 0 or (prev < 0) != (cur < 0):
            out.append((x, x2))
        x, prev = x2, cur
    return out


# ----------------------------------------------------------- FD + Richardson
def richardson(f: Callable, hs: Sequence[float], order: int = 2):
    """Extrapolate f(h) -> f(0) assuming an O(h^order) error series.

    NEVILLE TRIANGLE.  The naive in-place recurrence — combining the RAW
    f(h_k) with the already-extrapolated previous value at a fixed
    fac = (h_0/h_k)^order — reintroduces the leading error term for k >= 2
    and is *worse* than the 2-point result (Addendum BI-F1: 18x worse on
    1 + 2h^2 + 3h^4, and non-zero on a pure h^2 series where any correct
    scheme is exact).  The successive-ratio form below is exact on both.
    Fix supplied by the auditor; the bug was mine, and it shipped because
    an exported numerical primitive had no selftest check at all — the
    per-quantity rule applies to EVERY export, not just the delicate ones.
    """
    R = [[f(h)] for h in hs]
    for j in range(1, len(hs)):
        for k in range(j, len(hs)):
            fac = (hs[k - j] / hs[k]) ** order
            R[k].append(R[k][j - 1] + (R[k][j - 1] - R[k - 1][j - 1])
                        / (fac - 1))
    return R[-1][-1]


def _stencils(dev, h, hh):
    """The three monomial stencils.  NOTE the /4 on the last: a double
    second difference of a X^2 Y^2 term returns 4*c, not 2*c (BA)."""
    q2 = (dev(0, 0, h) + dev(0, 0, -h)) / (2 * h ** 2)
    xyz = sum(sx * sy * sz * dev(sx * h, sy * h, sz * h)
              for sx in (1, -1) for sy in (1, -1) for sz in (1, -1)) \
        / (8 * h ** 3)

    def d2(y):
        return (dev(hh, y, 0) - 2 * dev(0, y, 0) + dev(-hh, y, 0)) / hh ** 2
    x2y2 = (d2(hh) - 2 * d2(0) + d2(-hh)) / (4 * hh ** 2)
    return q2, xyz, x2y2


def _default_control():
    """Synthetic function with known monomial coefficients."""
    def g(X, Y, Z):
        return (3.0 * Z ** 2 - 7.0 * X * Y * Z + 11.0 * X ** 2 * Y ** 2
                + 0.5 * X ** 2 + 0.25 * Y ** 2 - 1.5 * X ** 2 * Z ** 2)
    return g, (3.0, -7.0, 11.0)


@dataclass
class Extraction:
    Z2: float
    XYZ: float
    X2Y2: float
    control_dev: Tuple[float, float, float]

    @property
    def total(self):
        return self.Z2 + self.XYZ + self.X2Y2


def extract(dev: Callable, h: float = 1e-5, hh: float = 1e-4,
            control=None, tol: float = 1e-6) -> Extraction:
    """The three monomial coefficients of `dev(X, Y, Z)` at the class
    X^a Y^b Z^c with a = b, a + c = 2 — i.e. Z^2, XYZ, X^2Y^2, which is
    exactly MONOMIAL_COUNT = 3 terms whenever the ratio identity xy = z
    holds (formal in i for every order of the hierarchy).

    EVERY STENCIL IS VALIDATED SEPARATELY against `control`, whose three
    coefficients are known.  `control` defaults to a synthetic function
    with known coefficients plus out-of-class contaminants; validation is
    never skipped (BI-F3: the safety property is the validation, not the
    argument being mandatory).  Passing one stencil is not passing the sum:
    the /4 regression lived in the third while the first was being checked
    against a template (PR-30 sec.3)."""
    if control is None:
        control, want = _default_control()
    else:
        control, want = control
    got = _stencils(control, h, hh)
    dev_rel = tuple(abs(g - w) / max(abs(w), 1.0) for g, w in zip(got, want))
    bad = [i for i, d in enumerate(dev_rel) if d > tol]
    if bad:
        names = ["Z^2", "XYZ", "X^2Y^2"]
        raise RuntimeError(
            "stencil control FAILED on " + ", ".join(names[i] for i in bad)
            + f" (rel dev {[f'{d:.2e}' for d in dev_rel]}) — "
              "extraction refused")
    z2, xyz, x2y2 = _stencils(dev, h, hh)
    return Extraction(z2, xyz, x2y2, dev_rel)


# ------------------------------------------------- completeness (Kimi K3)
def completeness(start: int, target: Fraction, max_len: int = 8,
                 levels: int = 18) -> List[Tuple[int, ...]]:
    """Exact-rational enumeration of closed walks from `start` whose bond
    ratio product equals `target`.  Bond (j,k) carries (j+2)/(k+2).

    Full kernel connectivity: no neighbour cutoff.  A +-3 cutoff agreed on
    every physical class tested but provably missed a jump-4 walk at target
    (6/10)^2 (BI-F4) -- and a completeness claim is a theorem, so it may not
    silently truncate.

    Donated by the auditor (Addenda AY/BC/BF).  This is what makes
    "no further terms exist" a theorem rather than an absence of
    counterexamples — and it is the only honest way to claim a monomial
    count is complete."""
    out: List[Tuple[int, ...]] = []
    tgt = Fraction(target)

    def walk(node, prod, path):
        if len(path) - 1 > max_len:
            return
        if node == start and len(path) > 1:
            if prod == tgt:
                out.append(tuple(path))
        if prod < tgt:            # every further bond only shrinks it
            return
        for nxt in range(0, levels):          # FULL connectivity (BI-F4):
            if nxt == node:                   # a +-3 cutoff silently misses
                continue                      # jump-4 walks; the prune below
                                              # already bounds the search.
            lo, hi = min(node, nxt) + 2, max(node, nxt) + 2
            walk(nxt, prod * Fraction(lo, hi), path + [nxt])

    walk(start, Fraction(1), [start])
    return out


# --------------------------------------- contamination discipline (Kimi K3)
def separation_hazard(p1: float, p2: float) -> bool:
    """True if two decay classes separate slower than HAZARD_RATIO per
    kappa — in which case a kappa-window fit is contaminated at ANY
    attainable kappa (AX.2, which cost a 1.3% phantom gap)."""
    r = abs(p1 / p2) if abs(p1) < abs(p2) else abs(p2 / p1)
    return r > HAZARD_RATIO


def power_fit(kappas: Sequence[float], devs: Sequence[float],
              classes: Sequence[float]) -> Dict:
    """Multi-component fit dev = sum_j c_j class_j^{2 kappa}, with hazard
    flags on every pair of classes that separates too slowly to trust."""
    import numpy as np
    Amat = np.array([[c ** (2 * k) for c in classes] for k in kappas],
                    dtype=float)
    coef, *_ = np.linalg.lstsq(Amat, np.array(devs, dtype=float), rcond=None)
    hazards = [(i, j) for i in range(len(classes))
               for j in range(i + 1, len(classes))
               if separation_hazard(classes[i] ** 2, classes[j] ** 2)]
    return {"coefficients": coef.tolist(), "hazards": hazards,
            "trustworthy": not hazards}


# ------------------------------------------------------------------ selftest

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
    """Targets are LEDGER LOOKUPS, not literals (BI-F2).

    The previous version hardcoded every target while its docstring
    claimed otherwise -- rebuilding exactly the stale-constant hazard the
    migration exists to kill.  Nothing here may certify against a number
    the ledger does not carry."""
    if mp is not None:
        mp.mp.dps = 40
    checks = []

    def chk(name, got, want, tol):
        checks.append((name, float(got), want, tol,
                       abs(float(got) - want) <= tol))

    for i in range(5):
        chk(f"A_{i}", A(i), _LF("A_FAMILY", i), 1e-7)
    for i in range(5):
        chk(f"C_{i}", C(i), _LF("C_FAMILY", i), 1e-5)
    chk("C3_3", C3(3), _L("C3_3"), 1e-5)
    chk("C3_4", C3(4), _L("C3_4"), 1e-5)

    chk("ratio A(3.54) dips", cancellation_ratio(channels_A(3.54)),
        0.00043617, 1e-5)
    chk("ratio C(3) > 1", cancellation_ratio(channels_C(3)), 1.47228, 1e-4)
    chk("ratio C3(3) > 1", cancellation_ratio(channels_C3(3)), 1.66393, 1e-4)

    br = scan_zeros(lambda x: float(A(x)), 1.3, 8.0, 0.25,
                    control=(lambda x: float(A(x)), 3.4, 3.7))
    chk("A has exactly one zero bracket", len(br), 1, 0)
    lo, hi = br[0]
    for _ in range(60):
        mid = (lo + hi) / 2
        if float(A(lo)) * float(A(mid)) > 0:
            lo = mid
        else:
            hi = mid
    chk("A_zero", (lo + hi) / 2, _L("A_zero"), 1e-6)
    chk("C has no zero", len(scan_zeros(lambda x: float(C(x)), 1.3, 8.0, 0.25,
                                        control=(lambda x: float(A(x)),
                                                 3.4, 3.7))), 0, 0)

    r2 = richardson(lambda h: 1.0 + 2.0 * h ** 2 + 3.0 * h ** 4,
                    [0.1, 0.05, 0.025])
    chk("richardson exact on 1+2h^2+3h^4", r2, 1.0, 1e-12)
    chk("richardson exact on pure h^2",
        richardson(lambda h: 5.0 + 7.0 * h ** 2, [0.1, 0.05, 0.025]),
        5.0, 1e-13)

    ex = extract(lambda X, Y, Z: 2.0 * Z ** 2 + 5.0 * X * Y * Z
                 - 3.0 * X ** 2 * Y ** 2)
    chk("extract Z^2", ex.Z2, 2.0, 1e-6)
    chk("extract XYZ", ex.XYZ, 5.0, 1e-5)
    chk("extract X^2Y^2", ex.X2Y2, -3.0, 1e-4)

    refused = False
    try:
        extract(lambda X, Y, Z: Z ** 2,
                control=(lambda X, Y, Z: Z ** 2, (99.0, 0.0, 0.0)))
    except RuntimeError:
        refused = True
    chk("extract refuses a failing control", 1.0 if refused else 0.0, 1.0, 0)

    w = completeness(4, Fraction(9, 16), max_len=6, levels=12)
    chk("completeness: 4 walks at the q^2 class", len(w), 4, 0)
    chk("hazard flag on r^2 vs q^2", 1.0 if separation_hazard(
        (6 / 7) ** 2, (6 / 8) ** 2) else 0.0, 0.0, 0.0)
    chk("hazard flag on (2/3)^2 vs 0.9604", 1.0 if separation_hazard(
        0.9604, 0.9795) else 0.0, 1.0, 0.0)

    ic = incident_couplings(3)
    chk("hierarchy order: 1st is up-from-upper", 1.0 if
        ic[0][2].startswith("A") else 0.0, 1.0, 0)
    chk("hierarchy order: 4th is (i+1,i+3)", 1.0 if
        ic[3][0] == "(4,6)" else 0.0, 1.0, 0)

    if verbose:
        for name, got, want, tol, ok in checks:
            print(f"[{'PASS' if ok else 'FAIL'}] {name:<34} "
                  f"{got:>18.10g} vs {want:<16.10g} tol {tol:g}")
        print(f"  {sum(c[4] for c in checks)}/{len(checks)} "
              f"expansion.py self-checks green")
    return all(c[4] for c in checks)
