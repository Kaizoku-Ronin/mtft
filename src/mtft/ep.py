"""
mtft.ep — exceptional points: extraction, winding, census, staircase
=====================================================================
MIT License — (c) 2026 Kimi K3 (MTFT independent auditor) — MTFT Research Program

Stage-3 build (Integration Plan v0.1, BH decision 1).  Absorbs the NINE
gsq() and EIGHT newton() hand copies counted in the migration census
(plan sec.2).  Built on mtft.chain: g and B come from the ONE kernel
diagonalisation in chain.internal — nothing here rebuilds T or
re-diagonalises the kernel.

THE OBJECT.  On a chain.Internal, H(u) = diag(g) - uB.  Crossings sit at
NEGATIVE u in the program convention; chain.crossing_limit(i) is |u| at
kappa -> inf and seeds every extraction.  At FINITE kappa the crossing
is AVOIDED (measured, stage-3 probe: min gap 1.019e-5 at kappa=30,
rung 0), so gsq has NO zero at the EP.  The rung's EP location is the
diabatic centre — the MINIMUM of the closest-pair gap^2 — located by
Newton on the DERIVATIVE of gsq from a shared 3-point stencil.  Never
Newton on gsq itself; never mp.diff (PR-29).

S3-1 (RESOLVED, Addenda BK/BL): the closest-approach observable
extracted here is the DIABATIC CENTRE; the sealed reduced-block
hierarchy describes the complex EP branch point.  The two differ by
exactly one copy of the within-pair repulsion constituent — the ledger
relation `dev_centre_vs_dev_EP`:

    dev_centre(i;k) = dev_EP(i;k) - rep_C(i) s^2k,
    dev_EP(i;k) = A_i r^2k - C_i s^2k + C3_i t^2k + K_i q^2k + K5_i q5^2k + ...

Derivation (BK, verified BL): in the 2x2 dressed block
H = [[a-up, -uv], [-uv, b-uq]] with D = p-q, m = b-a, the EP sits at
|u_EP| = m/sqrt(D^2+4v^2) -> u0 - 2 m v^2/D^3, while the diabatic
centre sits at u_min = mD/(D^2+4v^2) -> u0 - 4 m v^2/D^3; with
v^2 = e^{-(a+b)} s^2k (1 + O(1e-6), measured k=60) the difference is
rep_C(i) s^2k where rep_C(i) = 2 m e^{-(a+b)}/D^3 =
expansion.channels_C(i)[2].  Full-model residuals at k=60: relative
<= 1.2e-6 at rungs 1-4 — the hierarchy's own accuracy (rung 0 exempt:
its partner is the reference level lambda_0).  The selftest ASSERTS the
resolved relation; the strong-coupling LIMIT assertion stands as well.

STAIRCASE STRUCTURE (measured, stage 3): at rung i's EP the crossing
pair is ALWAYS sorted pair (0,1) — branches 0..i-1 have already crossed
above.  Closest-pair selection finds it automatically (which is why gsq
works at every rung); the rung's IDENTITY is carried by the characters
(the u=0 labels of the exchanging branches, (i, i+1)), certified by
pair_winding's diabatic-basis dominance labeling.  Sorted index !=
character index — the PR-22 rule in its final form.

LABELS COME FROM THE DIABATIC BASIS, NOT FROM FRAME-TO-FRAME TRACKING
(BK-F2, dispositioned BL): continuing labels by max-overlap between
successive u-frames ANTI-CONVERGES — finer stepping follows the
adiabatic eigenvalue branch through each veer and returns the wrong
characters (measured: rung 3 at kappa=30 flips [3,4,5,6,2,7] ->
[2,4,5,6,3,7] between 80 and 160 steps, then stays wrong).  'Character'
is defined against the FIXED strong-coupling basis, so labels are
argmax_a |V[a,j]| — parameter-free, no step knob.  A resolution-
stability gate on the old knob would catch the boundary but not the
converged-wrong regime; replacing the mechanism fixes the root.

CENSUS COUNTS EP PAIRS, NOT RUNGS (measured, stage 3): the pencil's
complex EPs attach to EVERY diabatic pair crossing, adjacent or not,
with imaginary offsets gamma spanning 1e-22 (non-adjacent) to ~0.64
(rung 4 at kappa=30).  On a wide strip (eps above every gamma) the count
equals the closed-form diabatic pair-crossing count — the robust,
assertable mode.  A tight strip is gamma-SELECTIVE (counts fewer); that
is documented behaviour, demonstrated in the selftest, not asserted
numerics.  For rung LOCATIONS use staircase().

RULES ENFORCED IN CODE (not documentation):
  1. CLOSEST-PAIR SELECTION ONLY.  No fixed-index gap is exported
     (index != identity, PR-22).
  2. PRECISION TIERS (PR-20): an f64 Internal cannot certify an EP below
     chain.F64_EP_FLOOR nor a winding bisection below
     chain.F64_WINDING_FLOOR — asking raises ValueError.
  3. NO mp.diff anywhere: explicit finite differences only (PR-29).
  4. nearest() takes levels= as a REQUIRED keyword — the caller names
     the continuation the EPs came from.  No default, no silent object.
  5. census() is argument principle on the moment discriminant
     det[tr H^{a+b}] — a count by theorem, not a scan (census != search).
     It counts EP PAIRS of the pencil in the strip (adjacent-rung EPs
     AND non-adjacent branch coalescences — see CENSUS above) and it
     RAISES unless the winding number is a clean non-negative even
     integer.
  6. selftest asserts against mtft.ledger, never literals (BI.F2).
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

from .chain import (Internal, internal, crossing_limit,
                    F64_EP_FLOOR, F64_WINDING_FLOOR)

__all__ = ["levels_of", "gsq", "closest_pair", "newton", "ep_of",
           "pair_winding", "Winding", "census", "staircase", "nearest",
           "selftest"]

NEWTON_H_MP = "1e-22"   # historical FD step (PR-28 onward)
NEWTON_H_F64 = 1e-7


# ------------------------------------------------------------- H(u) & levels
def _H(u, ic):
    """H(u) = diag(g) - uB on the Internal's own g, B (no re-diagonalisation)."""
    nb = ic.nb
    if ic.backend == "mp":
        H = mp.matrix(nb, nb)
        for a in range(nb):
            for b in range(nb):
                H[a, b] = -u * ic.B[a, b]
            H[a, a] += ic.g[a]
        return H
    return np.diag(np.asarray(ic.g, dtype=float)) - float(u) * np.asarray(ic.B)


def levels_of(u, ic, vectors: bool = False):
    """Sorted eigenvalues of H(u) — the continuation primitive.

    vectors=True also returns the sorted eigenvectors (columns), used by
    pair_winding's overlap certificate.  Uses the Internal's g and B;
    the kernel T is never touched here."""
    H = _H(u, ic)
    nb = ic.nb
    if ic.backend == "mp":
        if vectors:
            E, Z = mp.eigsy(H)
            order = sorted(range(nb), key=lambda j: E[j])
            V = mp.matrix(nb, nb)
            for c, j in enumerate(order):
                for r in range(nb):
                    V[r, c] = Z[r, j]
            return [E[j] for j in order], V
        return sorted(mp.eigsy(H, eigvals_only=True))
    if vectors:
        w, V = np.linalg.eigh(H)
        return list(w), V
    return list(np.linalg.eigvalsh(H))


# ------------------------------------------------------ gsq / closest pair
class _Shim:
    """Adapter so the historical gsq(u, g, B) signature still works."""
    def __init__(self, g, B):
        self.g = g
        self.B = B
        self.nb = len(g)
        self.backend = "mp" if (mp is not None and isinstance(g[0], mp.mpf)) \
            else "f64"


def gsq(u, g=None, B=None, ic=None):
    """Squared gap of the CLOSEST pair of H(u) levels.

    Closest-pair selection ONLY: there is deliberately no fixed-index
    variant to import (index != identity, PR-22).  Call as gsq(u, ic=ic)
    or in the historical form gsq(u, g, B)."""
    if ic is None:
        if g is None or B is None:
            raise TypeError("gsq needs either ic= or both g and B")
        ic = _Shim(g, B)
    lev = levels_of(u, ic)
    return min((lev[k + 1] - lev[k]) ** 2 for k in range(len(lev) - 1))


def closest_pair(u, ic) -> Tuple[int, Any]:
    """(index, gap) of the closest pair of H(u) levels — the identity
    check every EP claim must pass (index != identity, PR-22)."""
    lev = levels_of(u, ic)
    k = min(range(len(lev) - 1), key=lambda j: lev[j + 1] - lev[j])
    return k, lev[k + 1] - lev[k]


# ------------------------------------------------------------------ newton
def newton(f, x0, h=None, tol=None, maxiter: int = 100):
    """Newton iteration with an EXPLICIT finite-difference derivative.

    mp.diff is forbidden on implicit observables (PR-29: it returned
    exact zeros on a root-defined observable).  h defaults to the
    historical 1e-22 at mp precision (1e-7 in float64).  Absorbs the
    eight hand copies."""
    use_mp = mp is not None and isinstance(x0, mp.mpf)
    if h is None:
        h = mp.mpf(NEWTON_H_MP) if use_mp else NEWTON_H_F64
    if tol is None:
        tol = mp.mpf("1e-30") if use_mp else 1e-12
    x = x0
    fx = f(x)
    for _ in range(maxiter):
        dfx = (f(x + h) - fx) / h
        if dfx == 0:
            raise ArithmeticError("newton: zero derivative — no root in reach")
        dx = fx / dfx
        x = x - dx
        fx = f(x)
        if abs(dx) < tol:
            return x
    raise ArithmeticError(f"newton: no convergence in {maxiter} iterations")


# ------------------------------------------------------------------ ep_of
def ep_of(i: int, ic, tol=None, seed=None, h=None, maxiter: int = 80):
    """The rung-i EP at THIS kappa: the diabatic centre of the avoided
    crossing — minimum of the closest-pair gap^2 (see module docstring).

    OBSERVABLE (S3-1, resolved BK/BL): this returns the DIABATIC CENTRE,
    not the complex EP branch point; the two differ by exactly
    -rep_C(i) s^2k (ledger relation 'dev_centre_vs_dev_EP').  The sealed
    hierarchy (A, C, C3, K, K5) is the EP's; this observable's finite-k
    dev is the sealed value minus one copy of the repulsion constituent.

    Newton on gsq'(u) from a shared 3-point stencil (explicit FD only),
    seeded at -crossing_limit(i), with one h/10 refinement pass.  The
    seed sign is the program convention: crossings sit at NEGATIVE u
    (getting this wrong once cost the ledger's degen checks, stage 1).

    PR-20 tier rule: an f64 Internal cannot certify below F64_EP_FLOOR."""
    if tol is None:
        tol = mp.mpf("1e-22") if ic.backend == "mp" else 1e-10
    if ic.backend == "f64" and tol < F64_EP_FLOOR:
        raise ValueError(
            f"tol {tol} < F64_EP_FLOOR {F64_EP_FLOOR}: an f64 Internal "
            "caps any EP measurement at ~1e-15 regardless of solver "
            "precision (PR-20).  Rebuild with backend='mp'.")
    if seed is None:
        seed = -crossing_limit(i, backend=ic.backend)
    if ic.backend == "mp":
        u = mp.mpf(seed)
        h0 = mp.mpf(h) if h is not None else max(abs(u), 1) * mp.mpf("1e-6")
    else:
        u = float(seed)
        h0 = float(h) if h is not None else max(abs(u), 1.0) * 1e-6

    def G(x):
        return gsq(x, ic=ic)

    hh = h0
    for _pass in range(2):                      # h, then h/10 refinement
        for _ in range(maxiter):
            fm, f0, fp = G(u - hh), G(u), G(u + hh)
            d1 = (fp - fm) / (2 * hh)
            d2 = (fp - 2 * f0 + fm) / hh ** 2
            if d2 <= 0:
                raise ArithmeticError(
                    "ep_of: curvature went non-positive — the seed is "
                    "outside the rung's basin (wrong rung nearby?)")
            step = d1 / d2
            u = u - step
            if abs(step) < tol:
                break
        else:
            raise ArithmeticError(f"ep_of: no convergence in {maxiter} steps")
        hh = hh / 10
    return u


# ------------------------------------------------------------ pair winding
@dataclass
class Winding:
    """Certificate that one pair, and only that pair, exchanges character
    across a rung's EP.

    `pair` is the SORTED index of the crossing pair at the EP (found by
    closest-pair selection — at every rung measured this is (0,1); see
    STAIRCASE STRUCTURE in the module docstring).  `characters` is the
    rung's identity: the diabatic-basis labels of the exchanging
    branches — for rung i this must read (i, i+1)."""
    pair: Tuple[int, int]
    characters: Tuple[int, int]
    u_swap: Any
    exchanged: bool
    consistent: bool          # no label crosses the pair/spectator boundary
    spectator_swaps: int      # spectator slots that changed label (real
                              # pencil structure inside wide windows, BL)
    spectator_leak: float     # max |overlap| across the pair/spectator
                              # boundary ONLY (spectator-spectator events
                              # are not contamination, BL)
    overlap_diag: float
    overlap_offdiag: float


def _overlap(ic, Va, a, Vb, b):
    if ic.backend == "mp":
        return sum(Va[r, a] * Vb[r, b] for r in range(ic.nb))
    return float(np.dot(Va[:, a], Vb[:, b]))


def _dominance_labels(ic, u) -> Tuple[List[int], List[float]]:
    """Diabatic-basis character labels of the sorted levels at u:
    label[j] = argmax_a |V[a, j]| — the strong-coupling basis vector
    dominating instantaneous eigenvector j.  Returns (labels, dom) with
    dom[j] = w_best / w_second (character purity of that label).

    PARAMETER-FREE replacement for frame-to-frame continuation (BK-F2):
    marching labels from u=0 by max-overlap ANTI-CONVERGES — finer
    stepping resolves each veer's rotation and follows the adiabatic
    branch, flipping to the wrong characters (rung 3, kappa=30: correct
    at <=80 steps, wrong at >=160, converged-wrong thereafter).  The
    diabatic basis is the fixed reference that DEFINES character; at
    u=0 the eigenvectors ARE the basis vectors, so the two notions
    agree exactly where continuation starts."""
    _, V = levels_of(u, ic, vectors=True)
    labels: List[int] = []
    doms: List[float] = []
    for j in range(ic.nb):
        if ic.backend == "mp":
            w = [float(abs(V[a, j]) ** 2) for a in range(ic.nb)]
        else:
            col = np.asarray(V, dtype=float)[:, j]
            w = [float(x * x) for x in col]
        order = sorted(range(ic.nb), key=lambda a: -w[a])
        labels.append(order[0])
        doms.append(w[order[0]] / max(w[order[1]], 1e-300))
    return labels, doms


def pair_winding(ic, i: Optional[int] = None, u_ep=None,
                 half_width: float = 0.05, swap_tol: float = 1e-6,
                 max_bisect: int = 100, max_widen: int = 5) -> Winding:
    """Prove WHICH pair crosses at a rung, by eigenvector character.

    Two independent witnesses:
      (1) endpoint overlap — the crossing pair's eigenvectors at
          u_ep +/- half_width form the EXCHANGE matrix (off-diagonal ~ 1)
          while every spectator stays diagonal;
      (2) diabatic-basis labels — the characters (dominant strong-
          coupling basis vectors) of the sorted crossing pair SWAP
          between the endpoints while every spectator label is
          unchanged.  Labels are PARAMETER-FREE (BK-F2: frame-to-frame
          tracking anti-converges; see _dominance_labels).
    The swap point u_swap is bisected on the tracked character weight
    s(u) = |<v_{j+1}(u) | v_j(near)>|^2 crossing 1/2.  The sorted pair
    j itself is FOUND by closest-pair selection, never assumed.

    PURITY GATES AND REFUSAL.  An attempt is character-pure iff
    s_near <= 0.02, s_far >= 0.98 (the pair fully swaps across the
    window), the crossing pair's own labels dominate by >= 4:1, and the
    label sets are bijective.  If an attempt fails, half_width doubles
    while s_far < 0.9 (up to max_widen times), then refines by x1.25
    (up to 3 times) — doubling alone OVERSHOOTS crowded windows
    (measured: rung 3 at kappa=30 purifies to 0.985 only at hw=1.25,
    and hw=1.6 lands on the NEXT diabatic crossing).  Windows are
    capped at |u_ep|/2 so endpoints stay on their sides of the origin.
    If EVERY attempt fails the refusal distinguishes the two diagnoses
    by the s_far trajectory (BK-F1): s_far monotone rising -> a genuine
    half_width problem, message says increase half_width/max_widen;
    s_far NON-MONOTONE -> overlapping veers, NO character-pure endpoint
    exists at this coupling, a regime limit of the observable (small
    kappa), message says so.  It never mis-certifies, and never
    misdiagnoses which failure it is.

    SPECTATORS (measured structure, BL): wide windows REGULARLY contain
    spectator diabatic crossings (rung 3's clean window spans the
    (2,6) crossing 0.064 above u* and the (2,7) scramble near u*-1.2)
    — the census counts such crossings, they are real pencil structure,
    not contamination.  So `leak` measures only PAIR<->SPECTATOR
    mixing, `consistent` means no label crosses the pair/spectator
    BOUNDARY, and `spectator_swaps` reports how many spectator slots
    changed label (informational — expected nonzero at crowded rungs).

    PR-20 tier rule: an f64 Internal cannot bisect below
    F64_WINDING_FLOOR."""
    if ic.backend == "f64" and swap_tol < F64_WINDING_FLOOR:
        raise ValueError(
            f"swap_tol {swap_tol} < F64_WINDING_FLOOR "
            f"{F64_WINDING_FLOOR}: winding bisections floor at ~2e-4 in "
            "f64 (PR-20).  Rebuild with backend='mp'.")
    if u_ep is None:
        if i is None:
            raise ValueError("pair_winding needs i= or u_ep=")
        u_ep = ep_of(i, ic)
    to = mp.mpf if ic.backend == "mp" else float
    u_ep, w = to(u_ep), to(half_width)
    nb = ic.nb

    j, _ = closest_pair(u_ep, ic)            # FOUND, not assumed (PR-22)
    # `near` is the endpoint continuation from u=0 reaches first
    # (characters live there), `far` the one past the crossing.
    # Crossings sit at negative u in this program.
    accepted = None
    s_hist: List[Tuple[float, float, float]] = []
    n_double = 0
    n_refine = 0
    while True:
        near = u_ep + w if u_ep < 0 else u_ep - w
        far = u_ep - w if u_ep < 0 else u_ep + w
        _, V_near = levels_of(near, ic, vectors=True)
        _, V_far = levels_of(far, ic, vectors=True)

        def s(u, _Vn=V_near):
            """Tracked character weight |<v_{j+1}(u) | v_j(near)>|^2 —
            0 at near, 1 at far, through the veer (endpoint-referenced:
            no frame marching, BK-F2)."""
            _, V = levels_of(u, ic, vectors=True)
            o = _overlap(ic, V, j + 1, _Vn, j)
            return o * o

        s_near, s_far = s(near), s(far)
        s_hist.append((float(w), float(s_near), float(s_far)))
        lab_near, dom_near = _dominance_labels(ic, near)
        lab_far, _ = _dominance_labels(ic, far)
        bijective = (len(set(lab_near)) == nb and len(set(lab_far)) == nb)
        pair_dom = min(dom_near[j], dom_near[j + 1])
        if (s_near <= 0.02 and s_far >= 0.98 and bijective
                and pair_dom >= 4.0):
            accepted = (near, far, V_near, V_far, s, lab_near, lab_far, w)
            break
        # widen: double while far from purity, refine x1.25 near it —
        # pure doubling overshoots crowded windows (rung 3 purifies at
        # hw=1.25; hw=1.6 lands on the next diabatic crossing, BL).
        if float(s_far) < 0.9 and n_double < max_widen:
            w_next = 2 * w
        elif n_refine < 3:
            w_next = to(1.25) * w
            n_refine += 1
        else:
            break
        if w_next >= abs(u_ep) / 2:
            break
        if w_next == 2 * w:
            n_double += 1
        w = w_next
    if accepted is None:
        hist = ", ".join(f"hw={hw:g} -> s_near={sn:.4f}, s_far={sf:.4f}"
                         for hw, sn, sf in s_hist)
        sf_seq = [sf for _, _, sf in s_hist]
        monotone_up = all(b >= a - 0.02 for a, b in zip(sf_seq, sf_seq[1:]))
        if monotone_up:
            raise ValueError(
                "pair_winding: endpoints still not character-pure but "
                f"s_far IS converging ({hist}) — the window has not yet "
                "cleared the veer; increase half_width or max_widen "
                "(this IS a half_width problem).")
        raise ValueError(
            "pair_winding: NO CHARACTER-PURE ENDPOINT EXISTS at this "
            f"coupling ({hist}): s_far is NON-MONOTONE in half_width — "
            "the veers overlap and diabatic characters are ill-defined "
            "at every scale tried.  This is a REGIME limit (kappa too "
            "small for the winding certificate), not a half_width "
            "problem; refusing rather than mis-certifying.")
    near, far, V_near, V_far, s, lab_near, lab_far, w = accepted

    o_diag = (abs(_overlap(ic, V_near, j, V_far, j))
              + abs(_overlap(ic, V_near, j + 1, V_far, j + 1))) / 2
    o_off = (abs(_overlap(ic, V_near, j, V_far, j + 1))
             + abs(_overlap(ic, V_near, j + 1, V_far, j))) / 2
    pair_idx = {j, j + 1}
    leak = max(
        (abs(_overlap(ic, V_near, a, V_far, b))
         for a in range(nb) for b in range(nb)
         if (a in pair_idx) != (b in pair_idx)),
        default=0.0)

    characters = (lab_near[j], lab_near[j + 1])
    # decisive exchange: the exchange matrix is off-diagonal by a wide
    # margin AND the diabatic-basis labels swap exactly.
    exchanged = bool(o_off > 0.95 and o_off > 5 * o_diag
                     and lab_far[j] == lab_near[j + 1]
                     and lab_far[j + 1] == lab_near[j])
    # boundary semantics: the pair's characters stay inside the pair's
    # slots and no spectator label enters them.  Spectator-spectator
    # rearrangements are REAL pencil structure inside wide windows
    # (measured BL: rung 3's clean window spans the (2,6) diabatic
    # crossing) — reported as spectator_swaps, not gated as leak.
    consistent = bool(
        set(lab_far[a] for a in pair_idx) == set(characters)
        and all(lab_far[a] not in characters
                for a in range(nb) if a not in pair_idx))
    spec_swaps = sum(1 for a in range(nb)
                     if a not in pair_idx and lab_far[a] != lab_near[a])

    # bisect the character-swap point: t in [0,1], u = near + t(far-near)
    t_lo, t_hi = (mp.mpf(0), mp.mpf(1)) if ic.backend == "mp" else (0.0, 1.0)
    for _ in range(max_bisect):
        if (t_hi - t_lo) * w < swap_tol:
            break
        t_mid = (t_lo + t_hi) / 2
        if s(near + t_mid * (far - near)) < 0.5:
            t_lo = t_mid
        else:
            t_hi = t_mid
    t_swap = (t_lo + t_hi) / 2
    return Winding(pair=(j, j + 1), characters=characters,
                   u_swap=near + t_swap * (far - near), exchanged=exchanged,
                   consistent=consistent, spectator_swaps=int(spec_swaps),
                   spectator_leak=float(leak),
                   overlap_diag=float(o_diag), overlap_offdiag=float(o_off))


# --------------------------------------------------- argument-principle census
def _disc(u, ic):
    """Moment discriminant det[tr H^{a+b}]_{a,b=0..n-1} = prod (l_a-l_b)^2.

    Valid at COMPLEX u — no eigensolve, so no symmetry requirement.  The
    identity is det(V^T V) for the Vandermonde V of eigenvalues; the
    selftest verifies it against the direct product on the real axis.
    Census always runs at mp precision (it is a counting theorem)."""
    n = ic.nb
    if ic.backend == "mp":
        g, B = ic.g, ic.B
    else:
        g = [mp.mpf(float(x)) for x in ic.g]
        B = [[mp.mpf(float(ic.B[a, b])) for b in range(n)] for a in range(n)]
    H = mp.matrix(n, n)
    for a in range(n):
        for b in range(n):
            H[a, b] = -u * (B[a, b] if ic.backend == "mp" else B[a][b])
        H[a, a] += g[a]
    P = mp.eye(n)
    trs = [mp.mpf(n)]
    for _ in range(1, 2 * n - 1):
        P = P * H
        trs.append(sum(P[k, k] for k in range(n)))
    M = mp.matrix(n, n)
    for a in range(n):
        for b in range(n):
            M[a, b] = trs[a + b]
    return mp.det(M)


def census(ic, u_lo: float, u_hi: float, eps: float = 1.0,
           samples_per_unit: int = 20) -> Dict:
    """Count EP PAIRS of the pencil H(u) in the strip
    [u_lo, u_hi] x [-i eps, +i eps], by the argument principle on
    disc(u) — census, not search.

    This counts the EPs of EVERY branch pair whose coalescence sits in
    the strip, not only adjacent-rung EPs (see CENSUS in the module
    docstring): each EP contributes winding 2 (a double zero at an exact
    real crossing, a conjugate simple-zero pair at complex EPs).  On a
    strip wider than every gamma the count equals the closed-form
    diabatic pair-crossing count; a tighter strip is gamma-selective.
    Choose bounds that do not sit on a crossing.  RAISES unless the
    winding is a clean non-negative even integer — a non-integer means
    the contour was under-sampled or clipped a zero, and a silent count
    would be a lie.  For rung LOCATIONS use staircase().
    """
    if mp is None:                                            # pragma: no cover
        raise RuntimeError("census requires mpmath")
    if not u_lo < u_hi:
        raise ValueError("census needs u_lo < u_hi")
    n_h = max(16, int(samples_per_unit * (u_hi - u_lo)))
    n_v = max(8, int(samples_per_unit * 2 * eps))
    pts = ([mp.mpc(u_lo + (u_hi - u_lo) * k / n_h, -eps) for k in range(n_h)]
           + [mp.mpc(u_hi, -eps + 2 * eps * k / n_v) for k in range(n_v)]
           + [mp.mpc(u_hi - (u_hi - u_lo) * k / n_h, eps) for k in range(n_h)]
           + [mp.mpc(u_lo, eps - 2 * eps * k / n_v) for k in range(n_v)])

    total = mp.mpf(0)
    prev = _disc(pts[0], ic)
    for p in pts[1:] + [pts[0]]:
        cur = _disc(p, ic)
        if cur == 0 or prev == 0:
            raise RuntimeError("census: a zero lies ON the contour — "
                               "move the bounds")
        total += mp.arg(cur / prev)
        prev = cur
    w = total / (2 * mp.pi)
    W = int(round(float(w)))
    if abs(float(w) - W) > 0.1:
        raise RuntimeError(
            f"census: winding {mp.nstr(w, 6)} is not near-integer — "
            "contour under-sampled or clipped a zero; count refused")
    if W < 0 or W % 2:
        raise RuntimeError(
            f"census: winding {W} is negative or odd — a zero lies on or "
            "too near the contour; count refused")
    out = {"winding": W, "ep_pairs": W // 2, "raw": float(w)}
    if out["ep_pairs"] == 0:
        out["note"] = (
            "winding 0 — EITHER no EP pairs in this strip OR eps sits "
            "below every gamma (census is gamma-selective, BJ): widen "
            "eps or tighten [u_lo, u_hi] before concluding absence.  "
            "A zero count is not evidence of absence (BK-F3).")
    return out


# ----------------------------------------------------------- staircase
def staircase(ic, rungs: int = 5, tol=None) -> List:
    """The EP staircase: rungs 0..rungs-1 located by ep_of in rung order
    (continuation in rung index, each seeded from its own limit)."""
    return [ep_of(i, ic, tol=tol) for i in range(rungs)]


def nearest(u, *, levels: Sequence):
    """(index, value) of the EP nearest u.

    levels= is a REQUIRED keyword: the caller names the continuation the
    EPs came from.  A default would silently answer about the wrong
    object (index != identity; perturb the object the quantity is
    defined from)."""
    lv = list(levels)
    if not lv:
        raise ValueError("nearest: levels= is empty")
    k = min(range(len(lv)), key=lambda j: abs(lv[j] - u))
    return k, lv[k]


# ---------------------------------------------------------------- selftest
def selftest(verbose: bool = True):
    """Assert against mtft.ledger — never literals (BI.F2).

    Sharp checks assert the strong-coupling LIMITS (shared by every EP
    observable) AND the resolved S3-1 relation (dev_centre = sealed
    EP hierarchy - rep_C s^2k, k=60, rungs 1-4); the k=30 dev table is
    REPORTED and envelope-guarded as a wrong-basin check."""
    from . import ledger as L
    if mp is None:                                            # pragma: no cover
        raise RuntimeError("ep.selftest requires mpmath")
    mp.mp.dps = 50
    checks = []

    def chk(name, got, want, tol):
        ok = abs(got - want) <= tol
        checks.append((name, float(got), float(want), float(tol), bool(ok)))

    # newton (absorbed copies), both tiers
    chk("newton sqrt2 (mp)", abs(newton(lambda x: x * x - 2, mp.mpf("1.4"))
                                 - mp.sqrt(2)), 0, mp.mpf("1e-40"))
    chk("newton x^3-x-2 (f64)", abs(newton(lambda x: x ** 3 - x - 2, 1.5)
                                    - 1.5213797068045676), 0, 1e-9)

    ic30 = internal(30, nb=10, backend="mp", nsite=12)

    # continuation anchored at u=0 (positive control)
    dmax = max(abs(a - b) for a, b in zip(levels_of(mp.mpf(0), ic30), ic30.g))
    chk("levels_of(0) == ic.g", dmax, 0, mp.mpf("1e-35"))

    # historical gsq(u, g, B) signature preserved
    v1 = gsq(mp.mpf("-1.425"), ic=ic30)
    v2 = gsq(mp.mpf("-1.425"), ic30.g, ic30.B)
    chk("gsq(u,g,B) == gsq(u,ic)", abs(v1 - v2), 0, mp.mpf("1e-40"))

    # gsq dips at a crossing seed (positive control — a vacuous gsq fails this)
    chk("gsq dips at rung-0 seed",
        gsq(-crossing_limit(0), ic=ic30) / gsq(mp.mpf(0), ic=ic30), 0, 1e-3)

    # SHARP: rungs converge to the ledger limits (shared by all observables)
    ic100 = internal(100, nb=10, backend="mp", nsite=14)
    keys = ["u_ep_01", "u_ep_12", "u_ep_23", "u_ep_34", "u_ep_45"]
    for i, k in enumerate(keys):
        u = ep_of(i, ic100)
        chk(f"ep_of rung {i} (k=100) -> ledger {k}",
            abs(abs(u) - mp.mpf(float(L.LEDGER[k].value))), 0, mp.mpf("1e-7"))

    # staircase structure: at every rung the closest SORTED pair is (0,1),
    # found by closest-pair selection (index != identity, exercised)
    dev30 = []
    for i in range(5):
        u = ep_of(i, ic30)
        j, gp = closest_pair(u, ic30)
        chk(f"rung {i}: closest sorted pair is (0,1)", j, 0, 0)
        chk(f"rung {i}: closest gap << spacing", gp, 0, 5e-2)
        dev30.append((i, u, abs(u) - mp.mpf(float(L.LEDGER[keys[i]].value))))

    # loose wrong-basin guard on the k=30 dev (hierarchy envelope x 2.5;
    # the RESOLVED S3-1 relation itself is asserted at k=60 below)
    for i, u, dev in dev30:
        kap = 30
        r = (i + 3) / (i + 4)
        s_ = (i + 2) / (i + 3)
        t_ = (i + 1) / (i + 2)
        q = (i + 3) / (i + 5)
        q5 = (i + 2) / (i + 4)
        env = abs(float(L.A_of(i))) * r ** (2 * kap) \
            + abs(float(L.C_of(i))) * s_ ** (2 * kap)
        if i >= 1:
            env += abs(float(L.C3_of(i))) * t_ ** (2 * kap)
        for key, base in ((f"K_{i}", q), (f"K5_{i}", q5)):
            if key in L.LEDGER:
                env += abs(float(L.LEDGER[key].value)) * base ** (2 * kap)
        chk(f"dev rung {i} (k=30) within 2.5x envelope", abs(dev), 0,
            2.5 * env + 1e-9)

    # S3-1 RESOLVED (BK derivation, BL verification): the sealed
    # hierarchy describes the EP; ep_of's DIABATIC CENTRE differs by
    # exactly one copy of the repulsion constituent — ledger relation
    # 'dev_centre_vs_dev_EP'.  Assert the resolved relation on the full
    # model at k=60, rungs 1-4 (rung 0 exempt: its partner is the
    # reference level lambda_0).  All coefficients from the ledger and
    # expansion closed forms — never literals.
    from .expansion import channels_C
    ic60 = internal(60, nb=10, backend="mp", nsite=12)
    kap60 = 60
    for i in range(1, 5):
        u60 = ep_of(i, ic60, tol=mp.mpf("1e-10"))
        dev60 = abs(u60) - crossing_limit(i, backend="mp")
        r6 = mp.mpf(i + 3) / (i + 4)
        s6 = mp.mpf(i + 2) / (i + 3)
        t6 = mp.mpf(i + 1) / (i + 2)
        q6 = mp.mpf(i + 3) / (i + 5)
        q56 = mp.mpf(i + 2) / (i + 4)
        pred = L.A_of(i) * r6 ** (2 * kap60) - L.C_of(i) * s6 ** (2 * kap60) \
            + L.C3_of(i) * t6 ** (2 * kap60)
        for key, base in ((f"K_{i}", q6), (f"K5_{i}", q56)):
            if key in L.LEDGER:
                pred += mp.mpf(L.LEDGER[key].value) * base ** (2 * kap60)
        pred -= channels_C(i)[2] * s6 ** (2 * kap60)   # centre = EP - rep_C s^2k
        chk(f"S3-1 (resolved): dev rung {i} (k=60) == sealed - rep_C s^2k",
            abs(float((dev60 - pred) / pred)), 0, 1e-4)

    # winding certificate at rung 0 — both witnesses
    w0 = pair_winding(ic30, i=0)
    chk("winding: sorted pair is (0,1)", 1 if w0.pair == (0, 1) else 0, 1, 0)
    chk("winding: characters are (0,1)", 1 if w0.characters == (0, 1) else 0,
        1, 0)
    chk("winding: pair exchanged", 1 if w0.exchanged else 0, 1, 0)
    chk("winding: spectators keep labels", 1 if w0.consistent else 0, 1, 0)
    chk("winding: spectators stay diagonal", w0.spectator_leak, 0, 1e-2)
    chk("winding: u_swap vs ledger u_ep_01",
        abs(abs(w0.u_swap) - mp.mpf(float(L.LEDGER["u_ep_01"].value))), 0,
        mp.mpf("1e-5"))

    # index != identity, exercised at rung 2: sorted pair is STILL (0,1),
    # but the characters must read (2,3).  half_width=0.15 clears rung
    # 2's wider veer (avoided gap 2.4e-3 at kappa=30).
    w2 = pair_winding(ic30, i=2, half_width=0.25)
    chk("winding rung 2: sorted pair is (0,1)", 1 if w2.pair == (0, 1) else 0,
        1, 0)
    chk("winding rung 2: characters are (2,3)",
        1 if w2.characters == (2, 3) else 0, 1, 0)
    chk("winding rung 2: exchanged", 1 if w2.exchanged else 0, 1, 0)

    # rung 3 — the BK-F2 case.  Frame-to-frame tracking returned the
    # WRONG characters here at fine resolution (anti-convergence, flip
    # between 80 and 160 steps, converged-wrong thereafter); diabatic-
    # basis labeling is parameter-free and reads the character-faithful
    # answer.  The clean window is crowded: spectator_swaps reports the
    # (2,6) diabatic crossing etc. inside it — real pencil structure,
    # while the boundary witness (consistent) and the pair witnesses
    # stay decisive.
    w3 = pair_winding(ic30, i=3)
    chk("winding rung 3: sorted pair is (0,1)", 1 if w3.pair == (0, 1) else 0,
        1, 0)
    chk("winding rung 3: characters are (3,4) (BK-F2 fixed at root)",
        1 if w3.characters == (3, 4) else 0, 1, 0)
    chk("winding rung 3: exchanged (o_off/o_diag decisive)",
        1 if w3.exchanged else 0, 1, 0)
    chk("winding rung 3: no label crosses the pair boundary",
        1 if w3.consistent else 0, 1, 0)
    chk("winding rung 3: pair-spectator leak small",
        w3.spectator_leak, 0, 0.15)

    # BK-F1: two genuine instances of the overlapping-veers regime must
    # refuse with the REGIME diagnosis (s_far non-monotone), never the
    # misdiagnosing 'increase half_width'.
    #   (a) rung 4 at kappa=30: windows wide enough to purify swallow
    #       rung 3's veer (s_far collapses 0.82 -> 0.007, measured BL).
    #   (b) rung 1 at kappa*=5 (BK's probe): s_far 0.16/0.64/0.002/0.12
    #       across half_widths 0.05..1.0 — no character-pure endpoint
    #       exists at any scale.
    try:
        pair_winding(ic30, i=4, half_width=0.5)
        chk("BK-F1a: rung 4 (k=30) winding refuses", 0, 1, 0)
    except ValueError as e4:
        chk("BK-F1a: rung 4 (k=30) winding refuses", 1, 1, 0)
        chk("BK-F1a: refusal names the overlapping-veers regime",
            1 if "NO CHARACTER-PURE ENDPOINT" in str(e4) else 0, 1, 0)
    ic5 = internal(5, nb=14, backend="mp")
    try:
        pair_winding(ic5, i=1)
        chk("BK-F1b: kappa*=5 winding refuses (BK probe)", 0, 1, 0)
    except ValueError as e5:
        chk("BK-F1b: kappa*=5 winding refuses (BK probe)", 1, 1, 0)
        chk("BK-F1b: refusal names the regime, not half_width",
            1 if "NO CHARACTER-PURE ENDPOINT" in str(e5) else 0, 1, 0)

    # moment-discriminant identity (census core, per-quantity check)
    u0 = mp.mpf("-3.7")
    lev = levels_of(u0, ic30)
    d_eig = mp.mpf(1)
    for a in range(ic30.nb):
        for b in range(a + 1, ic30.nb):
            d_eig *= (lev[b] - lev[a]) ** 2
    chk("moment disc == prod (l_a-l_b)^2",
        abs((_disc(u0, ic30) - d_eig) / d_eig), 0, mp.mpf("1e-25"))

    # census: EP-pair counts asserted against the CLOSED-FORM diabatic
    # pair-crossing count (search-free expectation; wide-strip plateaus)
    ic8 = internal(30, nb=8, backend="mp", nsite=12)   # cheaper contour math

    def diabatic_pairs(ic, lo, hi):
        n = 0
        for a in range(ic.nb):
            for b in range(a + 1, ic.nb):
                u = -(ic.g[b] - ic.g[a]) / (ic.B[a, a] - ic.B[b, b])
                if lo < float(u) < hi:
                    n += 1
        return n

    c1 = census(ic8, -13.0, 0.0, eps=1.0)
    chk("census(-13,0,eps=1.0) == diabatic count", c1["ep_pairs"],
        diabatic_pairs(ic8, -13.0, 0.0), 0)
    c2 = census(ic8, -4.0, 0.0, eps=0.3)
    chk("census(-4,0,eps=0.3) == diabatic count", c2["ep_pairs"],
        diabatic_pairs(ic8, -4.0, 0.0), 0)
    c3 = census(ic8, -6.0, -4.0, eps=0.05)
    chk("census(-6,-4,eps=0.05) == diabatic count", c3["ep_pairs"],
        diabatic_pairs(ic8, -6.0, -4.0), 0)
    c4 = census(ic8, -4.0, 0.0, eps=0.05)
    chk("census: tight strip is gamma-selective",
        1 if c4["ep_pairs"] < c2["ep_pairs"] else 0, 1, 0)
    # BK-F3: a gamma-excluded zero count must SAY it is not evidence of
    # absence.  Probe: (-6,-4) at eps=0.005 sits below every gamma in
    # that window (BK probe P3), so the strip returns 0 pairs + note.
    c5 = census(ic8, -6.0, -4.0, eps=0.005)
    chk("BK-F3: gamma-excluded strip returns 0 pairs", c5["ep_pairs"], 0, 0)
    chk("BK-F3: zero count carries the not-absence note",
        1 if "note" in c5 else 0, 1, 0)

    # staircase + the mandatory levels= rule
    st = staircase(ic100, 3)
    for i in range(3):
        chk(f"staircase rung {i} -> ledger",
            abs(abs(st[i]) - mp.mpf(float(L.LEDGER[keys[i]].value))), 0, 1e-7)
    idx, val = nearest(mp.mpf("-5.0"), levels=st)
    chk("nearest(-5.0, levels=stair) -> rung 2", idx, 2, 0)
    try:
        nearest(-5.0)                                   # must raise TypeError
        chk("nearest() requires levels=", 0, 1, 0)
    except TypeError:
        chk("nearest() requires levels=", 1, 1, 0)

    # PR-20 tier rules raise
    icf = internal(30, nb=10, backend="f64", N=400)
    try:
        ep_of(0, icf, tol=1e-15)
        chk("f64 EP floor enforced (PR-20)", 0, 1, 0)
    except ValueError:
        chk("f64 EP floor enforced (PR-20)", 1, 1, 0)
    try:
        pair_winding(icf, i=0, swap_tol=1e-6)
        chk("f64 winding floor enforced (PR-20)", 0, 1, 0)
    except ValueError:
        chk("f64 winding floor enforced (PR-20)", 1, 1, 0)
    wf = pair_winding(icf, i=0, swap_tol=1e-3)          # f64 happy path
    chk("f64 winding smoke: exchanged", 1 if wf.exchanged else 0, 1, 0)

    if verbose:
        for name, got, want, tol, ok in checks:
            print(f"[{'PASS' if ok else 'FAIL'}] {name:<40} "
                  f"{got:>16.6g} vs {want:<12.6g} tol {tol:g}")
        n_ok = sum(c[4] for c in checks)
        print(f"  {n_ok}/{len(checks)} ep.py self-checks green")
        print("\n  k=30 finite-kappa dev report (S3-1 — REPORTED, not "
              "asserted; sealed hierarchy is the reduced-block observable):")
        for i, u, dev in dev30:
            print(f"    rung {i}: u* = {mp.nstr(u, 15)}  dev = {mp.nstr(dev, 4)}")
    return all(c[4] for c in checks)


if __name__ == "__main__":
    selftest()
