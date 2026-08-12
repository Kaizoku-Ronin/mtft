#!/usr/bin/env python3
"""
arithmetic_area_geometry.py — Fisher volume as arithmetic area
===============================================================

MIT License — Copyright (c) 2026 Roger Tano

Origin: an exploration proposed by an outside engine after reading
v0.14.0, audited here and promoted with four corrections.  The
exploration's own seven checks reproduce against the v0.14.0 tree;
what follows re-derives them, tests the central theorem OUT OF SAMPLE
on supports it was never fitted to, and files one correction against
the package.

THE STATISTIC PLANE.  Attach to each integer the sufficient-statistic
point X_n = (log n, w_n) with w_n = sum_{d|n} (log d)/d.  Since every
nontrivial proper divisor contributes a positive term,

    w_n >= (log n)/n,  with equality iff n = 1 or n is prime,

so y = x e^{-x} is an exact PRIME BOUNDARY in the plane and
compositeness is literally vertical displacement above it.  (The n = 1
case matters: X_1 = (0, 0) lies ON the curve and 1 is not prime.  The
exploration's boxed statement omitted it; A6 states the corrected
form.)

TWO IDENTITIES.  For any finite support with Gibbs weights p_i,
writing Delta_ijk for the oriented doubled area of the triangle
X_i X_j X_k and C_ij = det(X_i - mu, X_j - mu) about the Gibbs
centroid mu,

  A1   det g = sum_{i<j<k} p_i p_j p_k Delta_ijk^2
  A2   K = -N / (4 det(g)^2),
       N = sum_{i<j<k} p_i p_j p_k Delta_ijk C_ij C_ik C_jk

so the Fisher volume is a weighted sum of squared arithmetic triangle
areas, and the curvature is the normalized oriented-area chirality of
those triangles about the centroid.  A1 is the classical generalized-
variance/simplex-volume identity specialized to the arithmetic point
cloud; A2 is checked here against the INDEPENDENT Brioschi route in
mtft.curvature, which shares no step with the cumulant determinant.

THE RIGIDITY THEOREM (A3, proved; this replaces a mechanism sketch).
Let a two-dimensional discrete exponential family have support
consisting of any number (>= 2) of statistic points on a single affine
line, together with EXACTLY ONE point off that line, with arbitrary
positive base weights.  Then K == 1/4 identically on the whole
parameter manifold.

  Proof.  Gaussian curvature is a diffeomorphism invariant and the
  parameter change induced by an affine change of the statistic plane
  is linear, so normalize the line to y = 0 and the lone point to
  (0, 1).  Then Z(u, v) = A(u) + e^v with A the partition function of
  the line points.  Put a = log A, q = v - a(u), t = e^q/(1 + e^q).
  A direct computation gives psi_uv + a' psi_vv = 0, so in (u, q)
  coordinates the metric is DIAGONAL:
      ds^2 = (1 - t) V(u) du^2 + t(1 - t) dq^2,   V = a'' .
  Setting d(rho) = sqrt(V) du, t = sin^2(theta), y = 2 theta turns
  this into the warped form ds^2 = dy^2 + cos^2(y/2) d(rho)^2, and for
  ds^2 = dy^2 + f^2 d(rho)^2 one has K = -f''/f = 1/4.  []

  The hypothesis is SHARP in both directions, and A4 gates both: with
  ZERO off-line points the two statistics are affinely dependent,
  det g vanishes identically and there is no 2D metric at all; with
  TWO off-line points K is not 1/4 and not even constant.

  This is the mechanism behind the v0.14.0 lock K_four_atom = 1/4,
  because w_4 = (1/2) log 4 gives X_4 = 2 X_2 EXACTLY, so X_1, X_2,
  X_4 are collinear and X_3 is the lone off-line point (A5).  The
  theorem is tested OUT OF SAMPLE on {1,2,4,8} and {1,2,4,16}, which
  are in the rigidity class for the same reason and were never used
  to motivate it.

THE COLD LAW, CORRECTED (A7).  The cold curvature is a competition
between triangle products.  det g is dominated by the triangle (1,2,3)
with Boltzmann product 6, so det(g)^2 ~ Delta_123^4 36^{-beta}.  In N
the first candidate without the ground atom would be (2,3,4) with
product 24 — but its factor C_24 = cross(X_2, X_4) VANISHES by the
same collinearity, so the (3/2)^beta mode is killed identically.  The
first survivor is (2,3,5) with product 30, giving

    K(beta) ~ -c (36/30)^beta = -c (6/5)^beta ,

so 6/5 is not produced by atom 6 at all: the cold core is {1,2,3,5}
and the rate is the ratio of the squared metric triangle to the
leading curvature triangle.  The amplitude has a closed form in
log 2, log 3, log 5:

    c = Delta_235 C_23 C_25 C_35 / (4 C_23^4)
      = (9 L5^2)/(25 L2 L3) [1 - (9/5) L5/L3 + (4/5) L5/L2]
      = 0.27012646530542495706433719670365...

CORRECTION CC-02 (append-only).  v0.14.0 ledgers
curvature.CURVATURE["cold_amplitude"] = 0.270126465305424759517602.
That value is WRONG in the 16th digit.  It was extracted from the
SIX-atom model at beta = 200, and the six-atom model still carries an
atom-6 contamination of relative size (5/6)^beta; at beta = 200 that
is 1.46e-16, and the observed error is 1.98e-16.  A7 reproduces the
contaminated value from the six-atom model to the digit, confirming
the diagnosis, and replaces the constant by the closed form.  The old
value is preserved here as the discovery route.

CORRECTION CC-03 (API safety).  finite_atom_curvature accepts a dps
override that can LOWER precision below the safe floor.  The Brioschi
numerator cancels about 0.653*beta digits, so at beta = 80 with
dps = 90 the function returns ~1e+24 instead of 1/4 — silently.  A8
gates this failure and the accompanying patch floors the override.

Gates A1-A9.  Runtime ~3 min.  Writes
arithmetic_area_geometry_ledger.json next to itself.
"""

from __future__ import annotations

import json
import math
import os
import sys
import time
from itertools import combinations

import mpmath as mp

try:
    from mtft.curvature import brioschi, finite_atom_curvature
    from mtft.moments import weight_susceptibility
except ImportError:                                   # source tree
    for _m in [k for k in sys.modules if k.split(".")[0] == "mtft"]:
        del sys.modules[_m]
    for _p in ("/home/claude/v14/mtft-0.14.0/src",
               os.environ.get("MTFT_SRC", "")):
        if _p and _p not in sys.path:
            sys.path.insert(0, _p)
    from mtft.curvature import brioschi, finite_atom_curvature
    from mtft.moments import weight_susceptibility

HERE = os.path.dirname(os.path.abspath(__file__))
mp.mp.dps = 120


def weight(n):
    return sum(mp.log(d) / d for d in range(2, n + 1) if n % d == 0)


def point(n):
    return (mp.log(n), weight(n)) if n > 1 else (mp.mpf(0), mp.mpf(0))


def cross(a, b):
    return a[0] * b[1] - a[1] * b[0]


def sub(a, b):
    return (a[0] - b[0], a[1] - b[1])


def geometry(beta, atoms, custom_w=None):
    """Both routes to det g and K on a finite support."""
    b = mp.mpf(beta)
    X = []
    for n in atoms:
        wn = (mp.mpf(custom_w[n]) if custom_w and n in custom_w
              else (weight(n) if n > 1 else mp.mpf(0)))
        X.append((mp.log(n) if n > 1 else mp.mpf(0), wn))
    q = [mp.e ** (-b * x[0]) for x in X]
    Z = sum(q)
    p = [x / Z for x in q]
    m = len(X)
    mu = (sum(p[i] * X[i][0] for i in range(m)),
          sum(p[i] * X[i][1] for i in range(m)))
    C = [sub(x, mu) for x in X]
    E = sum(p[i] * C[i][0] ** 2 for i in range(m))
    F = sum(p[i] * C[i][0] * C[i][1] for i in range(m))
    G = sum(p[i] * C[i][1] ** 2 for i in range(m))
    detg = E * G - F * F
    Dtri = mp.mpf(0)
    Ngeo = mp.mpf(0)
    for i, j, k in combinations(range(m), 3):
        D = cross(sub(X[j], X[i]), sub(X[k], X[i]))
        wt = p[i] * p[j] * p[k]
        Dtri += wt * D ** 2
        Ngeo += wt * D * cross(C[i], C[j]) * cross(C[i], C[k]) \
            * cross(C[j], C[k])
    A = sum(p[i] * C[i][0] ** 3 for i in range(m))
    B = sum(p[i] * C[i][0] ** 2 * C[i][1] for i in range(m))
    Cc = sum(p[i] * C[i][0] * C[i][1] ** 2 for i in range(m))
    Dd = sum(p[i] * C[i][1] ** 3 for i in range(m))
    Ndir = mp.det(mp.matrix([[E, F, G], [A, B, Cc], [B, Cc, Dd]]))
    return dict(detg=detg, Dtri=Dtri, Ndir=Ndir, Ngeo=Ngeo,
                Kdir=-Ndir / (4 * detg ** 2),
                Kgeo=-Ngeo / (4 * Dtri ** 2))


def K_support(points, base, u, v):
    """Brioschi curvature for an arbitrary support and base weights."""
    def psi(a, b):
        return mp.log(sum(mp.mpf(base[i])
                          * mp.e ** (a * points[i][0] + b * points[i][1])
                          for i in range(len(points))))
    d = lambda i, j: mp.diff(psi, (mp.mpf(u), mp.mpf(v)), (i, j))
    return brioschi(d(2, 0), d(1, 1), d(0, 2), d(3, 0), d(2, 1),
                    d(2, 1), d(1, 2), d(1, 2), d(0, 3))


def area_zeta(s):
    """A(s) = sum Delta^2/(abc)^s, via zeta(s)^3 det g(s, 0)."""
    s = mp.mpf(s)
    E = mp.diff(lambda x: mp.log(mp.zeta(x)), s, 2)
    G = weight_susceptibility(s)
    F = mp.diff(mp.zeta, s + 1, 2)
    return mp.zeta(s) ** 3 * (E * G - F ** 2)


def cold_amplitude_closed():
    """The exact cold amplitude c, in log 2, log 3, log 5."""
    L2, L3, L5 = mp.log(2), mp.log(3), mp.log(5)
    return (9 * L5 ** 2) / (25 * L2 * L3) * (
        1 - mp.mpf(9) / 5 * L5 / L3 + mp.mpf(4) / 5 * L5 / L2)


def on_line(n, tol=mp.mpf("1e-60")):
    """Is X_n on the line y = x/2, i.e. X_n = k X_2 ?"""
    return abs(weight(n) - mp.log(n) / 2) < tol if n > 1 else True


def main() -> int:
    t0 = time.time()
    ledger = {"study": "arithmetic_area_geometry", "gates": {}}
    ok = True

    def gate(name, passed, **info):
        nonlocal ok
        ok &= bool(passed)
        ledger["gates"][name] = {"passed": bool(passed), **info}
        print(f"[{'PASS' if passed else 'FAIL'}] {name}  "
              + "  ".join(f"{k}={v}" for k, v in info.items()))

    # ── A1/A2: the two identities, against the Brioschi route ───────
    eD = eN = eK = eB = mp.mpf(0)
    for atoms in ((1, 2, 3, 4, 5, 6), (1, 2, 3, 5, 7, 10),
                  (1, 2, 3, 4, 6, 8, 9)):
        for beta in ("2.7", "5.3", "11"):
            g = geometry(beta, atoms)
            eD = max(eD, abs(g["detg"] - g["Dtri"]))
            eN = max(eN, abs(g["Ndir"] - g["Ngeo"]))
            eK = max(eK, abs(g["Kdir"] - g["Kgeo"]))
            pts = [point(n) for n in atoms]
            base = [mp.mpf(1)] * len(atoms)
            KB = K_support(pts, base, -mp.mpf(beta), 0)
            eB = max(eB, abs(g["Kgeo"] - KB))
    gate("A1_triangle_determinant_identity", eD < mp.mpf("1e-70"),
         max_error=mp.nstr(eD, 5),
         identity="det g = sum p_i p_j p_k Delta_ijk^2")
    gate("A2_chirality_identity_vs_brioschi",
         eN < mp.mpf("1e-70") and eK < mp.mpf("1e-60")
         and eB < mp.mpf("1e-55"),
         numerator_error=mp.nstr(eN, 5),
         vs_cumulant_determinant=mp.nstr(eK, 5),
         vs_independent_brioschi_route=mp.nstr(eB, 5))

    # ── A3: the rigidity theorem, including OUT OF SAMPLE ───────────
    line_pts = [(mp.mpf("-2.2"), 0), (0, 0), (mp.mpf("1.3"), 0),
                (mp.mpf("4.7"), 0), (mp.mpf("8.1"), 0), (0, 1)]
    base = [mp.mpf(x) for x in ("1.7", "2.1", "0.3", "5.0", "0.8",
                                "4.2")]
    synth = [K_support(line_pts, base, u, v)
             for u, v in ((mp.mpf(".2"), mp.mpf("-.3")),
                          (mp.mpf(1), mp.mpf(2)),
                          (mp.mpf(-2), mp.mpf(".5")))]
    oos = {}
    for atoms in ((1, 2, 4, 8), (1, 2, 4, 16)):
        oos[str(atoms)] = max(
            abs(finite_atom_curvature(b, atoms, dps=200) - mp.mpf(1) / 4)
            for b in (6, 17, 41))
    err = max([abs(x - mp.mpf(1) / 4) for x in synth]
              + list(oos.values()))
    gate("A3_rigidity_theorem", err < mp.mpf("1e-50"),
         max_error=mp.nstr(err, 5),
         synthetic="5 line points + 1 off-line, arbitrary base weights",
         out_of_sample={k: mp.nstr(v, 4) for k, v in oos.items()},
         theorem="K == 1/4")

    # ── A4: the hypothesis is sharp in both directions ──────────────
    degenerate = False
    try:
        finite_atom_curvature(6, (1, 2, 4), dps=120)
    except ZeroDivisionError:
        degenerate = True
    two_off = [max(abs(finite_atom_curvature(b, a, dps=200)
                       - mp.mpf(1) / 4) for b in (6, 17, 41))
               for a in ((1, 2, 4, 8, 16), (1, 2, 3, 4, 5))]
    gate("A4_hypothesis_is_sharp",
         degenerate and all(x > 1 for x in two_off),
         zero_off_line="det g = 0, no 2D metric",
         two_off_line=[mp.nstr(x, 4) for x in two_off])

    # ── A5: why {1,2,3,4} is in the class ───────────────────────────
    X2, X4 = point(2), point(4)
    exact = (abs(X4[0] - 2 * X2[0]) < mp.mpf("1e-100")
             and abs(X4[1] - 2 * X2[1]) < mp.mpf("1e-100"))
    gate("A5_dyadic_collinearity", exact
         and on_line(2) and on_line(4) and not on_line(3),
         identity="w_4 = (1/2) log 4  =>  X_4 = 2 X_2 exactly",
         consequence="{1,2,4} collinear, 3 is the lone off-line point")

    # ── A6: prime boundary, with the n = 1 correction ───────────────
    bnd_ok = True
    for n in range(1, 500):
        eq = abs(weight(n) - (mp.log(n) / n if n > 1 else mp.mpf(0))) \
            < mp.mpf("1e-70")
        isp = n > 1 and all(n % d for d in range(2, int(n ** .5) + 1))
        bnd_ok &= (eq == (isp or n == 1))
        bnd_ok &= weight(n) >= (mp.log(n) / n if n > 1 else mp.mpf(0)) \
            - mp.mpf("1e-70")
    gate("A6_prime_boundary_corrected", bnd_ok,
         statement="w_n >= (log n)/n, equality iff n = 1 OR n prime",
         correction="the n = 1 case was omitted upstream; X_1 = (0,0) "
                    "lies ON the curve and 1 is not prime",
         checked="1 <= n < 500")

    # ── A7: cold law, closed amplitude, and correction CC-02 ────────
    c_exact = cold_amplitude_closed()
    X3, X5 = point(3), point(5)
    C23, C25, C35 = cross(X2, X3), cross(X2, X5), cross(X3, X5)
    D235 = cross(sub(X3, X2), sub(X5, X2))
    c_geom = D235 * C23 * C25 * C35 / (4 * C23 ** 4)
    k150 = finite_atom_curvature(150, (1, 2, 3, 5), dps=250)
    k200 = finite_atom_curvature(200, (1, 2, 3, 5), dps=300)
    rate = (k200 / k150) ** (mp.mpf(1) / 50)
    c_core = -k200 / (mp.mpf(6) / 5) ** 200
    c_repo = mp.mpf("0.270126465305424759517602")
    k6 = finite_atom_curvature(200, (1, 2, 3, 4, 5, 6), dps=300)
    c_six = -k6 / (mp.mpf(6) / 5) ** 200
    reproduced = abs(c_six - c_repo) < mp.mpf("1e-22")
    contamination = abs(c_six - c_exact)
    predicted = (mp.mpf(5) / 6) ** 200
    gate("A7_cold_core_and_CC02",
         abs(c_geom - c_exact) < mp.mpf("1e-90")
         and abs(rate - mp.mpf(6) / 5) < mp.mpf("1e-12")
         and abs(c_core - c_exact) < mp.mpf("1e-16")
         and reproduced
         and contamination / predicted < 3,
         c_closed=mp.nstr(c_exact, 32),
         rate_from_1235=mp.nstr(rate, 16),
         repo_value_reproduced_by_six_atom=reproduced,
         contamination=mp.nstr(contamination, 4),
         predicted_scale_five_sixths_200=mp.nstr(predicted, 4),
         verdict="CC-02: the ledgered amplitude is contaminated by "
                 "atom 6; the closed form supersedes it")
    ledger["cold_amplitude_corrected"] = mp.nstr(c_exact, 32)
    ledger["cold_amplitude_retracted"] = str(c_repo)

    # ── A8: the (3/2) mode and the dps foot-gun (CC-03) ─────────────
    delta = mp.mpf("1e-4")
    cw = {1: mp.mpf(0), 2: mp.log(2) / 2, 3: mp.log(3) / 3,
          4: mp.log(2) + delta}
    cwm = dict(cw)
    cwm[4] = mp.log(2) - delta
    k40 = geometry(40, (1, 2, 3, 4), custom_w=cw)["Kdir"]
    k50 = geometry(50, (1, 2, 3, 4), custom_w=cw)["Kdir"]
    k40m = geometry(40, (1, 2, 3, 4), custom_w=cwm)["Kdir"]
    rate32 = (k50 / k40) ** (mp.mpf(1) / 10)
    amp = k50 / (mp.mpf(3) / 2) ** 50
    amp_lin = 3 * delta / mp.log(3)
    shield = (abs(rate32 - mp.mpf(3) / 2) < mp.mpf("2e-4")
              and k50 > 0 and k40m < 0
              and abs(amp / amp_lin - 1) < mp.mpf("0.01"))
    bad = finite_atom_curvature(80, (1, 2, 4, 8), dps=90)
    good = finite_atom_curvature(80, (1, 2, 4, 8), dps=200)
    footgun = (abs(bad - mp.mpf(1) / 4) > 1
               and abs(good - mp.mpf(1) / 4) < mp.mpf("1e-70"))
    gate("A8_dyadic_shield_and_CC03", shield and footgun,
         rate=mp.nstr(rate32, 14), amplitude=mp.nstr(amp, 12),
         linear_prediction=mp.nstr(amp_lin, 12),
         dps90_at_beta80=mp.nstr(bad, 6), dps200_at_beta80="1/4",
         verdict="CC-03: a dps override below the safe floor fails "
                 "silently; the patch floors it")

    # ── A9: the arithmetic area zeta ────────────────────────────────
    chi1 = weight_susceptibility(1)
    wall = [(e, mp.mpf(e) ** 5 * area_zeta(1 + mp.mpf(e)))
            for e in ("1e-2", "1e-3", "1e-4")]
    wall_ok = abs(wall[-1][1] / chi1 - 1) < mp.mpf("2e-4")
    cold = area_zeta(40) / (C23 ** 2 * mp.mpf(6) ** (-40))
    gate("A9_arithmetic_area_zeta",
         wall_ok and abs(cold - 1) < mp.mpf("2e-8"),
         identity="A(s) = zeta(s)^3 det g(s, 0)",
         wall_order=5, wall_coefficient=mp.nstr(chi1, 20),
         ladder=[mp.nstr(v, 12) for _, v in wall],
         cold_ratio_at_40=mp.nstr(cold, 12),
         cold_leading="Delta_123^2 6^-s")

    ledger["all_passed"] = ok
    ledger["runtime_s"] = round(time.time() - t0, 1)
    with open(os.path.join(HERE,
              "arithmetic_area_geometry_ledger.json"), "w") as f:
        json.dump(ledger, f, indent=2, default=str)
    print(f"\nledger written  [{ledger['runtime_s']}s]")
    print("ALL GATES PASS" if ok else "GATE FAILURE")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
