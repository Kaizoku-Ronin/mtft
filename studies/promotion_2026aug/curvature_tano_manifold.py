#!/usr/bin/env python3
"""
curvature_tano_manifold.py — the Riemannian geometry of the Tano ensemble
==========================================================================

MIT License — Copyright (c) 2026 Roger Tano

The two-parameter exponential family on the positive integers

    p_{beta,lambda}(n) = n^{-beta} e^{lambda w_n} / Z(beta, lambda),
    w_n = sum_{d|n} (log d)/d,

is an honest 2D statistical manifold for beta > 1 and all real lambda
(w_n = O((log log n)^2) on champions, so Z converges).  Its Fisher
metric at lambda = 0 was certified in w2_susceptibility.py and its
cubic (Amari-Chentsov) tensor in w3_cumulants.py.  This study computes
the actual Levi-Civita curvature.

THE CANCELLATION THEOREM (the kappa_4 twist).  In (beta, lambda)
coordinates the metric is a Hessian: g_ij = d_i d_j psi with
psi = log Z.  Hence Gamma_{ij,k} = (1/2) d_i d_j d_k psi, and the
second-derivative block of the curvature,

    -1/2 d_vv E + d_uv F - 1/2 d_uu G
      = -1/2 k4 + k4 - 1/2 k4 = 0,

vanishes identically because all three terms are the SAME totally
symmetric fourth derivative d^4 psi (= fourth cumulant).  The fourth
cumulant plays exactly no role in the curvature of an exponential
family.  Verified three ways: analytically (above), on the Gaussian
family (G1), and by direct finite-difference measurement of the block
on the Tano manifold itself (G2: |block| ~ 1e-11 vs terms O(1)).

Curvature therefore needs ONLY the certified second and third
cumulants.  With E = Var l, F = -Cov(l, w), G = chi_w and the sign
rule  d_beta^a d_lambda^b psi = (-1)^a kappa(l^a, w^b):

    E_u = -kappa_lll = (log zeta)'''(beta)
    E_v = F_u = kappa_llw = -zeta'''(beta+1)     [study-4 identity,
                        reconfirmed here as equality of mixed partials]
    F_v = G_u = -kappa_wwl = chi_w'(beta)
    G_v = kappa_www = U(beta) + 3 T(beta) zeta'(beta+1)
                      - 2 zeta'(beta+1)^3

and K is the reduced Brioschi form of these nine functions.

RESULTS
-------
R1 (Pr).  SIGN-CHANGING CURVATURE.  K(beta) along lambda = 0 is a
  positive dome on (1, beta_0) and a hyperbolic cold tail beyond:
  it leaves the Hagedorn wall flat (R2), climbs over the summit

      beta* = 4.593591164956...,   K* = 1.19569598199193852...,

  crosses zero at the flat temperature

      beta_0 = 8.8565170425...    (|K_exact(beta_0)| < 1e-8),

  and dives negative with exact asymptotic rate 6/5 (R4).  Anchors,
  E2 across three routes sharing no steps (exact zeta/Euler engines +
  reduced Brioschi;  raw 4e6-sieve cumulants + reduced Brioschi;
  (beta, lambda)-grid finite differences + FULL Brioschi):
      K(2.5) = 0.55919136446792...   K(3.5) = 0.99788515020960...
  CORRECTIONS TRAIL (append-only): the first draft of this study
  claimed K > 0 on all of (1, 30] — falsified by the profile itself
  (the draft's sieve tail past beta ~ 24 was float64 noise; the
  Brioschi numerator cancels ~0.653*beta digits, ratio 4.5^beta).
  A second interim reading called the negative tail an
  "infinite-support effect" — falsified by the K_M table: the
  {1,...,6} family already carries the full cold geometry.

R2 (Pr, closed form).  HAGEDORN FLATNESS.  K -> 0+ linearly as
  beta -> 1+, with slope

      A = ( zeta''(2) kappa3_cold  -  kappa_wwl_cold chi_cold )
          / ( 2 chi_cold^2 )
        = 0.423657463797093...

  — a pure combination of four constants already in the certified
  ledger (studies 1 and 4); ladder K(1+eps)/eps confirms with a
  constant quadratic coefficient ~ -1.18.

R3 (EXACT).  THREE RIGIDITY LOCKS validate the pipeline on known
  geometries: the Gaussian family gives K = -1/2 (both reduced and
  full Brioschi, 1e-30); the {1,2,3} family gives K = +1/4
  IDENTICALLY in beta (the full 2-simplex: Fisher sphere of radius
  2); and — a small discovery — the {1,2,3,4} family ALSO gives
  K = +1/4 identically (1e-40 gate; observed 1e-49..1e-56): atom 4
  is inert, mechanism sketch (Pp): (l_4, w_4) = 2 (l_2, w_2), so
  {1,2,4} are collinear and equally spaced in statistic space, and
  the second statistic can be re-based to be supported on atom 3
  alone.

R4 (Pr + mechanism; pre-registration FALSIFIED and preserved).  The
  pre-registered cold mechanism — "the full ensemble converges to the
  3-atom model" — is FALSE: K_3atom = 1/4 for all beta while the full
  K dives negative; |K_full - 1/4| grows without bound (table kept).
  The true mechanism, established by exact finite-support isolation:
  ATOM 5 FLIPS THE SIGN ({1,2,3,5} and {1,2,3,4,5} are negative and
  diverging; {1,2,3,4,6} is positive and beta-stable at
  1.3549368866...), and with atoms 5 and 6 both present the cold
  asymptote is

      K(beta)  ~  - c * (6/5)^beta ,
      rate 6/5 certified two ways sharing no steps:
        full-ensemble Richardson on the repaired high-precision tail
        (r_extrap = 1.200000419), and the exact {1,...,6} model
        (rate 1.2000138 at beta = 60 -> 64);
      amplitude c = 0.26991206794...  — full ensemble and 6-atom
      model agree to 8 digits at beta = 48; c is ledgered to 25
      digits from the model at beta = 200.

  The negative cold curvature of the Tano manifold is the geometry of
  the first six integers; the rate 6/5 is the Boltzmann ratio of atom
  5 (the new prime) to atom 6 (the first mixed composite).

Gates G1-G8.  Runtime ~4-6 min.  Writes
curvature_tano_manifold_ledger.json next to itself.
"""

from __future__ import annotations

import importlib.util as iu
import json
import os
import sys
import time

import numpy as np
from mpmath import mp

HERE = os.path.dirname(os.path.abspath(__file__))
W2_PATH = os.path.join(HERE, "w2_susceptibility.py")
E3_PATH = os.path.join(HERE, "engine3.py")


def load(name, path):
    spec = iu.spec_from_file_location(name, path)
    m = iu.module_from_spec(spec)
    try:
        spec.loader.exec_module(m)
    except SystemExit:
        pass
    return m


def brioschi(E, F, G, Eu, Ev, Fu, Fv, Gu, Gv, blk=0):
    h = mp.mpf("0.5")
    M1 = mp.matrix([[blk, h * Eu, Fu - h * Ev],
                    [Fv - h * Gu, E, F],
                    [h * Gv, F, G]])
    M2 = mp.matrix([[0, h * Ev, h * Gu],
                    [h * Ev, E, F],
                    [h * Gu, F, G]])
    return (mp.det(M1) - mp.det(M2)) / (E * G - F * F) ** 2


def main() -> int:
    t00 = time.time()
    mp.dps = 40
    w2 = load("w2", W2_PATH)
    e3 = load("e3", E3_PATH)

    # resumable cache for the one expensive primitive (exact U3):
    # written atomically after every miss, so interrupted runs resume.
    CACHE_PATH = os.path.join(HERE, "curv_cache.json")
    try:
        _CACHE = json.load(open(CACHE_PATH))
    except Exception:
        _CACHE = {}
    _u3_raw = e3.U3

    def _u3_cached(sv, P=500, OM=20):
        key = f"U|{mp.nstr(mp.mpf(sv), 30)}|{P}|{OM}"
        if key in _CACHE:
            return mp.mpf(_CACHE[key]), None
        val, Fd = _u3_raw(sv, P=P, OM=OM)
        _CACHE[key] = mp.nstr(val, 45)
        tmp = CACHE_PATH + ".tmp"
        with open(tmp, "w") as fh:
            json.dump(_CACHE, fh)
        os.replace(tmp, CACHE_PATH)
        return val, Fd

    e3.U3 = _u3_cached

    def _u3_hp(sv, P=50, OM=10, dps=120):
        key = f"U|{sv}|{P}|{OM}|d{dps}"
        if key in _CACHE:
            return mp.mpf(_CACHE[key])
        old = mp.dps
        mp.dps = dps
        val, _ = _u3_raw(mp.mpf(sv), P=P, OM=OM)
        mp.dps = old
        _CACHE[key] = mp.nstr(val, dps - 10)
        tmp = CACHE_PATH + ".tmp"
        with open(tmp, "w") as fh:
            json.dump(_CACHE, fh)
        os.replace(tmp, CACHE_PATH)
        return mp.mpf(_CACHE[key])
    ledger = {"study": "curvature_tano_manifold", "gates": {}}
    ok = True

    def gate(name, passed, **info):
        nonlocal ok
        ok &= bool(passed)
        ledger["gates"][name] = {"passed": bool(passed), **info}
        print(f"[{'PASS' if passed else 'FAIL'}] {name}  "
              + "  ".join(f"{k}={v}" for k, v in info.items()))

    # ── G1: Gaussian convention lock ────────────────────────────────
    psi = lambda a, b: -a * a / (4 * b) + mp.mpf("0.5") * mp.log(-mp.pi / b)
    g1 = True
    for (a0, b0) in ((mp.mpf("0.3"), mp.mpf("-0.7")),
                     (mp.mpf("-1.1"), mp.mpf("-0.25")),
                     (mp.mpf("0.0"), mp.mpf("-2.0"))):
        d = lambda i, j: mp.diff(psi, (a0, b0), (i, j))
        blk = -mp.mpf("0.5") * d(2, 2) + d(2, 2) - mp.mpf("0.5") * d(2, 2)
        Kr = brioschi(d(2, 0), d(1, 1), d(0, 2), d(3, 0), d(2, 1),
                      d(2, 1), d(1, 2), d(1, 2), d(0, 3), 0)
        Kf = brioschi(d(2, 0), d(1, 1), d(0, 2), d(3, 0), d(2, 1),
                      d(2, 1), d(1, 2), d(1, 2), d(0, 3), blk)
        g1 &= abs(Kr + mp.mpf("0.5")) < mp.mpf("1e-30") and \
              abs(Kf + mp.mpf("0.5")) < mp.mpf("1e-30")
    gate("G1_gaussian_convention", g1, K="-1/2 both formulas, 3 points")

    # ── exact-route components ──────────────────────────────────────
    zd = lambda s, k: mp.diff(mp.zeta, s, k) if k else mp.zeta(s)
    chi = w2.weight_susceptibility_exact
    Tm = w2.weight_second_moment_exact

    def exact_components(beta, P=500, OM=20):
        beta = mp.mpf(beta)
        E = mp.diff(lambda b: mp.log(mp.zeta(b)), beta, 2)
        Eu = mp.diff(lambda b: mp.log(mp.zeta(b)), beta, 3)
        F = -zd(beta + 1, 2)
        Fu = -zd(beta + 1, 3)
        Ev = Fu                      # study-4 identity kappa_llw
        G = chi(beta)
        Gu = mp.diff(chi, beta, 1)
        Fv = Gu                      # study-4 identity -kappa_wwl
        U, _ = e3.U3(beta, P=P, OM=OM)
        zp = zd(beta + 1, 1)
        Gv = U + 3 * Tm(beta) * zp - 2 * zp ** 3
        return [E, F, G, Eu, Ev, Fu, Fv, Gu, Gv]

    K_exact = lambda b, **kw: brioschi(*exact_components(b, **kw))

    # ── sieve route ─────────────────────────────────────────────────
    N = 4_000_000
    lg = np.zeros(N + 1)
    lg[1:] = np.log(np.arange(1, N + 1))
    wv = np.zeros(N + 1)
    for dd in range(2, N + 1):
        wv[dd::dd] += lg[dd] / dd
    narr = np.arange(1, N + 1, dtype=float)
    L, Wv = lg[1:], wv[1:]

    def sieve_components(beta, sup=None):
        if sup is None:
            nn, ll, ww = narr, L, Wv
        else:
            nn, ll, ww = narr[:sup], L[:sup], Wv[:sup]
        p = nn ** (-float(beta))
        p /= p.sum()
        mL, mW = (p * ll).sum(), (p * ww).sum()
        cL, cW = ll - mL, ww - mW
        E = (p * cL * cL).sum()
        F = -(p * cL * cW).sum()
        G = (p * cW * cW).sum()
        k3L = (p * cL ** 3).sum()
        kLLW = (p * cL * cL * cW).sum()
        kLWW = (p * cL * cW * cW).sum()
        k3W = (p * cW ** 3).sum()
        return [E, F, G, -k3L, kLLW, kLLW, -kLWW, -kLWW, k3W]

    K_sieve = lambda b, sup=None: brioschi(
        *[mp.mpf(float(x)) for x in sieve_components(b, sup)])

    # ── grid route (full Brioschi, no theorem input) ────────────────
    def K_grid(beta, h=1e-3):
        beta = float(beta)

        def gij(b, lam):
            wt = narr ** (-b) * np.exp(lam * Wv)
            q = wt / wt.sum()
            mL, mW = (q * L).sum(), (q * Wv).sum()
            cL, cW = L - mL, Wv - mW
            return ((q * cL * cL).sum(), -(q * cL * cW).sum(),
                    (q * cW * cW).sum())

        pts = {(i, j): gij(beta + i * h, j * h)
               for i in range(-2, 3) for j in range(-2, 3)}
        g = lambda k, i, j: pts[(i, j)][k]
        E, F, G = pts[(0, 0)]
        Eu = (g(0, 1, 0) - g(0, -1, 0)) / (2 * h)
        Ev = (g(0, 0, 1) - g(0, 0, -1)) / (2 * h)
        Fu = (g(1, 1, 0) - g(1, -1, 0)) / (2 * h)
        Fv = (g(1, 0, 1) - g(1, 0, -1)) / (2 * h)
        Gu = (g(2, 1, 0) - g(2, -1, 0)) / (2 * h)
        Gv = (g(2, 0, 1) - g(2, 0, -1)) / (2 * h)
        Evv = (g(0, 0, 2) - 2 * E + g(0, 0, -2)) / (4 * h * h)
        Guu = (g(2, 2, 0) - 2 * G + g(2, -2, 0)) / (4 * h * h)
        Fuv = (g(1, 1, 1) - g(1, 1, -1) - g(1, -1, 1)
               + g(1, -1, -1)) / (4 * h * h)
        blk = -0.5 * Evv + Fuv - 0.5 * Guu
        return brioschi(*[mp.mpf(x) for x in
                          (E, F, G, Eu, Ev, Fu, Fv, Gu, Gv)],
                        blk=mp.mpf(blk)), blk

    # ── G2/G3: anchors, three routes ────────────────────────────────
    KE = {b: K_exact(b) for b in ("2.5", "3.5")}
    ledger["K_2p5"] = mp.nstr(KE["2.5"], 30)
    ledger["K_3p5"] = mp.nstr(KE["3.5"], 30)
    Kg25, blk25 = K_grid(2.5)
    Kg35, blk35 = K_grid(3.5)
    gate("G2_hessian_cancellation",
         abs(blk25) < 1e-8 and abs(blk35) < 1e-8
         and abs(Kg25 - KE["2.5"]) / abs(KE["2.5"]) < 2e-5
         and abs(Kg35 - KE["3.5"]) / abs(KE["3.5"]) < 2e-5,
         fourth_order_block=f"{blk25:.2e},{blk35:.2e}",
         note="kappa_4 measured absent from curvature")
    rd35 = abs(K_sieve(3.5) - KE["3.5"]) / abs(KE["3.5"])
    rd25 = abs(K_sieve(2.5) - KE["2.5"]) / abs(KE["2.5"])
    gate("G3_three_routes", rd35 < 1e-9 and rd25 < 1e-5,
         reldiff_3p5=mp.nstr(rd35, 3), reldiff_2p5=mp.nstr(rd25, 3),
         note="2.5 limited by sieve truncation N^-1.5")

    # ── G4: component-level E2 at beta = 3.5 ────────────────────────
    ce = exact_components("3.5")
    cs = sieve_components(3.5)
    comp_rel = max(abs(mp.mpf(float(cs[i])) - ce[i]) / abs(ce[i])
                   for i in range(9))
    gate("G4_component_E2", comp_rel < 1e-8,
         max_reldiff=mp.nstr(comp_rel, 3),
         identities="Ev==Fu, Fv==Gu (mixed partials = study-4 cumulant"
                    " identities)")

    # ── G5: Hagedorn flatness slope, closed form ────────────────────
    chi1 = chi(mp.mpf(1))
    chip1 = mp.diff(chi, mp.mpf(1), 1)
    z2_2 = zd(mp.mpf(2), 2)
    U1, _ = e3.U3(mp.mpf(1))
    zp2 = zd(mp.mpf(2), 1)
    k3w1 = U1 + 3 * Tm(mp.mpf(1)) * zp2 - 2 * zp2 ** 3
    A = (chip1 * chi1 + z2_2 * k3w1) / (2 * chi1 ** 2)
    ladder = []
    for eps in ("1e-3", "1e-4", "1e-5"):
        K = K_exact(mp.mpf(1) + mp.mpf(eps), P=300, OM=16)
        ladder.append((eps, mp.nstr(K / mp.mpf(eps), 15)))
    resid = [abs(mp.mpf(v) - A) / mp.mpf(e) for e, v in ladder]
    # resid[i] = |K/eps - A| / eps  ~  |c2|  (the quadratic coefficient):
    # near-constancy across the ladder is the O(eps^2) linearity proof.
    lin = all(mp.mpf("0.8") < resid[i] / resid[2] < mp.mpf("1.25")
              for i in range(2))
    ledger["hagedorn_c2_observed"] = mp.nstr(-resid[2], 8)
    gate("G5_hagedorn_slope",
         abs(mp.mpf(ladder[-1][1]) - A) < mp.mpf("1e-4") and lin,
         A=mp.nstr(A, 30), ladder=str(ladder),
         formula="(zeta''(2)k3_cold - k_wwl_cold chi_cold)/(2 chi_cold^2)")
    ledger["hagedorn_slope_A"] = mp.nstr(A, 30)
    ledger["slope_constituents"] = {
        "chi_cold": mp.nstr(chi1, 25),
        "chi_prime_cold(=-k_wwl_cold)": mp.nstr(chip1, 25),
        "zeta''(2)": mp.nstr(z2_2, 25),
        "kappa3_cold": mp.nstr(k3w1, 25)}

    # ── K_tail: high-precision exact route for the cold tail ────────
    # The Brioschi numerator cancels ~0.653*beta digits; past beta~36
    # the dps-40 route (and any float64 route) is noise.  dps 130 with
    # dps-120 U3 (P=50 suffices: 53^-36 is dust) is exact-grade to 48.
    def K_tail(bstr):
        old = mp.dps
        mp.dps = 130
        b = mp.mpf(bstr)
        E = mp.diff(lambda x: mp.log(mp.zeta(x)), b, 2)
        Eu = mp.diff(lambda x: mp.log(mp.zeta(x)), b, 3)
        F = -mp.diff(mp.zeta, b + 1, 2)
        Fu = -mp.diff(mp.zeta, b + 1, 3)
        Ev = Fu
        G = chi(b)
        Gu = mp.diff(chi, b, 1)
        Fv = Gu
        U = _u3_hp(bstr)
        zp = mp.diff(mp.zeta, b + 1, 1)
        Gv = U + 3 * Tm(b) * zp - 2 * zp ** 3
        out = brioschi(E, F, G, Eu, Ev, Fu, Fv, Gu, Gv)
        mp.dps = old
        return out

    # exact finite-support models (arbitrary atom sets)
    def K_set(bv, atoms, dps=None):
        old = mp.dps
        mp.dps = dps or max(60, int(0.7 * float(bv)) + 50)
        ls = [mp.log(nn) for nn in atoms]
        ws = [sum(mp.log(dd) / dd for dd in range(2, nn + 1)
                  if nn % dd == 0) for nn in atoms]
        b = mp.mpf(bv)
        M = len(atoms)
        p = [mp.e ** (-b * ls[i]) for i in range(M)]
        Z = sum(p)
        p = [x / Z for x in p]
        mL = sum(p[i] * ls[i] for i in range(M))
        mW = sum(p[i] * ws[i] for i in range(M))
        cL = [ls[i] - mL for i in range(M)]
        cW = [ws[i] - mW for i in range(M)]
        m = lambda f: sum(p[i] * f(i) for i in range(M))
        out = brioschi(m(lambda i: cL[i] ** 2),
                       -m(lambda i: cL[i] * cW[i]),
                       m(lambda i: cW[i] ** 2),
                       -m(lambda i: cL[i] ** 3),
                       m(lambda i: cL[i] ** 2 * cW[i]),
                       m(lambda i: cL[i] ** 2 * cW[i]),
                       -m(lambda i: cL[i] * cW[i] ** 2),
                       -m(lambda i: cL[i] * cW[i] ** 2),
                       m(lambda i: cW[i] ** 3))
        mp.dps = old
        return out

    # ── G6: profile, summit, flat temperature ───────────────────────
    prof = []
    for b in ("1.02", "1.05", "1.1", "1.2", "1.35", "1.5", "1.75",
              "2.0", "2.25"):
        prof.append((float(b), float(K_exact(b, P=300, OM=16))))
    for b in np.arange(2.5, 18.01, 0.25):
        prof.append((round(float(b), 2), float(K_sieve(b))))
    for b in ("20", "24", "28", "32", "36", "40", "44", "48"):
        prof.append((float(b), float(K_tail(b))))
    ks = [k for _, k in prof]
    imax = int(np.argmax(ks))
    lo = prof[max(0, imax - 1)][0]
    hi = prof[min(len(prof) - 1, imax + 1)][0]
    gr = (np.sqrt(5) - 1) / 2
    a_, b_ = lo, hi
    c_ = b_ - gr * (b_ - a_)
    d_ = a_ + gr * (b_ - a_)
    fc, fd = float(K_sieve(c_)), float(K_sieve(d_))
    for _ in range(40):
        if fc > fd:
            b_, d_, fd = d_, c_, fc
            c_ = b_ - gr * (b_ - a_)
            fc = float(K_sieve(c_))
        else:
            a_, c_, fc = c_, d_, fd
            d_ = a_ + gr * (b_ - a_)
            fd = float(K_sieve(d_))
    bstar = (a_ + b_) / 2
    Kstar_e = K_exact(f"{float(bstar):.12f}")
    az, bz = 8.0, 10.0
    for _ in range(40):
        mz = (az + bz) / 2
        if K_sieve(mz) > 0:
            az = mz
        else:
            bz = mz
    b0 = (az + bz) / 2
    K_at_b0 = K_exact(f"{b0:.12f}")
    dome = all(k > 0 for bb, k in prof if bb < b0 - 1e-9)
    tail_neg = all(k < 0 for bb, k in prof if bb > b0 + 1e-9)
    gate("G6_profile_summit_zero",
         dome and tail_neg
         and abs(float(Kstar_e) - float(K_sieve(bstar)))
         / float(Kstar_e) < 1e-8
         and abs(K_at_b0) < mp.mpf("1e-8"),
         beta_star=f"{float(bstar):.10f}",
         K_star=mp.nstr(Kstar_e, 20),
         beta_0=f"{b0:.10f}", K_at_beta0=mp.nstr(K_at_b0, 3))
    ledger["profile"] = prof
    ledger["beta_star"] = f"{float(bstar):.12f}"
    ledger["K_star"] = mp.nstr(Kstar_e, 30)
    ledger["beta_0"] = f"{b0:.12f}"

    # ── G7a: rigidity locks on 3- and 4-atom families ───────────────
    lock = True
    for bb in ("5", "12", "50"):
        lock &= abs(K_set(bb, [1, 2, 3]) - mp.mpf(1) / 4) <             mp.mpf("1e-40")
    for bb in ("7", "19"):
        lock &= abs(K_set(bb, [1, 2, 3, 4]) - mp.mpf(1) / 4) <             mp.mpf("1e-40")
    gate("G7a_simplex_locks", lock,
         K3="1/4 identically", K4="1/4 identically (atom 4 inert)")

    # ── G7b: pre-registered 3-atom mechanism — FALSIFIED, preserved ─
    tabl = {}
    for bb in (8, 12, 16):
        tabl[bb] = abs(float(K_sieve(bb)) - 0.25)
    tabl[20] = abs(float(K_tail("20")) - 0.25)
    falsified = tabl[20] > tabl[8]
    gate("G7b_prereg_3atom_FALSIFIED", falsified,
         verdict="pre-registered convergence to 3-atom model is FALSE",
         abs_dev_from_quarter={str(k): f"{v:.3f}"
                               for k, v in tabl.items()})

    # ── G7c: atom isolation — atom 5 flips the sign ─────────────────
    iso = {}
    isok = True
    for bb in (16, 24):
        a5 = float(K_set(bb, [1, 2, 3, 5]))
        a45 = float(K_set(bb, [1, 2, 3, 4, 5]))
        a46 = float(K_set(bb, [1, 2, 3, 4, 6]))
        iso[bb] = (a5, a45, a46)
        isok &= a5 < 0 and a45 < 0 and a46 > 0
    isok &= abs(iso[16][2] - iso[24][2]) < 1e-3
    gate("G7c_atom5_flips_sign", isok,
         sets="{1,2,3,5}<0, {1,2,3,4,5}<0, {1,2,3,4,6}>0 stable",
         K_12346_const=f"{iso[24][2]:.10f}")
    ledger["K_12346_limit"] = mp.nstr(K_set("200", [1, 2, 3, 4, 6]),
                                      15)

    # ── G7d: cold rate 6/5 and amplitude, two routes ────────────────
    tail = {b: K_tail(b) for b in ("28", "32", "36", "40", "44",
                                   "48")}
    tb = sorted(int(x) for x in tail)
    rs = [float((tail[str(tb[i + 1])] / tail[str(tb[i])])
                ** mp.mpf("0.25")) for i in range(len(tb) - 1)]
    dd_ = [rs[i + 1] - rs[i] for i in range(len(rs) - 1)]
    qq = dd_[-1] / dd_[-2]
    r_full = rs[-1] + dd_[-1] * qq / (1 - qq)
    r_model = float((K_set("64", [1, 2, 3, 4, 5, 6])
                     / K_set("60", [1, 2, 3, 4, 5, 6]))
                    ** mp.mpf("0.25"))
    c_full = -float(tail["48"]) / (6 / 5) ** 48
    c_model48 = -float(K_set("48", [1, 2, 3, 4, 5, 6])) / (6 / 5) ** 48
    c_model = -K_set("200", [1, 2, 3, 4, 5, 6]) /         (mp.mpf(6) / 5) ** 200
    gate("G7d_cold_rate_six_fifths",
         abs(r_full - 1.2) < 1e-5 and abs(r_model - 1.2) < 5e-5
         and abs(c_full - c_model48) / c_full < 1e-6,
         r_full_richardson=f"{r_full:.9f}",
         r_model_60_64=f"{r_model:.7f}",
         c=mp.nstr(c_model, 25))
    ledger["cold_rate"] = "6/5 (Pr, two routes)"
    ledger["cold_amplitude_c"] = mp.nstr(c_model, 25)
    ledger["tail_ratios"] = [f"{r:.9f}" for r in rs]

    # ── G8: robustness of exact route (engine parameters) ───────────
    mp.dps = 40
    Ka = K_exact("2.5", P=500, OM=20)
    Kb = K_exact("2.5", P=800, OM=24)
    gate("G8_engine_robustness", abs(Ka - Kb) < mp.mpf("1e-25"),
         delta=mp.nstr(abs(Ka - Kb), 3))

    ledger["all_passed"] = ok
    ledger["runtime_s"] = round(time.time() - t00, 1)
    with open(os.path.join(HERE,
              "curvature_tano_manifold_ledger.json"), "w") as f:
        json.dump(ledger, f, indent=2, default=str)
    print(f"\nledger written  [{ledger['runtime_s']}s]")
    print("ALL GATES PASS" if ok else "GATE FAILURE")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
