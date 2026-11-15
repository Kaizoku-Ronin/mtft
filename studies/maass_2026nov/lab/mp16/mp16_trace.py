#!/usr/bin/env python3
"""MP16: the Selberg trace formula on X0(N), twisted by the Atkin-Lehner involutions W_Q and the mirror R: z -> -conj(z),
for odd Maass forms (which have no continuous spectrum).  Plan: ../MP16_H24_Mechanism_Plan.md (saved 16:26:25, 12 Nov 2026).

For squarefree N with k prime factors, G = <W_p (p | N), R> / Gamma0(N) = (Z/2)^(k+1); the sectors are its characters.
For Q || N the odd sectors satisfy (projecting with characters; identity, cusp and Eisenstein terms cancel)

    sum_{chi odd} eps_Q(chi) S_chi(h) = (1/2) [ T_{W_Q}(h) - T_{R W_Q}(h) ],      S_chi(h) = sum_{r_j in chi} h(r_j),

T_g(h) = sum over Gamma0(N)-classes in the coset g Gamma0(N) of their orbital integrals.  The classes are counted with
binary quadratic forms: an SL(2,Z)-class of integer matrices with det +-Q and trace Qs is a class of forms of
discriminant Delta = Q^2 s^2 -+ 4Q; its number of Gamma0(N)-classes (weighted by centralizer length) is
  (local factor 1 at p | Q) x (number of roots of the form mod p, at p | N/Q) x 2 log eps+(order),
and h+(D c^2) log eps+(D c^2) = sqrt(D c^2) L(1, chi_D) prod_{p | c} (1 - chi_D(p)/p).  Summed over classes:
  hyperbolic (det +Q, 2 cosh(l/2) = s sqrt Q):  2 sqrt(Q) g(l_s) W(Delta_s)
  glide      (det -Q, 2 sinh(l/2) = s sqrt Q):  2 sqrt(Q) g(l_s) W(Delta_s)
  W(Delta) = sum_{f | c0} f^-1 L(1, chi_D) prod_{p | c0/f} (1 - chi_D(p)/p) prod_{p | N/Q} rho_p,   Delta = D c0^2,
  rho_p = p + 1 if p | f, else 1 + (D (c0/f)^2 / p);
  elliptic (det +Q, trace 0): E_Q fixed points of W_Q, each contributing c_ell * int h(r) sech(pi r) dr;
  reflections (det -Q, trace 0): mirrors of total length L_Q = 2 sqrt(Q) W(4Q), contributing c_refl * L_Q * g(0).
The constants c_ell, c_refl and the hyperbolic and glide normalizations are measured by the V-checks below.

Usage:  python3 -I mp16_trace.py vchecks | arith | level143 | cross
Needs numpy, scipy, mpmath and cypari2 (PARI's quadclassunit gives h and the regulator; GRH-conditional in general,
checked against an independent L(1, chi) series on a sample)."""
import json
import math
import os
import sys
import time

import cypari2
import mpmath as mp
import numpy as np
from scipy.integrate import quad
from scipy.special import erfc, exp1

HERE = os.path.dirname(os.path.abspath(__file__))
pari = cypari2.Pari()
pari.set_real_precision(38)
T0 = time.time()


def dump(obj, name):
    with open(os.path.join(HERE, name), "w") as fh:
        json.dump(obj, fh, indent=1)


# ----------------------------------------------------------------------------------------------- V-checks (normalizations)
def vchecks():
    mp.mp.dps = 25
    out = {}
    for a in (1, 3):
        k = lambda u, a=a: mp.exp(-a * u)
        Qv = lambda v, a=a: mp.sqrt(mp.pi / a) * mp.exp(-a * v)
        g = lambda r: 2 * Qv(mp.sinh(r / 2) ** 2)                       # Iwaniec (1.62): g(r) = 2 Q(sinh^2(r/2))
        h = lambda t: 2 * mp.quad(lambda r: mp.cos(r * t) * g(r), [0, 1, 2, 4, 7])   # h(t) = int e^{irt} g(r) dr
        u_of = lambda z, w: abs(z - w) ** 2 / (4 * mp.im(z) * mp.im(w))
        rec = {}
        # V1: the Selberg transform, int_H k(u(i,w)) Im(w)^{1/2+it} dmu(w) = h(t)  (geodesic polar coordinates about i)
        for t in (0.37, 1.3):
            nu = -mp.mpf(1) / 2 + 1j * t
            val = 2 * mp.pi * mp.quad(lambda rho: k(mp.sinh(rho / 2) ** 2) * mp.sinh(rho) * mp.legenp(nu, 0, mp.cosh(rho), type=3), [0, 1, 3, 6, 9])
            rho0 = mp.mpf("1.1")                                           # spot check of the Legendre function against its angular integral
            ang = mp.quad(lambda al: (mp.cosh(rho0) - mp.sinh(rho0) * mp.cos(al)) ** (-nu - 1), [0, mp.pi, 2 * mp.pi]) / (2 * mp.pi)
            rec[f"V1 t={t}"] = {"integral": [float(mp.re(val)), float(mp.im(val))], "h": float(h(t)),
                                "legendre_spot": float(abs(ang - mp.legenp(nu, 0, mp.cosh(rho0), type=3)))}
        # V2-V4: strip integrals (centralizer length 1) for z -> e^l z, z -> -e^l conj z, z -> -conj z
        for l in (0.9, 2.3):
            el = mp.e ** l
            hyp = mp.quad(lambda ph: k(u_of(mp.expj(ph), el * mp.expj(ph))) / mp.sin(ph) ** 2, [0, mp.pi / 2, mp.pi])
            gli = mp.quad(lambda ph: k(u_of(mp.expj(ph), -el * mp.conj(mp.expj(ph)))) / mp.sin(ph) ** 2, [0, mp.pi / 2, mp.pi])
            rec[f"V2 hyperbolic l={l}"] = {"integral": float(hyp), "g/(2 sinh)": float(g(l) / (2 * mp.sinh(l / 2)))}
            rec[f"V3 glide l={l}"] = {"integral": float(gli), "g/(2 cosh)": float(g(l) / (2 * mp.cosh(l / 2)))}
        ref = mp.quad(lambda ph: k(u_of(mp.expj(ph), -mp.conj(mp.expj(ph)))) / mp.sin(ph) ** 2, [0, mp.pi / 2, mp.pi])
        rec["V4 reflection"] = {"integral": float(ref), "g(0)/2": float(g(0) / 2)}
        # V5: rotation by pi about i (z -> -1/z); the point at distance rho from i, direction al
        zpt = lambda rho, al: 1j * (1 + mp.tanh(rho / 2) * mp.expj(al)) / (1 - mp.tanh(rho / 2) * mp.expj(al))
        spread = max(abs(u_of(zpt(mp.mpf("0.8"), al), -1 / zpt(mp.mpf("0.8"), al)) - mp.sinh(mp.mpf("0.8")) ** 2) for al in (0.3, 1.7, 2.9, 4.4))
        ell = 2 * mp.pi * mp.quad(lambda rho: k(u_of(zpt(rho, 0.4), -1 / zpt(rho, 0.4))) * mp.sinh(rho), [0, 1, 2, 4, 6])
        Isech = 2 * mp.quad(lambda t: h(t) / mp.cosh(mp.pi * t), [0, 1, 2, 4, 8])
        Isech_via_g = mp.quad(lambda r: g(r) / mp.cosh(r / 2), [-7, -2, 0, 2, 7])   # the Fourier transform of sech(pi t) is sech(r/2)
        rec["V5 rotation by pi"] = {"integral": float(ell), "int h sech": float(Isech), "int g sech(r/2)": float(Isech_via_g),
                                    "c = integral / int h sech": float(ell / Isech), "u - sinh^2(rho) spread": float(spread)}
        I6 = 2 * mp.quad(lambda t: h(t) * t * mp.tanh(mp.pi * t), [0, 1, 2, 4, 8]) / (4 * mp.pi)
        rec["V6 identity"] = {"k(0)": 1.0, "(1/4pi) int h t tanh": float(I6)}
        out[f"k(u)=exp(-{a}u)"] = rec
        print(json.dumps(rec, indent=1), flush=True)
    dump(out, "mp16_vchecks.json")


# ----------------------------------------------------------------------------------------------- arithmetic
def primes_of(n):
    return [int(p) for p in pari.factor(n)[0]] if n > 1 else []


def divisors_from(ps_exp):
    divs = [1]
    for p, e in ps_exp:
        divs = [d * p ** j for d in divs for j in range(e + 1)]
    return divs


def kron(a, b):
    return int(pari.kronecker(a, b))


_L1 = {}


def L1(D):
    """L(1, chi_D) for a fundamental discriminant D > 0, from PARI: 2 h log(eps) / sqrt(D)."""
    if D not in _L1:
        q = pari.quadclassunit(D)
        _L1[D] = 2 * int(q[0]) * float(q[3]) / math.sqrt(D)
    return _L1[D]


def jacobi_vec(a, m):
    """Jacobi symbol (a/m) for a scalar a and an array of odd m > 0."""
    a = np.mod(np.int64(a), m).astype(np.int64)
    m = m.astype(np.int64).copy()
    t = np.ones(m.shape, dtype=np.int64)
    act = a != 0
    while act.any():
        while True:
            ev = act & (a % 2 == 0)
            if not ev.any():
                break
            a[ev] //= 2
            r = m[ev] % 8
            t[ev] = np.where((r == 3) | (r == 5), -t[ev], t[ev])
        idx = np.where(act)[0]
        aa, mm = a[idx].copy(), m[idx].copy()
        t[idx] = np.where((aa % 4 == 3) & (mm % 4 == 3), -t[idx], t[idx])
        a[idx], m[idx] = mm % aa, aa
        act = a != 0
    return np.where(m == 1, t, 0)


def kron_vec(D, n):
    n = n.astype(np.int64).copy()
    e2 = np.zeros(n.shape, dtype=np.int64)
    while True:
        ev = n % 2 == 0
        if not ev.any():
            break
        n[ev] //= 2
        e2[ev] += 1
    k2 = 0 if D % 2 == 0 else (1 if D % 8 in (1, 7) else -1)
    two = np.where(e2 == 0, 1, (k2 ** e2) if k2 != 0 else 0)
    return two * jacobi_vec(D, n)


def L1_series(D):
    """Second route: L(1, chi_D) = sum chi(n) [erfc(n sqrt(pi/D))/n + E1(pi n^2/D)/sqrt(D)]  (even primitive chi, D > 0)."""
    nmax = int(math.sqrt(D * 40 / math.pi)) + 10
    n = np.arange(1, nmax + 1, dtype=np.int64)
    chi = kron_vec(D, n).astype(float)
    x = math.pi * n.astype(float) ** 2 / D
    return float(np.sum(chi * (erfc(np.sqrt(x)) / n + exp1(x) / math.sqrt(D))))


def h_neg(D):
    """Number of reduced primitive positive definite forms of discriminant D < 0."""
    cnt, a = 0, 1
    while 3 * a * a <= -D:
        for b in range(-a + 1, a + 1):
            if (b * b - D) % (4 * a):
                continue
            c = (b * b - D) // (4 * a)
            if c < a or (b < 0 and a == c) or math.gcd(math.gcd(a, abs(b)), c) != 1:
                continue
            cnt += 1
        a += 1
    return cnt


def W_of(Delta, N, Q):
    """W(Delta) for the coset W_Q (or R W_Q) at level N; Delta > 0 and not a square."""
    cd = pari.coredisc(Delta, 1)
    D, c0 = int(cd[0]), int(cd[1])
    L = L1(D)
    fac = [(int(p), int(e)) for p, e in zip(*pari.factor(c0))] if c0 > 1 else []
    tot = 0.0
    for f in divisors_from(fac):
        c = c0 // f
        term = L / f
        for p, _ in fac:
            if c % p == 0:
                term *= 1 - kron(D, p) / p
        for p in primes_of(N // Q):
            term *= (p + 1) if f % p == 0 else (1 + (0 if c % p == 0 else kron(D, p)))
        tot += term
    return tot, D, c0


def elliptic_count(N, Q):
    """Fixed points of W_Q on X0(N): positive definite forms of discriminant -4Q/f^2 with the local factors at p | N/Q."""
    tot = 0.0
    for f in (1, 2):
        if (4 * Q) % (f * f):
            continue
        Dp = -4 * Q // (f * f)
        if Dp % 4 not in (0, 1):
            continue
        w = {-4: 2, -3: 3}.get(Dp, 1)
        rho = 1
        for p in primes_of(N // Q):
            rho *= (p + 1) if f % p == 0 else 1 + kron(Dp, p)
        tot += h_neg(Dp) * rho / w
    return tot


def mirror_length(N, Q):
    W, _, _ = W_of(4 * Q, N, Q)
    return 2 * math.sqrt(Q) * W


def geodesic_weights(N, Q, lmax):
    """Per s: (s, length, weight 2 sqrt(Q) W) for the hyperbolic (det +Q) and glide (det -Q) classes of the coset."""
    hyp, gli = [], []
    s = 1
    while True:
        x = s * math.sqrt(Q) / 2
        l_h = 2 * math.acosh(x) if x > 1 else None
        l_g = 2 * math.asinh(x)
        if l_g > lmax and (l_h is None or l_h > lmax):
            break
        if l_h is not None and l_h <= lmax:
            W, D, c0 = W_of(Q * Q * s * s - 4 * Q, N, Q)
            hyp.append((s, l_h, 2 * math.sqrt(Q) * W, D, c0))
        if l_g <= lmax:
            W, D, c0 = W_of(Q * Q * s * s + 4 * Q, N, Q)
            gli.append((s, l_g, 2 * math.sqrt(Q) * W, D, c0))
        s += 1
    return hyp, gli


# ----------------------------------------------------------------------------------------------- test functions, twisted traces
def testfun(r0, sig):
    h = lambda r: np.exp(-(r - r0) ** 2 / (2 * sig * sig)) + np.exp(-(r + r0) ** 2 / (2 * sig * sig))
    g = lambda u: (2 * sig / math.sqrt(2 * math.pi)) * np.cos(r0 * u) * np.exp(-sig * sig * u * u / 2)
    return h, g


def int_h_sech(h):
    return 2 * quad(lambda r: h(r) / math.cosh(math.pi * r), 0, 40, limit=400, epsabs=1e-14, epsrel=1e-13)[0]


class Coset:
    def __init__(self, N, Q, lmax, c_ell, c_refl, cache=None):
        self.N, self.Q = N, Q
        self.c_ell, self.c_refl = c_ell, c_refl
        if cache and os.path.exists(cache):                         # weights saved by an earlier run (same N, Q, lmax)
            z = np.load(cache)
            assert int(z["N"]) == N and int(z["Q"]) == Q and float(z["lmax"]) == lmax
            self.E, self.L = float(z["E"]), float(z["L"])
            self.hyp = [tuple(r) for r in z["hyp"].tolist()]
            self.gli = [tuple(r) for r in z["gli"].tolist()]
            return
        self.E = elliptic_count(N, Q)
        self.L = mirror_length(N, Q)
        self.hyp, self.gli = geodesic_weights(N, Q, lmax)
        if cache:
            np.savez(cache, N=N, Q=Q, lmax=lmax, E=self.E, L=self.L, hyp=np.array(self.hyp, dtype=float), gli=np.array(self.gli, dtype=float))

    def parts(self, h, g):
        """Returns the four pieces of T_{W_Q} - T_{RW_Q}: elliptic, -mirror, hyperbolic, -glide."""
        ell = self.E * self.c_ell * int_h_sech(h)
        mir = self.L * self.c_refl * float(g(0.0))
        hy = sum(float(g(l)) * w for _, l, w, _, _ in self.hyp)
        gl = sum(float(g(l)) * w for _, l, w, _, _ in self.gli)
        return {"elliptic": ell, "mirror": -mir, "hyperbolic": hy, "glide": -gl}


def constants():
    v = json.load(open(os.path.join(HERE, "mp16_vchecks.json")))
    rec = v["k(u)=exp(-1u)"]
    c_ell = rec["V5 rotation by pi"]["c = integral / int h sech"]
    c_refl = rec["V4 reflection"]["integral"] / (2 * rec["V4 reflection"]["g(0)/2"])
    return c_ell, c_refl


# ----------------------------------------------------------------------------------------------- arithmetic checks
def arith():
    out = {}
    # elliptic counts at 143 against Riemann-Hurwitz on the Atkin-Lehner dimensions (1,6,5,1) of S_2(143) by (eps11, eps13)
    dims = {(1, 1): 1, (1, -1): 6, (-1, 1): 5, (-1, -1): 1}
    g13 = 13
    rh = {}
    for Q, sel in ((11, lambda a, b: a), (13, lambda a, b: b), (143, lambda a, b: a * b)):
        gq = sum(d for (a, b), d in dims.items() if sel(a, b) == 1)
        rh[Q] = (2 * g13 - 2) - 2 * (2 * gq - 2)                      # Riemann-Hurwitz: 2g-2 = 2(2g'-2) + #fixed
    out["elliptic"] = {str(Q): {"class numbers": elliptic_count(143, Q), "Riemann-Hurwitz": rh[Q]} for Q in (11, 13, 143)}
    out["h(-143), h(-572), h(-52), h(-44), h(-11)"] = [h_neg(-143), h_neg(-572), h_neg(-52), h_neg(-44), h_neg(-11)]
    out["PARI qfbclassno same"] = [int(pari.qfbclassno(d)) for d in (-143, -572, -52, -44, -11)]
    # the R W_143 mirrors
    q = pari.quadclassunit(572)
    out["Q(sqrt 143)"] = {"h (wide)": int(q[0]), "log eps": float(q[3]), "2 log(12+sqrt143)": 2 * math.log(12 + math.sqrt(143)),
                          "mirror length L_143": mirror_length(143, 143), "4 log(12+sqrt143)": 4 * math.log(12 + math.sqrt(143))}
    out["mirror lengths at 143"] = {str(Q): mirror_length(143, Q) for Q in (11, 13, 143)}
    # second route for L(1, chi_D): the erfc/E1 series against PARI on a sample of the discriminants actually used
    rng = np.random.default_rng(16)
    Ds = set()
    for Q in (11, 13, 143):
        for s in rng.integers(1, 3000 if Q < 100 else 900, 25):
            for sgn in (-1, 1):
                Ds.add(int(pari.coredisc(Q * Q * int(s) * int(s) + sgn * 4 * Q)))
    Ds |= {5, 8, 13, 429, 572, 19877}
    worst, rows = 0.0, []
    for D in sorted(Ds):
        a, b = L1(D), L1_series(D)
        worst = max(worst, abs(a - b) / abs(a))
        rows.append((D, a, b))
    out["L(1,chi) two routes"] = {"n": len(rows), "max relative difference": worst, "largest D": max(Ds),
                                  "examples": rows[:3] + rows[-3:]}
    print(json.dumps(out, indent=1), flush=True)
    dump(out, "mp16_arith.json")


# ----------------------------------------------------------------------------------------------- the spectral side at 143
def odd_spectrum():
    d = json.load(open(os.path.join(HERE, "..", "mp15", "mp15_numbers.json")))
    L = d["level143"]
    rows = [(r["R"], sec) for sec, lst in L["sign_change_roots"].items() for r in lst] + [(r["R"], r["sector"]) for r in L["dip_roots"]]
    spec = {}
    for R, sec in rows:
        if not sec.startswith("odd"):
            continue
        a = 1 if "e11=+1" in sec else -1
        b = 1 if "e13=+1" in sec else -1
        spec.setdefault((a, b), []).append(R)
    return {k: sorted(v) for k, v in spec.items()}


def level143():
    c_ell, c_refl = constants()
    print(f"constants from the V-checks: c_ell = {c_ell:.12f}, c_refl = {c_refl:.12f}", flush=True)
    cos = {Q: Coset(143, Q, 23.5, c_ell, c_refl, cache=os.path.join(HERE, f"weights_143_{Q}.npz")) for Q in (11, 13, 143)}
    print(f"cosets built ({time.time() - T0:.0f} s): " + ", ".join(f"Q={Q}: E={c.E:g}, L={c.L:.6f}, {len(c.hyp)} hyp, {len(c.gli)} glide" for Q, c in cos.items()), flush=True)
    spec = odd_spectrum()
    eps = {11: lambda a, b: a, 13: lambda a, b: b, 143: lambda a, b: a * b}
    out = {"constants": {"c_ell": c_ell, "c_refl": c_refl}, "spectrum (odd, MP15)": {f"{k}": v for k, v in spec.items()},
           "cosets": {str(Q): {"fixed points": c.E, "mirror length": c.L, "hyperbolic classes (by s)": len(c.hyp), "glide classes (by s)": len(c.gli),
                               "shortest hyperbolic": c.hyp[0][:3], "shortest glide": c.gli[0][:3]} for Q, c in cos.items()},
           "identities": [], "fricke twist": {}}
    for sig in (0.35, 0.45, 0.6):
        lmax = 23.5 if sig < 0.4 else 17.0
        for r0 in (0.0, 0.374, 0.556, 0.845):
            h, g = testfun(r0, sig)
            tau = 0.25 * quad(lambda r: 14 * r * h(r), 1.6, 60, limit=200)[0]
            for Q in (11, 13, 143):
                c = cos[Q]
                keep = lambda lst: [t for t in lst if t[1] <= lmax]
                saved = (c.hyp, c.gli)
                c.hyp, c.gli = keep(c.hyp), keep(c.gli)
                parts = c.parts(h, g)
                c.hyp, c.gli = saved
                rhs = sum(parts.values()) / 8
                lhs = 0.25 * sum(eps[Q](a, b) * sum(float(h(R)) for R in Rs) for (a, b), Rs in spec.items())
                row = {"sigma": sig, "r0": r0, "Q": Q, "LHS": lhs, "RHS": rhs, "diff": lhs - rhs, "tau": tau,
                       "pass": abs(lhs - rhs) <= tau + 1e-6, "parts/8": {k: v / 8 for k, v in parts.items()}}
                out["identities"].append(row)
                print(f"sigma {sig} r0 {r0:5.3f} Q {Q:3d}: LHS {lhs:+.6f} RHS {rhs:+.6f} diff {lhs - rhs:+.2e} tau {tau:.1e} {'PASS' if row['pass'] else 'FAIL'}", flush=True)
            if r0 == 0.0 or (sig == 0.45 and r0 == 0.374):
                c = cos[143]
                parts = c.parts(h, g)
                terms = [("elliptic", parts["elliptic"] / 2), ("mirror", parts["mirror"] / 2)]
                terms += [(f"hyperbolic s={s} (l={l:.4f})", float(g(l)) * w / 2) for s, l, w, _, _ in c.hyp[:6]]
                terms += [(f"glide s={s} (l={l:.4f})", -float(g(l)) * w / 2) for s, l, w, _, _ in c.gli[:6]]
                total = sum(parts.values()) / 2
                osc = (parts["hyperbolic"] + parts["glide"]) / 2
                spec_side = sum(float(h(R)) for (a, b), Rs in spec.items() if a * b == 1 for R in Rs) - sum(float(h(R)) for (a, b), Rs in spec.items() if a * b == -1 for R in Rs)
                out["fricke twist"][f"sigma={sig} r0={r0}"] = {"geometric total": total, "spectral (MP15, R<1.6)": spec_side,
                                                               "smooth (elliptic + mirror)": (parts["elliptic"] + parts["mirror"]) / 2,
                                                               "oscillatory (geodesics + glides)": osc, "oscillatory share": abs(osc) / abs(total),
                                                               "largest terms": sorted(terms, key=lambda t: -abs(t[1]))[:8]}
                print(f"  Fricke twist sigma {sig} r0 {r0}: geometric {total:+.6f} (smooth {(parts['elliptic'] + parts['mirror']) / 2:+.6f}, oscillatory {osc:+.6f}); spectral {spec_side:+.6f}", flush=True)
    dump(out, "mp16_level143.json")


# ----------------------------------------------------------------------------------------------- other levels
def nu2(N):
    return 0 if N % 4 == 0 else math.prod(1 + kron(-4, p) for p in primes_of(N))


def nu3(N):
    return 0 if N % 9 == 0 else math.prod(1 + kron(-3, p) for p in primes_of(N))


def cross():
    c_ell, c_refl = constants()
    levels = [N for N in range(100, 201) if all(e == 1 for e in [int(x) for x in pari.factor(N)[1]])]
    rows = []
    for N in levels:
        c = Coset(N, N, 17.0, c_ell, c_refl)
        rec = {"N": N, "primes": primes_of(N), "nu2": nu2(N), "nu3": nu3(N), "fixed points": c.E, "mirror length": c.L,
               "shortest W_N geodesic": c.hyp[0][1] if c.hyp else None}
        for sig in (0.6, 0.45):
            h, g = testfun(0.0, sig)
            p = c.parts(h, g)
            rec[f"G sigma={sig}"] = sum(p.values()) / 2
            rec[f"S sigma={sig}"] = (p["elliptic"] + p["mirror"]) / 2
        rows.append(rec)
    out = {"levels": rows}
    for sig in (0.6, 0.45):
        G = np.array([r[f"G sigma={sig}"] for r in rows])
        S = np.array([r[f"S sigma={sig}"] for r in rows])
        order = sorted(rows, key=lambda r: -r[f"G sigma={sig}"])
        rank = [r["N"] for r in order].index(143) + 1
        out[f"sigma={sig}"] = {"pearson(S,G)": float(np.corrcoef(S, G)[0, 1]), "R^2 of G on S": float(np.corrcoef(S, G)[0, 1] ** 2),
                               "rank of 143": rank, "of": len(rows), "top 6": [(r["N"], round(r[f"G sigma={sig}"], 4)) for r in order[:6]],
                               "G_133": next(r[f"G sigma={sig}"] for r in rows if r["N"] == 133), "G_143": next(r[f"G sigma={sig}"] for r in rows if r["N"] == 143),
                               "G_145": next(r[f"G sigma={sig}"] for r in rows if r["N"] == 145)}
    cand = [115, 119, 123, 141, 155, 159, 161, 177, 187]
    assert all(nu2(N) == 0 and nu3(N) == 0 and len(primes_of(N)) == 2 for N in cand + [143])
    for sig in (0.45, 0.6):
        ordc = sorted(cand, key=lambda N: -next(r[f"G sigma={sig}"] for r in rows if r["N"] == N))
        out[f"C4 selection sigma={sig}"] = {"ordered": [(N, round(next(r[f"G sigma={sig}"] for r in rows if r["N"] == N), 4)) for N in ordc],
                                           "predict Fricke +1": ordc[:2], "predict Fricke -1": ordc[-2:]}
    print(json.dumps({k: v for k, v in out.items() if k != "levels"}, indent=1), flush=True)
    dump(out, "mp16_cross.json")


if __name__ == "__main__":
    {"vchecks": vchecks, "arith": arith, "level143": level143, "cross": cross}[sys.argv[1]]()
    print(f"done in {time.time() - T0:.0f} s")
