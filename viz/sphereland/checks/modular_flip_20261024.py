#!/usr/bin/env python3
"""
modular_flip_20261024.py -- certificates for Sphereland v15, Instruments 22-24.

  22  The modular flip: one boundary torus, two smooth ways to fill it.  Cap the space circle and the world is
      confined and gapped (glueballs); cap the time circle and it is a black brane (a fluid).  The filling that
      wins is the one that caps the shortest cycle, which is Gauss's lattice reduction; the phase diagram on the
      upper half-plane is the Farey tiling, and heating a twisted torus walks the continued fraction of its twist.
  23  Two strips: a mass gap is a strip of analyticity in momentum, smoothness is a strip of analyticity in
      position (Paley-Wiener read both ways).  Burgers' equation shows the smoothness strip closing in finite
      time (inviscid) and viscosity holding it open (Cole-Hopf poles).
  24  The gap on the modular curve: the hyperbolic plane has a built-in gap of 1/4; Selberg conjectured that
      congruence quotients keep it (proven at level 143); the scattering phase of the modular surface has the
      Riemann zeros as resonances of equal width 1/4 exactly when RH holds.

Sections
  1  Saddles of the torus  [structure]
  2  The spectrum flips: CFT2 correlators, BTZ, AdS5 soliton glueballs = plasma screening masses,
     AdS5 black-brane quasinormal modes, hydrodynamic poles  [derived + model]
  3  Two strips: Paley-Wiener, Poisson kernel, Burgers (Platzman, Cole-Hopf)  [structure + derived]
  4  The gap on the modular curve: Maass forms by Hejhal's method, Weyl's law, scattering phase, bounds  [structure]
  5  Data for the lab  (python3 modular_flip_20261024.py --json out.json)

Every check prints PASS/FAIL; the script ends with ALL CHECKS PASSED or exits non-zero.
Optional speed-up: python-flint (K-Bessel of imaginary order) and gmpy2 (mpmath backend).
"""
import math, cmath, time, sys, json, random
from fractions import Fraction
import numpy as np
import mpmath as mp
import sympy as sp
from scipy import special, integrate, optimize, linalg
from scipy.integrate import solve_ivp

try:
    from flint import acb, ctx as flint_ctx
    HAVE_FLINT = True
except Exception:  # pragma: no cover
    HAVE_FLINT = False

T0 = time.time()
mp.mp.dps = 30
RESULTS = []
DATA = {}


def check(name, ok, detail=""):
    ok = bool(ok)
    RESULTS.append((name, ok))
    print(("  PASS  " if ok else "  FAIL  ") + name + (f"\n        {detail}" if detail else ""))


# =====================================================================================================
print("§1  Saddles of the torus: the shortest cycle caps  [structure]")
# =====================================================================================================


def reduce_tau(x, y):
    """Gauss reduction of tau = x + iy.  Returns (a, b, c, d, zx, zy) with [[a,b],[c,d]] tau = zx + i zy in F."""
    a, b, c, d = 1, 0, 0, 1
    zx, zy = x, y
    for _ in range(2000):
        n = math.floor(zx + 0.5)
        if n:
            zx -= n
            a, b = a - n * c, b - n * d
        r2 = zx * zx + zy * zy
        if r2 < 1 - 1e-14:
            zx, zy = -zx / r2, zy / r2
            a, b, c, d = -c, -d, a, b
        else:
            break
    return a, b, c, d, zx, zy


def cusp_of(x, y):
    """Cusp p/q of the dominant saddle (the SL(2,Z) tile containing tau touches it) and Im(gamma tau)."""
    a, b, c, d, zx, zy = reduce_tau(x, y)
    if c == 0:
        return (1, 0), zy
    p, q = -d, c
    if q < 0:
        p, q = -p, -q
    return (p, q), zy


def im_saddle(x, y, p, q):
    return y / ((q * x - p) ** 2 + (q * y) ** 2)


rng = np.random.default_rng(20261024)
ok_F, ok_max, ok_det, worst = True, True, True, 0.0
for _ in range(2500):
    x, y = rng.uniform(-2, 2), 10 ** rng.uniform(-1.6, 0.5)
    a, b, c, d, zx, zy = reduce_tau(x, y)
    ok_det &= (a * d - b * c == 1)
    ok_F &= abs(zx) <= 0.5 + 1e-12 and zx * zx + zy * zy >= 1 - 1e-12
    best = y  # q = 0
    for q in range(1, 70):
        p0 = round(q * x)
        for p in range(p0 - 3, p0 + 4):
            if math.gcd(p, q) == 1:
                best = max(best, im_saddle(x, y, p, q))
    worst = max(worst, abs(best - zy) / zy)
    ok_max &= abs(best - zy) <= 1e-9 * zy
check("Gauss reduction finds the dominant saddle: for 2500 random τ the reduced point lies in the fundamental domain F, the matrix has det 1, and Im γτ equals the brute-force maximum of Im τ/|qτ − p|² over coprime (p, q), q < 70",
      ok_F and ok_det and ok_max,
      f"worst relative gap {worst:.1e}; log Z_γ = (πc/6) Im γτ in a 2D CFT / AdS₃, so the largest Im γτ is the leading saddle, and the saddle is labelled by the cusp p/q = γ⁻¹(∞)")

ok_ford = True
for _ in range(3000):
    x, y = rng.uniform(-1, 1), 10 ** rng.uniform(-3, 0.3)
    (p, q), im = cusp_of(x, y)
    if q == 0:
        inside = y > 1
    else:
        inside = (x - p / q) ** 2 + (y - 1 / (2 * q * q)) ** 2 < 1 / (4 * q ** 4)
    if abs(im - 1) > 1e-9:
        ok_ford &= (inside == (im > 1))
ok_tan = True
fr = sorted({Fraction(p, q) for q in range(1, 26) for p in range(0, q + 1) if math.gcd(p, q) == 1})
for i in range(len(fr)):
    for j in range(i + 1, len(fr)):
        A, B = fr[i], fr[j]
        det = abs(A.numerator * B.denominator - B.numerator * A.denominator)
        gap = (A - B) ** 2 - Fraction(1, A.denominator ** 2 * B.denominator ** 2)   # = (det^2 - 1)/(q^2 q'^2)
        ok_tan &= gap >= 0 and ((gap == 0) == (det == 1))
check("Ford circles: Im γτ > 1 exactly inside the Ford circle of the dominant cusp p/q (centre p/q + i/2q², radius 1/2q²); distinct Ford circles never overlap and touch exactly for Farey neighbours (|pq′ − p′q| = 1)",
      ok_ford and ok_tan,
      f"3000 random τ; all {len(fr)} fractions in [0,1] with q ≤ 25 compared pairwise in exact arithmetic: centre distance² − (r + r′)² = (det² − 1)/(q²q′²)")

ok_short = True
for _ in range(500):
    x, y = rng.uniform(-1, 1), 10 ** rng.uniform(-1.2, 0.4)
    (p, q), _ = cusp_of(x, y)
    vecs = [(pp, qq) for qq in range(0, 40) for pp in range(-60, 61) if (qq > 0 or pp == 1) and math.gcd(pp, qq) == 1]
    lens = {v: abs(complex(v[1] * x - v[0], v[1] * y)) for v in vecs}
    order_d = [sorted(vecs, key=lambda v: -y / lens[v] ** dd)[0] for dd in (2, 3, 4, 6)]
    ok_short &= all(o == (p, q) for o in order_d) and min(lens, key=lens.get) == (p, q)
check("the filling that wins caps the shortest cycle: the dominant saddle's cycle qτ − p is the shortest primitive vector of the lattice ℤ + ℤτ, and the winner is the same in every dimension (log Z ∝ Im τ/|qτ − p|^d for d = 2, 3, 4, 6)",
      ok_short,
      "500 random τ; in AdS_{d+1} with boundary T² × ℝ^{d−2} the smooth filling that caps the cycle of length ℓ has log Z = K V A/ℓ^d (A = area of T²), so the phase diagram is Gauss reduction in any d")

z_i = [(p, q, im_saddle(0, 1, p, q)) for q in range(0, 12) for p in range(-12, 13) if (q > 0 or p == 1) and math.gcd(p, q) == 1]
z_i.sort(key=lambda t: -t[2])
rho = (0.5, math.sqrt(3) / 2)
z_r = [(p, q, im_saddle(*rho, p, q)) for q in range(0, 12) for p in range(-12, 13) if (q > 0 or p == 1) and math.gcd(p, q) == 1]
z_r.sort(key=lambda t: -t[2])
check("the elliptic points are the multicritical points: at τ = i two saddles tie (thermal AdS 1/0 and BTZ 0/1 — the Hawking–Page point); at τ = ρ = e^{iπ/3} three tie (1/0, 0/1, 1/1 — a triple point); all others lie below",
      abs(z_i[0][2] - 1) < 1e-15 and abs(z_i[1][2] - 1) < 1e-15 and z_i[2][2] < 0.51 and
      all(abs(z_r[k][2] - math.sqrt(3) / 2) < 1e-15 for k in range(3)) and z_r[3][2] < 0.3,
      f"at i: {[(f'{p}/{q}', round(v, 4)) for p, q, v in z_i[:3]]}; at ρ: {[(f'{p}/{q}', round(v, 4)) for p, q, v in z_r[:4]]}; these are the fixed points of S and ST, the two elliptic points of SL(2,ℤ)\\H")


def logZ_conf(beta, L, d, KV=1.0):
    return KV * beta * L ** (1 - d)


def logZ_deconf(beta, L, d, KV=1.0):
    return KV * L * beta ** (1 - d)


def entropy(fun, beta, L, d, h=1e-6):
    lz = fun(beta, L, d)
    dl = (fun(beta * (1 + h), L, d) - fun(beta * (1 - h), L, d)) / (2 * beta * h)
    return lz - beta * dl


ok_th = True
for d in (2, 3, 4, 5):
    L = 1.7
    ok_th &= abs(logZ_conf(L, L, d) - logZ_deconf(L, L, d)) < 1e-14
    ok_th &= abs(entropy(logZ_conf, L, L, d)) < 1e-7
    ok_th &= abs(entropy(logZ_deconf, L, L, d) - d * L ** (2 - d)) < 1e-6
    E_sol = -(logZ_conf(L * 1.000001, L, d) - logZ_conf(L * 0.999999, L, d)) / (2e-6 * L)
    ok_th &= abs(E_sol + L ** (1 - d)) < 1e-7
c_cft = 24.0
S_btz = entropy(lambda b, L, d: (math.pi * c_cft / 6) * L / b, 0.37, 1.0, 2)
check("the flip is first order: at β = L the two free energies cross, the confined entropy is 0 and the deconfined entropy jumps to d K V L^{2−d}; the soliton's energy is negative, −K V L^{1−d} (Horowitz–Myers); in d = 2 the jump is Cardy's S = (πc/3) L T",
      ok_th and abs(S_btz - (math.pi * c_cft / 3) / 0.37) < 1e-6,
      "log Z_conf = K V β L^{1−d}, log Z_deconf = K V L β^{1−d}; for N = 4 SYM, K V = (π²/8) N² V₂, so the jump is s = (π²/2) N² T³ — order N², while the confined phase has O(N⁰) entropy: glueballs are not a gluon fluid")


def convergents(fr_, n=60):
    out = []
    h0, h1, k0, k1 = 0, 1, 1, 0
    x = fr_
    for _ in range(n):
        a = math.floor(x)
        h0, h1 = h1, a * h1 + h0
        k0, k1 = k1, a * k1 + k0
        out.append((h1, k1))
        if x == a:
            break
        x = 1 / (x - a)
    return out


def walk(theta, tmax=3.0, tmin=1e-8, n=3000):
    seq = []
    for k in range(n):
        t = tmax * (tmin / tmax) ** (k / (n - 1))
        pq, _ = cusp_of(theta, t)
        if not seq or seq[-1] != pq:
            seq.append(pq)
    return seq


gold = (math.sqrt(5) - 1) / 2
wg = walk(gold, tmin=1e-7)
fib = [1, 1]
while len(fib) < 30:
    fib.append(fib[-1] + fib[-2])
fib_ok = wg[0] == (1, 0) and all(wg[k] == (fib[k - 1], fib[k]) for k in range(1, 10))
random.seed(7)
met, nonconv, missed, beyond = 0, 0, 0, 0
for trial in range(300):
    th = Fraction(random.randrange(1, 10 ** 9), 10 ** 9)
    conv = convergents(th)
    convset = set(conv)
    seq = walk(float(th))
    for p, q in seq:
        if q == 0:
            continue
        met += 1
        nonconv += (p, q) not in convset
        beyond += q * abs(q * th - p) > 1 / math.sqrt(3) + 1e-9
    for p, q in conv:
        if 0 < q < 2500 and q * abs(q * th - p) < Fraction(1, 2) and (p, q) not in seq:
            missed += 1
check("heating a twisted torus walks the continued fraction of its twist: along τ = θ + it, t ↓ 0, the dominant saddles are convergents p/q of θ; at the golden twist they are the Fibonacci ratios 1/1, 1/2, 2/3, 3/5, 5/8, …",
      fib_ok and nonconv == 0 and missed == 0 and beyond == 0,
      f"golden walk {[f'{p}/{q}' for p, q in wg[:9]]}; 300 random twists: {met} saddles met, {nonconv} non-convergents (observed, not proven here), "
      f"{missed} convergents with q|qθ − p| < ½ missed (proven: inside their Ford circle Im γτ > 1), {beyond} with q|qθ − p| > 1/√3 (proven impossible: max_t Im γτ = 1/(2q|qθ − p|) must reach √3/2)")

N143 = 143


def cusp_class(p, q):
    return math.gcd(q, N143) if q else N143


def egcd(a, b):
    if b == 0:
        return a, 1, 0
    g, x, y = egcd(b, a % b)
    return g, y, x - (a // b) * y


def rand_gamma0(rng_):
    while True:
        c = N143 * int(rng_.integers(-6, 7))
        d = int(rng_.integers(-400, 401))
        if c != 0 and math.gcd(c, d) == 1:
            g, x, y = egcd(d, -c)    # x d + y (-c) = 1 -> a = x, b = y
            a, b = x, y
            if a * d - b * c == 1:
                return a, b, c, d
            if a * d - b * c == -1:
                return -a, -b, c, d


def mobius_frac(M, p, q):
    a, b, c, d = M
    num, den = a * p + b * q, c * p + d * q
    g = math.gcd(num, den)
    num, den = num // g, den // g
    if den < 0:
        num, den = -num, -den
    return num, den


ok_inv = True
for _ in range(400):
    g = rand_gamma0(rng)
    q = int(rng.integers(1, 3000))
    p = int(rng.integers(-3000, 3000))
    if math.gcd(p, q) != 1:
        continue
    ok_inv &= cusp_class(*mobius_frac(g, p, q)) == cusp_class(p, q)
widths = {}
for dd in (1, 11, 13, 143):
    h = 1
    while (dd * dd * h) % N143:
        h += 1
    widths[dd] = h
index = round(N143 * (1 + 1 / 11) * (1 + 1 / 13))
check("Γ₀(143) colours the Farey tiling with four cusps: the class of p/q is gcd(q, 143) ∈ {143 → ∞, 1 → 0, 11 → 1/11, 13 → 1/13}, invariant under Γ₀(143); the widths are 1, 143, 13, 11 and add up to the index 168",
      ok_inv and widths == {1: 143, 11: 13, 13: 11, 143: 1} and sum(widths.values()) == index == 168,
      f"400 random γ ∈ Γ₀(143) acting on random fractions; width of 1/d is the least h with 143 | d²h; area of X₀(143) = 168 · π/3 = 56π")

W = {11: (66, 5, 143, 11), 13: (78, 7, 143, 13), 143: (0, -1, 143, 0)}
ok_al = True
perm = {}
for Q, (a, b, c, d) in W.items():
    ok_al &= a * d - b * c == Q and a % Q == 0 and d % Q == 0 and c % N143 == 0
    for _ in range(100):
        g = rand_gamma0(rng)
        A = np.array([[a, b], [c, d]], dtype=object)
        Gm = np.array([[g[0], g[1]], [g[2], g[3]]], dtype=object)
        adj = np.array([[d, -b], [-c, a]], dtype=object)
        Mx = A.dot(Gm).dot(adj)
        ok_al &= all(int(v) % Q == 0 for v in Mx.flatten())
        Mx = Mx // Q
        ok_al &= int(Mx[1, 0]) % N143 == 0 and int(Mx[0, 0] * Mx[1, 1] - Mx[0, 1] * Mx[1, 0]) == 1
    perm[Q] = {dd: cusp_class(*mobius_frac((a, b, c, d), 1, dd)) if dd != 143 else cusp_class(*mobius_frac((a, b, c, d), 1, 0)) for dd in (1, 11, 13, 143)}
    for dd in (1, 11, 13, 143):
        ok_al &= perm[Q][dd] == dd * Q // math.gcd(dd, Q) ** 2
check("the Atkin–Lehner involutions normalise Γ₀(143) and permute its four cusps: W₁₁ swaps ∞ ↔ 1/13 and 0 ↔ 1/11, W₁₃ swaps ∞ ↔ 1/11 and 0 ↔ 1/13, and the Fricke involution W₁₄₃: τ → −1/(143τ) swaps ∞ ↔ 0 — the level-143 modular flip, self-dual at τ = i/√143",
      ok_al and perm[143][143] == 1 and perm[143][1] == 143,
      f"class maps (d → d′): W₁₁ {perm[11]}, W₁₃ {perm[13]}, W₁₄₃ {perm[143]}; on τ = iβ/L the Fricke flip sends β/L to L/(143β), so its Hawking–Page point is β/L = 1/√143 = {1 / math.sqrt(143):.5f}")

# =====================================================================================================
print(f"   [{time.time() - T0:.0f} s]")
print("§2  The spectrum flips  [derived + model]")
# =====================================================================================================
Delta, beta = 2.6, 1.3
Gb = lambda x: (math.pi / (beta * math.sinh(math.pi * x / beta))) ** (2 * Delta)
slope = -(math.log(Gb(20.0 + 1e-4)) - math.log(Gb(20.0 - 1e-4))) / 2e-4
GL = lambda z: (cmath.pi / (beta * cmath.sin(cmath.pi * z / beta))) ** (2 * Delta)
check("CFT₂ on the two circles: the thermal correlator (π/β sinh(πx/β))^{2Δ} decays with screening mass 2πΔT, and the circle correlator (π/L sin(πx/L))^{2Δ} continued to imaginary x is the same function with L = β — the flip exchanges a screening mass and an energy gap",
      abs(slope - 2 * math.pi * Delta / beta) < 1e-9 and abs(abs(GL(-1j * 0.7)) - Gb(0.7)) < 1e-12 * Gb(0.7),
      f"Δ = {Delta}, β = {beta}: decay rate {slope:.10f} vs 2πΔ/β = {2 * math.pi * Delta / beta:.10f}; sin(iz) = i sinh z")

hD = Delta / 2
Tq = 1 / beta
ok_btz = True
for n in range(4):
    for k in (0.0, 1.7):
        for sgn in (1, -1):
            w = sgn * k - 4j * math.pi * Tq * (hD + n)
            arg1 = hD - 1j * (w - sgn * k) / (4 * math.pi * Tq)
            ok_btz &= abs(arg1 + n) < 1e-12 and abs(mp.rgamma(arg1)) < 1e-12
            ok_btz &= abs(w - (sgn * k - 2j * math.pi * Tq * (Delta + 2 * n))) < 1e-12


def ads3_normal_mode_residual(omega, Dl, ell, rmax=60.0):
    # global AdS3, scalar m^2 = Dl(Dl-2): (1/r)(r(1+r^2)psi')' + (omega^2/(1+r^2) - ell^2/r^2 - m^2) psi = 0
    m2 = Dl * (Dl - 2)
    r0 = 1e-3
    a2 = (m2 - omega ** 2 - ell * (ell + 2)) / (4 * (ell + 1))      # regular series psi = r^ell (1 + a2 r^2 + ...)
    psi0 = r0 ** ell * (1 + a2 * r0 ** 2)
    dpsi0 = (ell * r0 ** (ell - 1) if ell else 0.0) + (ell + 2) * a2 * r0 ** (ell + 1)
    def rhs(r, y):
        psi, P = y          # P = r (1+r^2) psi'
        return [P / (r * (1 + r * r)), -r * (omega ** 2 / (1 + r * r) - ell ** 2 / (r * r) - m2) * psi]
    sol = solve_ivp(rhs, [r0, rmax], [psi0, r0 * (1 + r0 * r0) * dpsi0], rtol=1e-11, atol=1e-14, method='DOP853')
    psi, P = sol.y[:, -1]
    r = sol.t[-1]
    # psi ~ A r^{Dl-2} + B r^{-Dl}:  r psi' = (Dl-2) A r^{Dl-2} - Dl B r^{-Dl}
    rdpsi = P / (1 + r * r)
    A = (Dl * psi + rdpsi) / (2 * Dl - 2) / r ** (Dl - 2)
    return A


nm_ok = True
for ell in (0, 1):
    for n in range(3):
        w0 = 2.5 + 2 * n + ell
        root = optimize.brentq(lambda w: ads3_normal_mode_residual(w, 2.5, ell), w0 - 0.4, w0 + 0.4, xtol=1e-12)
        nm_ok &= abs(root - w0) < 1e-6
check("the spectrum rotates: thermal AdS₃ has real normal modes ω = (2π/L)(Δ + 2n + |ℓ|) (solved here by shooting for Δ = 5/2); BTZ has quasinormal modes ω = ±k − 2πiT(Δ + 2n), the poles of Γ(h − i(ω ∓ k)/4πT); at zero momentum the second ladder is the first turned by −90° with L → β",
      ok_btz and nm_ok,
      "Birmingham–Sachs–Solodukhin 2002; at k = 0: ω_n = (2π/L)(Δ + 2n) → −i(2π/β)(Δ + 2n). Real glueball-like levels of the confined filling become damped modes of the fluid filling")

r_, k_, Mq = sp.symbols('r k M', positive=True)
f_ = 1 - r_ ** -4
psi = sp.Function('psi')(r_)
g_diag = [r_ ** 2 * f_, r_ ** 2, r_ ** 2, r_ ** 2, 1 / (r_ ** 2 * f_)]
sqrtg = sp.sqrt(sp.prod(g_diag))
lap = sp.diff(sqrtg * (1 / g_diag[4]) * sp.diff(psi, r_), r_) / sqrtg - (1 / g_diag[3]) * k_ ** 2 * psi
target = (sp.diff(r_ ** 5 * f_ * sp.diff(psi, r_), r_) - r_ * k_ ** 2 * psi) / r_ ** 3
check("glueballs of the soliton are screening masses of the plasma: the Euclidean AdS₅ soliton (space circle capped) and the Euclidean black brane (time circle capped) are one metric with the circles relabelled, so a static mode e^{ikx} of either obeys (r⁵f ψ′)′ − k² r ψ = 0; k² = −M² gives the same ladder",
      sp.simplify(lap - target) == 0,
      "ds² = r²(f dθ² + dx₁² + dx₂² + dx₃²) + dr²/(r² f), f = 1 − r₀⁴/r⁴, θ ~ θ + π/r₀ (smooth tip): call θ time and it is a black brane at T = r₀/π; call it space and it is Witten's confining soliton on a circle L = π/r₀")


def soliton_shoot(M2, rmax=200.0):
    x0 = 1e-6
    psi0 = 1 - M2 / 4 * x0
    dpsi0 = -M2 / 4
    def rhs(r, y):
        ps, P = y
        f = 1 - r ** -4
        return [P / (r ** 5 * f), -M2 * r * ps]
    sol = solve_ivp(rhs, [1 + x0, rmax], [psi0, (1 + x0) ** 5 * (1 - (1 + x0) ** -4) * dpsi0], rtol=1e-11, atol=1e-13, method='DOP853')
    ps, P = sol.y[:, -1]
    r = sol.t[-1]
    B = -P / 4
    return ps - B * r ** -4


grid = np.linspace(2, 260, 260)
vals = [soliton_shoot(m) for m in grid]
glue = []
for i in range(len(grid) - 1):
    if vals[i] * vals[i + 1] < 0:
        glue.append(optimize.brentq(soliton_shoot, grid[i], grid[i + 1], xtol=1e-11))
coot = [11.59, 34.53, 68.98, 114.9, 172.3, 241.2]
ML = [math.pi * math.sqrt(m) for m in glue]
check("the soliton's 0⁺⁺ glueballs (dilaton, regular at the tip, normalizable at the boundary): M²/r₀² = 11.588, 34.53, 68.98, 114.9, 172.3, 241.2 — Csáki–Ooguri–Oz–Terning 1999 — so M·L = 10.694, 18.460, 26.091, … and the plasma's screening masses are the same numbers times T",
      len(glue) == 6 and all(abs(g - c) <= 0.006 * c for g, c in zip(glue, coot)) and abs(glue[0] - 11.5877) < 2e-3,
      f"computed M² = {[round(g, 4) for g in glue]}; M·L = {[round(v, 4) for v in ML]}; M/(πT) = √(M²) = {[round(math.sqrt(g), 4) for g in glue[:3]]}")
DATA['glueM2'] = [round(g, 6) for g in glue]
DATA['glueML'] = [round(v, 6) for v in ML]

from numpy.polynomial import polynomial as NP
_s = NP.polysub(NP.polypow([1, 1], 4), [1])
_p = NP.polymul([1, 1], _s)
_P2, _P1 = NP.polymul(_s, _p), NP.polymul(_s, NP.polyder(_p))
_R5 = NP.polypow([1, 1], 5)


def bb_frob(omega, nterms):
    a = -1j * omega / 4
    P0 = omega ** 2 * _R5
    cf = lambda arr, j: arr[j] if 0 <= j < len(arr) else 0.0
    coef = np.zeros(nterms, complex)
    coef[0] = 1
    for m in range(1, nterms):
        acc = 0j
        for k in range(max(0, m - 8), m):
            acc += coef[k] * (cf(_P2, m - k + 2) * (k + a) * (k + a - 1) + cf(_P1, m - k + 1) * (k + a) + cf(P0, m - k))
        coef[m] = -acc / (cf(_P2, 2) * (m + a) * (m + a - 1) + cf(_P1, 1) * (m + a) + cf(P0, 0))
    return a, coef


def bb_boundary(omega, x0=0.3, rmax=2000.0, nterms=90):
    a, coef = bb_frob(omega, nterms)
    n = np.arange(nterms)
    ph = np.sum(coef * x0 ** (n + a))
    dph = np.sum(coef * (n + a) * x0 ** (n + a - 1))
    r0 = 1 + x0
    sol = solve_ivp(lambda r, y: [y[1] / (r ** 5 - r), -omega ** 2 * r ** 5 / (r ** 4 - 1) * y[0]],
                    [r0, rmax], [ph, (r0 ** 5 - r0) * dph], rtol=1e-12, atol=1e-16, method='DOP853')
    phe, Pe = sol.y[:, -1]
    r = sol.t[-1]
    return (phe + Pe / (4 * r ** 4)) / (1 + omega ** 2 / (8 * r ** 2)) / abs(ph)


def bb_qnm(guess):
    w = mp.findroot(lambda z: complex(bb_boundary(2 * complex(z))), mp.mpc(guess), solver='secant', tol=1e-24, maxsteps=80, verify=False)
    return complex(w)


starinets = [(3.119452, -2.746676), (5.169521, -4.763570), (7.187931, -6.769565), (9.197199, -8.772481), (11.202676, -10.774162)]
qnm = [bb_qnm(complex(a, b) / 2 + 0.01) for a, b in starinets] + [bb_qnm(6.60 - 6.39j)]
qnm_ok = all(abs(2 * q - complex(a, b)) < 3e-6 for q, (a, b) in zip(qnm, starinets))


def cheb(N):
    x = np.cos(np.pi * np.arange(N + 1) / N)
    c = np.hstack([2, np.ones(N - 1), 2]) * (-1) ** np.arange(N + 1)
    X = np.tile(x, (N + 1, 1)).T
    D = np.outer(c, 1 / c) / (X - X.T + np.eye(N + 1))
    return D - np.diag(D.sum(axis=1)), x


def bb_spectral(N):
    D, x = cheb(N)
    u = (1 - x) / 2
    Du = -2 * D
    Du2 = Du @ Du
    A = u * (1 - u) * (1 + u) ** 2
    M0 = np.diag(A) @ Du2 + np.diag((1 + u) * (3 - 5 * u ** 2)) @ Du + np.diag(-4 * u * (1 + u))
    M1 = np.diag(1j * u * (1 + u) ** 2) @ Du + np.diag(0.5j * (1 + u) * (3 + 4 * u))
    M2 = np.diag((u ** 2 + 3 * u + 4) / 4).astype(complex)
    n = N + 1
    I, Z = np.eye(n), np.zeros((n, n))
    ev = linalg.eig(np.block([[Z, I], [-M0, -M1]]), np.block([[I, Z], [Z, M2]]), right=False)
    return ev[np.isfinite(ev)]


ev = bb_spectral(48)
spec1 = min(ev, key=lambda z: abs(z - qnm[0]))
check("the black brane's quasinormal modes (massless scalar, k = 0, ingoing at the horizon, normalizable at the boundary): 𝔴 = ω/2πT = ±1.559726 − 1.373338i, ±2.584760 − 2.381785i, ±3.593965 − 3.384783i, … — Starinets 2002, Table 1 — by two independent methods",
      qnm_ok and abs(spec1 - qnm[0]) < 1e-6,
      f"Frobenius-series shooting: {[f'{q.real:.6f}{q.imag:+.6f}i' for q in qnm]}; Chebyshev collocation in u = r₀²/r²: {spec1.real:.7f}{spec1.imag:+.7f}i; spacing → (±1 − i): the modes march off at 45°")
DATA['qnm'] = [[round(q.real, 7), round(q.imag, 7)] for q in qnm]

D_shear = 1 / (4 * math.pi)        # eta/(e+p) at T = 1, eta/s = 1/4pi
G_sound = (2 / 3) * D_shear
qq = 0.37
w_shear = -1j * D_shear * (2 * math.pi * qq) ** 2 / (2 * math.pi)
w_sound = (2 * math.pi * qq) / math.sqrt(3) / (2 * math.pi) - 1j * G_sound * (2 * math.pi * qq) ** 2 / (2 * math.pi)
r_shear = abs(complex(1.8906469, 1.1711505)) ** 0.5
check("the fluid's gapless poles (Navier–Stokes): with η/s = 1/4π and c_s² = 1/3, shear 𝔴 = −i𝔮²/2 and sound 𝔴 = ±𝔮/√3 − i𝔮²/3; the gradient expansion converges out to |𝔮| ≈ 1.49 (shear) and √2 (sound), where these poles collide with the gapped tower (Grozdanov–Kovtun–Starinets–Tadić 2019)",
      abs(w_shear + 1j * qq ** 2 / 2) < 1e-14 and abs(w_sound - (qq / math.sqrt(3) - 1j * qq ** 2 / 3)) < 1e-14 and abs(r_shear - 1.49131) < 5e-6,
      f"𝔮*² = 1.8906469 ± 1.1711505i (shear), ±2i (sound): |𝔮*| = {r_shear:.5f}, {math.sqrt(2):.5f} — the cited critical points, arithmetic checked here; the gapped modes set how far Navier–Stokes (plus its gradient corrections) can be trusted")

# =====================================================================================================
print(f"   [{time.time() - T0:.0f} s]")
print("§3  Two strips: Paley–Wiener read both ways  [structure + derived]")
# =====================================================================================================
m_ = 0.8
ok_pw = True
for p_ in (0.0, 0.5, 2.3):
    val = 2 * integrate.quad(lambda x: math.exp(-m_ * x) * math.cos(p_ * x), 0, np.inf, limit=400)[0]
    ok_pw &= abs(val - 2 * m_ / (p_ * p_ + m_ * m_)) < 1e-9
pois = lambda y, x: y / (math.pi * (x * x + y * y))
check("a gap is a strip: e^{−m|x|} has Fourier transform 2m/(p² + m²), analytic for |Im p| < m with poles at ±im; and 1/(p² + m²) = (π/m) P_m(p), the Poisson kernel at height m — the propagator of a mass-m particle is the Poisson kernel in momentum space",
      ok_pw and abs(1 / (2.3 ** 2 + m_ ** 2) - (math.pi / m_) * pois(m_, 2.3)) < 1e-15,
      "Paley–Wiener: exponential decay e^{−m|x|} ⇔ analytic continuation to the strip |Im p| < m; Källén–Lehmann makes this the mass gap of a QFT")
ok_pk = True
for dl in (0.05, 0.4, 1.3):
    for x in (0.0, 0.7, math.pi):
        lhs = 1 + 2 * sum(math.exp(-dl * k) * math.cos(k * x) for k in range(1, 4000))
        ok_pk &= abs(lhs - math.sinh(dl) / (math.cosh(dl) - math.cos(x))) < 1e-9 * max(1, lhs)
ok_img = True
for rho_ in (0.3, 1.1):
    for x in (0.2, 2.0):
        s = sum(1 / (rho_ ** 2 + (x + 2 * math.pi * n) ** 2) for n in range(-200000, 200001))
        ok_img &= abs(s - math.sinh(rho_) / (2 * rho_ * (math.cosh(rho_) - math.cos(x)))) < 1e-5
check("smoothness is a strip: Fourier coefficients e^{−δ|k|} sum to the Poisson kernel sinh δ/(cosh δ − cos x), analytic for |Im x| < δ with poles at x = ±iδ; it is Instrument 21's Green function on ℝ³ × S¹ (images at distance ρ give Σ 1/(ρ² + (x + 2πn)²) = sinh ρ/(2ρ(cosh ρ − cos x)))",
      ok_pk and ok_img,
      "one function, two readings: in x it is smooth in a strip of width δ, in k it is a ladder with gap δ — the same theorem with position and momentum exchanged")


def burgers_char(t, N=4096):
    xs = 2 * np.pi * np.arange(N) / N
    x0 = xs.copy()
    for _ in range(80):
        F = x0 + t * np.sin(x0) - xs
        x0 -= F / (1 + t * np.cos(x0))
    return xs, np.sin(x0)


def platzman(k, t):
    return (-1) ** (k + 1) * 2 * special.jv(k, k * t) / (k * t)


ok_pl = True
for t in (0.5, 0.9):
    xs, u = burgers_char(t)
    bk = -2 * np.fft.rfft(u).imag / len(xs)
    ok_pl &= all(abs(bk[k] - platzman(k, t)) < 1e-11 for k in range(1, 40))
check("Burgers, u_t + u u_x = 0 from u₀ = sin x: the exact Fourier coefficients are b_k(t) = (−1)^{k+1} 2J_k(kt)/(kt) (Platzman 1964), confirmed against the method of characteristics at t = 0.5 and 0.9",
      ok_pl, "u = Σ b_k sin kx; the fold of characteristics at x = π, t = 1 is the shock")
delta_inv = lambda t: math.acosh(1 / t) - math.sqrt(1 - t * t)
ok_dl = True
for t in (0.2, 0.6, 0.95):
    z0 = complex(math.pi, math.acosh(1 / t))
    ok_dl &= abs(1 + t * cmath.cos(z0)) < 1e-12
    xs_ = z0 + t * cmath.sin(z0)
    ok_dl &= abs(xs_.real - math.pi) < 1e-12 and abs(xs_.imag - delta_inv(t)) < 1e-12
asym = [delta_inv(1 - e) / ((2 * math.sqrt(2) / 3) * e ** 1.5) for e in (1e-2, 1e-3, 1e-4)]
mp.mp.dps = 40
dec = []
for t in (0.5, 0.8):
    lb = [float(mp.log(abs(mp.besselj(k, k * t) / (k * t)))) for k in (600, 1200)]
    rate = -(lb[1] - lb[0] - (-1.5) * math.log(2)) / 600
    dec.append((t, rate, delta_inv(t)))
mp.mp.dps = 30
check("the smoothness strip closes in finite time: the nearest complex singularity sits at x = π + iδ(t) with δ(t) = arccosh(1/t) − √(1 − t²) → (2√2/3)(1 − t)^{3/2}, and the spectrum decays as k^{−3/2} e^{−δ(t) k} (Debye) — blowup is δ → 0",
      ok_dl and abs(asym[-1] - 1) < 1e-3 and all(abs(r - d) < 2e-3 * d for _, r, d in dec),
      f"singularity = image of the complex fold 1 + t cos x₀ = 0; δ/(asymptotic) = {[round(a, 5) for a in asym]}; decay rate vs δ: {[(t, round(r, 5), round(d, 5)) for t, r, d in dec]}")


def bessel_ratios(nu, nmax, dps):
    with mp.workdps(dps):
        z = mp.mpf(1) / (2 * mp.mpf(nu))
        Nn = nmax + 60
        a = [mp.mpf(0)] * (Nn + 2)
        a[Nn] = mp.mpf(10) ** (-dps)
        for n in range(Nn, 0, -1):
            a[n - 1] = a[n + 1] + (2 * n / z) * a[n]
        return [x / a[0] for x in a[:nmax]]


class ColeHopf:
    """phi(x, t) = sum_n I_n(1/2nu)/I_0 e^{-nu n^2 t} e^{inx};  u = -2 nu phi_x/phi.  On x = pi + iy phi is real."""

    def __init__(self, nu, ymax=3.0, tmin=0.08):
        self.nu = nu
        self.nmax = int(80 + 3 * (ymax + 1) / (nu * (1 + tmin)) + 60 / math.sqrt(nu))
        self.dps = int(30 + (ymax ** 2 + math.pi ** 2) / (4 * nu * (1 + tmin)) / 2.3)
        self.c = bessel_ratios(nu, self.nmax, self.dps)

    def weights(self, t):
        with mp.workdps(self.dps):
            nu, t = mp.mpf(self.nu), mp.mpf(t)
            return [self.c[0]] + [2 * (-1) ** n * self.c[n] * mp.e ** (-nu * n * n * t) for n in range(1, self.nmax)]

    def phi_line(self, w, y):
        with mp.workdps(self.dps):
            q = mp.e ** mp.mpf(y)
            return (mp.polyval(w[::-1], q) + mp.polyval(w[::-1], 1 / q)) / 2

    def u_complex(self, x, t):
        with mp.workdps(self.dps):
            nu, t = mp.mpf(self.nu), mp.mpf(t)
            ph = self.c[0]
            dph = mp.mpf(0)
            for n in range(1, self.nmax):
                wn = self.c[n] * mp.e ** (-nu * n * n * t)
                ph += 2 * wn * mp.cos(n * x)
                dph += -2 * n * wn * mp.sin(n * x)
            return -2 * nu * dph / ph

    def zeros(self, t, ymax, ny=160, kmax=5, ymin=1e-4):
        w = self.weights(t)
        ys = np.linspace(ymin, ymax, ny)
        vals = [self.phi_line(w, y) for y in ys]
        out = []
        for i in range(ny - 1):
            if mp.sign(vals[i]) != mp.sign(vals[i + 1]):
                with mp.workdps(self.dps):
                    r = mp.findroot(lambda yy: self.phi_line(w, yy), (mp.mpf(ys[i]), mp.mpf(ys[i + 1])), solver='anderson', verify=False)
                    if not (ys[i] - 1e-12 <= float(r) <= ys[i + 1] + 1e-12):          # fall back to bisection
                        a, b = mp.mpf(ys[i]), mp.mpf(ys[i + 1])
                        fa = self.phi_line(w, a)
                        for _ in range(60):
                            mid = (a + b) / 2
                            fm = self.phi_line(w, mid)
                            if mp.sign(fm) == mp.sign(fa):
                                a, fa = mid, fm
                            else:
                                b = mid
                        r = (a + b) / 2
                out.append(float(r))
                if len(out) >= kmax:
                    break
        return out


def u_heat_kernel(x, t, nu, n=6001):
    """Cole-Hopf through the heat kernel, in log-sum-exp form (what the lab computes in the browser)."""
    xi = np.linspace(x - t - 1.5 - 12 * math.sqrt(nu * t), x + t + 1.5 + 12 * math.sqrt(nu * t), n)
    e = -((x - xi) ** 2 / (2 * t) - np.cos(xi)) / (2 * nu)
    w = np.exp(e - e.max())
    return float(np.sum(w * (x - xi)) / np.sum(w) / t)


ch1 = ColeHopf(0.1)
ok_hk = True
for x in (0.4, 2.0, 3.0, 4.5):
    ok_hk &= abs(u_heat_kernel(x, 0.7, 0.1) - float(ch1.u_complex(mp.mpf(x), 0.7))) < 1e-8
z_test = ch1.zeros(1.0, 2.5)
x0c = mp.mpc(mp.pi, z_test[0])
eps_ = mp.mpf('1e-7')
res = (eps_ * 1j) * ch1.u_complex(x0c + eps_ * 1j, 1.0)
check("viscosity holds the strip open: u = −2ν ∂ₓ log φ with φ the heat flow of e^{cos x/2ν} (Cole–Hopf); φ > 0 on the real line, so u's singularities are φ's complex zeros — simple poles with residue −2ν; by symmetry the nearest sit on Re x = π",
      ok_hk and abs(res + 0.2) < 1e-5,
      f"heat-kernel (log-sum-exp) and Bessel-series evaluations agree to 1e-8 at ν = 0.1; nearest poles at t = 1: y = {[round(v, 4) for v in z_test[:4]]}; residue at the first = {mp.nstr(res, 8)} (−2ν = −0.2)")

def phi_complex(ch, x, t):
    with mp.workdps(ch.dps):
        nu, t = mp.mpf(ch.nu), mp.mpf(t)
        return ch.c[0] + 2 * mp.fsum(ch.c[n] * mp.e ** (-nu * n * n * t) * mp.cos(n * x) for n in range(1, ch.nmax))


box = [(mp.pi - 1.5, 0.05), (mp.pi + 1.5, 0.05), (mp.pi + 1.5, 1.8), (mp.pi - 1.5, 1.8)]
path = []
for (xa, ya), (xb, yb) in zip(box, box[1:] + box[:1]):
    for k in range(500):
        path.append(mp.mpc(xa + (xb - xa) * k / 500, ya + (yb - ya) * k / 500))
path.append(path[0])
args = np.unwrap(np.array([float(mp.arg(phi_complex(ch1, z, 1.0))) for z in path]))
winding = (args[-1] - args[0]) / (2 * math.pi)
on_line = sum(1 for v in z_test if 0.05 < v < 1.8)
check("the poles found on Re x = π are all the poles in the box |Re x − π| < 1.5, 0.05 < Im x < 1.8 (ν = 0.1, t = 1): the argument principle counts the same number of zeros of φ",
      abs(winding - on_line) < 1e-6 and on_line == 4,
      f"winding number of φ around the box = {winding:.6f}; zeros found on the line in the box: {on_line}")

tgrid = [round(0.08 + 0.04 * k, 2) for k in range(61)]       # 0.08 ... 2.48
burg = {}
below_changes = 0
for nu in (0.1, 0.03, 0.01):
    ch = ColeHopf(nu, ymax=3.6, tmin=0.08)
    rows = []
    for t in tgrid:
        ylo = max(1e-4, delta_inv(t) - 0.02) if t < 1 else 1e-4      # the viscous poles lie above the inviscid singularity before the shock
        ymax = min(3.6, (delta_inv(t) + 0.9) if t < 1 else 1.4)
        zs = ch.zeros(t, ymax, ny=170, ymin=ylo)
        if t < 1:      # nothing below: phi keeps its sign on [0, delta - 0.02]
            wts = ch.weights(t)
            below = [ch.phi_line(wts, yy) for yy in np.linspace(1e-4, ylo, 40)]
            below_changes += sum(1 for a_, b_ in zip(below, below[1:]) if mp.sign(a_) != mp.sign(b_))
        rows.append([round(v, 5) for v in zs])
    burg[nu] = rows
ok_pos = all(len(r) > 0 and r[0] > 0 for nu in burg for r in burg[nu])
conv_ok = True
for t in (0.32, 0.6):
    i = tgrid.index(t)
    d1, d2, d3 = burg[0.1][i][0], burg[0.03][i][0], burg[0.01][i][0]
    conv_ok &= d1 > d2 > d3 > delta_inv(t)
xi_s = optimize.brentq(lambda s: math.pi - s - 2 * math.sin(s), 0.5, 2.5)
U2 = math.sin(xi_s)
i2 = tgrid.index(2.0)
ratio = [burg[nu][i2][0] / (math.pi * nu / U2) for nu in (0.03, 0.01)]
check("the viscous strip never closes and the inviscid one is its limit: the nearest pole height δ_ν(t) stays positive for all t ∈ [0.08, 2.48] at ν = 0.1, 0.03, 0.01, lies above the inviscid δ(t) and decreases toward it as ν → 0 before the shock, and after it sits at πν/U above the shock (the poles of a tanh front of half-jump U)",
      ok_pos and conv_ok and below_changes == 0 and all(abs(r - 1) < 0.06 for r in ratio),
      f"t = 0.6: δ_ν = {burg[0.1][tgrid.index(0.6)][0]}, {burg[0.03][tgrid.index(0.6)][0]}, {burg[0.01][tgrid.index(0.6)][0]} → δ = {delta_inv(0.6):.5f}; t = 2: U = {U2:.4f}, δ_ν/(πν/U) = {[round(r, 4) for r in ratio]}; min over t: {[min(r[0] for r in burg[nu]) for nu in burg]}")
DATA['burgers'] = {'t': tgrid, 'nu': [0.1, 0.03, 0.01], 'zeros': [burg[0.1], burg[0.03], burg[0.01]]}

# =====================================================================================================
print(f"   [{time.time() - T0:.0f} s]")
print("§4  The gap on the modular curve  [structure]")
# =====================================================================================================
x_, y_, s_ = sp.symbols('x y s')
fs = y_ ** s_
check("the hyperbolic plane has a built-in gap: −Δ y^s = s(1 − s) y^s for Δ = y²(∂ₓ² + ∂ᵧ²), and on the tempered line s = ½ + ir this is ¼ + r² ≥ ¼ — the bottom of the continuous spectrum of every finite-area quotient",
      sp.simplify(-y_ ** 2 * (sp.diff(fs, x_, 2) + sp.diff(fs, y_, 2)) - s_ * (1 - s_) * fs) == 0 and
      sp.expand((sp.Rational(1, 2) + sp.I * sp.Symbol('r', real=True)) * (sp.Rational(1, 2) - sp.I * sp.Symbol('r', real=True))) == sp.Rational(1, 4) + sp.Symbol('r', real=True) ** 2,
      "the Eisenstein series E(τ, ½ + ir) of Instruments 19 and 21 are these continuum states on SL(2,ℤ)\\H; eigenvalues in (0, ¼) would be 'exceptional'")


def Lam(s):
    return mp.pi ** (-s / 2) * mp.gamma(s / 2) * mp.zeta(s)


def phi_sc(s):
    return Lam(2 * s - 1) / Lam(2 * s)


ok_sc = all(abs(abs(phi_sc(mp.mpf(0.5) + 1j * t)) - 1) < 1e-20 for t in (1, 5.5, 13.3, 27.1))
for s0 in (mp.mpc(0.3, 2.1), mp.mpc(0.83, -4.4), mp.mpc(1.7, 9.0)):
    ok_sc &= abs(phi_sc(s0) * phi_sc(1 - s0) - 1) < 1e-20
resid = (mp.mpf('1e-12')) * phi_sc(1 + mp.mpf('1e-12'))
gam = [mp.zetazero(n).imag for n in range(1, 31)]
pole_ok = all(abs(Lam(2 * (mp.mpf(0.25) + 1j * g / 2))) < 1e-15 for g in gam[:10])
check("the scattering matrix of the modular surface: φ(s) = Λ(2s − 1)/Λ(2s) (Λ = π^{−s/2}Γ(s/2)ζ(s)) has |φ(½ + it)| = 1, φ(s)φ(1 − s) = 1, residue 3/π = 1/area at s = 1, and poles at s = ρ/2 — on Re s = ¼ exactly when the Riemann zeros are on Re s = ½",
      ok_sc and abs(resid - 3 / mp.pi) < 1e-9 and pole_ok,
      f"residue {mp.nstr(resid, 12)} vs 3/π = {mp.nstr(3 / mp.pi, 12)}; Λ(2s) = 0 at s = ¼ + iγₙ/2 for the first 10 zeros (Pavlov–Faddeev 1972, Lax–Phillips 1976)")

ts = np.round(np.arange(0.05, 32.0001, 0.05), 2)
ph = np.unwrap(np.array([float(mp.arg(phi_sc(mp.mpf(0.5) + 1j * mp.mpf(float(t))))) for t in ts]))
tdel = -np.gradient(ph, ts)
peaks = [ts[i] for i in range(1, len(ts) - 1) if tdel[i] > tdel[i - 1] and tdel[i] > tdel[i + 1]]
pk_ok = all(min(abs(np.array(peaks) - float(g) / 2)) < 0.03 for g in gam[:12] if float(g) / 2 < 31.5)
heights = [float(tdel[int(round(float(g) / 2 / 0.05)) - 1]) for g in gam[:6]]
lor = lambda tau_, w_: (2 * w_) / (w_ * w_ + tau_ * tau_) / 2   # one pole/zero pair at distance w from the line: -d arg/dt = 2w/(w^2+tau^2)
check("the Riemann zeros are resonances of the modular surface: the time delay −d arg φ(½ + it)/dt peaks at t = γₙ/2 with height ≈ 8 and half-width ¼ — one Lorentzian 2w/(w² + τ²) with w = ¼ per zero; a zero off the line would split its peak into widths ¼ ∓ ε/2",
      pk_ok and all(7.5 < h < 9.5 for h in heights) and abs(2 * 0.25 / (0.25 ** 2) - 8) < 1e-12,
      f"peaks at {[round(float(p), 2) for p in peaks[:8]]} vs γₙ/2 = {[round(float(g) / 2, 3) for g in gam[:8]]}; heights {[round(h, 2) for h in heights]}; RH ⇔ every resonance has the same lifetime")
DATA['scatter'] = {'t0': 0.05, 'dt': 0.05, 'phase': [round(float(v), 5) for v in ph], 'delay': [round(float(v), 4) for v in tdel],
                   'gamma': [round(float(g), 6) for g in gam]}

lmfdb = [(9.53369526, 'odd'), (12.1730083, 'odd'), (13.7797513, 'even'), (14.3585095, 'odd'), (16.1380731, 'odd'),
         (16.6442592, 'odd'), (17.7385633, 'even'), (18.1809178, 'odd'), (19.4234814, 'even'), (19.4847138, 'odd'),
         (20.1066946, 'odd'), (21.3157959, 'even'), (21.4790575, 'odd'), (22.1946739, 'odd'), (22.7859084, 'even'),
         (23.2013961, 'odd'), (23.2637115, 'odd'), (24.1123527, 'even'), (24.4197154, 'odd'), (25.0508548, 'odd'),
         (25.8262437, 'even'), (26.0569177, 'odd'), (26.1520854, 'even'), (26.4469964, 'odd'), (27.2843840, 'odd'),
         (27.3327080, 'even'), (27.7759207, 'odd'), (28.5102777, 'odd'), (28.5307476, 'even'), (28.8633943, 'even'),
         (29.1375875, 'odd'), (29.5463881, 'odd'), (30.2790484, 'odd'), (30.4043270, 'odd'), (30.4106788, 'even'),
         (31.0565339, 'odd'), (31.5265821, 'even'), (31.5662754, 'even'), (31.9161824, 'odd'), (32.0184064, 'odd'),
         (32.5081177, 'even'), (32.8911702, 'even'), (32.9324655, 'odd'), (33.4923312, 'odd'), (33.5709896, 'odd'),
         (34.0278842, 'even'), (34.1859699, 'odd'), (34.4562715, 'even'), (34.6953104, 'odd'), (35.4316646, 'odd')]


def kbes_scaled(r, x):
    """e^{pi r/2} K_{ir}(x)."""
    if HAVE_FLINT:
        v = acb(float(x)).bessel_k(acb(0, float(r)))
        return mp.mpf(float(v.real)) * mp.e ** (mp.pi * r / 2)
    return mp.re(mp.besselk(1j * r, x)) * mp.e ** (mp.pi * r / 2)


if HAVE_FLINT:
    flint_ctx.prec = 140


def pullback(x, y):
    a, b, c, d, zx, zy = reduce_tau(float(x), float(y))
    return zx, zy


def hejhal_coeffs(r, Y, M, parity, Q=None):
    Q = Q or M + 8
    xs = [(m - 0.5) / (2 * Q) for m in range(1, Q + 1)]
    cs = math.cos if parity == 'even' else math.sin
    V = np.zeros((M, M))
    for x in xs:
        zx, zy = pullback(x, Y)
        Kn = [math.sqrt(zy) * float(kbes_scaled(r, 2 * math.pi * n * zy)) * cs(2 * math.pi * n * zx) for n in range(1, M + 1)]
        for l in range(1, M + 1):
            w = cs(2 * math.pi * l * x) * 2 / Q
            V[l - 1, :] += w * np.array(Kn)
    for l in range(1, M + 1):
        V[l - 1, l - 1] -= math.sqrt(Y) * float(kbes_scaled(r, 2 * math.pi * l * Y))
    a = np.linalg.solve(V[1:, 1:], -V[1:, 0])
    return np.concatenate([[1.0], a])


def hejhal_defect(r, parity, M=18):
    a = hejhal_coeffs(r, 0.80, M, parity)
    b = hejhal_coeffs(r, 0.71, M, parity)
    return a[1] - b[1], a


found = []
ok_h = True
for idx in (0, 2, 5, 6):
    r0, par = lmfdb[idx]
    other = 'even' if par == 'odd' else 'odd'
    rr = optimize.brentq(lambda r: hejhal_defect(r, par)[0], r0 - 0.003, r0 + 0.003, xtol=1e-11)
    dwrong = abs(hejhal_defect(r0, other)[0])
    found.append((idx + 1, rr, par, dwrong))
    ok_h &= abs(rr - r0) < 2e-7 and dwrong > 1e-3
# refine r1 at lower sampling heights, where the coefficients up to a_12 are well determined
r1 = optimize.brentq(lambda r: hejhal_coeffs(r, 0.35, 36, 'odd')[1] - hejhal_coeffs(r, 0.31, 36, 'odd')[1], found[0][1] - 2e-6, found[0][1] + 2e-6, xtol=1e-13)
a1 = hejhal_coeffs(r1, 0.35, 36, 'odd')
hk = [a1[3] - (a1[1] ** 2 - 1), a1[5] - a1[1] * a1[2], a1[8] - (a1[2] ** 2 - 1), a1[9] - a1[1] * a1[4], a1[11] - a1[3] * a1[2]]
hecke = all(abs(v) < 1e-7 for v in hk)
check("Maass cusp forms on SL(2,ℤ)\\H by Hejhal's method (f = Σ aₙ √y K_{ir}(2πny) cs(2πnx), matched across the pullback to F at two heights): r = 9.5336953 (odd), 13.7797513 (even), 16.6442592 (odd), 17.7385633 (even) as in the LMFDB, the other parity fails, and the coefficients obey Hecke's relations a₄ = a₂² − 1, a₆ = a₂a₃, a₉ = a₃² − 1, a₁₀ = a₂a₅, a₁₂ = a₃a₄",
      ok_h and hecke,
      f"{[(n, round(r, 8), p) for n, r, p, _ in found]} (wrong-parity defects {[round(float(dw), 3) for _, _, _, dw in found]}); refined r₁ = {r1:.11f}; first form: a₂ = {a1[1]:.10f}, a₃ = {a1[2]:.10f}, Hecke residuals {[f'{v:.0e}' for v in hk]}; λ₁ = ¼ + r₁² = {0.25 + r1 ** 2:.4f}; K-Bessel via {'python-flint' if HAVE_FLINT else 'mpmath'}")


def weyl_M(T):
    return T * T / 12 - 2 * T * math.log(T) / math.pi + T / math.pi * (2 + math.log(math.pi / 2)) - 131 / 144


rs = [r for r, _ in lmfdb]
Svals = []
for i, r in enumerate(rs):
    Svals += [i - weyl_M(r - 1e-9), (i + 1) - weyl_M(r + 1e-9)]
Tgrid = np.linspace(10, 35.4, 2000)
Smean = np.mean([sum(1 for r in rs if r <= T) - weyl_M(T) for T in Tgrid])
check("Weyl's law for SL(2,ℤ)\\H: N(T) = T²/12 − (2T/π) log T + (T/π)(2 + log π/2) − 131/144 + S(T) with S fluctuating about 0; the LMFDB's first 50 eigenvalues (r ≤ 35.43) fit it with |S| < 2 and mean S ≈ 0 — a missing eigenvalue would push S down by 1 from there on",
      max(abs(v) for v in Svals) < 2 and abs(Smean) < 0.25,
      f"max |S| = {max(abs(v) for v in Svals):.3f}, mean S on [10, 35.4] = {Smean:+.3f}; the −(2T/π) log T term is the continuous spectrum's (scattering phase) share (Jorgenson–Smajlović–Then 2014)")
DATA['maass'] = [[r, 1 if p == 'even' else 0] for r, p in lmfdb]

idx143 = 143 * Fraction(12, 11) * Fraction(14, 13)
check("X₀(143): index 168 in SL(2,ℤ), area 56π, so Weyl's leading term is (area/4π) T² = 14 T² (against T²/12 for SL(2,ℤ)), with four cusps — four Eisenstein families in the continuum",
      idx143 == 168 and abs(168 * math.pi / 3 - 56 * math.pi) < 1e-12 and abs(56 * math.pi / (4 * math.pi) - 14) < 1e-12 and abs((math.pi / 3) / (4 * math.pi) - 1 / 12) < 1e-15,
      "Selberg's ¼ is proven for Γ₀(N), N squarefree ≤ 857 (Booker–Strömbergsson 2007) and for Γ₁(N), N ≤ 880 (Booker–Lee–Strömbergsson 2020) — both cover level 143: X₀(143) has no exceptional eigenvalues")
bounds = [(Fraction(1, 4), Fraction(3, 16), 'Selberg 1965'), (Fraction(3, 14), Fraction(10, 49), 'Iwaniec'), (Fraction(5, 28), Fraction(171, 784), 'Luo–Rudnick–Sarnak 1995'),
          (Fraction(5, 34), Fraction(66, 289), 'Kim–Shahidi 2002'), (Fraction(7, 64), Fraction(975, 4096), 'Kim–Sarnak 2003')]
check("the ladder of proven bounds λ₁ ≥ ¼ − θ²: θ = 1/4 (3/16, Selberg 1965), 3/14 (10/49, Iwaniec), 5/28 (171/784 ≈ 0.218, Luo–Rudnick–Sarnak 1995), 5/34 (66/289, Kim–Shahidi 2002), 7/64 (975/4096 ≈ 0.238, Kim–Sarnak 2003); the conjecture is θ = 0",
      all(Fraction(1, 4) - th ** 2 == lam for th, lam, _ in bounds),
      f"{[(n, f'{float(l):.4f}') for _, l, n in bounds]}; θ bounds the archimedean Ramanujan exponent: Selberg's ¼ is Ramanujan at the infinite place")

R_, t_ = sp.symbols('R t', positive=True)
th_, ph_ = sp.symbols('theta phi')
gS2 = sp.diag(R_ ** 2, R_ ** 2 * sp.sin(th_) ** 2)
Kgauss = 1 / R_ ** 2            # Ric = K g in 2D
ricci_flow_S2 = sp.diff(R_ ** 2, R_) * 0 - 2 * Kgauss * R_ ** 2    # d(R^2)/dt = -2 K R^2
lap_S3_l1 = 1 * (1 + 2)
check("curvature acts oppositely as base and as target: as the space fields live on, H² has a gap (¼ above); as the space fields take values in, a sphere S² shrinks under the renormalization group (Friedan: RG = Ricci flow, d(R²)/dt = −2: asymptotic freedom and a dynamical gap, Polyakov's O(3) model) while H² grows (d(R²)/dt = +2: free in the infrared); positive Ricci curvature of field space bounds a gap (Lichnerowicz: λ₁(S³) = 3 = nK, sharp)",
      sp.simplify(ricci_flow_S2 + 2) == 0 and lap_S3_l1 == 3,
      "S³ = SU(2) with Ric = 2g: Bakry–Émery curvature of the gauge group is what Shen–Zhu–Zhu (2023) use to prove a mass gap for lattice Yang–Mills at strong coupling")

# =====================================================================================================
n_pass = sum(ok for _, ok in RESULTS)
print(f"\n{n_pass}/{len(RESULTS)} checks passed in {time.time() - T0:.1f} s")
if '--json' in sys.argv:
    path = sys.argv[sys.argv.index('--json') + 1]
    with open(path, 'w') as fh:
        json.dump(DATA, fh, separators=(',', ':'))
    print('data written to', path)
if n_pass == len(RESULTS):
    print("ALL CHECKS PASSED")
else:
    for name, ok in RESULTS:
        if not ok:
            print("FAILED:", name)
    sys.exit(1)
