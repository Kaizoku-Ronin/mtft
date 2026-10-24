"""Gravity's knot: the dimensional formula of G is Kepler's cusp, whose link is the trefoil; the Kepler problem is the Hopf map's shadow of a
4D oscillator; orbital resonances are torus knots; in 2+1 gravity G is the knot parameter (Chern-Simons level) and the modular knots are
BTZ horizons; MTFT's own G (Paper 31) is the modulus of the Dedekind eta whose phase is the knot's linking number; fractal knots by
substitution are cables.  2026-10-21.  EXACT unless marked (Paper 31's formula is the programme's model, reproduced, not endorsed)."""
import math, cmath, itertools, bisect
import numpy as np
import sympy as sp
import mpmath as mp
mp.mp.dps = 30
rng = np.random.default_rng(11)
t = sp.symbols('t')

# ---------------------------------------------------------------- 0. the Alexander polynomial of a polygonal link from a projection (as in the earlier scripts)
def alexander(components, tries=6):
    for _ in range(tries):
        Q, _ = np.linalg.qr(rng.normal(size=(3, 3))); comps = [np.asarray(c)[:, :3] @ Q.T for c in components]; cross = []; ok = True
        for ci, P in enumerate(comps):
            for cj, Qp in enumerate(comps):
                if cj < ci: continue
                N, M = len(P), len(Qp); p = P[:, :2]; rr = np.roll(P, -1, 0)[:, :2] - p; q = Qp[:, :2]; s = np.roll(Qp, -1, 0)[:, :2] - q
                dp = P[:, 2]; dq = Qp[:, 2]; dpn = np.roll(dp, -1); dqn = np.roll(dq, -1)
                for i in range(N):
                    js = np.arange(M)
                    if ci == cj: js = js[(js > i + 1) & ~((i == 0) & (js == N - 1))]
                    if len(js) == 0: continue
                    den = rr[i, 0] * s[js, 1] - rr[i, 1] * s[js, 0]; good = np.abs(den) > 1e-14; qp = q[js] - p[i]
                    ti = np.where(good, (qp[:, 0] * s[js, 1] - qp[:, 1] * s[js, 0]) / np.where(good, den, 1), -1); tj = np.where(good, (qp[:, 0] * rr[i, 1] - qp[:, 1] * rr[i, 0]) / np.where(good, den, 1), -1)
                    hit = good & (ti >= 0) & (ti < 1) & (tj >= 0) & (tj < 1)
                    for k in np.nonzero(hit)[0]:
                        j = js[k]
                        if min(ti[k], 1 - ti[k], tj[k], 1 - tj[k]) < 1e-6: ok = False
                        di = dp[i] + ti[k] * (dpn[i] - dp[i]); dj = dq[j] + tj[k] * (dqn[j] - dq[j]); over_i = di > dj
                        d_over, d_under = (rr[i], s[j]) if over_i else (s[j], rr[i]); sign = np.sign(d_over[0] * d_under[1] - d_over[1] * d_under[0])
                        cross.append((ci, i + ti[k], cj, j + tj[k], over_i, sign))
        if not ok: continue
        unders = [[] for _ in comps]
        for ci, si, cj, sj, over_i, sign in cross:
            if over_i: unders[cj].append(sj)
            else: unders[ci].append(si)
        for uu in unders: uu.sort()
        offset = []; narcs = 0
        for uu in unders: offset.append(narcs); narcs += max(1, len(uu))
        def arc_of(c, s_):
            uu = unders[c]
            if not uu: return offset[c]
            return offset[c] + (bisect.bisect_right(uu, s_) - 1) % len(uu)
        n = len(cross)
        if n == 0: return sp.Integer(1), 0
        A = sp.zeros(n, narcs)
        for row, (ci, si, cj, sj, over_i, sign) in enumerate(cross):
            if over_i: co, so, cu, su = ci, si, cj, sj
            else: co, so, cu, su = cj, sj, ci, si
            kk = arc_of(co, so); uu = unders[cu]; idx = bisect.bisect_left(uu, su); out_arc = offset[cu] + idx; in_arc = offset[cu] + (idx - 1) % len(uu)
            A[row, kk] += 1 - t
            if sign > 0: A[row, in_arc] += t; A[row, out_arc] += -1
            else: A[row, in_arc] += -1; A[row, out_arc] += t
        D = sp.expand(A[:-1, :-1].det(method='berkowitz'))
        if D == 0: return sp.Integer(0), n
        P_ = sp.Poly(D, t); low = min(mon[0] for mon in P_.monoms()); D = sp.expand(sp.cancel(D / t ** low))
        if sp.Poly(D, t).LC() < 0: D = -D
        return D, n
    raise RuntimeError('no generic projection')
def same_poly(D, E): return sp.expand(D - E) == 0 or sp.expand(D - sp.expand(sp.cancel(E.subs(t, 1 / t) * t ** sp.degree(E, t)))) == 0
torus_poly = lambda p, q: sp.expand(sp.cancel((t ** (p * q) - 1) * (t - 1) / ((t ** p - 1) * (t ** q - 1)))) if math.gcd(p, q) == 1 else None

# ---------------------------------------------------------------- 1. the dimensional formula of G is Kepler's cusp, and the cusp's link is the trefoil
M_, L_, T_, G_, Msun, r_, v_, d_ = sp.symbols('M L T G M_sun r v d', positive=True)
# F = G M m / r^{d-2} (Gauss's law in d-1 spatial dimensions): [G] = [F] L^{d-2} / M^2 = M^{-1} L^{d-1} T^{-2}
for d in (3, 4, 5, 6, 7, 8):
    dimG = {'M': -1, 'L': d - 1, 'T': -2}
    assert dimG == {'M': 1 - 2, 'L': 1 + (d - 2), 'T': -2}                                       # [F] = M L T^-2, times L^{d-2}, over M^2
    # circular orbit: v^2 / r = G M / r^{d-2}  =>  T^2 = (2 pi r / v)^2 = 4 pi^2 r^{d-1} / (G M)
    v2 = G_ * Msun / r_ ** (d - 3); T2 = sp.simplify((2 * sp.pi * r_) ** 2 / v2)
    assert sp.simplify(T2 - 4 * sp.pi ** 2 * r_ ** (d - 1) / (G_ * Msun)) == 0
print("[G] = M^-1 L^{d-1} T^-2 in d spacetime dimensions (Gauss's law), and Kepler's law for circular orbits reads T^2 = 4 pi^2 r^{d-1}/(GM): the exponents of L and T in [G] are the exponents of Kepler's cusp r^{d-1} = c T^2")
for d in (3, 4, 5, 6, 7, 8):
    g = math.gcd(2, d - 1); name = 'the Hopf link' if d == 3 else 'the trefoil' if d == 4 else ('the torus link T(2,%d), %d components' % (d - 1, g) if g > 1 else 'the torus knot T(2,%d)' % (d - 1))
    print("   d = %d: [G] = M^-1 L^%d T^-2, cusp y^2 = x^%d, link %s%s" % (d, d - 1, d - 1, name, ', Alexander ' + str(torus_poly(2, d - 1)) if g == 1 and d > 3 else ''))
# the cusp y^2 = x^3 meets the sphere |x|^2 + |y|^2 = 1 in the curve (x, y) = (r^2 e^{2i th}, r^3 e^{3i th}), r^4 + r^6 = 1: the trefoil
r0 = sp.nsolve(sp.Symbol('r') ** 4 + sp.Symbol('r') ** 6 - 1, sp.Symbol('r'), 0.8); r0 = float(r0)
cusp = []
for th in np.linspace(0, 2 * math.pi, 600, endpoint=False):
    x = r0 ** 2 * cmath.exp(2j * th); y = r0 ** 3 * cmath.exp(3j * th); assert abs(y * y - x ** 3) < 1e-12 and abs(abs(x) ** 2 + abs(y) ** 2 - 1) < 1e-12
    v = np.array([x.real, x.imag, y.real, y.imag]); cusp.append(v[:3] / (1 - v[3]))                 # stereographic projection from (0,0,0,1)
D, nc = alexander([cusp]); assert same_poly(D, t ** 2 - t + 1), D
print("   the link of y^2 = x^3 (its intersection with the unit 3-sphere, projected): Alexander polynomial %s = the trefoil; the value of G rescales the cusp a^3 = (GM/4pi^2) T^2 without changing its knot" % D)
print("   Bertrand 1873: among power-law forces only Newton's 1/r^2 (d = 4) and Hooke's r have all bounded orbits closed; d = 3 is Chern-Simons gravity, where the observable is linking (the Hopf link is the (2,2) cusp)")
# the planets: a^3 / T^2 = 1 AU^3 / yr^2, slope 3/2, and GM_sun = 4 pi^2 AU^3/yr^2
planets = {'Mercury': (0.38709927, 0.240846), 'Venus': (0.72333566, 0.615198), 'Earth': (1.00000261, 1.0000174), 'Mars': (1.52371034, 1.880848), 'Jupiter': (5.20288700, 11.8618), 'Saturn': (9.53667594, 29.4571), 'Uranus': (19.18916464, 84.0205), 'Neptune': (30.06992276, 164.79), 'Pluto': (39.48211675, 247.94)}   # Standish's mean elements (a in AU), sidereal periods (yr)
ratios = {n: a ** 3 / T ** 2 for n, (a, T) in planets.items()}; assert max(abs(x - 1) for x in ratios.values()) < 2.5e-3, ratios
assert abs(ratios['Jupiter'] - 1 - 1 / 1047.35) < 4e-4                                              # a^3/T^2 = 1 + m/M_sun: Jupiter's mass shows
la = np.log([a for a, T in planets.values()]); lT = np.log([T for a, T in planets.values()]); slope = np.polyfit(la, lT, 1)[0]; assert abs(slope - 1.5) < 1e-3
AU, yr = 1.495978707e11, 365.25 * 86400; GM = 4 * math.pi ** 2 * AU ** 3 / yr ** 2; assert abs(GM / 1.32712440018e20 - 1) < 2e-4
print("   the planets: a^3/T^2 = 1 AU^3/yr^2 to %.1e (Mercury to Pluto; Jupiter's excess %.5f is its mass, 1/1047), log-log slope %.4f; GM_sun = 4 pi^2 AU^3/yr^2 = %.4e m^3/s^2 (IAU 1.3271244e20); the Gaussian constant k = sqrt(GM) in AU, day, M_sun units = %.11f (Gauss: 0.01720209895)" % (max(abs(x - 1) for x in ratios.values()), ratios['Jupiter'] - 1, slope, GM, math.sqrt(GM / AU ** 3 * 86400 ** 2)))

# ---------------------------------------------------------------- 2. Kustaanheimo-Stiefel: the Kepler problem is the Hopf map's shadow of a 4D harmonic oscillator
def KS(u):                                                                                        # R^4 -> R^3, the Hopf map in KS's coordinates
    u1, u2, u3, u4 = u; return np.array([u1 * u1 - u2 * u2 - u3 * u3 + u4 * u4, 2 * (u1 * u2 - u3 * u4), 2 * (u1 * u3 + u2 * u4)])
for _ in range(50):
    u = rng.normal(size=4); x = KS(u); assert abs(np.linalg.norm(x) - np.dot(u, u)) < 1e-12
    ph = rng.uniform(0, 2 * math.pi); u1, u2, u3, u4 = u                                            # the fibre: a circle in R^4 mapping to the same point
    uf = np.array([u1 * math.cos(ph) - u4 * math.sin(ph), u2 * math.cos(ph) + u3 * math.sin(ph), u3 * math.cos(ph) - u2 * math.sin(ph), u4 * math.cos(ph) + u1 * math.sin(ph)])
    assert np.allclose(KS(uf), x) and abs(np.linalg.norm(uf) - np.linalg.norm(u)) < 1e-12
# a 4D oscillator orbit u'' = -w^2 u obeying the bilinear constraint u4 u1' - u1 u4' + u2 u3' - u3 u2' = 0, mapped down, is a Kepler ellipse with the focus at the origin
w = 0.7; Acoef = np.array([1.0, 0.0, 0.0, 0.0]); Bcoef = np.array([0.0, 0.6, 0.0, 0.0])             # u = A cos ws + B sin ws, A.B = 0 -> constraint holds
s_grid = np.linspace(0, 2 * math.pi / w, 4001); U = np.outer(np.cos(w * s_grid), Acoef) + np.outer(np.sin(w * s_grid), Bcoef); Up = w * (np.outer(-np.sin(w * s_grid), Acoef) + np.outer(np.cos(w * s_grid), Bcoef))
bil = U[:, 3] * Up[:, 0] - U[:, 0] * Up[:, 3] + U[:, 1] * Up[:, 2] - U[:, 2] * Up[:, 1]; assert np.max(np.abs(bil)) < 1e-12
X = np.array([KS(u) for u in U]); rad = np.linalg.norm(X, axis=1); tt = np.concatenate([[0], np.cumsum(0.5 * (rad[1:] + rad[:-1]) * np.diff(s_grid))])   # dt = |u|^2 ds = r ds
# Kepler: x'' = -mu x / r^3 with mu = 2|u'|^2 + 2 w^2 |u|^2 (Stiefel-Scheifele), energy h = -2 w^2
mu = 2 * np.dot(Up[0], Up[0]) + 2 * w * w * np.dot(U[0], U[0])
Xt = np.gradient(X, tt, axis=0); Xtt = np.gradient(Xt, tt, axis=0); mid = slice(200, -200)
resid = np.linalg.norm(Xtt[mid] + mu * X[mid] / rad[mid, None] ** 3, axis=1) / np.linalg.norm(mu * X[mid] / rad[mid, None] ** 3, axis=1)
assert np.max(resid) < 2e-3, np.max(resid)
h = 0.5 * np.einsum('ij,ij->i', Xt[mid], Xt[mid]) - mu / rad[mid]; assert np.std(h) < 1e-3 * abs(np.mean(h)) and abs(np.mean(h) + 2 * w * w) < 2e-3 * 2 * w * w
Lvec = np.cross(X[mid], Xt[mid]); ecc = np.cross(Xt[mid], Lvec) / mu - X[mid] / rad[mid, None]; e_ = np.linalg.norm(np.mean(ecc, axis=0)); assert np.max(np.linalg.norm(ecc - np.mean(ecc, axis=0), axis=1)) < 2e-3
print("Kustaanheimo-Stiefel: the map u -> (u1^2 - u2^2 - u3^2 + u4^2, 2(u1 u2 - u3 u4), 2(u1 u3 + u2 u4)) is the Hopf map (|x| = |u|^2, circle fibres); a 4D oscillator orbit u'' = -w^2 u with dt = |u|^2 ds maps to a Kepler orbit: x'' = -mu x/r^3 to %.1e, energy -2w^2 conserved, eccentricity vector constant, e = %.3f" % (np.max(resid), e_))
print("   so Newton's gravity, regularised, lives one dimension up as a harmonic oscillator with the Hopf fibre as its gauge (Instrument 11); the bound Kepler problem's hidden SO(4) is Instrument 17's rotation group, and hydrogen's n^2 degeneracy is its irreducible representations")

# ---------------------------------------------------------------- 3. resonances are torus knots
res = {'Neptune : Pluto': (164.8, 247.94, (3, 2)), 'Io : Europa': (1.769138, 3.551181, (2, 1)), 'Europa : Ganymede': (3.551181, 7.154553, (2, 1)), 'Hilda asteroids : Jupiter': (7.91, 11.862, (3, 2))}
for name, (T1, T2, (p, q)) in res.items():
    ratio = T2 / T1; assert abs(ratio - p / q) < 0.02 * (p / q), (name, ratio)
print("resonances: a pair of periods in the ratio p:q traces a closed curve on the torus of the two mean longitudes, the torus knot T(p,q): " + "; ".join("%s = %.4f ~ %d:%d" % (n, v[1] / v[0], *v[2]) for n, v in res.items()))
tor = [((2 + 0.6 * math.cos(3 * s)) * math.cos(2 * s), (2 + 0.6 * math.cos(3 * s)) * math.sin(2 * s), 0.6 * math.sin(3 * s)) for s in np.linspace(0, 2 * math.pi, 600, endpoint=False)]
D, nc = alexander([tor]); assert same_poly(D, t ** 2 - t + 1)
print("   Neptune-Pluto's 3:2 on the longitude torus: Alexander polynomial %s, the trefoil -- gravity's period clock; the Kirkwood gaps sit at 4:1, 3:1, 5:2, 7:3, 2:1: the unknot, the unknot, T(2,5), T(3,7), the unknot" % D)

# ---------------------------------------------------------------- 4. 2+1 gravity: G is the knot parameter, and the modular knots are BTZ horizons
ell, Gn, k_, rp, rm = sp.symbols('ell G k r_+ r_-', positive=True)
cBH = 3 * ell / (2 * Gn); kCS = ell / (4 * Gn); assert sp.simplify(cBH - 6 * kCS) == 0
print("2+1 gravity with Lambda = -1/ell^2 is Chern-Simons theory at level k = ell/(4G) (Achucarro-Townsend 1986, Witten 1988), Brown-Henneaux central charge c = 3 ell/(2G) = 6k; Wilson loops are knot invariants at q = e^{2 pi i/(k+2)} (Witten 1989): Newton's constant is the quantum parameter of the Jones polynomial")
# BTZ with J = 0: M = r+^2/(8 G ell^2), T = r+/(2 pi ell^2), S = 2 pi r+/(4G); first law dM = T dS
Mb = rp ** 2 / (8 * Gn * ell ** 2); Tb = rp / (2 * sp.pi * ell ** 2); Sb = 2 * sp.pi * rp / (4 * Gn)
assert sp.simplify(sp.diff(Mb, rp) - Tb * sp.diff(Sb, rp)) == 0
# Euclidean BTZ = H^3 / <dilation by e^{2 pi r+/ell}> (Carlip-Teitelboim 1995); a hyperbolic A in SL(2,R) acts on H^3 with translation length 2 arccosh(tr A/2)
def h3_dist(p, q):                                                                                # points (z, h) of the upper half-space, z complex, h > 0
    return math.acosh(1 + (abs(p[0] - q[0]) ** 2 + (p[1] - q[1]) ** 2) / (2 * p[1] * q[1]))
def act_h3(A, p):                                                                                 # PSL(2,C) on the upper half-space (quaternion form)
    (a, b), (c, d) = A; z, h = p; den = abs(c * z + d) ** 2 + abs(c) ** 2 * h ** 2
    return ((a * z + b) * np.conj(c * z + d) + a * np.conj(c) * h ** 2) / den, h / den
cases = {'LR': [[1, 1], [1, 2]], 'L²R': [[1, 1], [2, 3]], 'L³R': [[1, 1], [3, 4]], '[84 -47; 143 -80]': [[84, -47], [143, -80]], '[12 1; 143 12]': [[12, 1], [143, 12]]}
print("   Euclidean BTZ is H^3 modulo a hyperbolic element; A in SL(2,R) translates its axis by ell_A = 2 arccosh(tr A/2), so the horizon has circumference 2 pi r+ = ell * ell_A: r+/ell = ell_A/2pi, G M = ell_A^2/(32 pi^2), G S/ell = ell_A/4, ell T = ell_A/(4 pi^2)")
for name, A in cases.items():
    (a, b), (c, d) = A; tr = a + d; lA = 2 * math.acosh(tr / 2)
    # a point on the axis: the fixed points of A on the real line are the roots of c x^2 + (d - a) x - b = 0; the axis is the semicircle over them; take its top point
    x1, x2 = sorted(np.roots([c, d - a, -b]).real); p = (complex((x1 + x2) / 2), (x2 - x1) / 2)
    q = act_h3(A, p); disp = h3_dist(p, q); assert abs(disp - lA) < 1e-9, (name, disp, lA)
    print("   A = %-18s tr %2d: ell_A = %.4f (translation on H^3 checked %.4f); r+/ell = %.4f, G M = %.5f, G S/ell = %.4f, ell T = %.5f" % (name, tr, lA, disp, lA / (2 * math.pi), lA ** 2 / (32 * math.pi ** 2), lA / 4, lA / (4 * math.pi ** 2)))
print("   and SL(2,Z)\\SL(2,R) with the Killing metric is a locally AdS_3 vacuum solution on S^3 minus the trefoil (every quotient of AdS_3 = SL(2,R) by a discrete isometry group acting freely is one), with the Seifert fibres as its closed timelike curves")
# the Jones polynomial of the trefoil at the gravitational q, and the volume conjecture for the figure-eight knot (a theorem for 4_1)
V_tref = lambda q: q + q ** 3 - q ** 4
assert abs(V_tref(1) - 1) < 1e-15
for kk in (1, 2, 3, 6, 12):
    q = cmath.exp(2j * math.pi / (kk + 2)); print("   k = %2d (c = %3d): q = e^{2 pi i/%d}, V_trefoil(q) = %.4f %+.4fi, |V| = %.4f" % (kk, 6 * kk, kk + 2, V_tref(q).real, V_tref(q).imag, abs(V_tref(q))))
def J_fig8(N):                                                                                    # Kashaev's invariant <4_1>_N = sum_k prod_{j<=k} |1 - e^{2 pi i j/N}|^2
    total = mp.mpf(1); prod = mp.mpf(1)
    for k in range(1, N):
        prod *= (2 * mp.sin(mp.pi * k / N)) ** 2; total += prod
    return total
VOL41 = 2.0298832128193072
seq = []
for N in (50, 100, 200, 400, 800, 1600):
    val = float(2 * mp.pi / N * mp.log(J_fig8(N))); corr = val - 3 * math.pi * math.log(N) / N; seq.append((N, val, corr))
assert all(seq[i][2] < seq[i + 1][2] for i in range(len(seq) - 1)) and abs(seq[-1][2] - VOL41) < 0.03 and abs(seq[-1][1] - VOL41) < 0.06, seq
print("   volume conjecture (a theorem for 4_1): (2 pi/N) log|J_N(4_1; e^{2 pi i/N})| -> Vol(S^3 - 4_1) = %.6f: " % VOL41 + ", ".join("N = %d: %.4f (%.4f with the 3/2 log N term removed)" % s for s in seq))
print("   the Euclidean gravity action of a hyperbolic knot complement is (k/pi) Vol = ell Vol/(4 pi G): the Jones polynomial's growth rate IS the gravitational action; torus knots (the trefoil) have Vol = 0 -- Seifert fibred, the AdS_3 quotient instead")

# ---------------------------------------------------------------- 5. MTFT's G (Paper 31): the modulus of the eta function whose phase is the knot's linking number  [MODEL, reproduced]
hbar, c_light, m_e, alpha_inv = 1.054571817e-34, 2.99792458e8, 9.1093837015e-31, 137.035999084
G_meas = 6.67430e-11
def eta_mod(y, N=400):                                                                            # |eta(i y)| = e^{-pi y/12} prod (1 - e^{-2 pi n y})
    return math.exp(-math.pi * y / 12) * math.prod(1 - math.exp(-2 * math.pi * n * y) for n in range(1, N + 1))
assert abs(eta_mod(1.0) - float(abs(mp.eta(1j) if hasattr(mp, 'eta') else mp.mpf(0.7682254223260567)))) < 1e-12 or abs(eta_mod(1.0) - 0.7682254223260567) < 1e-12
assert abs(eta_mod(0.4) - eta_mod(1 / 0.4) * math.sqrt(1 / 0.4)) < 1e-12                          # eta(i/y) = sqrt(y) eta(iy)
kexp = 2 * alpha_inv / math.pi
MM = 0.2614972128476428; yc_MM = MM * math.log(2); yc_lab = 0.18174
for label, yc in (('y_c = M ln 2 (Paper 31)', yc_MM), ('y = 0.18174 (the lab\'s confinement depth)', yc_lab)):
    Gm = (hbar * c_light / m_e ** 2) * eta_mod(yc) ** (2 * kexp)
    print("   Paper 31: G = (hbar c/m_e^2) |eta(i y)|^{2k}, k = 2/(pi alpha) = %.3f, %s = %.5f: |eta| = %.5f, |eta|^{2k} = %.3e, G = %.4e (measured %.4e, %+.2f%%)" % (kexp, label, yc, eta_mod(yc), eta_mod(yc) ** (2 * kexp), Gm, G_meas, 100 * (Gm / G_meas - 1)))
# sensitivity: d ln G / dy = 2k d ln|eta(iy)|/dy = -(pi k/6) E_2(iy)
def E2_iy(y, N=400): return 1 - 24 * sum(sum(d for d in range(1, n + 1) if n % d == 0) * math.exp(-2 * math.pi * n * y) for n in range(1, N + 1))
sens = -(math.pi * kexp / 6) * E2_iy(yc_MM); num = (math.log(eta_mod(yc_MM + 1e-6) ** (2 * kexp)) - math.log(eta_mod(yc_MM - 1e-6) ** (2 * kexp))) / 2e-6
assert abs(sens - num) < 1e-4 * abs(sens)
print("   sensitivity d ln G/dy = -(pi k/6) E_2(i y_c) = %.1f: a shift of 1e-5 in y_c moves G by %.2f%%; the 10^45 weakness is |eta|^{%.0f} -- the modulus of the function whose phase round a modular knot is the linking number (Instrument 17)" % (sens, 100 * (math.exp(sens * 1e-5) - 1), 2 * kexp))
print("   status: the bridge m_e/M_P = |eta(i y_c)|^k is the programme's model (Paper 22); given it, G follows by definition of the Planck mass -- reproduced here, not derived")

# ---------------------------------------------------------------- 6. fractal knots by substitution are cables
# the (p,q)-cable of K (p strands winding round K): Delta(t) = Delta_K(t^p) Delta_{T(p,q)}(t); it is the link of a plane-curve singularity (two Puiseux pairs) iff q > p * (crossings' pq) : for the trefoil and p = 2, q > 12
cable = lambda DK, p, q: sp.expand(DK.subs(t, t ** p) * torus_poly(p, q))
D13 = cable(t ** 2 - t + 1, 2, 13); assert sp.degree(D13, t) == 4 + 12 and D13.subs(t, 1) == 1
D15 = cable(D13, 2, 55); assert sp.degree(D15, t) == 2 * 16 + 54
print("fractal knots: Fathauer's iteration substitutes copies of the knot into itself (wild in the limit); the algebraic version is iterated cabling: the (2,q)-cable of the trefoil has Alexander polynomial Phi_6(t^2) Delta_{T(2,q)}(t), degree %d for q = 13, and is the link of a plane-curve singularity (Puiseux pairs (2,3),(2,13)) exactly when q > 2*2*3 = 12; iterating again, (2,55): degree %d; crossing numbers multiply, Delta(1) stays 1" % (sp.degree(D13, t), sp.degree(D15, t)))
print("   the three names of the figure -- hyperbolic, toroidal, elliptic -- are where the substitution accumulates: at the boundary, periodically, at the centre; the same three words name the elements of SL(2,R) by where their fixed points sit, which is a pun until a map is written")
print("ALL CHECKS PASSED")
