"""Bessel knots and the knot's zeta: nodal knots in the Helmholtz solution space, the Alexander polynomial as a zeta function, and the
ladder from the trefoil to a knotted S^5 in S^7 (eight real dimensions).  2026-10-19.  EXACT unless marked.

(1) The Helmholtz solution space is Bessel's: J_m(b R) e^{i m phi} e^{i z sqrt(1 - b^2)} in R^3, r^{1 - d/2} J_{l + d/2 - 1}(r) Y_l in R^d.
(2) Berry-Dennis (Proc. R. Soc. A 457, 2001): three J_m beams tuned so that the loop of zeros is doubly degenerate, plus a weak J_0 beam,
    have the (m, 2) torus knot for a nodal line, threaded by an m-stranded helix.  We reproduce their trefoil numbers, track the nodal
    line of the field and compute its Alexander polynomial from a projection: t^2 - t + 1 (and the cinquefoil for m = 5, the Hopf link
    for m = 2).  Enciso-Peralta-Salas (Adv. Math. 2013, Thm 1.3): every compact codimension-2 submanifold with trivial normal bundle --
    every n-knot in R^{n+2} -- is, after a C^p-small diffeomorphism, a union of components of {u1 = u2 = 0}, u_r solutions of any
    second-order elliptic equation with analytic coefficients, Helmholtz included.
(3) The Alexander polynomial is a zeta function three ways: Burau/Lin-Wang det(I - B(t)) (1 - t)/(1 - t^n); Milnor's
    1/Delta(t) = exp(sum_k tr(M^k) t^k/k) for the monodromy M; Fox's |H_1(Sigma_k)| = |det(I - M^k)|, the knot's 'point count', with the
    Weil zeta of 143a1 beside it (the analogy is the hypothesis, both sides are theorems).  Ihara's zeta of a graph is the same
    identity on the non-backtracking matrix (Bass).
(4) Sebastiani-Thom: adding a square to the polynomial multiplies the monodromy by -1, so z1^2 + z2^3 (trefoil, S^3) -> + z3^2 (L(3,1) in
    S^5, Delta = Phi_3) -> + z4^2 (a knotted S^5 in S^7 of R^8, Delta = Phi_6, Milnor number 2).  Levine (1970): simple (2q-1)-knots,
    q >= 3, are classified by the Seifert form, so in R^8 the algebra is complete.
(5) Bessel with zeta, exactly: Rayleigh's sums are the zeta function of the Bessel zeros (the ball's Dirichlet spectrum); the Mellin
    transform of K_nu is the Gamma factor of an L-function; the Eisenstein series on the modular surface (whose unit tangent bundle is
    the trefoil complement) has K-Bessel Fourier modes with coefficients sigma_{1-2s}(n)/zeta(2s), and -(1/2) d/ds sigma_{1-2s}(n) at
    s = 1 is the MTFT weight w_n = sum_{d|n} log(d)/d."""
import math, itertools, bisect
import numpy as np
import sympy as sp
from scipy.special import jv, jvp, jn_zeros, kv, gamma as Gamma
from scipy.optimize import fsolve
import mpmath as mp
rng = np.random.default_rng(3)
t = sp.symbols('t')

# ---------------------------------------------------------------- 1. the Helmholtz solution space is Bessel's
def beam(m, b, R, phi, z): return jv(m, b * R) * np.exp(1j * m * phi) * np.exp(1j * z * math.sqrt(1 - b * b))
def helmholtz_residual(f, x, y, z, h=1e-3):                               # (Laplacian + 1) f at (x, y, z), central differences
    F = lambda X, Y, Z: f(np.hypot(X, Y), np.arctan2(Y, X), Z)
    lap = (F(x + h, y, z) + F(x - h, y, z) + F(x, y + h, z) + F(x, y - h, z) + F(x, y, z + h) + F(x, y, z - h) - 6 * F(x, y, z)) / h ** 2
    return abs(lap + F(x, y, z)) / max(1e-12, abs(F(x, y, z)))
for m, b in ((3, 1.0), (3, 1 / 3), (0, 0.25), (5, 2 / 3)):
    assert helmholtz_residual(lambda R, p, z: beam(m, b, R, p, z), 2.3, 1.1, 0.7) < 1e-5
def radial_residual(d, l, r):                                             # u = r^{1-d/2} J_{l+d/2-1}(r) solves u'' + (d-1)/r u' + (1 - l(l+d-2)/r^2) u = 0
    nu = l + d / 2 - 1; h = 1e-4; u = lambda s: s ** (1 - d / 2) * jv(nu, s)
    return abs((u(r + h) - 2 * u(r) + u(r - h)) / h ** 2 + (d - 1) / r * (u(r + h) - u(r - h)) / (2 * h) + (1 - l * (l + d - 2) / r ** 2) * u(r)) / abs(u(r))
for d, l in ((3, 0), (3, 2), (8, 0), (8, 2), (8, 5)): assert radial_residual(d, l, 3.7) < 1e-5
print("Helmholtz: J_m(bR) e^{im phi} e^{iz sqrt(1-b^2)} solves (Laplacian + 1) psi = 0 in R^3 (1e-5); in R^8 the radial solutions are r^{-3} J_{l+3}(r): the Bessel solution space is the Helmholtz solution space in every dimension")
print("   compact nodal sets exist for Helmholtz (the spheres r = j_{3,k} = %s in R^8), never for harmonic functions (maximum principle)" % ", ".join("%.3f" % x for x in jn_zeros(3, 3)))

# ---------------------------------------------------------------- 2. Berry-Dennis: a (m,2) torus knot as a nodal line of four Bessel beams
B = (1.0, 1 / 3, 2 / 3); KZ = tuple(math.sqrt(1 - b * b) for b in B)
def coalesce(m, a2, a3, R):                                               # F(R,0) = F_R(R,0) = F_z(R,0)/i = 0 with a1 = 1: the two smallest zeros of F(.,0) coalesce
    a = (1.0, a2, a3)
    return [sum(a[l] * jv(m, B[l] * R) for l in range(3)), sum(a[l] * B[l] * jvp(m, B[l] * R) for l in range(3)), sum(a[l] * KZ[l] * jv(m, B[l] * R) for l in range(3))]
def solve_coalescence(m, guess):
    sol, info, ok, msg = fsolve(lambda v: coalesce(m, *v), guess, full_output=True, xtol=1e-13)
    assert ok == 1 and max(abs(x) for x in coalesce(m, *sol)) < 1e-10, (m, msg)
    return sol
paper = (10.0302, -3.18960, 5.44992)                                       # Berry-Dennis's trefoil: a2, a3, R*
a2, a3, Rs = solve_coalescence(3, paper)
assert abs(a2 - paper[0]) < 2e-4 and abs(a3 - paper[1]) < 2e-5 and abs(Rs - paper[2]) < 2e-5, (a2, a3, Rs)
print("Berry-Dennis trefoil: b = (1, 1/3, 2/3), a = (1, %.5f, %.5f), R* = %.5f -- the paper's (1, 10.0302, -3.18960), 5.44992 reproduced" % (a2, a3, Rs))
class Field:
    def __init__(self, m, a, eps=0.02, bp=0.25):
        self.m, self.a, self.eps, self.bp, self.kp = m, a, eps, bp, math.sqrt(1 - bp * bp)
    def F(self, R, z):   return sum(self.a[l] * jv(self.m, B[l] * R) * np.exp(1j * KZ[l] * z) for l in range(3))
    def FR(self, R, z):  return sum(self.a[l] * B[l] * jvp(self.m, B[l] * R) * np.exp(1j * KZ[l] * z) for l in range(3))
    def Fz(self, R, z):  return sum(self.a[l] * 1j * KZ[l] * jv(self.m, B[l] * R) * np.exp(1j * KZ[l] * z) for l in range(3))
    def P(self, R, z):   return jv(0, self.bp * R) * np.exp(1j * self.kp * z)
    def psi(self, R, phi, z): return np.exp(1j * self.m * phi) * self.F(R, z) - self.eps * self.P(R, z)
    def grad(self, R, phi, z):                                            # d psi / dR, d psi / dz (complex)
        e = np.exp(1j * self.m * phi)
        return e * self.FR(R, z) + self.eps * self.bp * jv(1, self.bp * R) * np.exp(1j * self.kp * z), e * self.Fz(R, z) - self.eps * 1j * self.kp * self.P(R, z)
    def newton(self, R, phi, z, iters=30):
        for _ in range(iters):
            f = self.psi(R, phi, z); dR, dz = self.grad(R, phi, z)
            J = np.array([[dR.real, dz.real], [dR.imag, dz.imag]]); step = np.linalg.solve(J, -np.array([f.real, f.imag]))
            R += step[0]; z += step[1]
            if np.hypot(*step) < 1e-13: break
        return R, z, abs(self.psi(R, phi, z))
    def roots_at(self, phi, Rc, r0):                                      # the roots of psi(., phi, .) near the ring
        found = []
        for th in np.linspace(0, 2 * np.pi, 16, endpoint=False):
            R, z, res = self.newton(Rc + r0 * math.cos(th), phi, r0 * math.sin(th))
            if res < 1e-9 and abs(R - Rc) < 6 * r0 and abs(z) < 6 * r0 and all(np.hypot(R - p[0], z - p[1]) > 1e-6 for p in found): found.append((R, z))
        return found
    def track(self, phi0, R, z, turns, steps=720):                        # continue a root in phi; returns the curve in (x, y, z)
        pts = []
        for i in range(steps * turns):
            phi = phi0 + 2 * np.pi * i / steps
            R, z, res = self.newton(R, phi, z); assert res < 1e-9
            pts.append((R * math.cos(phi), R * math.sin(phi), z, R, phi))
        return pts
def local_winding(fld, Rc, r=1e-2, n=720):                                 # degree of w -> F(R* + w) on a small circle (2 = a doubly degenerate loop of the knotting kind)
    vals = np.array([fld.F(Rc + r * math.cos(th), r * math.sin(th)) for th in np.linspace(0, 2 * np.pi, n, endpoint=False)])
    return int(round(np.sum(np.angle(np.roll(vals, -1) / vals)) / (2 * np.pi)))

# ---- the Alexander polynomial of a polygonal link from a generic projection (Wirtinger presentation, Fox calculus)
def alexander(components, tries=6):
    for _ in range(tries):
        Q, _ = np.linalg.qr(rng.normal(size=(3, 3)))
        comps = [np.asarray(c)[:, :3] @ Q.T for c in components]
        cross = []                                                            # (comp_i, s_i, comp_j, s_j, over_is_i, sign)
        ok = True
        for ci, P in enumerate(comps):
            for cj, Qp in enumerate(comps):
                if cj < ci: continue
                N, M = len(P), len(Qp)
                p = P[:, :2]; r = np.roll(P, -1, 0)[:, :2] - p; q = Qp[:, :2]; s = np.roll(Qp, -1, 0)[:, :2] - q
                dp = P[:, 2]; dq = Qp[:, 2]; dpn = np.roll(dp, -1); dqn = np.roll(dq, -1)
                for i in range(N):
                    js = np.arange(M)
                    if ci == cj: js = js[(js > i + 1) & ~((i == 0) & (js == N - 1))]
                    if len(js) == 0: continue
                    den = r[i, 0] * s[js, 1] - r[i, 1] * s[js, 0]
                    good = np.abs(den) > 1e-14
                    qp = q[js] - p[i]
                    ti = np.where(good, (qp[:, 0] * s[js, 1] - qp[:, 1] * s[js, 0]) / np.where(good, den, 1), -1)
                    tj = np.where(good, (qp[:, 0] * r[i, 1] - qp[:, 1] * r[i, 0]) / np.where(good, den, 1), -1)
                    hit = good & (ti >= 0) & (ti < 1) & (tj >= 0) & (tj < 1)
                    for k in np.nonzero(hit)[0]:
                        j = js[k]
                        if min(ti[k], 1 - ti[k], tj[k], 1 - tj[k]) < 1e-6: ok = False
                        di = dp[i] + ti[k] * (dpn[i] - dp[i]); dj = dq[j] + tj[k] * (dqn[j] - dq[j])
                        over_i = di > dj
                        d_over, d_under = (r[i], s[j]) if over_i else (s[j], r[i])
                        sign = np.sign(d_over[0] * d_under[1] - d_over[1] * d_under[0])
                        cross.append((ci, i + ti[k], cj, j + tj[k], over_i, sign))
        if not ok: continue
        # arcs: per component, sorted under-crossing positions
        unders = [[] for _ in comps]
        for ci, si, cj, sj, over_i, sign in cross:
            if over_i: unders[cj].append(sj)
            else: unders[ci].append(si)
        for u in unders: u.sort()
        offset = []; narcs = 0
        for u in unders: offset.append(narcs); narcs += max(1, len(u))
        def arc_of(c, s):
            u = unders[c]
            if not u: return offset[c]
            return offset[c] + (bisect.bisect_right(u, s) - 1) % len(u)
        n = len(cross)
        if n == 0: return sp.Integer(1)
        A = sp.zeros(n, narcs)
        for row, (ci, si, cj, sj, over_i, sign) in enumerate(cross):
            if over_i: co, so, cu, su = ci, si, cj, sj
            else: co, so, cu, su = cj, sj, ci, si
            k = arc_of(co, so); u = unders[cu]; idx = bisect.bisect_left(u, su); assert abs(u[idx] - su) < 1e-12
            out_arc = offset[cu] + idx; in_arc = offset[cu] + (idx - 1) % len(u)
            A[row, k] += 1 - t
            if sign > 0: A[row, in_arc] += t; A[row, out_arc] += -1
            else: A[row, in_arc] += -1; A[row, out_arc] += t
        assert narcs == n, (narcs, n)
        D = sp.expand(A[:-1, :-1].det(method='berkowitz'))
        if D == 0: return sp.Integer(0)
        P = sp.Poly(D, t); low = min(mon[0] for mon in P.monoms()); D = sp.expand(sp.cancel(D / t ** low))
        if sp.Poly(D, t).LC() < 0: D = -D
        return D, n
    raise RuntimeError('no generic projection found')
def same_poly(D, E): return sp.expand(D - E) == 0 or sp.expand(D - sp.expand(sp.cancel(E.subs(t, 1 / t) * t ** sp.degree(E, t)))) == 0
def torus_curve(p, q, n=600): return [((2 + math.cos(q * s)) * math.cos(p * s), (2 + math.cos(q * s)) * math.sin(p * s), math.sin(q * s)) for s in np.linspace(0, 2 * np.pi, n, endpoint=False)]
D3, c3 = alexander([torus_curve(2, 3)]); assert same_poly(D3, t ** 2 - t + 1), D3
D5, c5 = alexander([torus_curve(2, 5)]); assert same_poly(D5, t ** 4 - t ** 3 + t ** 2 - t + 1), D5
D8, c8 = alexander([[((2 + math.cos(2 * s)) * math.cos(3 * s), (2 + math.cos(2 * s)) * math.sin(3 * s), math.sin(4 * s)) for s in np.linspace(0, 2 * np.pi, 600, endpoint=False)]]); assert same_poly(D8, t ** 2 - 3 * t + 1), D8
circ = [(math.cos(s), math.sin(s), 0.0) for s in np.linspace(0, 2 * np.pi, 200, endpoint=False)]; circ2 = [(1 + math.cos(s), 0.0, math.sin(s)) for s in np.linspace(0, 2 * np.pi, 200, endpoint=False)]
DH, cH = alexander([circ, circ2]); assert same_poly(DH, t - 1), DH
assert alexander([circ]) == 1
print("Alexander polynomial from a projection (Wirtinger + Fox): T(2,3) -> %s, T(2,5) -> %s, figure-eight -> %s, Hopf link -> %s, unknot -> 1 (test curves)" % (D3, D5, D8, DH))

results = {}
for m, guess in ((3, (a2, a3, Rs)), (5, (10.0, -3.0, 7.5)), (2, (4.73341, -2.70176, 4.32636))):
    sol = solve_coalescence(m, guess); fld = Field(m, (1.0, sol[0], sol[1]))
    Rc = sol[2]; wind = local_winding(fld, Rc)
    FRR = (fld.F(Rc + 1e-4, 0) - 2 * fld.F(Rc, 0) + fld.F(Rc - 1e-4, 0)).real / 1e-8; Fzz = (fld.F(Rc, 1e-4) - 2 * fld.F(Rc, 0) + fld.F(Rc, -1e-4)).real / 1e-8
    roots = fld.roots_at(0.0, Rc, 0.12)
    assert len(roots) == 2, (m, roots)
    curve = fld.track(0.0, *roots[0], turns=1)
    end = fld.newton(curve[-1][3], 2 * np.pi, curve[-1][2])                   # where the root arrives after one circuit
    swapped = np.hypot(end[0] - roots[1][0], end[1] - roots[1][1]) < 1e-6
    if swapped:
        curve = fld.track(0.0, *roots[0], turns=2); comps = [[p[:3] for p in curve]]
    else:
        comps = [[p[:3] for p in curve], [p[:3] for p in fld.track(0.0, *roots[1], turns=1)]]
    # windings: round the axis (longitude) and round the ring (meridian), from the tracked curve
    merid = sum(np.angle(((c[3] - Rc) + 1j * c[2]) / ((curve[i - 1][3] - Rc) + 1j * curve[i - 1][2])) for i, c in enumerate(curve) if i > 0) / (2 * np.pi)
    longi = (curve[-1][4] - curve[0][4] + 2 * np.pi / 720) / (2 * np.pi)
    # embedded: the closed curve has no near self-contacts
    P = np.array(comps[0]); dmin = min(np.min(np.linalg.norm(P[i] - P[i + 3:], axis=1)) for i in range(0, len(P) - 3, 7))
    D, nc = alexander(comps)
    expected = sp.expand(sp.cancel((1 - (-t) ** m) / (1 + t)))               # Burau: det(I - B(sigma_1^m)) (1 - t)/(1 - t^2) for the braid sigma_1^m
    assert same_poly(D, expected), (m, D, expected)
    assert wind == 2 and FRR * Fzz < 0 and abs(FRR + Fzz) < 1e-5 * abs(FRR) and round(abs(merid)) == (m if swapped else m // 2) and round(longi) == (2 if swapped else 1) and dmin > 1e-3, (m, wind, FRR, Fzz, merid, longi, dmin)
    # F_RR = -F_zz exactly at the double zero: the Helmholtz equation in cylindrical coordinates, F_RR + F_R/R - m^2 F/R^2 + F_zz + F = 0, with F = F_R = 0 there
    # the m-stranded helix near the axis: m roots of psi(x, y, z) = 0 at each height
    strands = []
    for z in np.linspace(-3, 3, 7):
        Rg, Pg = np.meshgrid(np.linspace(0.2, 3.0, 57), np.linspace(0, 2 * np.pi, 72, endpoint=False), indexing='ij')
        V = np.abs(fld.psi(Rg, Pg, z)); cand = sorted(zip(V.ravel(), Rg.ravel(), Pg.ravel()))[:40]
        found = []
        for v, R, ph2 in cand:
            for _ in range(40):                                               # Newton in (x, y) at fixed z, via (R, phi)
                f = fld.psi(R, ph2, z); dR, dz = fld.grad(R, ph2, z); dphi = 1j * m * np.exp(1j * m * ph2) * fld.F(R, z)
                J = np.array([[dR.real, dphi.real], [dR.imag, dphi.imag]]); st = np.linalg.solve(J, -np.array([f.real, f.imag])); R += st[0]; ph2 += st[1]
                if abs(st[0]) + R * abs(st[1]) < 1e-13: break
            if abs(fld.psi(R, ph2, z)) < 1e-10 and 0 < R < 3.0:
                p = (R * math.cos(ph2), R * math.sin(ph2))
                if all(np.hypot(p[0] - q[0], p[1] - q[1]) > 1e-6 for q in found): found.append(p)
        assert len(found) == m, (m, z, found)
        strands.append(found)
    rad = [np.hypot(*p) for s_ in strands for p in s_]
    results[m] = dict(a=(1.0, sol[0], sol[1]), R=Rc, D=D, ncross=nc, comps=len(comps), merid=merid, longi=longi, npts=len(curve), eps=fld.eps)
    name = {3: 'trefoil T(2,3)', 5: 'cinquefoil T(2,5)', 2: 'Hopf link T(2,2)'}[m]
    print("   m = %d: a = (1, %.5f, %.5f), R* = %.5f, local degree %d (F_RR F_zz < 0); nodal line tracked through %d x 2pi: %d component(s), %d turns round the ring, %d round the axis; Alexander polynomial from %d crossings: %s = %s" % (m, sol[0], sol[1], Rc, wind, 2 if swapped else 1, len(comps), round(abs(merid)), round(longi), nc, D, name))
    print("      threaded by a %d-stranded helix (%d roots at radius %.2f-%.2f at each of 7 heights); eps = %.2f, J_0(R/4) e^{i sqrt(15) z/4} perturbation as in the paper; F_RR = -F_zz at the loop (Helmholtz with F = F_R = 0), so the loop is always of the knotting kind" % (m, m, min(rad), max(rad), fld.eps))
print("   Berry-Dennis: the loop conditions cannot be met paraxially (no dislocation loops of strength > 1 in paraxial waves): Helmholtz, i.e. Bessel beams, is needed")
print("Enciso-Peralta-Salas 2013, Thm 1.3: compact codimension-m (m >= 2) submanifolds of R^n with trivial normal bundle are, after a C^p-small diffeomorphism, unions of components of joint level sets of solutions to second-order elliptic equations with analytic coefficients -- every n-knot in R^{n+2}, in the Helmholtz solution space, in every dimension")

# ---------------------------------------------------------------- 3. the Alexander polynomial as a zeta function
# (a) Burau / Lin-Tian-Wang: the reduced Burau matrix of sigma_1 in B_2 is (-t); Delta(closure of beta) = det(I - B(beta)) (1 - t)/(1 - t^n)
for m in (3, 5, 7):
    D = sp.expand(sp.cancel((1 - (-t) ** m) * (1 - t) / (1 - t ** 2))); assert same_poly(D, sp.cyclotomic_poly(2 * m, t)), (m, D)
assert sp.expand(sp.cancel((1 - (-t) ** 9) * (1 - t) / (1 - t ** 2)) - sp.cyclotomic_poly(6, t) * sp.cyclotomic_poly(18, t)) == 0
print("Burau (Lin-Tian-Wang): Delta(T(2,m)) = det(I - B(sigma_1^m)) (1-t)/(1-t^2) = (1 + t^m)/(1 + t) = Phi_{2m} for m = 3, 5, 7 (Phi_6 Phi_18 for m = 9); Lin-Wang (1998): 1/Delta is the Ihara-Selberg zeta of a weighted random walk on the diagram")
# (b) Milnor: the monodromy M of the trefoil fibration has charpoly t^2 - t + 1; 1/Delta(t) = 1/det(I - tM) = exp(sum tr(M^k) t^k / k)
M = sp.Matrix([[0, -1], [1, 1]]); assert sp.expand((sp.eye(2) - t * M).det()) == t ** 2 - t + 1 and M ** 6 == sp.eye(2) and M ** 3 == -sp.eye(2)
traces = [int((M ** k).trace()) for k in range(1, 25)]; assert traces == [1, -1, -2, -1, 1, 2] * 4
from fractions import Fraction
def exp_series(g, order):                                                   # f = exp(g), g(0) = 0: n f_n = sum_k k g_k f_{n-k}  (exact)
    f = [Fraction(1)] + [Fraction(0)] * order
    for n in range(1, order + 1): f[n] = sum(k * g[k] * f[n - k] for k in range(1, n + 1)) / n
    return f
def inverse_series(p, order):                                               # 1/p for p(0) = 1
    f = [Fraction(1)] + [Fraction(0)] * order
    for n in range(1, order + 1): f[n] = -sum(p[k] * f[n - k] for k in range(1, min(n, len(p) - 1) + 1))
    return f
g = [Fraction(0)] + [Fraction(traces[k - 1], k) for k in range(1, 25)]
assert exp_series(g, 24) == inverse_series([1, -1, 1], 24) and inverse_series([1, -1, 1], 24)[:12] == [1, 1, 0, -1, -1, 0] * 2
covers = [abs(int((sp.eye(2) - M ** k).det())) for k in range(1, 13)]; assert covers == [1, 3, 4, 3, 1, 0] * 2
import cmath
fox = [abs(math.prod((lambda z: z * z - z + 1)(cmath.exp(2j * math.pi * j / k)) for j in range(k))) for k in range(1, 8)]
assert all(abs(fox[k - 1] - covers[k - 1]) < 1e-9 for k in range(1, 8))
print("Milnor: 1/Delta(t) = exp(sum_k tr(M^k) t^k/k) for the order-6 monodromy M = ST (series agree to t^24: 1/Delta = 1 + t - t^3 - t^4 + t^6 + ...); tr M^k = %s, period 6 = the period clock" % traces[:6])
print("Fox: |H_1(Sigma_k)| = prod_j |Delta(zeta_k^j)| = |det(I - M^k)| = %s for k = 1..7 -- the knot's point counts; Zeeman: the k-twist-spun trefoil is a fibred 2-knot in S^4 with fibre the punctured Sigma_k (k = 5: the Poincaré sphere, Delta = 1; k = 6: the Nil manifold, infinite)" % covers[:7])
# (c) the Weil zeta of 143a1: y^2 + y = x^3 - x^2 - x - 2, counted by hand, beside the knot's
def count_Fp(p):
    n = 1                                                                     # the point at infinity
    for x in range(p):
        rhs = (x ** 3 - x * x - x - 2) % p
        n += sum(1 for y in range(p) if (y * y + y - rhs) % p == 0)
    return n
table = {2: 0, 3: -1, 5: -1, 7: -2, 11: -1, 13: -1, 17: -4, 19: 2, 23: 7, 29: -2, 31: -3, 37: -11, 41: 10, 43: -4, 47: -4, 53: 2, 59: -1, 61: -2}   # the package's a_p (x0_143.CURVE_143A1)
ap = {p: p + 1 - count_Fp(p) for p in table}
assert all(ap[p] == table[p] for p in table if p not in (11, 13)), ap    # at the bad primes 11, 13 the count gives the a_p of the singular fibre, which also matches here
assert ap[11] == table[11] and ap[13] == table[13]
def count_Fp2(p):                                                             # F_{p^2} = F_p[s]/(s^2 - d), d a non-residue; brute force
    d = next(x for x in range(2, p) if pow(x, (p - 1) // 2, p) == p - 1)
    els = [(u, v) for u in range(p) for v in range(p)]
    mul = lambda a, b: ((a[0] * b[0] + d * a[1] * b[1]) % p, (a[0] * b[1] + a[1] * b[0]) % p)
    add = lambda a, b: ((a[0] + b[0]) % p, (a[1] + b[1]) % p)
    n = 1
    for x in els:
        x2 = mul(x, x); x3 = mul(x2, x); rhs = add(add(x3, ((-x2[0]) % p, (-x2[1]) % p)), ((-x[0] - 2) % p, (-x[1]) % p))
        for y in els:
            if add(mul(y, y), y) == rhs: n += 1
    return n
for p in (3, 5):
    al = (ap[p] + 1j * math.sqrt(4 * p - ap[p] ** 2)) / 2; N2 = round(p * p + 1 - 2 * (al ** 2).real)
    assert count_Fp2(p) == N2, (p, count_Fp2(p), N2)
print("143a1 (y^2 + y = x^3 - x^2 - x - 2) counted by hand: a_p = %s; #E(F_{p^2}) = p^2 + 1 - (alpha^2 + conj) checked for p = 3, 5 -- the Weil zeta's local factor 1 - a_p T + p T^2 is det(1 - F T | H^1)" % [ap[p] for p in (2, 3, 5, 7, 11, 13, 17, 19, 23)])
print("   the dictionary (hypothesis): Delta(t) = 1 - t + t^2 reads as 1 - a t + p t^2 with a = 1, p = 1; |H_1(Sigma_k)| = 1 + 1 - tr M^k beside #E(F_{p^k}) = p^k + 1 - tr F^k; both are exact, the analogy (Mazur, Morishita: Alexander <-> Iwasawa) is the conjecture layer")
# (d) Ihara's zeta of a graph: the same identity on the non-backtracking matrix (Bass)
u = sp.symbols('u'); V = 4; E = [(i, j) for i in range(V) for j in range(i + 1, V)]                       # K4
darts = [(i, j) for i, j in E] + [(j, i) for i, j in E]; nd = len(darts)
T = sp.zeros(nd, nd)
for a, (i, j) in enumerate(darts):
    for b, (k, l) in enumerate(darts):
        if j == k and not (l == i): T[a, b] = 1
x = sp.symbols('x'); cp = sp.Poly(T.charpoly(x).as_expr(), x); zeta_inv = sp.factor(sp.expand(sum(c * u ** (nd - k) for k, c in enumerate(cp.all_coeffs()[::-1]))))   # det(I - uT) = u^nd charpoly_T(1/u)
A = sp.zeros(V, V)
for i, j in E: A[i, j] = A[j, i] = 1
bass = sp.factor((1 - u ** 2) ** (len(E) - V) * (sp.eye(V) - u * A + 2 * u ** 2 * sp.eye(V)).det())
assert sp.expand(zeta_inv - bass) == 0
walks = [int((T ** k).trace()) for k in range(1, 13)]
zcoef = [Fraction(int(c)) for c in sp.Poly(sp.expand(zeta_inv), u).all_coeffs()[::-1]]
assert exp_series([Fraction(0)] + [Fraction(walks[k - 1], k) for k in range(1, 13)], 12) == inverse_series(zcoef, 12)
print("Ihara zeta of K4: 1/zeta(u) = det(I - uT) on the 12 darts = %s = Bass's (1-u^2)^{|E|-|V|} det(I - Au + 2u^2); zeta = exp(sum_k N_k u^k/k) with N_k = closed non-backtracking walks %s (series to u^12) -- the same engine as Milnor's and Weil's" % (bass, walks[:6]))

# ---------------------------------------------------------------- 4. the ladder to eight dimensions: Milnor-Orlik and Sebastiani-Thom
def brieskorn_alexander(exps):                                              # eigenvalues of the monodromy: products of non-trivial a_i-th roots of unity (Pham, Milnor-Orlik), as exact fractions of a turn
    from fractions import Fraction
    from collections import Counter
    qs = Counter()
    for combo in itertools.product(*[range(1, a) for a in exps]): qs[sum(Fraction(j, a) for j, a in zip(combo, exps)) % 1] += 1
    D = sp.Integer(1); byden = {}
    for q, n in qs.items(): byden.setdefault(q.denominator, []).append((q.numerator, n))
    for d, lst in byden.items():                                              # Galois-stable: every primitive residue mod d appears equally often -> Phi_d^(count/phi(d))
        prim = [r for r in range(1, d + 1) if math.gcd(r, d) == 1]; cnt = dict(lst)
        assert set(cnt) == set(prim) if d > 1 else True; ns = set(cnt.values()); assert len(ns) == 1
        D *= sp.cyclotomic_poly(d, t) ** ns.pop()
    return sp.expand(D)
mu = lambda exps: sp.prod(a - 1 for a in exps)
ladder = []
for exps, dim, link in (((2, 3), 'S^1 in S^3 (R^4)', 'the trefoil'), ((2, 3, 2), 'a 3-manifold in S^5 (R^6)', 'the lens space L(3,1)'), ((2, 3, 2, 2), 'a 5-manifold in S^7 (R^8)', 'a homotopy 5-sphere = S^5 (Theta_5 = 0), knotted')):
    D = brieskorn_alexander(exps); ladder.append((exps, D, mu(exps), dim, link))
assert sp.expand(ladder[0][1] - (t ** 2 - t + 1)) == 0 and sp.expand(ladder[1][1] - (t ** 2 + t + 1)) == 0 and sp.expand(ladder[2][1] - (t ** 2 - t + 1)) == 0
assert all(l[2] == 2 for l in ladder)
assert ladder[0][1].subs(t, 1) == 1 and ladder[1][1].subs(t, 1) == 3 and ladder[2][1].subs(t, 1) == 1
# Sebastiani-Thom: adding z^2 joins with the A_1 fibre (two points, monodromy -1): eigenvalues flip sign
for k in range(2):
    assert sp.expand(ladder[k + 1][1] - sp.expand(ladder[k][1].subs(t, -t) * (1 if k == 1 else 1))) == 0 or sp.expand(ladder[k + 1][1] + sp.expand(ladder[k][1].subs(t, -t))) == 0
for m in (5, 7):
    D4 = brieskorn_alexander((2, m, 2, 2)); assert sp.expand(D4 - sp.cyclotomic_poly(2 * m, t)) == 0 and mu((2, m, 2, 2)) == m - 1
print("Milnor-Orlik / Sebastiani-Thom ladder (exact): z1^2 + z2^3 -> Delta = t^2 - t + 1, mu = 2, the trefoil in S^3; + z3^2 -> Delta(-t) = t^2 + t + 1, |H_1| = 3: L(3,1) in S^5; + z4^2 -> t^2 - t + 1 again, Delta(1) = 1: a knotted S^5 in S^7 of R^8 with the trefoil's monodromy (and T(2,m) -> Phi_{2m} in S^7, mu = m - 1)")
print("   Levine 1970 / Trotter / Kearton: simple (2q-1)-knots with q >= 3 are classified by the S-equivalence class of the Seifert matrix (the Blanchfield form on the Alexander module) -- in S^7 the algebra is a complete invariant, which fails in S^3; Levine 1966: Delta(1) = +-1 for every n-knot, Kinoshita 1961: every such Delta is a 2-knot's")

# ---------------------------------------------------------------- 5. Bessel with zeta, exactly
# (a) Rayleigh: sum_k j_{nu,k}^{-2} = 1/(4(nu+1)), sum_k j_{nu,k}^{-4} = 1/(16 (nu+1)^2 (nu+2)): the zeta function of the Bessel zeros = the Dirichlet spectrum of the disc / the 8-ball
for nu in (0, 1, 3):
    z = jn_zeros(nu, 4000); K = len(z); tail = 1 / (math.pi ** 2 * (K + nu / 2 - 0.25 + 0.5))
    assert abs(np.sum(z ** -2.0) + tail - 1 / (4 * (nu + 1))) < 3e-8, nu
    assert abs(np.sum(z ** -4.0) - 1 / (16 * (nu + 1) ** 2 * (nu + 2))) < 1e-10, nu
print("Rayleigh: zeta_{J_nu}(2) = sum j_{nu,k}^{-2} = 1/(4(nu+1)) and zeta_{J_nu}(4) = 1/(16(nu+1)^2(nu+2)) for nu = 0, 1, 3 (4000 zeros + tail): the spectral zeta of the ball's Dirichlet Laplacian is a Bessel zeta; nu = 3 is the radial index of R^8")
# (b) the Mellin transform of K_nu is a pair of Gamma factors: the archimedean factor of a Maass form's L-function has nu = i r
mp.mp.dps = 20
for nu, s in ((0.5, 2.0), (0.3j, 1.7), (1.0, 3.2)):
    lhs = mp.quad(lambda y: mp.besselk(nu, y) * y ** (s - 1), [0, 1, mp.inf]); rhs = 2 ** (s - 2) * mp.gamma((s - nu) / 2) * mp.gamma((s + nu) / 2)
    assert abs(lhs - rhs) < 1e-12, (nu, s, lhs, rhs)
print("Mellin: int_0^inf K_nu(y) y^{s-1} dy = 2^{s-2} Gamma((s-nu)/2) Gamma((s+nu)/2) (1e-12, nu real and imaginary) -- the Gamma factor of a completed L-function is a Bessel transform")
# (c) the Eisenstein series of SL(2,Z): lattice sum against the K-Bessel Fourier expansion; its eigen-equation; the functional equation; the weights
def sigma(n, k): return sum(sp.Integer(d) ** k for d in range(1, n + 1) if n % d == 0)
def E_lattice(z, s, Nmax=400):
    x, y = z.real, z.imag; tot = 0.0
    for c in range(-Nmax, Nmax + 1):
        for d in range(-Nmax, Nmax + 1):
            if (c, d) == (0, 0) or math.gcd(c, d) != 1: continue
            tot += y ** s / abs(c * z + d) ** (2 * s)
    return tot / 2
def E_fourier(z, s, nmax=12):
    x, y = z.real, z.imag
    phi = math.sqrt(math.pi) * Gamma(s - 0.5) * float(sp.zeta(2 * s - 1)) / (Gamma(s) * float(sp.zeta(2 * s)))
    tot = y ** s + phi * y ** (1 - s)
    for n in range(1, nmax + 1):
        tot += 2 * (2 * math.pi ** s * math.sqrt(y) / (Gamma(s) * float(sp.zeta(2 * s)))) * n ** (s - 0.5) * float(sigma(n, 1 - 2 * s)) * kv(s - 0.5, 2 * math.pi * n * y) * math.cos(2 * math.pi * n * x)
    return tot
z0 = 0.3 + 0.8j
for s, tol in ((2.0, 2e-5), (2.5, 3e-7)):                                   # the lattice sum is truncated at |c|,|d| <= 400: error ~ 400^{2-2s}
    lat, fou = E_lattice(z0, s), E_fourier(z0, s)
    assert abs(lat - fou) < tol, (s, lat, fou)
s = 2.5; h = 1e-3
lap = -(z0.imag ** 2) * ((E_fourier(z0 + h, s) + E_fourier(z0 - h, s) + E_fourier(z0 + 1j * h, s) + E_fourier(z0 - 1j * h, s) - 4 * E_fourier(z0, s)) / h ** 2)
assert abs(lap - s * (1 - s) * E_fourier(z0, s)) < 1e-5 * abs(E_fourier(z0, s)), (lap, s * (1 - s) * E_fourier(z0, s))
phi = lambda s: mp.sqrt(mp.pi) * mp.gamma(s - 0.5) * mp.zeta(2 * s - 1) / (mp.gamma(s) * mp.zeta(2 * s))
for s_ in (0.3 + 1.1j, 0.7, 2.2): assert abs(phi(mp.mpc(s_)) * phi(1 - mp.mpc(s_)) - 1) < 1e-12
w = lambda n: sum(math.log(d) / d for d in range(1, n + 1) if n % d == 0)
for n in (2, 6, 12, 143):
    ds = 1e-6; der = (float(sigma(n, 1 - 2 * (1 + ds))) - float(sigma(n, 1 - 2 * (1 - ds)))) / (2 * ds)
    assert abs(-der / 2 - w(n)) < 1e-7, (n, der, w(n))
print("Eisenstein series E(z,s) of SL(2,Z) at z = 0.3 + 0.8i: lattice sum = y^s + phi(s) y^{1-s} + (2 pi^s sqrt(y)/(Gamma(s) zeta(2s))) sum n^{s-1/2} sigma_{1-2s}(n) K_{s-1/2}(2 pi n y) e(nx) (2e-5 at s = 2, 3e-7 at s = 2.5, the lattice truncation); -y^2 Laplacian E = s(1-s) E; phi(s) phi(1-s) = 1 (the functional equation of zeta)")
print("   the Fourier modes sqrt(y) K_{s-1/2}(2 pi n y) are the Bessel solution space of the hyperbolic Helmholtz equation on the modular surface, whose unit tangent bundle is the trefoil complement; and -(1/2) d/ds sigma_{1-2s}(n)|_{s=1} = sum_{d|n} log(d)/d = w_n (checked n = 2, 6, 12, 143): the MTFT weights are the s-derivative at s = 1 of the Eisenstein coefficients")
print("ALL CHECKS PASSED")
