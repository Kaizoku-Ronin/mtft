#!/usr/bin/env python3
"""The Modular Plane, version 8, Route A: Room 12 (gravity: the frozen field).

    python3 -I mp8_verify.py        ->  mp8_numbers.json next to this file

Pre-registered in MP08_Version8_Plan.md (written before this script was first run). sympy for the exact identities,
numpy and scipy for the rest; seeded 143. A failed check prints MISS and is counted; the script does not stop at it.
Conventions: a = x - iy, b = t - i - z, d = r^2 - (t - i)^2; labels alpha = 1 - 2ib/d, beta = 2a/d; U = [[alpha, beta], [-beta*, alpha*]];
the chart at t = 0 is (x, y, z) = (Re beta, -Im beta, Im alpha)/(1 - Re alpha); F = n/d^3 = E + iB."""
import bisect, json, math, pathlib
import numpy as np, sympy as sp
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

rng = np.random.default_rng(143)
t = sp.Symbol('t')                       # the variable of the Alexander polynomial (times are tm, tv below)
NPASS = 0; NFAIL = 0; NREC = 0; NUM = {}
def ok(name, cond, detail='', recorded_miss=False):
    """recorded_miss: a pre-registered prediction that the first run missed; it is kept, re-evaluated and reported, never edited away"""
    global NPASS, NFAIL, NREC
    cond = bool(cond)
    if recorded_miss:
        NREC += 1; print(('  MISS (recorded) ' if not cond else '  ?? a recorded miss now passes: ') + name + (('   [' + detail + ']') if detail else ''))
        if cond: NFAIL += 1
        return
    if cond: NPASS += 1
    else: NFAIL += 1
    print(('  ok   ' if cond else '  MISS ') + name + (('   [' + detail + ']') if detail else ''))
def cj(zv): return [float(np.real(zv)), float(np.imag(zv))]            # a complex number as the page stores it
def fl(v): return [float(u) for u in np.asarray(v).ravel()]

# ---------------------------------------------------------------- the hopfion's two numbers, its labels, chart, field and flow
def ab(tv, X): X = np.asarray(X, float); return X[..., 0] - 1j * X[..., 1], tv - 1j - X[..., 2]
def dd(tv, X): X = np.asarray(X, float); return np.sum(X * X, -1) - (tv - 1j) ** 2
def labels(tv, X): a, b = ab(tv, X); d = dd(tv, X); return 1 - 2j * b / d, 2 * a / d
def nvec(tv, X): a, b = ab(tv, X); return np.stack([4 * (b * b - a * a), -4j * (a * a + b * b), 8 * a * b], -1)
def Fpq(tv, X, p=1, q=1): al, be = labels(tv, X); return (p * q * al ** (p - 1) * be ** (q - 1))[..., None] * nvec(tv, X) / dd(tv, X)[..., None] ** 3
def vflow(tv, X): n = nvec(tv, X); return np.real(1j * np.cross(n, np.conj(n))) / np.real(np.sum(n * np.conj(n), -1))[..., None]
def chart(tv, al, be):
    w = 1 / (1 - al); tau = (tv - w.imag) / w.real; xy = be * (1 + 1j * tau) * w
    return np.stack([xy.real, -xy.imag, tv - tau], -1)
PAULI = {1: np.array([[0, 1], [1, 0]], complex), 2: np.array([[0, -1j], [1j, 0]]), 3: np.array([[1, 0], [0, -1]], complex)}
def flow_orbit(side, m, al, be, N):
    """first row of exp(i phi sigma_m) U (side 'L') or of U exp(i phi sigma_m) (side 'R'), phi over one turn"""
    phi = np.linspace(0, 2 * np.pi, N, endpoint=False); c, s = np.cos(phi), np.sin(phi); S = PAULI[m]
    e11, e12, e21, e22 = c + 1j * s * S[0, 0], 1j * s * S[0, 1], 1j * s * S[1, 0], c + 1j * s * S[1, 1]
    if side == 'L': return e11 * al + e12 * (-np.conj(be)), e11 * be + e12 * np.conj(al)
    return al * e11 + be * e21, al * e12 + be * e22
def hopf_point(eta, xi1, xi2): return np.cos(eta) * np.exp(1j * xi1), np.sin(eta) * np.exp(1j * xi2)
def core(p, q, psi0, N):
    phi = np.linspace(0, 2 * np.pi, N, endpoint=False)
    return math.sqrt(p / (p + q)) * np.exp(1j * (q * phi + psi0 / p)), math.sqrt(q / (p + q)) * np.exp(-1j * p * phi)

# ---------------------------------------------------------------- closed polygons: Gauss integral, writhe, crossings, twist
def kl_omega(P1, P2, Q1, Q2):
    """signed solid angle of segment P1->P2 against Q1->Q2 (Klenin and Langowski 2000, method 1a, with arcsin(n.n') written as pi/2 - atan2(|n x n'|, n.n'):
    the first run of this script used arcsin and lost half the digits on a planar curve, writhe -1.4e-7 instead of 0)"""
    r13 = Q1 - P1; r14 = Q2 - P1; r23 = Q1 - P2; r24 = Q2 - P2
    n1 = np.cross(r13, r14); n2 = np.cross(r14, r24); n3 = np.cross(r24, r23); n4 = np.cross(r23, r13)
    ang = lambda u, v: np.arctan2(np.linalg.norm(np.cross(u, v), axis=-1), np.sum(u * v, -1))
    om = 2 * np.pi - ang(n1, n2) - ang(n2, n3) - ang(n3, n4) - ang(n4, n1)
    return om * np.sign(np.sum(np.cross(Q2 - Q1, P2 - P1) * r13, -1))
def gauss_poly(C1, C2): return float(kl_omega(C1[:, None], np.roll(C1, -1, 0)[:, None], C2[None], np.roll(C2, -1, 0)[None]).sum() / (4 * np.pi))
def writhe_poly(C):
    N = len(C); om = kl_omega(C[:, None], np.roll(C, -1, 0)[:, None], C[None], np.roll(C, -1, 0)[None]); i, j = np.indices((N, N))
    return float(np.where((np.abs(i - j) > 1) & (np.abs(i - j) < N - 1), om, 0).sum() / (4 * np.pi))
VIEW = np.array([[0.36, 0.48, 0.8], [-0.8, 0.6, 0.0], [-0.48, -0.64, 0.6]])     # a fixed proper rotation; its first row is the line of sight
def lk_cross(C1, C2, view=VIEW):
    """half the signed crossings of two closed polygons seen along view[0]: an integer"""
    A = C1 @ view[[1, 2, 0]].T; B = C2 @ view[[1, 2, 0]].T
    p = A[:, None, :]; r = (np.roll(A, -1, 0) - A)[:, None, :]; q = B[None, :, :]; s = (np.roll(B, -1, 0) - B)[None, :, :]
    den = r[..., 0] * s[..., 1] - r[..., 1] * s[..., 0]; qp = q - p
    with np.errstate(divide='ignore', invalid='ignore'):
        ti = (qp[..., 0] * s[..., 1] - qp[..., 1] * s[..., 0]) / den; tj = (qp[..., 0] * r[..., 1] - qp[..., 1] * r[..., 0]) / den
        hit = (np.abs(den) > 1e-14) & (ti >= 0) & (ti < 1) & (tj >= 0) & (tj < 1)
        h1 = p[..., 2] + ti * r[..., 2]; h2 = q[..., 2] + tj * s[..., 2]
        sg = np.where(hit, np.sign(den) * np.sign(h1 - h2), 0)
    return float(0.5 * sg.sum()), int(hit.sum())
def dper(C):
    N = len(C); k = np.fft.fftfreq(N, 1.0 / N)
    if N % 2 == 0: k[N // 2] = 0
    return np.real(np.fft.ifft(1j * k[:, None] * np.fft.fft(C, axis=0), axis=0))
def gauss_smooth(C1, C2):
    T1, T2 = dper(C1), dper(C2); Rr = C1[:, None, :] - C2[None, :, :]; cr = np.cross(T1[:, None, :], T2[None, :, :])
    return float(np.sum(np.einsum('ijk,ijk->ij', Rr, cr) / np.linalg.norm(Rr, axis=2) ** 3) * (2 * np.pi / len(C1)) * (2 * np.pi / len(C2)) / (4 * np.pi))
def twist(C, Cn):
    T = dper(C); T /= np.linalg.norm(T, axis=1)[:, None]; U = Cn - C; U -= np.sum(U * T, 1)[:, None] * T; U /= np.linalg.norm(U, axis=1)[:, None]
    return float(np.sum(np.einsum('ij,ij->i', np.cross(T, U), dper(U))) / len(C))
sinang = lambda A, B: np.linalg.norm(np.cross(A, B), axis=-1) / (np.linalg.norm(A, axis=-1) * np.linalg.norm(B, axis=-1))

# ---------------------------------------------------------------- the Alexander polynomial of a polygonal link from a projection (from gravity_knot_20261021.py)
def alexander(components, tries=6):
    for _ in range(tries):
        Qm, _ = np.linalg.qr(rng.normal(size=(3, 3))); comps = [np.asarray(c)[:, :3] @ Qm.T for c in components]; cross = []; good_proj = True
        for ci, Pp in enumerate(comps):
            for cj, Qp in enumerate(comps):
                if cj < ci: continue
                Nn, Mm = len(Pp), len(Qp); p = Pp[:, :2]; rr = np.roll(Pp, -1, 0)[:, :2] - p; q = Qp[:, :2]; s = np.roll(Qp, -1, 0)[:, :2] - q
                dp = Pp[:, 2]; dq = Qp[:, 2]; dpn = np.roll(dp, -1); dqn = np.roll(dq, -1)
                for i in range(Nn):
                    js = np.arange(Mm)
                    if ci == cj: js = js[(js > i + 1) & ~((i == 0) & (js == Nn - 1))]
                    if len(js) == 0: continue
                    den = rr[i, 0] * s[js, 1] - rr[i, 1] * s[js, 0]; good = np.abs(den) > 1e-14; qp = q[js] - p[i]
                    ti = np.where(good, (qp[:, 0] * s[js, 1] - qp[:, 1] * s[js, 0]) / np.where(good, den, 1), -1); tj = np.where(good, (qp[:, 0] * rr[i, 1] - qp[:, 1] * rr[i, 0]) / np.where(good, den, 1), -1)
                    hit = good & (ti >= 0) & (ti < 1) & (tj >= 0) & (tj < 1)
                    for k in np.nonzero(hit)[0]:
                        j = js[k]
                        if min(ti[k], 1 - ti[k], tj[k], 1 - tj[k]) < 1e-6: good_proj = False
                        di = dp[i] + ti[k] * (dpn[i] - dp[i]); dj = dq[j] + tj[k] * (dqn[j] - dq[j]); over_i = di > dj
                        d_over, d_under = (rr[i], s[j]) if over_i else (s[j], rr[i]); sign = np.sign(d_over[0] * d_under[1] - d_over[1] * d_under[0])
                        cross.append((ci, i + ti[k], cj, j + tj[k], over_i, sign))
        if not good_proj: continue
        unders = [[] for _ in comps]
        for ci, si, cj, sj, over_i, sign in cross:
            if over_i: unders[cj].append(sj)
            else: unders[ci].append(si)
        for uu_ in unders: uu_.sort()
        offset = []; narcs = 0
        for uu_ in unders: offset.append(narcs); narcs += max(1, len(uu_))
        def arc_of(c, s__):
            uu_ = unders[c]
            if not uu_: return offset[c]
            return offset[c] + (bisect.bisect_right(uu_, s__) - 1) % len(uu_)
        nc = len(cross)
        if nc == 0: return sp.Integer(1), 0
        Am = sp.zeros(nc, narcs)
        for row, (ci, si, cj, sj, over_i, sign) in enumerate(cross):
            if over_i: co, so, cu, su = ci, si, cj, sj
            else: co, so, cu, su = cj, sj, ci, si
            kk = arc_of(co, so); uu_ = unders[cu]; idx = bisect.bisect_left(uu_, su); out_arc = offset[cu] + idx; in_arc = offset[cu] + (idx - 1) % len(uu_)
            Am[row, kk] += 1 - t
            if sign > 0: Am[row, in_arc] += t; Am[row, out_arc] += -1
            else: Am[row, in_arc] += -1; Am[row, out_arc] += t
        D = sp.expand(Am[:-1, :-1].det(method='berkowitz'))
        if D == 0: return sp.Integer(0), nc
        P_ = sp.Poly(D, t); low = min(mon[0] for mon in P_.monoms()); D = sp.expand(sp.cancel(D / t ** low))
        if sp.Poly(D, t).LC() < 0: D = -D
        return D, nc
    raise RuntimeError('no generic projection')


# ================================================================ 12.1 three flows on one matrix
print('12.1 three flows on one matrix')
REF_PTS = [hopf_point(0.9, 0.4, 1.1), hopf_point(0.5, 2.0, -0.7), hopf_point(1.2, -1.3, 0.2)]
def circle_fit(C):
    """centre, radius and worst departure of a closed curve from the circle through three of its points"""
    A, B, Cc = C[0], C[len(C) // 3], C[2 * len(C) // 3]; u = B - A; v = Cc - A; nrm = np.cross(u, v)
    M = np.array([u, v, nrm]); rhs = np.array([u @ (A + B) / 2, v @ (A + Cc) / 2, nrm @ A]); cen = np.linalg.solve(M, rhs); R = np.linalg.norm(A - cen)
    return cen, R, max(np.abs(np.linalg.norm(C - cen, axis=1) - R).max(), np.abs((C - cen) @ nrm / np.linalg.norm(nrm)).max())
worst = 0; N1 = 240
for (al, be) in REF_PTS:
    for side in 'LR':
        for m in (1, 2, 3): worst = max(worst, circle_fit(chart(0.0, *flow_orbit(side, m, al, be, N1)))[2])
ok('A1 every orbit of the three left flows and the three right flows is a round circle in the chart (NUMERICAL)', worst < 1e-9, 'worst departure %.1e' % worst)
# tangents of the three left flows at a point are mutually orthogonal (the chart is conformal)
def tangent(side, m, al, be, h=1e-6):
    S = PAULI[m]; U = np.array([[al, be], [-np.conj(be), np.conj(al)]])
    def at(ph):
        E = np.cos(ph) * np.eye(2) + 1j * np.sin(ph) * S; V = E @ U if side == 'L' else U @ E; return chart(0.0, V[0, 0], V[0, 1])
    return (at(h) - at(-h)) / (2 * h)
worst = 0
for (al, be) in REF_PTS:
    T = [tangent('L', m, al, be) for m in (1, 2, 3)]; worst = max(worst, max(abs(T[i] @ T[j]) / (np.linalg.norm(T[i]) * np.linalg.norm(T[j])) for i in range(3) for j in range(i)))
ok('A2 the three left flows through a point leave it in three mutually orthogonal directions (NUMERICAL)', worst < 1e-8, 'worst cosine %.1e' % worst)
# the left flows are the lines of v, E, B of the hopfion at t = 0
worst = 0
for (al, be) in REF_PTS:
    X0 = chart(0.0, al, be); F0 = Fpq(0.0, X0); fields = {3: vflow(0.0, X0), 2: F0.real, 1: F0.imag}
    worst = max(worst, max(sinang(tangent('L', m, al, be), fields[m]) for m in (1, 2, 3)))
ok('A3 exp(i phi sigma_3) U follows the flow v, exp(i phi sigma_2) U the electric lines, exp(i phi sigma_1) U the magnetic lines of F = n/d^3 at t = 0 (NUMERICAL)', worst < 1e-7, 'worst sine %.1e' % worst)
lks = {}
for side in 'LR':
    for m in (1, 2, 3):
        vals = []
        for i in range(3):
            for j in range(i):
                C1 = chart(0.0, *flow_orbit(side, m, *REF_PTS[i], N1)); C2 = chart(0.0, *flow_orbit(side, m, *REF_PTS[j], N1))
                vals.append((gauss_poly(C1, C2), lk_cross(C1, C2)[0], gauss_smooth(C1, C2)))
        lks[side + str(m)] = vals
flat = [v for vals in lks.values() for v in vals]
ok('A4 the polygon Gauss integral is an integer and equals the signed crossing count; the smooth Gauss integral agrees (NUMERICAL)',
   max(abs(g - c) for g, c, s in flat) < 1e-9 and max(abs(s - c) for g, c, s in flat) < 1e-6, 'worst |Gauss - crossings| %.1e, smooth %.1e' % (max(abs(g - c) for g, c, s in flat), max(abs(s - c) for g, c, s in flat)))
left = sorted({round(c) for k, vals in lks.items() if k[0] == 'L' for g, c, s in vals}); right = sorted({round(c) for k, vals in lks.items() if k[0] == 'R' for g, c, s in vals})
ok('A5 PREDICTION any two orbits of one left flow link +1, any two orbits of one right flow link -1 (the mirror family)', left == [1] and right == [-1], 'left %s, right %s' % (left, right))
# the right flow with sigma_1 is the Kustaanheimo-Stiefel fibre of gravity_knot_20261021.py, u = (Re alpha, Im alpha, Re beta, Im beta)
KS = lambda u: np.array([u[0] ** 2 - u[1] ** 2 - u[2] ** 2 + u[3] ** 2, 2 * (u[0] * u[1] - u[2] * u[3]), 2 * (u[0] * u[2] + u[1] * u[3])])
al, be = REF_PTS[0]; u0 = np.array([al.real, al.imag, be.real, be.imag]); worst = 0
for ph in (0.3, 1.7, 4.0):
    c, s = math.cos(ph), math.sin(ph); fib = np.array([u0[0] * c - u0[3] * s, u0[1] * c + u0[2] * s, u0[2] * c - u0[1] * s, u0[3] * c + u0[0] * s])
    E = np.cos(ph) * np.eye(2) + 1j * np.sin(ph) * PAULI[1]; V = np.array([[al, be], [-np.conj(be), np.conj(al)]]) @ E
    worst = max(worst, np.abs(fib - np.array([V[0, 0].real, V[0, 0].imag, V[0, 1].real, V[0, 1].imag])).max(), np.abs(KS(fib) - KS(u0)).max())
ok('A6 U exp(i phi sigma_1) is the Kustaanheimo-Stiefel fibre u -> u e^{k phi}, and the map KS(u) does not see it (NUMERICAL)', worst < 1e-12, 'worst %.1e' % worst)
C1 = chart(0.0, *flow_orbit('L', 1, *REF_PTS[0], N1)); C2 = chart(0.0, *flow_orbit('L', 1, *REF_PTS[1], N1)); cen, Rc, _ = circle_fit(C1)
NUM['flows'] = {'points': [[cj(a), cj(b)] for a, b in REF_PTS], 'hopf': [[0.9, 0.4, 1.1], [0.5, 2.0, -0.7], [1.2, -1.3, 0.2]], 'N': N1, 'chart0': fl(chart(0.0, *REF_PTS[0])),
                'circle_L1_p0': {'centre': fl(cen), 'radius': float(Rc)}, 'lk_left': 1, 'lk_right': -1, 'gauss_L1_01': gauss_poly(C1, C2), 'crossings_L1_01': lk_cross(C1, C2)[1],
                'orbit_L1_p0_at_1': fl(chart(0.0, *[v[40] for v in flow_orbit('L', 1, *REF_PTS[0], N1)])), 'orbit_R1_p0_at_1': fl(chart(0.0, *[v[40] for v in flow_orbit('R', 1, *REF_PTS[0], N1)]))}

# ================================================================ 12.2 Kepler's orbit is a square
print('12.2 Kepler\'s orbit is a square')
def kepler(aw, bw, psi):
    Rm = np.array([[math.cos(psi), -math.sin(psi)], [math.sin(psi), math.cos(psi)]]); M = Rm @ np.diag([aw, bw]); G = M @ M.T
    a = (aw * aw + bw * bw) / 2; b = aw * bw; c = (aw * aw - bw * bw) / 2
    return dict(M=M, G=G, a=a, b=b, c=c, e=c / a, d=math.acosh(a / b), psi=psi, aw=aw, bw=bw)
K = kepler(1.2, 0.6, 0.35); s = np.linspace(0, 2 * np.pi, 400, endpoint=False); w = K['M'] @ np.array([np.cos(s), np.sin(s)]); wz = w[0] + 1j * w[1]; zz = wz ** 2
focus2 = 2 * K['c'] * np.exp(2j * K['psi'])
ok('B1 z = w^2 is an ellipse with one focus at 0 and the other at 2c e^{2i psi}: |z| + |z - 2c e^{2i psi}| = 2a (NUMERICAL)', np.abs(np.abs(zz) + np.abs(zz - focus2) - 2 * K['a']).max() < 1e-12, 'a = %.4f, b = %.4f, c = %.4f, e = %.4f' % (K['a'], K['b'], K['c'], K['e']))
ev = np.linalg.eigvalsh(K['G']); trG = np.trace(K['G']); detG = np.linalg.det(K['G'])
ok('B2 G = M M^T: tr G = 2a, det G = b^2, tr^2 - 4 det = (2c)^2, eigenvalues = perihelion and aphelion distances a -+ c (NUMERICAL)',
   max(abs(trG - 2 * K['a']), abs(detG - K['b'] ** 2), abs(trG ** 2 - 4 * detG - 4 * K['c'] ** 2), abs(ev[0] - (K['a'] - K['c'])), abs(ev[1] - (K['a'] + K['c'])), abs(np.abs(zz).min() - ev[0]), abs(np.abs(zz).max() - ev[1])) < 1e-12)
xh, th = K['G'][0, 1] / K['G'][1, 1], math.sqrt(detG) / K['G'][1, 1]; dh = math.acosh(1 + (xh ** 2 + (th - 1) ** 2) / (2 * th))
ok('B3 e = tanh d, d the hyperbolic distance from i to the point of G/sqrt(det G); the oscillator\'s axis ratio is e^d and a^2 = b^2 + c^2 reads E^2 = m^2 + p^2 with rapidity d (NUMERICAL)',
   max(abs(math.tanh(dh) - K['e']), abs(dh - K['d']), abs(math.exp(dh) - K['aw'] / K['bw']), abs(K['a'] - K['b'] * math.cosh(dh)), abs(K['c'] - K['b'] * math.sinh(dh))) < 1e-12, 'd = %.6f = ln 2: %s' % (dh, abs(dh - math.log(2)) < 1e-12))
# dynamics with GM = 1: start at the aphelion, compare Newton's orbit with z(s) at t(s) = a^{3/2}(2s + e sin 2s)
a_, e_, b_ = K['a'], K['e'], K['b']; z0 = (a_ + K['c']) * np.exp(2j * K['psi']); v0 = 1j * np.exp(2j * K['psi']) * (b_ / math.sqrt(a_)) / (a_ + K['c'])
ss = np.array([0.4, 1.1, 2.0, 3.0]); ts = a_ ** 1.5 * (2 * ss + e_ * np.sin(2 * ss))
sol = solve_ivp(lambda tt, yv: [yv[2], yv[3], -yv[0] / (yv[0] ** 2 + yv[1] ** 2) ** 1.5, -yv[1] / (yv[0] ** 2 + yv[1] ** 2) ** 1.5], (0, ts[-1]), [z0.real, z0.imag, v0.real, v0.imag], t_eval=ts, method='DOP853', rtol=1e-12, atol=1e-13)
zs = (K['M'] @ np.array([np.cos(ss), np.sin(ss)])); zs = (zs[0] + 1j * zs[1]) ** 2
ok('B4 the square of the oscillator is Newton\'s orbit: with dt = 2 sqrt(a) r ds, t(s) = a^{3/2}(2s + e sin 2s), Kepler\'s equation from the aphelion; energy -1/tr G, L^2 = 2 det G/tr G (NUMERICAL, GM = 1)',
   np.abs(sol.y[0] + 1j * sol.y[1] - zs).max() < 1e-9 and abs(abs(v0) ** 2 / 2 - 1 / abs(z0) + 1 / trG) < 1e-12 and abs((z0.conjugate() * v0).imag ** 2 - 2 * detG / trG) < 1e-12, 'worst position error %.1e; period %.6f' % (np.abs(sol.y[0] + 1j * sol.y[1] - zs).max(), 2 * math.pi * a_ ** 1.5))
J = K['M'] @ np.array([1, 1j]); S = np.array([abs(J[0]) ** 2 + abs(J[1]) ** 2, abs(J[0]) ** 2 - abs(J[1]) ** 2, 2 * (J[0] * np.conj(J[1])).real, 2 * (np.conj(J[0]) * J[1]).imag])
chi = math.atan2(K['bw'], K['aw']); Mj = np.array([[math.cos(K['psi']) * math.cos(chi), -math.sin(K['psi']) * math.sin(chi)], [math.sin(K['psi']) * math.cos(chi), math.cos(K['psi']) * math.sin(chi)]])
ok('B5 the Stokes vector of the oscillator\'s ellipse (10.5, Jones vector M(1, i)^T) is S/S_0 = (e cos 2psi, e sin 2psi, sqrt(1 - e^2)): the eccentricity vector, pointing to the aphelion, and L/L_circ; M/sqrt(tr G) is 10.5\'s jonesOf(psi, chi) with cos 2chi = e (NUMERICAL)',
   max(abs(S[0] - trG), abs(S[1] / S[0] - e_ * math.cos(2 * K['psi'])), abs(S[2] / S[0] - e_ * math.sin(2 * K['psi'])), abs(S[3] / S[0] - math.sqrt(1 - e_ ** 2)), abs(math.cos(2 * chi) - e_), np.abs(K['M'] / math.sqrt(trG) - Mj).max(), abs(np.angle(z0) - 2 * K['psi'])) < 1e-12, 'S/S0 = (%.4f, %.4f, %.4f)' % tuple(S[1:] / S[0]))
uu = lambda sv: np.array([(K['M'] @ [math.cos(sv), math.sin(sv)])[0], (K['M'] @ [math.cos(sv), math.sin(sv)])[1], 0.0, 0.0])
ok('B6 one revolution of the planet is s -> s + pi and sends the oscillator to its negative: u(s + pi) = -u(s) = u(s) e^{k pi}, half a turn of the hidden circle, with KS(u) unchanged (NUMERICAL)',
   max(np.abs(uu(0.7 + math.pi) + uu(0.7)).max(), np.abs(KS(uu(0.7 + math.pi)) - KS(uu(0.7))).max(), abs(KS(uu(0.7))[0] + 1j * KS(uu(0.7))[1] - (uu(0.7)[0] + 1j * uu(0.7)[1]) ** 2)) < 1e-12)
NUM['kepler'] = {'in': [1.2, 0.6, 0.35], 'M': fl(K['M']), 'G': fl(K['G']), 'a': K['a'], 'b': K['b'], 'c': K['c'], 'e': K['e'], 'd': K['d'], 'rApo': float(ev[1]), 'rPeri': float(ev[0]), 'energy': -1 / trG, 'L': b_ / math.sqrt(a_), 'period': 2 * math.pi * a_ ** 1.5,
                 'stokes': fl(S / S[0]), 'chi': chi, 'point': [xh, th], 'times': fl(ts), 's': fl(ss), 'z': [cj(v) for v in zs], 'aphelion': cj(z0)}

