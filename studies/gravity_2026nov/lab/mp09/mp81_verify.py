#!/usr/bin/env python3
"""The Modular Plane, version 8.1, Route A: T-G5, the 2 x 2 matrix of field lines round a core line (instrument 12.8).

    python3 -I mp81_verify.py        ->  mp81_numbers.json next to this file

Needs mp8_verify_part1.py (its helpers are reused, nothing of it is run) and mp8_numbers.json (MP08's writhes) beside it.
Pre-registered in MP09_Version81_Plan.md (written before this script was first run). Sections P1 to P7 work in space, on the field
B = Im F_pq itself; the closed forms they are compared with are the plan's. A failed check prints MISS and is counted."""
import json, math, pathlib
import numpy as np, sympy as sp
from scipy.integrate import solve_ivp, quad
from scipy.optimize import brentq
HERE = pathlib.Path(__file__).resolve().parent
exec(compile((HERE / 'mp8_verify_part1.py').read_text().split('# ================================================================ 12.1')[0], 'mp8_helpers', 'exec'))
V8 = json.loads((HERE / 'mp8_numbers.json').read_text())
PQ = [(1, 1), (2, 3), (3, 2), (2, 5)]; TV = 0.8
is_sq = lambda p, q: math.isqrt(2 * p * q * (p + q)) ** 2 == 2 * p * q * (p + q)
tol = lambda p, q: 1e-4 if is_sq(p, q) else 1e-6          # near the identity the angle is the square root of the trace's error (plan, amendment before the run)
nu_of = lambda p, q: math.sqrt(2 * p * q * (p + q))
hmax_of = lambda p, q: (p / (p + q)) ** (p / 2) * (q / (p + q)) ** (q / 2)

# ---------------------------------------------------------------- the plan's closed forms
def lw(p, q, u): return p * u + q * np.log1p(-np.exp(u))                      # ln(x^p (1 - x)^q) at x = e^u
def levels(p, q, h):
    """u = ln x at the two roots of x^p (1 - x)^q = h^2"""
    l2h = 2 * math.log(h); u1 = brentq(lambda u: lw(p, q, u) - l2h, -700.0, math.log(p / (p + q)))
    v2 = brentq(lambda v: p * math.log1p(-math.exp(v)) + q * v - l2h, -700.0, math.log(q / (p + q)))     # v = ln(1 - x)
    return u1, math.log1p(-math.exp(v2))
def iota(p, q, h):
    u1, u2 = levels(p, q, h); mid, half = (u1 + u2) / 2, (u2 - u1) / 2; l2h = 2 * math.log(h)
    dmax = lw(p, q, math.log(p / (p + q))) - l2h; slope = lambda ue: abs(p + q * math.exp(ue) / math.expm1(ue))          # |d lambda/du| at an end
    def f(th):
        d = lw(p, q, mid + half * math.sin(th)) - l2h
        if d < 1e-7 * dmax: d = slope(u2 if th > 0 else u1) * half * 2 * math.sin((math.pi / 2 - abs(th)) / 2) ** 2      # the last sliver by its tangent, free of cancellation
        return half * math.cos(th) / math.sqrt(math.expm1(d))
    return 2 * math.pi * q / quad(f, -math.pi / 2, math.pi / 2, limit=800, epsabs=1e-13, epsrel=1e-13)[0]
def monodromy(p, q):
    nu = nu_of(p, q); c, s, k = math.cos(2 * math.pi * nu), math.sin(2 * math.pi * nu), math.sqrt(2 * (p + q)); return np.array([[c, -s / k], [k * s, c]])

# ---------------------------------------------------------------- E: the derivation, exactly
print('E the derivation')
eta, x1, x2 = sp.symbols('eta xi_1 xi_2', real=True); Ps, Qs = sp.symbols('p q', positive=True); crd = (eta, x1, x2)
Th = Ps * x1 + Qs * x2; Xm = sp.Matrix([-sp.sin(Th) * sp.sin(eta) * sp.cos(eta), sp.cos(Th) * sp.sin(eta) ** 2, -sp.cos(Th) * sp.cos(eta) ** 2])
Hs = sp.cos(eta) ** Ps * sp.sin(eta) ** Qs * sp.cos(Th)
ok('E1 H = cos^p(eta) sin^q(eta) cos(Theta) is constant along the plan\'s line field, for symbolic p, q (EXACT)', sp.simplify(sum(Xm[i] * sp.diff(Hs, crd[i]) for i in range(3)) / Hs) == 0)
worst = 0
for (p, q) in PQ + [(1, 8)]:
    f = sp.cos(eta) ** p * sp.exp(sp.I * p * x1); g = sp.sin(eta) ** q * sp.exp(sp.I * q * x2)
    Om = lambda a, b: sp.im(sp.diff(f, crd[a]) * sp.diff(g, crd[b]) - sp.diff(f, crd[b]) * sp.diff(g, crd[a]))
    ker = [Om(1, 2), -Om(0, 2), Om(0, 1)]; fac = p * q * sp.cos(eta) ** (p - 1) * sp.sin(eta) ** (q - 1)
    for pt in ((0.7, 0.4, 1.9), (1.1, -2.2, 0.3), (0.3, 2.9, -1.4)):
        sub = {eta: sp.Float(pt[0], 30), x1: sp.Float(pt[1], 30), x2: sp.Float(pt[2], 30), Ps: p, Qs: q}
        worst = max(worst, max(abs(sp.N(ker[i].subs(sub), 30) - sp.N((fac * Xm[i]).subs(sub), 30)) for i in range(3)))
ok('E2 the kernel of Im d(alpha^p) ^ d(beta^q) on the label sphere is pq cos^{p-1} sin^{q-1} times the plan\'s line field, five pairs, three points, 30 digits', worst < 1e-25, 'worst %.1e' % worst)
red = sp.Matrix([-sp.sin(Th) * sp.sin(eta) * sp.cos(eta), sp.cos(Th) * (Ps * sp.sin(eta) ** 2 - Qs * sp.cos(eta) ** 2)]); Ths = sp.Symbol('Theta', real=True)
Jr = red.subs(Th, Ths).jacobian([eta, Ths]).subs({Ths: 0, eta: sp.atan(sp.sqrt(Qs / Ps))}).applyfunc(sp.simplify); ss = sp.sqrt(Ps * Qs) / (Ps + Qs)
ok('E3 linearised at the core line: L = [[0, -s], [2(p+q)s, 0]], s = sqrt(pq)/(p+q), det L = 2pq/(p+q), and (2 pi (p+q))^2 det L = (2 pi)^2 2pq(p+q) (EXACT)',
   (Jr - sp.Matrix([[0, -ss], [2 * (Ps + Qs) * ss, 0]])).applyfunc(sp.simplify) == sp.zeros(2) and sp.simplify(Jr.det() - 2 * Ps * Qs / (Ps + Qs)) == 0)

# ---------------------------------------------------------------- the space side
def bun(tv, X, p, q):
    B = Fpq(tv, X, p, q).imag; return B / np.linalg.norm(B, axis=-1)[..., None]
def rhs_var(tv, p, q, hh=1e-5):
    E = np.eye(3) * hh
    def f(s, y):
        X = y[:3]; b = bun(tv, np.vstack([X, X + E, X - E]), p, q); J = ((b[1:4] - b[4:7]) / (2 * hh)).T
        return np.concatenate([b[0], (J @ y[3:].reshape(3, 3)).ravel()])
    return f
def circuit(p, q, tv, psi0=0.0):
    C = chart(tv, *core(p, q, psi0, 4000)); L0 = np.linalg.norm(np.roll(C, -1, 0) - C, axis=1).sum(); X0 = C[0]; b0 = bun(tv, X0, p, q)
    sol = solve_ivp(rhs_var(tv, p, q), (0, 1.03 * L0), np.concatenate([X0, np.eye(3).ravel()]), method='DOP853', rtol=1e-10, atol=1e-12, dense_output=True)
    s1 = brentq(lambda s: (sol.sol(s)[:3] - X0) @ b0, 0.97 * L0, 1.03 * L0, xtol=1e-13)
    return sol, s1, X0, (1.0 if (C[1] - C[0]) @ b0 > 0 else -1.0)
def torus_frame_count(p, q, tv, sol, s1, K=6000, J=720, eps=1e-6, N=20000):
    """the solutions of the variational equation in the closed torus frame (U toward the neighbour on the same torus, V = T x U): the return map and its rotation number"""
    ss_ = np.linspace(0, s1, K + 1); Y = sol.sol(ss_); X = Y[:3].T; PHI = Y[3:].T.reshape(K + 1, 3, 3); T = bun(tv, X, p, q)
    al, be = labels(tv, X); U = chart(tv, al * np.exp(1j * eps), be) - X; U -= np.sum(U * T, 1)[:, None] * T; U /= np.linalg.norm(U, axis=1)[:, None]; V = np.cross(T, U)
    A1 = PHI @ U[0]; A2 = PHI @ V[0]; c = np.array([[np.sum(A1 * U, 1), np.sum(A1 * V, 1)], [np.sum(A2 * U, 1), np.sum(A2 * V, 1)]])
    a = np.linspace(0, 2 * np.pi, J, endpoint=False)[:, None]
    ang = np.unwrap(np.arctan2(np.cos(a) * c[0, 1] + np.sin(a) * c[1, 1], np.cos(a) * c[0, 0] + np.sin(a) * c[1, 0]), axis=1); th1 = ang[:, -1] - ang[:, 0]
    x = 0.3; tot = 0.0
    for _ in range(N):
        dlt = float(np.interp(x % (2 * np.pi), a[:, 0], th1, period=2 * np.pi)); tot += dlt; x += dlt
    return tot / (2 * np.pi * N), np.array([[c[0, 0, -1], c[1, 0, -1]], [c[0, 1, -1], c[1, 1, -1]]])

print('P1 the conserved function')
worst = 0
for (p, q) in PQ:
    for tv in (0.0, TV):
        so = solve_ivp(lambda s, X: bun(tv, X, p, q), (0, 12.0), [0.6, -0.3, 0.5], method='DOP853', rtol=1e-11, atol=1e-13, dense_output=True)
        al, be = labels(tv, so.sol(np.linspace(0, 12, 400)).T); Hh = (al ** p * be ** q).real; worst = max(worst, np.abs(Hh - Hh[0]).max() / hmax_of(p, q))
ok('P1 PREDICTION H = Re(alpha^p beta^q) is constant along a magnetic line integrated in space, t = 0 and t = 0.8, four pairs (1e-8 of h_max)', worst < 1e-8, 'worst drift %.1e' % worst)

print('P2 to P4 round the core line')
rows = {}; okP2 = okP3 = okP4 = True; det_w = clo_w = 0
for (p, q) in PQ + [(1, 8)]:
    for tv in (0.0, TV):
        sol, s1, X0, sgn = circuit(p, q, tv); Phi = sol.sol(s1)[3:].reshape(3, 3); clo = float(np.linalg.norm(sol.sol(s1)[:3] - X0)); trp = float(np.trace(Phi) - 1)
        nu_avg, M2 = torus_frame_count(p, q, tv, sol, s1); the = math.acos(max(-1.0, min(1.0, trp / 2))) / (2 * math.pi); m0 = round(nu_avg)
        nu = min((m + sg * the for m in (m0 - 1, m0, m0 + 1) for sg in (1, -1)), key=lambda v: abs(v - nu_avg))
        lkt = lk_cross(chart(tv, *core(p, q, 0.0, 480)), chart(tv, *core(p, q, 0.3, 480)))[0]
        wr = (math.sqrt(p * q) - p * q if tv == 0 else V8['carried']['writhe_smooth'].get('%d,%d' % (p, q))) if (p, q) in PQ else None
        rows[(p, q, tv)] = dict(p=p, q=q, t=tv, length=float(s1), closure=clo, det=float(np.linalg.det(Phi)), tr=trp, tr_frame=float(np.trace(M2)), det_frame=float(np.linalg.det(M2)), nu_avg=float(nu_avg), nu=float(nu), lk_torus=float(lkt),
                                lk_phys=float(nu + lkt), twist=(float(nu + lkt - wr) if wr is not None else None), sense=sgn, M=fl(M2))
        det_w = max(det_w, abs(np.linalg.det(Phi) - 1), abs(np.linalg.det(M2) - 1)); clo_w = max(clo_w, clo)
        okP2 &= abs(trp - 2 * math.cos(2 * math.pi * nu_of(p, q))) < 1e-6 and abs(np.trace(M2) - trp) < 1e-6
        okP3 &= abs(nu - nu_of(p, q)) < tol(p, q) and abs(nu_avg - nu_of(p, q)) < 2e-3
        okP4 &= lkt == -p * q
        print('     (%d,%d) t = %.1f: length %.4f, closure %.1e, tr %.7f (plan %.7f), nu %.6f (mean %.4f; plan %.6f), torus frame %+d, against Seifert %.6f%s' % (
            p, q, tv, s1, clo, trp, 2 * math.cos(2 * math.pi * nu_of(p, q)), nu, nu_avg, nu_of(p, q), lkt, nu + lkt, '' if wr is None else ', twist %.6f' % (nu + lkt - wr)))
ok('P2 PREDICTION the return map round the core line has det 1 and trace 2 cos(2 pi sqrt(2pq(p+q))): the identity for (1,1), elliptic for (2,3), (3,2), (2,5); same at t = 0 and 0.8 (1e-6)', okP2 and det_w < 1e-6, 'worst |det - 1| %.1e, worst closure %.1e' % (det_w, clo_w))
ok('P3 PREDICTION the rotation number in the closed torus frame is +sqrt(2pq(p+q)): 2, 7.745967, 7.745967, 11.832160 (counted continuously; fraction from the eigenvalue, 1e-6)', okP3)
tw_plan = {(1, 1, 0.0): 1.0, (2, 3, 0.0): 5.296477, (3, 2, 0.0): 5.296477, (2, 5, 0.0): 8.669882, (1, 1, TV): 1.176192, (2, 3, TV): 5.751232, (3, 2, TV): 5.703596, (2, 5, TV): 9.289322}
ok('P4 PREDICTION the torus frame links -pq; winding against Seifert nu - pq = 1, 1.745967, 1.745967, 1.832160 at both times; twist nu - pq - Wr as pre-registered (2e-6)',
   okP4 and all(abs(rows[k]['lk_phys'] - (nu_of(k[0], k[1]) - k[0] * k[1])) < tol(k[0], k[1]) for k in tw_plan) and all(abs(rows[k]['twist'] - v) < 2 * tol(k[0], k[1]) for k, v in tw_plan.items()),
   'twist: ' + ', '.join('(%d,%d) t=%.1f %.6f' % (k[0], k[1], k[2], rows[k]['twist']) for k in tw_plan))

print('P5 the profile')
def line_returns(p, q, h, nret, tau_max=400.0):
    """a magnetic line of the torus H = h at t = 0, integrated in space at the speed (1 + |X|^2)/2 (arc length on the label sphere); its returns to Theta = 0 on the outer side"""
    u1, u2 = levels(p, q, h); xs = math.exp(u2); X0 = chart(0.0, np.array(math.sqrt(xs) + 0j), np.array(math.sqrt(1 - xs) + 0j))
    f = lambda s, X: bun(0.0, X, p, q) * (1 + X @ X) / 2
    def ev(s, X): al, be = labels(0.0, X); return float((al ** p * be ** q).imag)
    ev.direction = -1 if ev(0, X0 + 1e-6 * f(0, X0)) < 0 else 1; got = []; t0 = 0.0; Xc = X0; sols = []
    while len(got) < nret and t0 < tau_max:
        so = solve_ivp(f, (t0, t0 + 4.0), Xc, method='DOP853', rtol=1e-11, atol=1e-13, events=ev, dense_output=True)
        got += [(te, ye) for te, ye in zip(so.t_events[0], so.y_events[0]) if te > 1e-6]; sols.append(so); t0 = so.t[-1]; Xc = so.y[:, -1]
    return X0, got[:nret], sols
def dphi(p, q, X0, got, sols):
    """the advance of phi between successive returns, from the unwrapped phase of beta (phi = -xi_2/p), exact at the returns"""
    ts = np.concatenate([np.linspace(so.t[0], so.t[-1], 800, endpoint=False) for so in sols]); Xs = np.concatenate([so.sol(np.linspace(so.t[0], so.t[-1], 800, endpoint=False)).T for so in sols])
    x2u = np.unwrap(np.angle(labels(0.0, Xs)[1])); out = [float(np.angle(labels(0.0, X0)[1]))]
    for te, ye in got:
        w = float(np.angle(labels(0.0, ye)[1])); out.append(w + 2 * np.pi * round((np.interp(te, ts, x2u) - w) / (2 * np.pi)))
    return -np.diff(out) / p
prof = {}; worst = 0; hop = 0
for (p, q) in PQ:
    for lam in (0.9, 0.5, 0.1):
        h = lam * hmax_of(p, q); X0, got, sols = line_returns(p, q, h, 3); io = 2 * np.pi / np.abs(dphi(p, q, X0, got, sols)); prof['%d,%d,%.1f' % (p, q, lam)] = dict(space=float(io.mean()), formula=iota(p, q, h))
        worst = max(worst, np.abs(io - iota(p, q, h)).max());  hop = max(hop, abs(io.mean() - 2)) if (p, q) == (1, 1) else hop
edge = {'%d,%d' % (p, q): iota(p, q, 1e-6 * hmax_of(p, q)) for p, q in PQ}; axis = {'%d,%d' % (p, q): iota(p, q, (1 - 1e-6) * hmax_of(p, q)) for p, q in PQ}
mono = all(all(np.diff([iota(p, q, lam * hmax_of(p, q)) for lam in np.linspace(0.02, 0.999, 40)]) < 1e-12) for p, q in PQ[1:])
ok('P5 PREDICTION iota(h) from lines integrated in space at h/h_max = 0.9, 0.5, 0.1 equals the plan\'s integral (1e-6); 2 on every torus of the hopfion; within 2 % of 2pq at 1e-6 h_max; -> nu on the axis',
   worst < 1e-6 and hop < 1e-6 and all(abs(edge['%d,%d' % (p, q)] / (2 * p * q) - 1) < 0.02 for p, q in PQ) and all(abs(axis['%d,%d' % (p, q)] - nu_of(p, q)) < 1e-5 for p, q in PQ),
   'worst %.1e; edge %s; falls from the edge to the axis without a turn: %s' % (worst, {k: round(v, 4) for k, v in edge.items()}, mono))

print('P6 when the neighbours close')
sq = [(p, q) for p in range(1, 61) for q in range(p, 61) if math.gcd(p, q) == 1 and math.isqrt(2 * p * q * (p + q)) ** 2 == 2 * p * q * (p + q)]
pell = [(xx - 1) // 2 for xx, mm in ((3, 2), (17, 12), (99, 70), (577, 408)) if xx * xx - 2 * mm * mm == 1 and 2 * ((xx - 1) // 2) * ((xx - 1) // 2 + 1) == mm * mm]
r18 = rows[(1, 8, 0.0)]
ok('P6 PREDICTION 2pq(p+q) is a perfect square, coprime p <= q <= 60, only for (1,1), (1,8), (1,49), (32,49); p = 1 is Pell\'s (2q+1)^2 - 2m^2 = 1: q = 1, 8, 49, 288; for (1,8) the return map from the field is the identity with nu = 12',
   sq == [(1, 1), (1, 8), (1, 49), (32, 49)] and pell == [1, 8, 49, 288] and abs(r18['tr'] - 2) < 1e-6 and abs(r18['nu_avg'] - 12) < 2e-3 and np.abs(np.array(r18['M']).reshape(2, 2) - np.eye(2)).max() < 1e-4, 'found %s; (1,8): tr %.8f, M - I %.1e, against Seifert %.4f' % (sq, r18['tr'], np.abs(np.array(r18['M']).reshape(2, 2) - np.eye(2)).max(), r18['nu_avg'] - 8))

print('P7 a closed line off the axis')
p, q = 2, 3; h8 = brentq(lambda h: iota(p, q, h) - 8, 1e-4 * hmax_of(p, q), 0.9999 * hmax_of(p, q)); X0, got, sols = line_returns(p, q, h8, 8); te, ye = got[-1]; clo8 = float(np.linalg.norm(ye - X0))
tt = np.linspace(0, te, 1600, endpoint=False); bounds = [so.t[0] for so in sols] + [sols[-1].t[-1]]; idx = np.clip(np.searchsorted(bounds, tt, side='right') - 1, 0, len(sols) - 1)
Cl = np.array([sols[i].sol(tv_) for i, tv_ in zip(idx, tt)]); Cc = chart(0.0, *core(p, q, 0.0, 800)); sense = rows[(2, 3, 0.0)]['sense']; lk8, nc8 = lk_cross(Cl, Cc); g8 = gauss_poly(Cl[::2], Cc[::2])
ok('P7 PREDICTION on the torus iota = 8 of the (2,3) field a magnetic line closes after one circuit (1e-6) and links the core line +2, both along B (integer by crossings = Gauss integral)',
   clo8 < 1e-6 and lk8 * sense == 2 and abs(g8 - lk8) < 1e-6, 'h/h_max = %.6f, closure %.1e, linking %+d (%d crossings), Gauss %.8f' % (h8 / hmax_of(p, q), clo8, int(lk8 * sense), nc8, g8 * sense))

NUM = {'pairs': [[p, q] for p, q in PQ + [(1, 8)]], 'nu': {'%d,%d' % (p, q): nu_of(p, q) for p, q in PQ + [(1, 8)]}, 'hmax': {'%d,%d' % (p, q): hmax_of(p, q) for p, q in PQ + [(1, 8)]},
       'M': {'%d,%d' % (p, q): fl(monodromy(p, q)) for p, q in PQ + [(1, 8)]}, 'rows': list(rows.values()), 'profile': prof, 'edge': edge,
       'iota': {'%d,%d' % (p, q): [[lam, iota(p, q, lam * hmax_of(p, q))] for lam in (0.05, 0.2, 0.5, 0.8, 0.95, 0.99)] for p, q in PQ + [(1, 8)]},
       'levels': {'%d,%d' % (p, q): [math.exp(v) for v in levels(p, q, 0.5 * hmax_of(p, q))] for p, q in PQ + [(1, 8)]}, 'squares': [list(v) for v in sq], 'pell': pell,
       'closed': {'p': 2, 'q': 3, 'iota': 8, 'level': h8 / hmax_of(2, 3), 'closure': clo8, 'linking': int(lk8 * sense)}, 'twist_plan': {'%d,%d,%.1f' % k: v for k, v in tw_plan.items()}}
out = HERE / 'mp81_numbers.json'; out.write_text(json.dumps(NUM, separators=(',', ':'), default=lambda o: o.item() if hasattr(o, 'item') else str(o)))
print('\n%d checks passed, %d failed; numbers written to %s' % (NPASS, NFAIL, out.name)); raise SystemExit(1 if NFAIL else 0)
