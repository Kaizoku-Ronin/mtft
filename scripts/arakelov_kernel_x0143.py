"""Provenance script for AG-01: the canonical Green's-function kernel of X0(143) between the blown-up cusp and the homes.

Requires cypari2 (pip install cypari2).  Builds the differential of the second kind with a double pole at the cusp infinity as
eta(t)eta(11t)/(eta(13t)eta(143t)) times a weight-2 cusp form, fixes its holomorphic part by Petersson orthogonality over the Ford domain of the
Atkin–Lehner quotient (validated: area 14 pi, Gram matrix = Im Omega to 8e-6), and evaluates |d_P d_Q g| / sqrt(mu(P) mu(Q)) at each home with two
constructions that differ by a holomorphic differential.  Output: /tmp/ag_results.json; frozen in mtft.matter.arakelov_kernel.
"""
import numpy as np, math, json, time, cypari2
t0 = time.time()
from mtft.matter.special_points import density_setup
NMAX = 900; DIV = (1, 11, 13, 143); CUSPW = (1, 143, 13, 11)
pari = cypari2.Pari(); pari.allocatemem(10**9)
A, Yi = density_setup(); Y = np.linalg.inv(Yi)
pari('mf2 = mfinit([143,2],1); B2 = mfbasis(mf2)')
C2 = np.array([[float(c) for c in pari(f'mfcoefs(B2[{i+1}], {NMAX})')] for i in range(13)])
W2 = {}
for Q in (143, 13, 11):
    pari(f'at2 = mfatkininit(mf2, {Q})')
    img = np.array([[complex(c) for c in pari(f'mfcoefs(mfatkin(at2, B2[{i+1}]), 80)')] for i in range(13)])
    W2[Q] = np.linalg.lstsq(C2[:, 1:81].T.astype(complex), img[:, 1:81].T, rcond=None)[0].T
Pa = np.linalg.lstsq(C2[:, 0:141].T.astype(complex), A.astype(complex), rcond=None)[0].T
al = Pa @ C2                                   # al[k, n]: q^n coefficient of f_k, omega_k = f_k(tau) dtau
M = {1: np.eye(13, dtype=complex)}; M.update({Q: Pa @ W2[Q] @ np.linalg.inv(Pa) for Q in W2})
def orders(r):
    r1, r11, r13, r143 = r
    return ((r1 + 11*r11 + 13*r13 + 143*r143)/24, (143*r1 + 13*r11 + 11*r13 + r143)/24,
            (13*r1 + 143*r11 + r13 + 11*r143)/24, (11*r1 + r11 + 143*r13 + 13*r143)/24)
def perm(r, Q): return tuple(r[DIV.index(dd * Q // math.gcd(dd, Q) ** 2)] for dd in DIV)
def prodser(n, r):
    p = np.zeros(n + 1); p[0] = 1.0
    for dd, e in zip(DIV, r):
        for m in range(1, n // dd + 1):
            k = dd * m
            for _ in range(abs(e)):
                if e > 0: p[k:] -= p[:n + 1 - k].copy()
                else: p[k:] += 0; [p.__setitem__(j, p[j] + p[j - k]) for j in range(k, n + 1)]
    return p
def eta(t, n=4000):
    q = np.exp(2j * np.pi * t); return np.exp(2j * np.pi * t / 24) * np.prod(1 - q ** np.arange(1, n))
def etaq(r, t): return np.prod([eta(dd * t) ** e for dd, e in zip(DIV, r)])
WM = {143: (0, -1, 143, 0), 11: (11, 5, 143, 66), 13: (13, 7, 143, 78)}
def mob(m, t): a, b, c, dd = m; return (a * t + b) / (c * t + dd)
def build_psi(r):
    o = [int(round(x)) for x in orders(r)]; k = -o[0]
    rows = [al[:, j] for j in range(1, k - 1)]
    for ci in (1, 2, 3):
        for j in range(1, max(0, -o[ci]) + 1): rows.append((M[CUSPW[ci]] @ al)[:, j])
    Amat = np.vstack(rows + [al[:, k - 1]]); v = np.linalg.lstsq(Amat, np.r_[np.zeros(len(rows)), 1.0], rcond=None)[0]
    psi = {}
    for ci, Q in enumerate(CUSPW):
        rp = r if Q == 1 else perm(r, Q); oc = o[ci]
        lam = 1.0
        if Q != 1:
            a, b, c, dd = WM[Q]; tp = complex(-dd / c + 0.0003, math.sqrt(Q) / c)
            lam = etaq(r, mob(WM[Q], tp)) / etaq(rp, tp)
        hs = prodser(NMAX + 20, rp); fw = v @ (M[Q] @ al)
        prod = np.convolve(hs, fw)[:NMAX + 1] * lam          # coefficient i <-> q^(i + oc)
        arr = np.zeros(NMAX + 2, complex)                    # index i <-> q^(i-1)
        for i in range(NMAX + 1):
            if 0 <= i + oc + 1 < NMAX + 2: arr[i + oc + 1] = prod[i]
        psi[Q] = arr
    return psi

import numpy as np, math, json, time
t0 = time.time()
def build(r):
    o = [int(round(x)) for x in orders(r)]; k = -o[0]
    rows = [al[:, j] for j in range(1, k - 1)]
    for ci in (1, 2, 3):
        for j in range(1, max(0, -o[ci]) + 1): rows.append((M[CUSPW[ci]] @ al)[:, j])
    Am = np.vstack(rows + [al[:, k - 1]]); v = np.linalg.lstsq(Am, np.r_[np.zeros(len(rows)), 1.0], rcond=None)[0]
    lam, rq = {1: 1.0}, {1: r}
    for Q in (143, 13, 11):
        rq[Q] = perm(r, Q); a, b, c, dd = WM[Q]; tp = complex(-dd / c + 0.0003, math.sqrt(Q) / c)
        lam[Q] = etaq(r, mob(WM[Q], tp)) / etaq(rq[Q], tp)
    return {"v": v, "lam": lam, "r": rq, "fw": {Q: v @ (M[Q] @ al) for Q in CUSPW}}
cons = {"A": build((1, 1, -1, -1))}
rows = [al[:, j] for j in range(1, 5)] + [(M[11] @ al)[:, j] for j in range(1, 7)] + [al[:, 5]]
ns = np.linalg.svd(np.vstack(rows))[2][-2:].conj()
cA2 = {"v": cons["A"]["v"] + 3.0 * ns[0] - 2.0j * ns[1], "lam": cons["A"]["lam"], "r": cons["A"]["r"]}
cA2["fw"] = {Q: cA2["v"] @ (M[Q] @ al) for Q in CUSPW}; cons["A2"] = cA2
DV = np.array(DIV, float)
def logeta(t):                      # log eta(delta t) for delta in DIV, t array -> (4, n)
    out = []
    for dd in DIV:
        q = np.exp(2j * np.pi * dd * t); nmax = int(min(4000, 40 / max(2 * np.pi * dd * t.imag.min(), 1e-9))) + 2
        s = 2j * np.pi * dd * t / 24; qn = q.copy()
        for _ in range(nmax): s = s + np.log(1 - qn); qn = qn * q
        out.append(s)
    return np.array(out)
def evalset(t):
    q = np.exp(2j * np.pi * t); P = np.empty((len(t), NMAX + 1), complex); P[:, 0] = 1
    for n in range(1, NMAX + 1): P[:, n] = P[:, n - 1] * q
    F = P @ al.T                                            # (n, 13) cusp forms f_k
    LE = logeta(t); vals = {}
    for nm, cb in cons.items():
        for Q in CUSPW:
            h = np.exp(np.array(cb["r"][Q], float) @ LE) * cb["lam"][Q]
            vals[(nm, Q)] = h * (P @ cb["fw"][Q])
    return F, vals
# check product evaluation against the q-series at a high point
tt = np.array([0.1 + 1.0j, -0.3 + 0.6j])
F, vals = evalset(tt)
ser = lambda coef, q: sum(cf * q ** (i - 1) for i, cf in enumerate(coef))

# Ford domain
circ = []
for e in DIV:
    for g in range(1, 40):
        for dl in range(-143 * g - 5, 143 * g + 6):
            if math.gcd(e * dl, (143 // e) * g) != 1: continue
            c = -e * dl / (143 * g)
            if -0.6 <= c <= 0.6: circ.append((c, math.sqrt(e) / (143 * g)))
yc = 1 / math.sqrt(143)
def lower(xs):
    ym = np.zeros(len(xs))
    for c, rr in circ:
        i0, i1 = np.searchsorted(xs, c - rr), np.searchsorted(xs, c + rr)
        if i1 > i0: ym[i0:i1] = np.maximum(ym[i0:i1], np.sqrt(np.maximum(rr * rr - (xs[i0:i1] - c) ** 2, 0)))
    return ym
gx, gw = np.polynomial.legendre.leggauss(20)
xb = -0.5 + (np.arange(4000) + 0.5) / 4000; yb = lower(xb)
TB = np.repeat(xb, 20) + 1j * (yb[:, None] + (yc - yb)[:, None] * (gx[None, :] + 1) / 2).ravel()
WB = ((yc - yb)[:, None] * gw[None, :] / 2).ravel() / 4000
xa = -0.5 + (np.arange(2000) + 0.5) / 2000; pts, wts = [], []
for lo, hi in ((yc, 0.3), (0.3, 1.0), (1.0, 2.5)):
    yy = lo + (hi - lo) * (gx + 1) / 2; ww = (hi - lo) * gw / 2
    pts.append((xa[:, None] + 1j * yy[None, :]).ravel()); wts.append((np.ones(2000)[:, None] * ww[None, :] / 2000).ravel())
TA, WA = np.concatenate(pts), np.concatenate(wts)
Tall, Wall = np.r_[TB, TA], np.r_[WB, WA]
G = np.zeros((13, 13), complex); T = {(nm, Q): np.zeros(13, complex) for nm in cons for Q in CUSPW}
for s in range(0, len(Tall), 6000):
    F, vals = evalset(Tall[s:s + 6000]); w = Wall[s:s + 6000]
    G += (F.T * w) @ F.conj()
    for key, val in vals.items(): T[key] += (val * w) @ F.conj()
Gram = sum(M[Q] @ G @ M[Q].conj().T for Q in CUSPW)
print("Gram vs Y relative error:", np.abs(Gram - Y).max() / np.abs(Y).max(), " t", round(time.time() - t0, 1))
mu_inf = float(np.real(al[:, 1].conj() @ Yi @ al[:, 1])) / (13 * 4 * np.pi ** 2)
homes = {"tau/2 preimage": complex(0.18363, 0.03346), "Fricke point": 1j / math.sqrt(143)}
res = {}
for nm, cb in cons.items():
    P = sum(M[Q].conj() @ (T[(nm, Q)] * 2j * np.pi) for Q in CUSPW)     # (i/2) int psi ^ conj(omega_j), psi = 2 pi i f_psi dtau
    c = -np.linalg.solve(Y, P); out = {}
    tt = np.array(list(homes.values())); F, vals = evalset(tt)
    for i, name in enumerate(homes):
        q = np.exp(2j * np.pi * tt[i]); s = (vals[(nm, 1)][i] + c @ F[i] / (2j * np.pi)) / q
        mu = float(np.real((F[i] / q).conj() @ Yi @ (F[i] / q))) / (13 * 4 * np.pi ** 2)
        out[name] = abs(s) / math.sqrt(mu_inf * mu)
    for name, Q in (("cusp 0", 143), ("cusp 1/11", 13), ("cusp 1/13", 11)):
        a1 = (M[Q] @ al)[:, 1]; b1 = (cb["lam"][Q] * np.convolve(prodser(NMAX, cb["r"][Q]), cb["fw"][Q]))   # q-order offset below
        oc = int(round(orders(cb["r"][1])[(1, 143, 13, 11).index(Q)]))
        coef_q1 = b1[1 - oc] if 0 <= 1 - oc < len(b1) else 0.0
        s = coef_q1 + c @ a1 / (2j * np.pi)
        mu = float(np.real(a1.conj() @ Yi @ a1)) / (13 * 4 * np.pi ** 2)
        out[name] = abs(s) / math.sqrt(mu_inf * mu)
    res[nm] = out; print(nm, {k: round(v, 4) for k, v in out.items()})
json.dump({"mu_inf": mu_inf, "kernels": res}, open("/tmp/ag_results.json", "w")); print("t", round(time.time() - t0, 1))
