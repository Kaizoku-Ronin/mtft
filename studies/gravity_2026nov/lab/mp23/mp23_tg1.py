#!/usr/bin/env python3
"""MP23, T-G1 and T-G2: a tetrad for the squared hopfion, frozen into its null flow; the two helicities.
Plan: ../MP23_TG1_TG2_Frozen_Tetrad_Plan.md (written 14 November 2026, 09:11 New York time, and filed, before this file existed).

    python3 -I mp23_tg1.py [Q0 Q1 Q2 Q3 Q4]      ->  mp23_numbers.json next to this file (sections named, or all)

Conventions: x^mu = (t, x, y, z), signature (-,+,+,+), gradients carry a lower index.  a = x - iy, b = t - i - z,
d = r^2 - (t - i)^2, alpha = 1 - 2ib/d, beta = 2a/d, n = (4(b^2 - a^2), -4i(a^2 + b^2), 8ab), W = n (x) n / d^5 = E + iB,
v = Re(i n x conj n)/(n . conj n), k = (1, v).
The frozen class of the plan: B_{a nu} = Re(Q_a(alpha, beta) d_nu beta), so that F^a = dB^a = Re(G^a d alpha ^ d beta) with
G = dQ/d alpha, and h = B + B^T.  Linearised curvature by central second differences at 40 digits.
A failed check prints MISS and is counted; forecasts print held / missed."""
import itertools
import json
import math
import pathlib
import sys
import time

import mpmath as mp
import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
T0 = time.time()
mp.mp.dps = 40
NPASS = NFAIL = 0
OUTF = HERE / "mp23_numbers.json"
OUT = json.loads(OUTF.read_text()) if OUTF.exists() else {}
ETA = np.diag([-1.0, 1.0, 1.0, 1.0])
EPS3 = np.zeros((3, 3, 3))
for i_, j_, k_ in itertools.permutations(range(3)):
    EPS3[i_, j_, k_] = np.linalg.det(np.eye(3)[[i_, j_, k_]])
I = mp.mpc(0, 1)
HST = mp.mpf("1e-10")


def ok(name, cond, detail=""):
    global NPASS, NFAIL
    cond = bool(cond)
    NPASS += cond
    NFAIL += not cond
    print(("  ok   " if cond else "  MISS ") + name + (("   [" + detail + "]") if detail else ""), flush=True)
    return cond


def forecast(name, cond, detail=""):
    print(("  held   " if cond else "  missed ") + name + (("   [" + detail + "]") if detail else ""), flush=True)
    return bool(cond)


# ---------------------------------------------------------------- the fields, at mp precision
def fields(x):
    t, X, Y, Z = x
    a, b = X - I * Y, t - I - Z
    d = X * X + Y * Y + Z * Z - (t - I) ** 2
    al, be = 1 - 2 * I * b / d, 2 * a / d
    da, db, dd_ = [0, 1, -I, 0], [1, 0, 0, -1], [-2 * (t - I), 2 * X, 2 * Y, 2 * Z]
    dal = [-2 * I * db[m] / d + 2 * I * b * dd_[m] / d ** 2 for m in range(4)]
    dbe = [2 * da[m] / d - 2 * a * dd_[m] / d ** 2 for m in range(4)]
    n = [4 * (b * b - a * a), -4 * I * (a * a + b * b), 8 * a * b]
    return {"al": al, "be": be, "d": d, "dal": dal, "dbe": dbe, "n": n}


def flow(f):
    n = f["n"]
    nc = [mp.conj(v) for v in n]
    cr = [n[1] * nc[2] - n[2] * nc[1], n[2] * nc[0] - n[0] * nc[2], n[0] * nc[1] - n[1] * nc[0]]
    nn = sum(n[i] * nc[i] for i in range(3))
    v = [mp.re(I * cr[i]) / mp.re(nn) for i in range(3)]
    return [mp.mpf(1)] + v                                         # k^mu


def W_of(f):
    n, d = f["n"], f["d"]
    return np.array([[complex(n[i] * n[j] / d ** 5) for j in range(3)] for i in range(3)])


def second(func, x):
    """all second derivatives of a list-valued function: D[rho][sigma] is a list, by central differences"""
    def at(shift):
        return func([x[m] + shift[m] for m in range(4)])
    f0 = at([0, 0, 0, 0])
    nf = len(f0)
    D = [[None] * 4 for _ in range(4)]
    plus, minus = [], []
    for r in range(4):
        e = [HST if m == r else 0 for m in range(4)]
        plus.append(at(e))
        minus.append(at([-v for v in e]))
        D[r][r] = [(plus[r][q] - 2 * f0[q] + minus[r][q]) / HST ** 2 for q in range(nf)]
    for r in range(4):
        for s in range(r + 1, 4):
            sh = lambda sr, ss: [(sr * HST if m == r else 0) + (ss * HST if m == s else 0) for m in range(4)]
            pp, pm, mp_, mm = at(sh(1, 1)), at(sh(1, -1)), at(sh(-1, 1)), at(sh(-1, -1))
            D[r][s] = D[s][r] = [(pp[q] - pm[q] - mp_[q] + mm[q]) / (4 * HST ** 2) for q in range(nf)]
    return D


def first(func, x):
    def at(shift):
        return func([x[m] + shift[m] for m in range(4)])
    out = []
    for r in range(4):
        e = [HST if m == r else 0 for m in range(4)]
        fp, fm = at(e), at([-v for v in e])
        out.append([(fp[q] - fm[q]) / (2 * HST) for q in range(len(fp))])
    return out


def curvature(D2h):
    """D2h[rho, sigma, mu, nu] = d_rho d_sigma h_{mu nu} (complex128) -> E_ij = R_0i0j, B_ij = (1/2) eps_jkl R_0ikl, Ricci"""
    R = 0.5 * (np.einsum("nrms->mnrs", D2h) + np.einsum("msnr->mnrs", D2h) - np.einsum("mrns->mnrs", D2h) - np.einsum("nsmr->mnrs", D2h))
    E = R[0, 1:, 0, 1:]
    B = 0.5 * np.einsum("jkl,ikl->ij", EPS3, R[0, 1:, 1:, 1:])
    Ric = np.einsum("rs,rmsn->mn", ETA, R)
    return E, B, Ric, R


def to_np(D, shape):
    return np.array([[[complex(v) for v in D[r][s]] for s in range(4)] for r in range(4)]).reshape((4, 4) + shape)


def monos(deg):
    return [(i, tot - i) for tot in range(deg + 1) for i in range(tot, -1, -1)]         # G contains alpha^i beta^j


def T_list(x, MON):
    """for each monomial of G, the vector q d_nu beta with q = alpha^(i+1) beta^j / (i+1); flattened"""
    f = fields(x)
    out = []
    for i, j in MON:
        q = f["al"] ** (i + 1) * f["be"] ** j / (i + 1)
        out += [q * f["dbe"][nu] for nu in range(4)]
    return out


def basis_curvatures(x, MON):
    """E, B, Ricci of the complex potentials h_{mu nu} = delta_{mu c} T_nu + delta_{nu c} T_mu, for every monomial and c"""
    D = to_np(second(lambda y: T_list(y, MON), x), (len(MON), 4))             # [rho, sigma, mono, nu]
    out = []
    for m_ in range(len(MON)):
        for c in range(4):
            D2h = np.zeros((4, 4, 4, 4), complex)
            D2h[:, :, c, :] += D[:, :, m_, :]
            D2h[:, :, :, c] += D[:, :, m_, :]
            out.append(curvature(D2h)[:3])
    return out


def rand_points(n, seed):
    r = np.random.default_rng(seed)
    return [[mp.mpf(float(r.uniform(-1, 2)))] + [mp.mpf(float(v)) for v in r.normal(0, 1.2, 3)] for _ in range(n)]


def fit(MON, pts, sign):
    """real least squares for the complex coefficients c_j: Re(sum c_j E_j) = Re W, Re(sum c_j B_j) = sign Im W, Re(sum c_j Ric_j) = 0"""
    rows, rhs = [], []
    for x in pts:
        W = W_of(fields(x))
        sc = np.linalg.norm(W)
        bc = basis_curvatures(x, MON)
        for which, tgt in ((0, W.real), (1, sign * W.imag), (2, np.zeros((4, 4)))):
            A = np.array([b_[which] for b_ in bc])                              # [basis, i, j]
            A = A.reshape(len(bc), -1)
            rows.append(np.concatenate([A.real, -A.imag], axis=0).T / sc)
            rhs.append(tgt.reshape(-1) / sc)
    M, y = np.concatenate(rows), np.concatenate(rhs)
    u, res, rank, sv = np.linalg.lstsq(M, y, rcond=1e-11)
    return u[:len(u) // 2] + 1j * u[len(u) // 2:], float(np.abs(M @ u - y).max()), int(rank), sv, M, y


SECT = set(sys.argv[1:]) or {"Q0", "Q1", "Q2", "Q3", "Q4"}

# ================================================================ Q0: the setting
if "Q0" in SECT:
    print("Q0 the setting")
    w = [mp.mpf(0)] * 6
    for x in rand_points(12, 5):
        f = fields(x)
        k = flow(f)
        g = lambda u_, v_: -u_[0] * v_[0] + sum(u_[i] * v_[i] for i in (1, 2, 3))        # the Minkowski product of two covectors
        sc = sum(abs(v) for v in f["dal"]) + sum(abs(v) for v in f["dbe"])
        w[0] = max(w[0], abs(sum(k[m] * f["dal"][m] for m in range(4))) / sc, abs(sum(k[m] * f["dbe"][m] for m in range(4))) / sc)
        w[1] = max(w[1], abs(g(f["dal"], f["dal"])) / sc ** 2, abs(g(f["dbe"], f["dbe"])) / sc ** 2, abs(g(f["dal"], f["dbe"])) / sc ** 2)
        sig = [f["al"] * f["dbe"][m] - f["be"] * f["dal"][m] for m in range(4)]
        N = [0, 2 * f["al"], -2 * I * f["al"], -2 * I * f["be"]]
        xl = [-x[0], x[1], x[2], x[3]]
        w[2] = max(w[2], max(abs(f["d"] * sig[m] - (N[m] - 2 * f["be"] * xl[m])) for m in range(4)) / (abs(f["d"]) * sum(abs(v) for v in sig)))
        n, d = f["n"], f["d"]
        w[3] = max(w[3], abs(sum(n[i] * n[i] for i in range(3))) / sum(abs(n[i]) ** 2 for i in range(3)))
        w[4] = max(w[4], abs(sum(n[i] * k[i + 1] for i in range(3))) / mp.sqrt(sum(abs(n[i]) ** 2 for i in range(3))))
        w[5] = max(w[5], abs(-k[0] ** 2 + sum(k[i] ** 2 for i in (1, 2, 3))))
        # the labels are constant along the ray through x
        f2 = fields([x[m] + mp.mpf("0.7") * k[m] for m in range(4)])
        w[0] = max(w[0], abs(f2["al"] - f["al"]), abs(f2["be"] - f["be"]))
    OUT["Q0"] = {"k.d(alpha), k.d(beta), and the labels along a ray": float(w[0]), "the gradients are null and orthogonal": float(w[1]), "d sigma_mu = N_mu - 2 beta x_mu": float(w[2]),
                 "n.n = 0": float(w[3]), "n.v = 0": float(w[4]), "k is null": float(w[5])}
    ok("Q0 the labels are carried by k; their gradients are null and orthogonal; d sigma_mu = N_mu - 2 beta x_mu; W is trace-free and across v (1e-30)", max(w) < mp.mpf("1e-30"), "worst %.1e" % float(max(w)))
    # the curvature routine on a plane wave: h_xx = -h_yy = cos(t - z) has E_xx = -E_yy = -cos(t - z)/2 ... checked as vacuum, type N
    pw = lambda y: [mp.cos(y[0] - y[3]) if (m, n_) == (1, 1) else -mp.cos(y[0] - y[3]) if (m, n_) == (2, 2) else mp.mpf(0) for m in range(4) for n_ in range(4)]
    x = rand_points(1, 9)[0]
    E, B, Ric, R = curvature(to_np(second(pw, x), (4, 4)))
    cph = float(mp.cos(x[0] - x[3]))
    okpw = abs(E[0, 0] - 0.5 * cph) < 1e-15 and abs(E[1, 1] + 0.5 * cph) < 1e-15 and np.abs(Ric).max() < 1e-15 and abs(abs(B[0, 1]) - 0.5 * abs(cph)) < 1e-15 and abs(B[0, 1] - B[1, 0]) < 1e-15
    gz = lambda y: [(y[m] ** 2 * y[n_] + y[n_] ** 2 * y[m]) * mp.sin(y[0] + 2 * y[1]) + 0 * I for m in range(4) for n_ in range(4)]     # not a gauge term: only a non-zero control
    xi = lambda y: [mp.sin(y[0] * y[2]) * (m + 1) + mp.cos(y[1] + y[3] * y[3]) * (m * m - 2) for m in range(4)]
    dxi = first(xi, x)                                                                # d_mu xi_nu, then its second derivatives by nesting
    hg = lambda y: (lambda dv: [dv[m][n_] + dv[n_][m] for m in range(4) for n_ in range(4)])(first(xi, y))
    Eg, Bg, Ricg, Rg = curvature(to_np(second(hg, x), (4, 4)))
    okg = np.abs(Rg).max() < 1e-7
    OUT["Q0"]["the curvature routine: plane wave, and a pure gauge term"] = [bool(okpw), float(np.abs(Rg).max())]
    ok("Q0' the curvature routine: a plane wave has E_xx = cos/2, |B_xy| = |cos|/2, no Ricci (1e-15); a gauge term d_mu xi_nu + d_nu xi_mu has no curvature (1e-7, nested differences)", okpw and okg,
       "gauge term's curvature %.1e" % np.abs(Rg).max())

# ================================================================ Q1: the square of the potential
if "Q1" in SECT:
    print("Q1 the square of the potential")
    def H1(y):
        f = fields(y)
        sig = [f["al"] * f["dbe"][m] - f["be"] * f["dal"][m] for m in range(4)]
        return [f["d"] * sig[m] * sig[n_] for m in range(4) for n_ in range(4)]
    best = None
    pts = rand_points(16, 11)
    for sign in (1, -1):
        rows, rhs = [], []
        for x in pts:
            W = W_of(fields(x))
            sc = np.linalg.norm(W)
            E, B, Ric, R = curvature(to_np(second(H1, x), (4, 4)))
            for A, tgt in ((E, W.real), (B, sign * W.imag), (Ric, np.zeros((4, 4)))):
                A = A.reshape(-1)
                rows.append(np.stack([A.real, -A.imag], axis=1) / sc)
                rhs.append(tgt.reshape(-1) / sc)
        M, y = np.concatenate(rows), np.concatenate(rhs)
        u = np.linalg.lstsq(M, y, rcond=None)[0]
        r = float(np.abs(M @ u - y).max())
        if best is None or r < best[0]:
            best = (r, sign, complex(u[0], u[1]))
    # how large is the Ricci tensor of H1 itself, and its E, B against W, without any fit
    x = pts[0]
    E, B, Ric, R = curvature(to_np(second(H1, x), (4, 4)))
    W = W_of(fields(x))
    OUT["Q1"] = {"best residual / |W|": best[0], "sign of B": best[1], "c": [best[2].real, best[2].imag], "at one point: |Ricci|/|W|": float(np.abs(Ric).max() / np.linalg.norm(W)),
                 "|E + iB| / |W|": float(np.linalg.norm(E + 1j * B) / np.linalg.norm(W)), "|E - iB| / |W|": float(np.linalg.norm(E - 1j * B) / np.linalg.norm(W))}
    OUT["Q1"]["forecast held"] = forecast("Q1 FORECAST (0.75) the curvature of Re(c d sigma (x) sigma) is the squared hopfion for one constant c", best[0] < 1e-9,
                                          "residual %.1e of |W|; c = %.6f%+.6fi; sign %+d; Ricci of d sigma sigma itself %.1e of |W|" % (best[0], best[2].real, best[2].imag, best[1], OUT["Q1"]["at one point: |Ricci|/|W|"]))

# ================================================================ Q2: the search in the frozen class
if "Q2" in SECT:
    print("Q2 the search in the frozen class")
    Q2 = {}
    found = None
    for deg in (2, 4):
        MON = monos(deg)
        pts = rand_points(40, 23)
        for sign in (1, -1):
            c, res, rank, sv, M, y = fit(MON, pts, sign)
            Q2["degree <= %d, sign %+d" % (deg, sign)] = {"unknowns (complex)": len(c), "residual / |W|": res, "rank": rank, "real unknowns": 2 * len(c),
                                                        "singular values, largest and five smallest": [float(sv[0])] + [float(v) for v in sv[-5:]]}
            print("    degree <= %d, sign %+d: residual %.2e of |W|, rank %d of %d (%.0f s)" % (deg, sign, res, rank, 2 * len(c), time.time() - T0), flush=True)
            if res < 1e-10 and found is None:
                found = {"degree": deg, "sign": sign, "monomials": MON, "c": [[v.real, v.imag] for v in c], "rank": rank, "kernel dimension (real)": 2 * len(c) - rank}
        if found:
            break
    Q2["found"] = found
    OUT["Q2"] = Q2
    if found:
        c = np.array([complex(*v) for v in found["c"]]).reshape(len(found["monomials"]), 4)
        print("    the coefficients of G^a (rows: monomials %s; columns: a = t, x, y, z), least-norm solution:" % found["monomials"])
        for mrow, cr in zip(found["monomials"], c):
            print("      alpha^%d beta^%d: " % mrow + "  ".join("%+.6f%+.6fi" % (v.real, v.imag) for v in cr))
    Q2["forecast held"] = forecast("Q2 FORECAST (0.65) a frozen tetrad of the holomorphic class with G of degree <= 2 reproduces the squared hopfion (1e-10 of |W|)", bool(found and found["degree"] == 2),
                                   "found at degree %s" % (found["degree"] if found else "none"))

# ================================================================ Q3: the exhibit (the closed form of Amendment 1)
def B_closed(y):
    """B_{mu nu} = Re[(i a^2 (d_mu b)(d_nu gamma) - b^2 (d_mu a)(d_nu beta)) / d^2], gamma = 1 - alpha = 2ib/d; flattened [mu][nu]"""
    t, X, Y, Z = y
    a, b = X - I * Y, t - I - Z
    d = X * X + Y * Y + Z * Z - (t - I) ** 2
    da, db, dd_ = [0, 1, -I, 0], [1, 0, 0, -1], [-2 * (t - I), 2 * X, 2 * Y, 2 * Z]
    dgam = [2 * I * db[m] / d - 2 * I * b * dd_[m] / d ** 2 for m in range(4)]
    dbe = [2 * da[m] / d - 2 * a * dd_[m] / d ** 2 for m in range(4)]
    return [mp.re((I * a * a * db[m] * dgam[n_] - b * b * da[m] * dbe[n_]) / d ** 2) for m in range(4) for n_ in range(4)]


if "Q3" in SECT:
    print("Q3 the exhibit")
    found = OUT["Q2"]["found"]
    cfit = np.array([complex(*v) for v in found["c"]]).reshape(len(found["monomials"]), 4)
    MONF = [tuple(v) for v in found["monomials"]]

    def B_fit(y):
        f = fields(y)
        Q = [sum(mp.mpc(cfit[m_, c]) * f["al"] ** (i + 1) * f["be"] ** j / (i + 1) for m_, (i, j) in enumerate(MONF)) for c in range(4)]
        return [mp.re(Q[m] * f["dbe"][n_]) for m in range(4) for n_ in range(4)]
    w_fro = w_E = w_B = w_Ric = w_Rk = w_same = 0.0
    for x in rand_points(60, 77):
        f = fields(x)
        k = [float(v) for v in flow(f)]
        W = W_of(f)
        sc = np.linalg.norm(W)
        dB = np.array([[float(v) for v in row] for row in first(B_closed, x)]).reshape(4, 4, 4)       # [rho, mu, nu]
        F = dB - np.transpose(dB, (2, 1, 0))                                                    # F[rho, mu, nu] = d_rho B_{mu nu} - d_nu B_{mu rho}
        w_fro = max(w_fro, np.abs(np.einsum("r,rmn->mn", k, F)).max() / np.abs(F).max())
        dBf = np.array([[float(v) for v in row] for row in first(B_fit, x)]).reshape(4, 4, 4)
        Ff = dBf - np.transpose(dBf, (2, 1, 0))
        w_same = max(w_same, np.abs(F - Ff).max() / np.abs(F).max())
        hfun = lambda y: (lambda Bv: [Bv[4 * m + n_] + Bv[4 * n_ + m] for m in range(4) for n_ in range(4)])(B_closed(y))
        E, Bm, Ric, R = curvature(to_np(second(hfun, x), (4, 4)))
        w_E = max(w_E, np.abs(E - W.real).max() / sc)
        w_B = max(w_B, np.abs(Bm + W.imag).max() / sc)
        w_Ric = max(w_Ric, np.abs(Ric).max() / sc)
        w_Rk = max(w_Rk, np.abs(np.einsum("mnrs,s->mnr", R, k)).max() / np.abs(R).max())
    OUT["Q3"] = {"|i_k F| / |F|": w_fro, "E - Re W": w_E, "B + Im W (B_ij = (1/2) eps_jkl R_0ikl)": w_B, "Ricci": w_Ric, "R.k / |R|": w_Rk, "F of the fitted potential against F of the closed form": w_same}
    ok("Q3 PREDICTION (T-G1's criterion) the tetrad e = dx + B is frozen into k (1e-12 of |F|), the curvature of h = B + B^T is the squared hopfion (1e-10 of |W|), and its Ricci tensor vanishes (1e-10): 60 points, t from -1 to 2",
       w_fro < 1e-12 and w_E < 1e-10 and w_B < 1e-10 and w_Ric < 1e-10, "i_k F %.1e; E %.1e; B %.1e; Ricci %.1e; R.k %.1e; fitted and closed forms differ by %.1e" % (w_fro, w_E, w_B, w_Ric, w_Rk, w_same))

# ================================================================ Q4: the two helicities
if "Q4" in SECT:
    print("Q4 the two helicities")
    DA, DB = np.array([0, 1, -1j, 0]), np.array([1, 0, 0, -1], complex)
    SIGN = np.array([-1.0, 1.0, 1.0, 1.0])

    def three_forms(al, be, dal, dbe):
        """B_a ^ dB^a and A ^ dA on three tangent vectors, given the derivatives of the labels along them (arrays [3, ...])"""
        gam, dgam = 1 - al, -dal
        Bv = 0.25 * np.real(gam ** 2 * DA[:, None, None] * dbe[None] + 1j * be ** 2 * DB[:, None, None] * dgam[None])            # [mu, i, ...]
        pair = lambda u, v, j, k: u[j] * v[k] - u[k] * v[j]
        dBv = lambda j, k: 0.5 * np.real((gam * DA[:, None] - 1j * be * DB[:, None]) * pair(dgam, dbe, j, k)[None])            # [mu, ...]
        Kg = sum(SIGN[m] * (Bv[m, 0] * dBv(1, 2)[m] - Bv[m, 1] * dBv(0, 2)[m] + Bv[m, 2] * dBv(0, 1)[m]) for m in range(4))
        Av = 0.5 * np.real(al * dbe - be * dal)
        dAv = lambda j, k: np.real(pair(dal, dbe, j, k))
        Ke = Av[0] * dAv(1, 2) - Av[1] * dAv(0, 2) + Av[2] * dAv(0, 1)
        return Kg, Ke
    # (a) on the label sphere: alpha = cos(eta) e^{i xi_1}, beta = sin(eta) e^{i xi_2}
    ge, we = np.polynomial.legendre.leggauss(40)
    eta, we = (ge + 1) * math.pi / 4, we * math.pi / 4
    xi = np.arange(48) * 2 * math.pi / 48
    Et, X1, X2 = np.meshgrid(eta, xi, xi, indexing="ij")
    al, be = np.cos(Et) * np.exp(1j * X1), np.sin(Et) * np.exp(1j * X2)
    dal = np.array([-np.sin(Et) * np.exp(1j * X1), 1j * al, 0 * al])
    dbe = np.array([np.cos(Et) * np.exp(1j * X2), 0 * be, 1j * be])
    Kg, Ke = three_forms(al.reshape(-1), be.reshape(-1), dal.reshape(3, -1), dbe.reshape(3, -1))
    wgt = (we[:, None, None] * np.ones_like(Et) * (2 * math.pi / 48) ** 2).reshape(-1)
    Hg_s, He_s = float(np.sum(Kg * wgt)), float(np.sum(Ke * wgt))
    # (b) in space, on the slices t = 0 and t = 0.8
    def slice_integrals(tv, n=88, L=1.2):
        g, w = np.polynomial.legendre.leggauss(n)
        xs, ws = L * g / (1 - g * g), w * L * (1 + g * g) / (1 - g * g) ** 2
        X, Y, Z = np.meshgrid(xs, xs, xs, indexing="ij")
        Wt = (ws[:, None, None] * ws[None, :, None] * ws[None, None, :]).reshape(-1)
        X, Y, Z = X.reshape(-1), Y.reshape(-1), Z.reshape(-1)
        a, b = X - 1j * Y, tv - 1j - Z
        d = X * X + Y * Y + Z * Z - (tv - 1j) ** 2
        al, be = 1 - 2j * b / d, 2 * a / d
        da, db, dd_ = [1, -1j, 0], [0, 0, -1], [2 * X, 2 * Y, 2 * Z]                                     # the three space derivatives
        dal = np.array([-2j * db[i] / d + 2j * b * dd_[i] / d ** 2 for i in range(3)])
        dbe = np.array([2 * da[i] / d - 2 * a * dd_[i] / d ** 2 for i in range(3)])
        Kg, Ke = three_forms(al, be, dal, dbe)
        return float(np.sum(Kg * Wt)), float(np.sum(Ke * Wt))
    sl = {tv: slice_integrals(tv) for tv in (0.0, 0.8)}
    # the normalisation: grad(alpha) x grad(beta) against n / d^3, and the electric vector of the 2-form against i times it
    wn = 0.0
    for x in rand_points(8, 31):
        f = fields(x)
        A_, B_ = [complex(v) for v in f["dal"]], [complex(v) for v in f["dbe"]]
        cr = np.cross(A_[1:], B_[1:])
        Fv = np.array([complex(v / f["d"] ** 3) for v in f["n"]])
        Ev = np.array([A_[i] * B_[0] - A_[0] * B_[i] for i in (1, 2, 3)])                          # (d alpha ^ d beta)_{i0}
        wn = max(wn, np.abs(cr - Fv).max() / np.abs(Fv).max(), np.abs(Ev - 1j * Fv).max() / np.abs(Fv).max())
    Q4 = {"on the label sphere": {"H_g": Hg_s, "H_em": He_s, "ratio": Hg_s / He_s}, "pi^2": math.pi ** 2, "3 pi^2 / 4": 0.75 * math.pi ** 2,
          "in space, t = 0": {"H_g": sl[0.0][0], "H_em": sl[0.0][1], "ratio": sl[0.0][0] / sl[0.0][1]}, "in space, t = 0.8": {"H_g": sl[0.8][0], "H_em": sl[0.8][1], "ratio": sl[0.8][0] / sl[0.8][1]},
          "grad alpha x grad beta = n/d^3, and the electric vector of d alpha ^ d beta is i n/d^3": wn}
    OUT["Q4"] = Q4
    worst = max(abs(abs(v["H_em"]) / math.pi ** 2 - 1) for v in (Q4["on the label sphere"], Q4["in space, t = 0"], Q4["in space, t = 0.8"]))
    worstr = max(abs(v["ratio"] - 0.75) for v in (Q4["on the label sphere"], Q4["in space, t = 0"], Q4["in space, t = 0.8"]))
    ok("Q4 PREDICTION (T-G2) |H_em| = pi^2 and H_g / H_em = 3/4, on the label sphere and on the slices t = 0 and t = 0.8 in space (1e-6); the light field is n/d^3 up to a duality turn",
       worst < 1e-6 and worstr < 1e-6 and wn < 1e-12, "H_em/pi^2 off by %.1e; ratio off 3/4 by %.1e; sphere: H_g = %.8f, H_em = %.8f; space t=0: %.8f, %.8f; t=0.8: %.8f, %.8f; normalisation %.1e"
       % (worst, worstr, Hg_s, He_s, sl[0.0][0], sl[0.0][1], sl[0.8][0], sl[0.8][1], wn))

OUT.setdefault("checks", {})["run of sections %s" % " ".join(sorted(SECT))] = {"passed": NPASS, "failed": NFAIL, "seconds": time.time() - T0}
OUTF.write_text(json.dumps(OUT, indent=1))
print("%d passed, %d failed (%.0f s)" % (NPASS, NFAIL, time.time() - T0))
