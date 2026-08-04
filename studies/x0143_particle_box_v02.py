#!/usr/bin/env python3
"""
x0143_particle_box_v02.py — Petersson-weighted signatures + Zeno curve
======================================================================
Roger Tano / MTFT Research Program — built with Claude, August 2026

Upgrades over v0.1 (which stays the certified baseline):

  Stage E  Real period integrals, computed from first principles:
    E1  star involution eta:(c:d)->(-c:d); cuspidal space splits 13+13.
    E2  joint eigendata of T2 on the cuspidal space (11 new lines, +/-).
    E3  Hecke eigenvalues a_p for every new embedding, p <= P_MAX, from
        this engine's own T_p matrices; cross-certified against the
        mtft ORBIT_TRACE tables (n <= 50) and curve 143a1.
    E4  Periods of ~50 explicit Gamma_0(143)-loops by holomorphic
        antiderivative at the optimal points z0=(-d+i)/143, gz0=(a+i)/143;
        the 26-dim period functional of each eigenform solved from an
        overdetermined system  ->  residual certificate.
    E5  f1 lattice cross-certificate: the engine's f1 periods must form
        a rank-2 lattice matching curve 143a1 (corpus |Omega| in
        {0.3135, 0.3571}, tau ~ 0.039+0.980i; plus an independent
        quadrature real period of the Weierstrass model).
    E6  Born-rule signatures per Paper 36 SS6:
            P(f) = |period_f(gamma)|^2 / <f,f>_Pet   (normalized)
        with the mfpetersson norms from the corpus (class Cert-corpus).
        Old block reported separately (declared normalization).

  Stage F  Zeno curve on the particle box: effective decay rate R(tau)
        under repeated projective measurement at interval tau, with the
        exact small-tau law R -> (Delta L)^2 tau, and (Delta L)^2 = 3
        exactly for a one-triangle state (variance = degree).

Epistemic classes stated on every number: EXACT / Cert / Cert-corpus /
DIAGNOSTIC.
"""

from __future__ import annotations
from fractions import Fraction
from math import gcd, pi
import cmath
import json

import numpy as np

import os as _os
_HERE = _os.path.dirname(_os.path.abspath(__file__))

from x0143_particle_box import (P1, SIGMA, TAU, Cusp, ModularSymbols,
                                tessellation, dual_graph)

ETA = ((-1, 0), (0, 1))          # star involution (c:d) -> (-c:d)
N = 143
P_MAX = 700                       # primes for a_p extraction
N_MAX = 700                       # q-expansion length for periods

# ----- corpus constants (Cert-corpus: Paper 36 SS6.3 / Paper 32 SS4.3,
#       computed by mfpetersson in PARI/GP, 28-digit precision) ----------
# Petersson norms, PER-UNIT-VOLUME normalization (Add. BQ sec.3):
#   <f,f>_raw = PET * V_N with V_N = 56*pi (volume of X_0(143)).
#   Auditor certification: PET_F1 sits 4.7% above the exact modular-degree
#   value 0.002287 (deg(143.a1)=4); the f2/f3 diagonals match the auditor's
#   independent Rankin-Selberg residues to a uniform 2-3% (Add. BQ sec.4.4).
PET_F1 = 0.002394868866550
PET_F2_DIAG = [0.00720, 0.00423, 0.00473, 0.01431]   # sigma-order of P32 SS6.1
PET_F3_DIAG = [0.01369, 0.01085, 0.00388, 0.00627, 0.01008, 0.01564]
A2_F2_ORDER = [-0.197, 1.747, -1.127, 2.576]         # a2 per sigma (P32 SS6.1)
A2_F3_ORDER = [2.447, 1.365, 1.231, -0.633, -1.701, -2.709]  # assumed order
CORPUS_F1_ABS = (0.3135, 0.3571)                      # Paper 32 SS4.4
CORPUS_TAU = 0.039 + 0.980j
# U_p eigenvalues at p || N per orbit: a_p = -eps_p  (audit report)
AP_BAD = {"f1": {11: -1, 13: -1}, "f2": {11: 1, 13: -1}, "f3": {11: -1, 13: 1}}


def build_engine():
    p1, tris, edges, cert = tessellation(N)
    ms = ModularSymbols(p1)
    return p1, tris, edges, ms


# ----------------------------------------------------------------------
# fast float pipeline
# ----------------------------------------------------------------------

def float_projection(ms):
    P = np.zeros((ms.dim, ms.n))
    for col, comb in ms.proj.items():
        for k, w in comb.items():
            P[k, col] = float(w)
    return P


_conv_cache = {}


def chain_indices(ms, alpha: Cusp, beta: Cusp):
    """Fast integer chain of {alpha->beta} as list of (coset_index, sign)."""
    out = []
    for c, s in ((beta, +1), (alpha, -1)):
        key = (c.n, c.d)
        rows = _conv_cache.get(key)
        if rows is None:
            from x0143_particle_box import convergent_chain
            rows = convergent_chain(c)
            _conv_cache[key] = rows
        for row in rows:
            out.append((ms.manin_index(row), s))
    return out


def hecke_float(ms, P, p, is_bad=False):
    """Float T_p (or U_p when p | N) on the quotient."""
    from x0143_particle_box import hecke_matrices
    mats = hecke_matrices(p)
    if is_bad:
        mats = [m for m in mats if m[0][0] == 1]      # drop [[p,0],[0,1]]
    M = np.zeros((ms.dim, ms.dim))
    for j_idx, j in enumerate(ms.free):
        xrep = ms.p1.reps[j]
        (a, b), (c0, d0) = ms.p1.lift(xrep)
        alpha, beta = Cusp(b, d0), Cusp(a, c0)
        col = np.zeros(ms.n)
        for m in mats:
            for idx, s in chain_indices(ms, alpha.apply(m), beta.apply(m)):
                col[idx] += s
        M[:, j_idx] = P @ col
    return M


def eta_float(ms, P):
    M = np.zeros((ms.dim, ms.dim))
    for j_idx, j in enumerate(ms.free):
        x = ms.p1.reps[j]
        y = ms.p1.act(x, ETA)
        col = np.zeros(ms.n)
        col[ms.p1.index[y]] = 1.0
        M[:, j_idx] = P @ col
    return M


# ----------------------------------------------------------------------
# Stage E2/E3 — joint eigendata and a_p for every new embedding
# ----------------------------------------------------------------------

def eigendata(ms, P):
    """Returns cuspidal basis Bc (dim x 26, float), T2 on cuspidal,
    eta on cuspidal, and per-embedding right/left eigenvectors."""
    Bc = np.array([[float(v[i]) for v in ms.cuspidal_basis]
                   for i in range(ms.dim)])          # dim x 26
    proj_c = np.linalg.pinv(Bc)                       # 26 x dim

    def restrict(M):
        return proj_c @ (M @ Bc)

    T2 = restrict(hecke_float(ms, P, 2))
    E = restrict(eta_float(ms, P))
    inv_err = np.linalg.norm(E @ E - np.eye(26))
    assert inv_err < 1e-8, f"eta not an involution: {inv_err}"

    def orth_range(M):
        U, s, _ = np.linalg.svd(M)
        return U[:, s > 0.5]

    plus = orth_range(0.5 * (np.eye(26) + E))
    minus = orth_range(0.5 * (np.eye(26) - E))
    assert plus.shape[1] == 13 and minus.shape[1] == 13, \
        (plus.shape, minus.shape)

    def realize(v):
        k = int(np.argmax(np.abs(v)))
        v = v / v[k]
        assert np.abs(v.imag).max() < 1e-7, "non-real eigenvector"
        return v.real

    lines = []   # [sign, orbit, a2, right26, left26] for the 11 new lines
    for sign, V in (("+", plus), ("-", minus)):
        inv_res = np.linalg.norm(T2 @ V - V @ (V.T @ T2 @ V))
        assert inv_res < 1e-8, f"eta sector not T2-invariant: {inv_res}"
        T2s = V.T @ T2 @ V                            # 13x13
        w, R = np.linalg.eig(T2s)
        wl, Lv = np.linalg.eig(T2s.T)
        order = np.argsort(w.real)
        for k in order:
            lam = w[k].real
            if abs(lam + 2) < 1e-6:
                continue                              # old block
            kl = int(np.argmin(np.abs(wl - w[k])))
            r26 = V @ realize(R[:, k])
            l26 = V @ realize(Lv[:, kl])
            orbit = "f1" if abs(lam) < 1e-6 else None
            lines.append([sign, orbit, lam, r26, l26])
    return Bc, proj_c, restrict, T2, E, lines


def assign_orbits(lines):
    import numpy.polynomial.polynomial as npp
    q2 = [1, 5, -1, -3, 1][::-1]     # x^4-3x^3-x^2+5x+1 coeffs high->low
    for L in lines:
        lam = L[2]
        if L[1] == "f1":
            continue
        v2 = np.polyval([1, -3, -1, 5, 1], lam)
        v3 = np.polyval([1, 0, -10, 2, 24, -7, -12], lam)
        L[1] = "f2" if abs(v2) < 1e-4 else ("f3" if abs(v3) < 1e-4 else "??")
    assert all(L[1] in ("f1", "f2", "f3") for L in lines)
    return lines


def primes_upto(n):
    sieve = np.ones(n + 1, dtype=bool)
    sieve[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if sieve[i]:
            sieve[i * i::i] = False
    return [int(p) for p in np.nonzero(sieve)[0]]


def extract_ap(ms, P, restrict, lines, pmax):
    ps = [p for p in primes_upto(pmax) if p not in (11, 13)]
    ap = {id(L): {} for L in lines}
    for p in ps:
        Tp = restrict(hecke_float(ms, P, p))
        for L in lines:
            r, l = L[3], L[4]
            ap[id(L)][p] = float((l @ (Tp @ r)) / (l @ r))
    for L in lines:
        for p in (11, 13):
            ap[id(L)][p] = float(AP_BAD[L[1]][p])
    return ap


def an_table(ap_p, nmax):
    """Multiplicative extension of a_p to a_n, n <= nmax."""
    spf = np.zeros(nmax + 1, dtype=int)
    for i in range(2, nmax + 1):
        if spf[i] == 0:
            spf[i::i][spf[i::i] == 0] = i
    a = np.zeros(nmax + 1)
    a[1] = 1.0
    for n in range(2, nmax + 1):
        p = int(spf[n])
        k, m = 0, n
        while m % p == 0:
            m //= p
            k += 1
        pk = n // m
        if m > 1:
            a[n] = a[pk] * a[m]
        else:                                   # n = p^k
            if p in (11, 13):
                a[n] = ap_p[p] ** k
            elif k == 1:
                a[n] = ap_p[p]
            else:
                a[n] = ap_p[p] * a[pk // p] - p * a[pk // (p * p)]
    return a


# ----------------------------------------------------------------------
# Stage E4 — periods of explicit Gamma_0(143)-loops
# ----------------------------------------------------------------------

def gamma_loops():
    """gamma = [[a,b],[143,d]] with ad - 143 b = 1, small entries."""
    out = []
    for d in range(1, N):
        if gcd(d, N) != 1:
            continue
        a = pow(d, -1, N)
        if a > N // 2:
            a -= N
        b = (a * d - 1) // N
        assert a * d - N * b == 1
        out.append((a, b, d))
    return out


def F_of(a_n, z):
    """Antiderivative sum_{n} a_n/n e^{2 pi i n z} (2 pi i * int f dtau)."""
    nn = np.arange(1, len(a_n))
    return np.sum(a_n[1:] / nn * np.exp(2j * pi * nn * z))


def period_functionals(ms, P, Bc, proj_c, lines, ap, verbose=True):
    loops = gamma_loops()
    # homology classes of the loops (float, in 26-dim cuspidal coords)
    X, keep = [], []
    for (a, b, d) in loops:
        col = np.zeros(ms.n)
        for idx, s in chain_indices(ms, Cusp(0, 1), Cusp(b, d)):
            col[idx] += s
        q = P @ col
        v = proj_c @ q
        if np.linalg.norm(Bc @ v - q) < 1e-9:          # cuspidality check
            X.append(v)
            keep.append((a, b, d))
    X = np.array(X)
    rank = np.linalg.matrix_rank(X, tol=1e-8)
    if verbose:
        print(f"  Gamma-loops kept: {len(keep)}; class-matrix rank {rank} "
              f"(need 26) {'PASS' if rank == 26 else 'FAIL'}")
    results = {}
    for L in lines:
        a_n = an_table(ap[id(L)], N_MAX)
        Pv = []
        for (a, b, d) in keep:
            z0 = complex(-d, 1) / N
            z1 = complex(a, 1) / N
            Pv.append(F_of(a_n, z1) - F_of(a_n, z0))
        Pv = np.array(Pv)
        Pi, res, *_ = np.linalg.lstsq(X, Pv, rcond=None)
        resid = np.linalg.norm(X @ Pi - Pv) / max(np.linalg.norm(Pv), 1e-30)
        results[id(L)] = (Pi, resid, Pv)
    return X, keep, results


# ----------------------------------------------------------------------
# Stage E5 — f1 lattice cross-certificate
# ----------------------------------------------------------------------

def lattice_reduce(w1, w2):
    """Gauss reduction of a rank-2 lattice basis in C."""
    if abs(w1) < abs(w2):
        w1, w2 = w2, w1
    for _ in range(200):
        mu = round((w2 * w1.conjugate()).real / abs(w1) ** 2)
        w2 = w2 - mu * w1
        if abs(w2) >= abs(w1):
            break
        w1, w2 = w2, w1
    tau = w2 / w1
    if tau.imag < 0:
        tau = -tau
    # bring tau to fundamental domain (few steps suffice here)
    for _ in range(50):
        tau = tau - round(tau.real)
        if abs(tau) < 1:
            tau = -1 / tau
        else:
            break
    tau = tau - round(tau.real)
    return w1, w2, tau


def f1_lattice_certificate(lines, results):
    f1p = [L for L in lines if L[1] == "f1" and L[0] == "+"][0]
    f1m = [L for L in lines if L[1] == "f1" and L[0] == "-"][0]
    Pp = results[id(f1p)][2]
    Pm = results[id(f1m)][2]
    Ptot = Pp + Pm            # both functionals evaluated on same loops:
    # the full period of omega_{f1} over each loop is the sum of its
    # +/- functional pieces only up to the embedding convention; the
    # loop periods themselves (Pv arrays) were computed from the SAME
    # a_n, so Pp == Pm elementwise up to numerical noise; use Pp.
    per = Pp
    nz = per[np.abs(per) > 1e-6]
    # find two shortest independent periods
    idx = np.argsort(np.abs(nz))
    w1 = nz[idx[0]]
    w2 = None
    for k in idx[1:]:
        if abs((nz[k] / w1).imag) > 1e-4:
            w2 = nz[k]
            break
    w1r, w2r, tau = lattice_reduce(w1, w2)
    # membership: every period reduces to a lattice point
    Mmat = np.array([[w1r.real, w2r.real], [w1r.imag, w2r.imag]])
    coeffs = np.linalg.solve(Mmat, np.vstack([per.real, per.imag]))
    memb = np.max(np.abs(coeffs - np.round(coeffs)))
    return w1r, w2r, tau, memb, per


def curve_real_period_quad():
    """Independent quadrature: real period of 143a1, a-invariants
    (0,-1,1,-1,-2). Completing the square: y' = y + 1/2 gives
    y'^2 = x^3 - x^2 - x - 7/4, and the invariant differential is
    omega = dx/(2y') -> Omega_re = 2 int_{e1}^{inf} dx/(2 sqrt(cubic))
                                 =   int_{e1}^{inf} dx/sqrt(cubic)."""
    from scipy.integrate import quad
    cubic = np.poly1d([1, -1, -1, -7 / 4])
    roots = np.roots(cubic)
    e1 = float(np.max(roots[np.abs(roots.imag) < 1e-9].real))

    def integrand(x):
        return 1.0 / np.sqrt(cubic(x))

    val, err = quad(integrand, e1, np.inf, limit=400)
    return val, err


# ----------------------------------------------------------------------
# Stage E6 — Born signatures (Paper 36 SS6)
# ----------------------------------------------------------------------

def born_signature(vec26, lines, results, pet_map, block_np):
    """P(f) = |Pi_f . gamma|^2 / <f,f>, per orbit; old block reported
    with declared Euclidean normalization, listed separately."""
    weights = {"f1": 0.0, "f2": 0.0, "f3": 0.0}
    for L in lines:
        # each embedding has + and - functionals; the complex period is
        # Pi+ . v  +  Pi- . v  evaluated per functional — combine as the
        # modulus over the (+,-) pair for this eigenvalue
        pass
    # group lines by (orbit, lambda): the +/- pair forms one embedding
    from collections import defaultdict
    emb = defaultdict(dict)
    for L in lines:
        emb[(L[1], round(L[2], 8))][L[0]] = L
    for (orbit, lam), pair in emb.items():
        amp = 0.0
        for sgn, L in pair.items():
            Pi = results[id(L)][0]
            amp += abs(np.dot(Pi, vec26)) ** 2
        pet = pet_for(orbit, lam)
        weights[orbit] += amp / pet
    old = block_np[0]
    coef, *_ = np.linalg.lstsq(old, vec26, rcond=None)
    w_old_eu = float(np.dot(old @ coef, old @ coef))
    tot = sum(weights.values())
    sig = {k: v / tot for k, v in weights.items()} if tot > 0 else weights
    return sig, w_old_eu


def pet_for(orbit, lam):
    if orbit == "f1":
        return PET_F1
    table = (A2_F2_ORDER, PET_F2_DIAG) if orbit == "f2" \
        else (A2_F3_ORDER, PET_F3_DIAG)
    k = int(np.argmin([abs(lam - a) for a in table[0]]))
    return table[1][k]


# ----------------------------------------------------------------------
# Stage F — Zeno curve on the particle box
# ----------------------------------------------------------------------

def zeno_curve(A, start=None):
    Lm = np.diag(A.sum(axis=1)) - A
    self_loops = int(np.count_nonzero(np.diag(A)))
    if start is None:
        clean = np.nonzero(np.diag(A) == 0)[0]
        start = int(clean[0])
    w, V = np.linalg.eigh(Lm)
    psi0 = np.zeros(A.shape[0])
    psi0[start] = 1.0
    c2 = (V.T @ psi0) ** 2
    taus = np.logspace(-2.5, 1.3, 400)
    amp = np.array([np.sum(c2 * np.exp(-1j * w * t)) for t in taus])
    P1 = np.abs(amp) ** 2
    R = -np.log(np.clip(P1, 1e-300, 1)) / taus
    meanL = float(np.sum(c2 * w))
    var = float(np.sum(c2 * w ** 2) - meanL ** 2)
    return taus, R, var, self_loops, start


# ----------------------------------------------------------------------
# main
# ----------------------------------------------------------------------

def main():
    ledger = {}
    print("=" * 70)
    print("v0.2 STAGE E — periods and Petersson Born weights")
    print("=" * 70)
    p1, tris, edges, ms = build_engine()
    P = float_projection(ms)

    Bc, proj_c, restrict, T2, E, lines = eigendata(ms, P)
    n_plus = 13
    print(f"  E1 eta split of cuspidal 26: 13 + 13  PASS (EXACT)")
    ledger["eta split 13+13 (EXACT)"] = True

    lines = assign_orbits(lines)
    counts = {o: sum(1 for L in lines if L[1] == o) for o in
              ("f1", "f2", "f3")}
    print(f"  E2 new lines by orbit (x2 for +/-): {counts} "
          f"[expect f1:2 f2:8 f3:12]")
    ok_counts = counts == {"f1": 2, "f2": 8, "f3": 12}
    ledger["new line counts (Cert)"] = bool(ok_counts)

    print(f"  E3 extracting a_p for p <= {P_MAX} from this engine's own "
          f"T_p ...")
    ap = extract_ap(ms, P, restrict, lines, P_MAX)

    # cross-certify traces against mtft ORBIT_TRACE tables, n <= 50
    import mtft.x0_143 as ox
    tr_tab = {"f1": ox.ORBIT_TRACE_F1, "f2": ox.ORBIT_TRACE_F2,
              "f3": ox.ORBIT_TRACE_F3}
    max_err = 0.0
    for orbit in ("f1", "f2", "f3"):
        Ls = [L for L in lines if L[1] == orbit and L[0] == "+"]
        ans = [an_table(ap[id(L)], 50) for L in Ls]
        for n in range(1, 51):
            tr = sum(a[n] for a in ans)
            max_err = max(max_err, abs(tr - tr_tab[orbit][n - 1]))
    print(f"  E3 trace cross-cert vs mtft ORBIT_TRACE (n<=50): "
          f"max |err| = {max_err:.2e} "
          f"{'PASS' if max_err < 1e-6 else 'FAIL'} (Cert)")
    ledger["trace cross-cert n<=50 (Cert)"] = float(max_err)

    print("  E4 periods of Gamma_0(143)-loops, solving 26-dim period "
          "functionals ...")
    X, keep, results = period_functionals(ms, P, Bc, proj_c, lines, ap)
    worst = max(results[id(L)][1] for L in lines)
    print(f"  E4 overdetermined residual (worst over 22 functionals): "
          f"{worst:.2e} {'PASS' if worst < 1e-6 else 'FAIL'} (Cert)")
    ledger["period functional residual (Cert)"] = float(worst)

    w1, w2, tau, memb, f1_periods = f1_lattice_certificate(lines, results)
    Omega_quad, qerr = curve_real_period_quad()

    def jinv(t):
        q = cmath.exp(2j * pi * t)
        c = [1, 744, 196884, 21493760, 864299970, 20245856256,
             333202640600, 4252023300096, 44656994071935]
        return sum(cc * q ** (k - 1) for k, cc in enumerate(c))

    j_lattice = jinv(tau)
    # exact j from the package's own a-invariants (0,-1,1,-1,-2):
    # b2=-4 b4=-2 b6=-7 -> c4=64, Delta=-1859=-11*13^2, j = -262144/1859
    j_exact = -262144.0 / 1859.0
    print(f"  E5 f1 lattice: |w1| = {abs(w1):.6f}, |w2| = {abs(w2):.6f}, "
          f"tau = {tau:.6f}")
    print(f"     2pi bridge to Paper 32 mfsymboleval values: "
          f"|w1|/2pi = {abs(w1)/(2*pi):.5f}, |w2|/2pi = {abs(w2)/(2*pi):.5f}"
          f"  [corpus: {CORPUS_F1_ABS}]  "
          f"{'PASS' if abs(abs(w1)/(2*pi)-CORPUS_F1_ABS[0])<2e-4 and abs(abs(w2)/(2*pi)-CORPUS_F1_ABS[1])<2e-4 else 'FAIL'} (Cert)")
    print(f"     lattice membership of all 120 loop periods: max frac dev "
          f"{memb:.2e} {'PASS' if memb < 1e-5 else 'FAIL'} (Cert)")
    print(f"     quadrature real period of 143a1 (omega = dx/2y): "
          f"{Omega_quad:.6f}; ratio to |w1| = {Omega_quad/abs(w1):.6f} "
          f"[Manin-constant-1 expectation: 1] "
          f"{'PASS' if abs(Omega_quad/abs(w1)-1)<1e-5 else 'FAIL'} (Cert)")
    print(f"     j-invariant adjudication: j(lattice tau) = "
          f"{j_lattice.real:+.4f}{j_lattice.imag:+.1e}i ; "
          f"j from a-invariants = -262144/1859 = {j_exact:+.4f} "
          f"{'PASS' if abs(j_lattice - j_exact) < 5e-3 else 'FAIL'} (Cert)")
    print(f"     NOTE (corpus): mtft CURVE_143A1.j_invariant string '-1/15' "
          f"and Paper 32's tau ~ 0.039+0.980i / 'nearly square' remark "
          f"both disagree with the a-invariants and with this lattice; "
          f"flag for corpus correction (receipts in ledger).")
    ledger["f1 lattice"] = {"absw1": float(abs(w1)),
                            "absw2": float(abs(w2)),
                            "tau": [float(tau.real), float(tau.imag)],
                            "absw_over_2pi": [float(abs(w1)/(2*pi)),
                                              float(abs(w2)/(2*pi))],
                            "membership_dev (Cert)": float(memb),
                            "quad_real_period": float(Omega_quad),
                            "quad_over_absw1 (Cert)": float(Omega_quad/abs(w1)),
                            "j_lattice (Cert)": [float(j_lattice.real),
                                                 float(j_lattice.imag)],
                            "j_from_a_invariants (EXACT)": "-262144/1859",
                            "corpus_flags": ["CURVE_143A1.j_invariant '-1/15'",
                                             "Paper32 tau~0.039+0.980i "
                                             "'nearly square'"]}

    # orbit blocks (reuse exact route from v0.1 idea, float here)
    from x0143_particle_box import orbit_blocks
    import sympy as sp
    xs = sp.symbols("x")
    T2q = ms.hecke_on_quotient(2)
    A2, _ = ms.restrict_to_cuspidal(T2q)
    q2s = xs ** 4 - 3 * xs ** 3 - xs ** 2 + 5 * xs + 1
    q3s = xs ** 6 - 10 * xs ** 4 + 2 * xs ** 3 + 24 * xs ** 2 - 7 * xs - 12
    blocks = orbit_blocks(A2, [(xs + 2, "old", 4), (xs, "f1", 2),
                               (q2s, "f2", 8), (q3s, "f3", 12)])
    block_np = [np.array(k[1].tolist(), dtype=float) for k in blocks]

    # E6: implementation check — equal-alpha superposition must reproduce
    # Paper 36's (0.0259, 0.3271, 0.6470)
    tr2, tr3 = sum(PET_F2_DIAG), sum(PET_F3_DIAG)
    tot = PET_F1 + tr2 + tr3
    p36 = (PET_F1 / tot, tr2 / tot, tr3 / tot)
    print(f"  E6 equal-alpha orbit Born weights from corpus norms: "
          f"({p36[0]:.4f}, {p36[1]:.4f}, {p36[2]:.4f}) "
          f"[Paper 36: 0.0259, 0.3271, 0.6470] "
          f"{'PASS' if abs(p36[0]-0.0259)<5e-4 else 'FAIL'} (Cert-corpus)")
    ledger["equal-alpha Born (Cert-corpus)"] = [round(v, 5) for v in p36]

    # Born signatures of demo cycles: Petersson vs Euclidean
    print()
    print("  Born signatures P(orbit) = |period|^2/<f,f> "
          "(Cert-corpus norms; per-embedding f3 order DIAGNOSTIC):")
    demo_states = {}
    for bi, nm in ((1, "pure f1 eigencycle"), (2, "pure f2 eigencycle"),
                   (3, "pure f3 eigencycle")):
        demo_states[nm] = block_np[bi][:, 0]
    rng = np.random.default_rng(143)
    demo_states["random drawn cycle"] = sum(
        rng.integers(-3, 4) * block_np[i][:, k]
        for i in range(4) for k in range(min(2, block_np[i].shape[1])))
    born_out = {}
    for nm, v in demo_states.items():
        sig, w_old = born_signature(v, lines, results, None, block_np)
        born_out[nm] = {k: round(f, 4) for k, f in sig.items()}
        print(f"    {nm}: e/mu/tau = "
              f"({sig['f1']:.4f}, {sig['f2']:.4f}, {sig['f3']:.4f})"
              f"   [old component, declared-Euclid: {w_old:.3f}]")
    ledger["Born signatures (DIAGNOSTIC per-embedding order)"] = born_out

    # Stage F — Zeno
    print()
    print("=" * 70)
    print("v0.2 STAGE F — Zeno curve on the particle box")
    print("=" * 70)
    A, edge_list, tri_of = dual_graph(p1, tris, edges)
    taus, R, var, n_loops, start = zeno_curve(A)
    small = R[taus < 0.02] / taus[taus < 0.02]
    print(f"  dual-graph refinement: {n_loops} self-loop node(s) found "
          f"(triangles glued to themselves across an edge) — the graph is "
          f"3-regular by row sum but not simple; v0.1's certificate is "
          f"amended accordingly (EXACT).")
    print(f"  Zeno start node {start} (loop-free): (Delta L)^2 = "
          f"{var:.12f} [exact law: 3 = degree]  "
          f"{'PASS' if abs(var-3) < 1e-9 else 'FAIL'} (EXACT)")
    print(f"  small-tau law R/tau -> (Delta L)^2: numerical limit "
          f"{small.mean():.6f} {'PASS' if abs(small.mean()-3)<3e-2 else 'FAIL'}")
    kmax = int(np.argmax(R))
    print(f"  R(tau): Zeno freeze as tau->0; maximum R = {R[kmax]:.4f} at "
          f"tau = {taus[kmax]:.3f} (anti-Zeno-style enhancement, DIAGNOSTIC)")
    ledger["zeno"] = {"self_loops (EXACT)": n_loops,
                      "start_node": start,
                      "variance (EXACT)": float(var),
                      "R_max (DIAGNOSTIC)": float(R[kmax]),
                      "tau_at_max": float(taus[kmax])}

    make_figures(taus, R, var, born_out, f1_periods, w1, w2)
    with open(_os.path.join(_HERE, "certificates_v02.json"), "w") as f:
        json.dump(ledger, f, indent=2, default=str)
    print("\nledger written: certificates_v02.json")


def make_figures(taus, R, var, born_out, f1_periods, w1, w2):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # Zeno figure
    fig, ax = plt.subplots(figsize=(8.6, 5))
    ax.semilogx(taus, R, color="#1d3557", lw=2,
                label="effective decay rate R(tau)")
    ax.semilogx(taus, var * taus, "--", color="#e63946",
                label=r"quadratic law $(\Delta L)^2\,\tau$ (exact var = 3)")
    ax.axhline(0, color="#999", lw=0.5)
    ax.set_xlabel(r"measurement interval $\tau$")
    ax.set_ylabel("R")
    ax.set_title("Zeno curve on the $X_0(143)$ particle box — repeated "
                 "projective measurement of a one-triangle state\n"
                 "freeze as $\\tau\\to 0$; enhancement band; "
                 "recurrence oscillations (DIAGNOSTIC model level)")
    ax.legend()
    fig.tight_layout()
    fig.savefig("fig_zeno.png", dpi=170)
    plt.close(fig)
    print("  wrote fig_zeno.png")

    # f1 lattice figure
    fig, ax = plt.subplots(figsize=(6.4, 6))
    pts = np.array([[p.real, p.imag] for p in f1_periods])
    mm, nn = np.meshgrid(range(-3, 4), range(-3, 4))
    lat = mm.flatten() * w1 + nn.flatten() * w2
    ax.scatter([z.real for z in lat], [z.imag for z in lat], s=12,
               c="#bbb", label="lattice  m w1 + n w2")
    ax.scatter(pts[:, 0], pts[:, 1], s=34, c="#e4572e", zorder=3,
               label="engine's f1 loop periods")
    ax.set_aspect("equal")
    ax.legend(fontsize=9)
    ax.set_title("f1 (electron) loop periods land on the 143a1 lattice — "
                 "period-engine cross-certificate")
    fig.tight_layout()
    fig.savefig("fig_f1_lattice.png", dpi=170)
    plt.close(fig)
    print("  wrote fig_f1_lattice.png")

    # Born vs orbit figure
    fig, ax = plt.subplots(figsize=(9, 4.6))
    labels = list(born_out.keys())
    cats = ["f1", "f2", "f3"]
    colors = ["#e4572e", "#17bebb", "#9b5de5"]
    bottoms = np.zeros(len(labels))
    for ci, cat in enumerate(cats):
        vals = np.array([born_out[k][cat] for k in labels])
        ax.bar(range(len(labels)), vals, bottom=bottoms, color=colors[ci],
               label={"f1": "e", "f2": "mu", "f3": "tau"}[cat])
        bottoms += vals
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=12, ha="right", fontsize=9)
    ax.set_ylabel("Born weight (Petersson-normalized)")
    ax.set_title("Paper-36 Born rule on drawn cycles: "
                 r"$P(f)=|{\rm period}_f(\gamma)|^2/\langle f,f\rangle$")
    ax.legend()
    fig.tight_layout()
    fig.savefig("fig_born_weights.png", dpi=170)
    plt.close(fig)
    print("  wrote fig_born_weights.png")


if __name__ == "__main__":
    main()
