# --- stage 5 (Integration Plan v0.1): re-pointed to the mtft package. ---
# internal() delegates to mtft.chain.internal (certified spectrally
# identical: g to 9e-15, full H(u) spectrum to 2.5e-14; B differs only
# by gauge inside near-degenerate blocks).  Where a local gsq() remains
# it is retained VERBATIM: it evaluates complex u in the f64 backend,
# which mtft.ep deliberately does not offer (PR-20 floors) -- S5-1.
# ------------------------------------------------------------------------
#!/usr/bin/env python3
"""
rung5b_tree.py — Space Becomes Arithmetic: the Bruhat-Tits Tree
================================================================
MIT License — (c) 2026 Roger Tano — MTFT Research Program

Rung 5 coupled the marked gas to Z.  Rung 5b replaces Z by the
(p+1)-regular tree = the Bruhat-Tits building of PGL_2(Q_p), whose
adjacency operator IS the Hecke operator T_p.  Space stops being a
stand-in and becomes arithmetic.

MODEL.  H = h(kappa) (x) 1  -  tau A (x) B(kappa),  A = adjacency.
A is not translation-generated (the tree is non-amenable, no Bloch
momentum), but it has an exact spherical-Plancherel decomposition:
the l^2-spectrum of A is the KESTEN band [-2 sqrt p, 2 sqrt p] with

    dmu_q(x) = q sqrt(4(q-1) - x^2) / (2 pi (q^2 - x^2)) dx,  q = p+1.

So the whole rung-5 machinery ports with  2 cos k  ->  x ~ dmu_q:
    H(x) = h(kappa) - tau x B(kappa).

WHY IT IS NOT COSMETIC.  The 1-D measure 1/(pi sqrt(4 - x^2)) DIVERGES
at the band edges (van Hove); the Kesten measure VANISHES there like a
square root.  Every edge-sensitive statement therefore flips:

  Pr Q  purely a.c. spectrum (min-max, B strictly PD) — ports.
  Pr R  SOFT EDGES: rho ~ (x_max - x)^{+1/2} on the tree vs
        (x_max - x)^{-1/2} on Z.  No van Hove singularities at all.
  Pr S  BINDING THRESHOLD (answers PR-7(a) in both geometries):
        on Z the binding integral diverges => ANY attraction binds,
        V_b = 0; on the tree it converges => V_b > 0, finite and
        computable.  The classic low-dimension/high-dimension
        dichotomy, produced here by arithmetic rather than by d.
  Pr T  ARITHMETIC: finite quotients of the tree by the discriminant-143
        quaternion order are the Brandt / supersingular isogeny graphs.
        Their adjacency eigenvalues are {p+1} U {a_p(f)} over the
        newforms of X_0(143) (Eichler / Jacquet-Langlands), and
        Ramanujan (|a_p| <= 2 sqrt p) says exactly that the arithmetic
        spectrum lies INSIDE the Kesten band: the arithmetic spatial
        sector is an optimal approximation to the tree.

Gates: TG0 Kesten measure vs exact closed-walk counts; TG1 a.c.;
TG2 soft-edge exponent two-leg; TG3 binding threshold tree vs Z;
TG4 band-merging tau_c on the tree; TG5 arithmetic (corpus Hecke data
in the Kesten band); TG6 mixing-rate comparison (DIAGNOSTIC).

Run:  py rung5b_tree.py
"""
from __future__ import annotations
import math
import numpy as np
import mpmath as mp
from mtft.chain import internal as _chain_internal


BETA, KSTAR, P = 2.0, 5.0, 2
Q = P + 1
XMAX = 2.0 * math.sqrt(P)
NB = 60
REPORT = []


def rec(name, gtype, value, cls, ok, note=""):
    REPORT.append((name, gtype, value, cls, bool(ok), note))
    print(f"[{'PASS' if ok else 'FAIL'}] {name:<28} {gtype:<12} "
          f"{value:<32} {cls:<20} {note}", flush=True)


def internal(N=1600, kappa=KSTAR, nb=NB, gcap=200.0):
    ic = _chain_internal(kappa, nb=nb, backend="f64", N=N, gcap=gcap)
    return np.asarray(ic.g), np.asarray(ic.B)


G, B = internal()


def kesten(x):
    return Q * np.sqrt(np.maximum(4.0 * (Q - 1) - x ** 2, 0.0)) / (
        2 * math.pi * (Q ** 2 - x ** 2))


def gauss_x(n=400):
    """Kesten nodes via x = XMAX sin(theta): the substitution cancels the
    sqrt endpoint factor, so the integrand is analytic in theta and
    Gauss-Legendre converges spectrally (the raw-x grid stalls at ~1e-6,
    the sqrt being non-polynomial at the edges)."""
    t, w = np.polynomial.legendre.leggauss(n)
    th = 0.5 * math.pi * t
    x = XMAX * np.sin(th)
    dens = Q * XMAX ** 2 * np.cos(th) ** 2 / (
        2 * math.pi * (Q ** 2 - x ** 2))
    return x, w * 0.5 * math.pi * dens


def gauss_z(n=400):
    """1-D chain: x = 2 cos k, uniform in k (smooth parametrization)."""
    t, w = np.polynomial.legendre.leggauss(n)
    k = 0.5 * math.pi * (t + 1.0)
    return 2.0 * np.cos(k), w * 0.5 * math.pi / math.pi


# ------------------------------------------------------------------ TG0
def tg0(nquad=600):
    x, w = gauss_x(nquad)
    norm = float(np.sum(w))
    m = {j: float(np.sum(w * x ** j)) for j in (2, 4, 6)}
    # exact closed-walk counts on the q-regular tree, depth-limited BFS
    depth = 8
    adj = {0: []}
    frontier = [(0, None)]
    nxt = 1
    for d in range(depth):
        newf = []
        for v, par in frontier:
            deg = Q if v == 0 else Q - 1
            for _ in range(deg):
                adj[nxt] = []
                adj[v].append(nxt); adj[nxt].append(v)
                newf.append((nxt, v)); nxt += 1
        frontier = newf
    idx = {v: i for i, v in enumerate(adj)}
    A = np.zeros((len(adj), len(adj)))
    for v, ns in adj.items():
        for u in ns:
            A[idx[v], idx[u]] = 1.0
    e0 = np.zeros(len(adj)); e0[idx[0]] = 1.0
    walks = {}
    vec = e0.copy()
    for step in range(1, 7):
        vec = A @ vec
        if step in (2, 4, 6):
            walks[step] = float(vec[idx[0]])
    devs = [abs(m[j] - walks[j]) for j in (2, 4, 6)]
    rel = max(devs[i] / m[j] for i, j in enumerate((2, 4, 6)))
    ok = abs(norm - 1) < 1e-13 and rel < 1e-13
    rec("TG0 Kesten measure", "Identity",
        f"norm={norm:.12f}; m2,m4,m6={m[2]:.4f},{m[4]:.4f},{m[6]:.4f}",
        "EXACT(walks)", ok,
        f"closed-walk counts on the {Q}-regular tree "
        f"{[walks[j] for j in (2,4,6)]}; max rel dev {rel:.1e}")


# ------------------------------------------------------------------ TG1
def tg1(tau=0.05, nb=12):
    worst = np.inf
    for xv in np.linspace(-XMAX, XMAX, 9):
        w, U = np.linalg.eigh(np.diag(G) - tau * xv * B)
        for i in range(nb):
            worst = min(worst, float(U[:, i] @ (B @ U[:, i])))
    ok = worst > 1e-8
    rec("TG1 a.c. on the tree", "Theorem",
        f"min <phi|B|phi> = {worst:.4f}", "EXACT(HF)", ok,
        f"branches strictly monotone in x over the Kesten support "
        f"[-{XMAX:.4f}, {XMAX:.4f}] => no flat bands; Plancherel is "
        f"dmu_q, not Lebesgue-in-k")


# ------------------------------------------------------------------ TG2
def tg2():
    eps = np.array([10.0 ** e for e in (-5, -5.5, -6, -6.5, -7)])
    rho_t = kesten(XMAX - eps)
    a_t = float(np.polyfit(np.log(eps), np.log(rho_t), 1)[0])
    rho_z = 1.0 / (math.pi * np.sqrt(np.maximum(4.0 - (2.0 - eps) ** 2, 1e-300)))
    a_z = float(np.polyfit(np.log(eps), np.log(rho_z), 1)[0])
    ok = abs(a_t - 0.5) < 1e-3 and abs(a_z + 0.5) < 1e-3
    rec("TG2 soft edge (tree) vs van Hove (Z)", "Theorem",
        f"tree {a_t:+.5f} vs Z {a_z:+.5f}", "CERTIFIED(1e-3)", ok,
        "rho ~ (x_max - x)^(+1/2) on the tree, ^(-1/2) on Z: the edge "
        "singularity is REMOVED by arithmetic geometry")


# ------------------------------------------------------------------ TG3
def tg3(tau=0.05, d=0, nquad=800):
    """Binding threshold: 1 = V G_dd(E), E below band 0's bottom.  The
    defect sits on the band-0 orbital, so the near-edge weight |c_0|^2
    is O(1) and the geometry — not the weight — decides."""
    def Gdd(E, nodes):
        xs, ws = nodes
        tot = 0.0
        for xv, wv in zip(xs, ws):
            ev, U = np.linalg.eigh(np.diag(G) - tau * xv * B)
            tot += wv * float(np.sum(U[d, :] ** 2 / (E - ev)))
        return tot

    nt = gauss_x(nquad)
    nz = gauss_z(4 * nquad)
    e_t = min(np.linalg.eigvalsh(np.diag(G) - tau * XMAX * B)[0],
              np.linalg.eigvalsh(np.diag(G) + tau * XMAX * B)[0])
    e_z = min(np.linalg.eigvalsh(np.diag(G) - tau * 2.0 * B)[0],
              np.linalg.eigvalsh(np.diag(G) + tau * 2.0 * B)[0])
    ds = np.array([1e-3, 3e-4, 1e-4, 3e-5, 1e-5])
    seq_t = np.array([abs(Gdd(e_t - dd, nt)) for dd in ds])
    seq_z = np.array([abs(Gdd(e_z - dd, nz)) for dd in ds])
    # tree: rho ~ sqrt(edge) => G(delta) = G(0) - c sqrt(delta), finite
    ct = np.polyfit(np.sqrt(ds), seq_t, 2)   # G0 - c u + d u^2, u=sqrt(delta)
    G0 = float(ct[2])
    fit_res = float(np.max(np.abs(np.polyval(ct, np.sqrt(ds)) - seq_t))) / G0
    # Z: rho ~ 1/sqrt(edge) => G ~ delta^{-1/2}, no finite limit
    pz = float(np.polyfit(np.log(ds), np.log(seq_z), 1)[0])
    Vb_tree = 1.0 / G0
    ok = fit_res < 5e-3 and abs(pz + 0.5) < 0.02 and G0 > 0
    rec("TG3 binding threshold", "Theorem",
        f"|V_b| tree = {Vb_tree:.5f}, Z = 0", "CERTIFIED(extrap)", ok,
        f"tree G(delta)=G(0)-c sqrt(delta)+O(delta) -> G(0)={G0:.4f} (resid "
        f"{fit_res:.1e}); Z exponent {pz:+.4f} vs -1/2 exactly (no finite "
        f"limit): PR-7(a) answered — ANY attraction binds on Z, a FINITE "
        f"threshold |V_b|={Vb_tree:.5f} on the tree")


# ------------------------------------------------------------------ TG4
def tg4():
    mu0, mu1, m = B[0, 0], B[1, 1], G[1]
    pred = m / (XMAX * (mu0 + mu1))

    def ind_gap(tau):
        lo2 = np.linalg.eigvalsh(np.diag(G) - tau * XMAX * B)[1]
        hi1 = np.linalg.eigvalsh(np.diag(G) + tau * XMAX * B)[0]
        return lo2 - hi1

    lo, hi = 0.25 * pred, 3.0 * pred
    assert ind_gap(lo) > 0 > ind_gap(hi)
    for _ in range(50):
        mid = 0.5 * (lo + hi)
        if ind_gap(mid) > 0: lo = mid
        else: hi = mid
    tc = 0.5 * (lo + hi)
    ratio = tc / 0.230032
    ok = abs(tc / pred - 1) < 0.05
    rec("TG4 tau_c on the tree", "Theorem",
        f"tau_c = {tc:.5f} vs pred {pred:.5f}", "MEASURED vs 1st-order",
        ok, f"Z value 0.230032; ratio {ratio:.4f} vs the geometric "
            f"prediction 2/x_max = {2.0/XMAX:.4f} (band merging happens "
            f"EARLIER on the tree: wider spectral support)")


# ------------------------------------------------------------------ TG5
def tg5():
    import mtft
    OT = mtft.ORBIT_TRACES_VERIFIED
    rows, ok_all = [], True
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23):
        a_p = OT[p][0]                      # dim-1 Galois orbit: trace = a_p
        band = 2.0 * math.sqrt(p)
        inside = abs(a_p) <= band
        ok_all &= inside
        rows.append((p, a_p, band, abs(a_p) / band))
    sat = max(r[3] for r in rows)
    from mtft.modular_curve import X0
    gen = X0(143).genus
    ok = ok_all and gen == 13
    rec("TG5 arithmetic in the band", "Theorem",
        f"9 primes, max |a_p|/(2 sqrt p) = {sat:.3f}", "EXACT(corpus)",
        ok, f"X_0(143) genus {gen}, index 168; dim-1 orbit a_p = "
            f"{[r[1] for r in rows]}; Ramanujan <=> the Brandt/isogeny "
            f"graph spectrum lies inside the Kesten band")


# ------------------------------------------------------------------ TG6
def tg6(tau=0.05, nquad=400):
    xt, wt = gauss_x(nquad)
    xz, wz = gauss_z(nquad)
    def spread(xs, ws):
        e = np.array([np.linalg.eigvalsh(np.diag(G) - tau * xv * B)[0]
                      for xv in xs])
        mean = float(np.sum(ws * e))
        return math.sqrt(abs(float(np.sum(ws * (e - mean) ** 2))))
    st, sz = spread(xt, wt), spread(xz, wz)
    rec("TG6 band-0 dispersion width", "Diagnostic",
        f"tree {st:.5f} vs Z {sz:.5f}", "DIAGNOSTIC", True,
        f"ratio {st/sz:.4f}; tree support is sqrt(p) wider but its "
        f"measure is edge-depleted — the two effects partly cancel")


if __name__ == "__main__":
    print("=" * 104)
    print(f"  RUNG 5b — THE BRUHAT-TITS TREE   p={P}, q={Q}-regular, "
          f"Kesten band [-{XMAX:.4f}, {XMAX:.4f}]   [space becomes "
          f"arithmetic]")
    print("=" * 104)
    tg0(); tg1(); tg2(); tg3(); tg4(); tg5(); tg6()
    print("-" * 104)
    n = sum(1 for r in REPORT if r[4])
    print(f"  {n}/{len(REPORT)} gates green")
    print("=" * 104)
