"""PB-01 (2026-10-02): messenger-exchange Yukawas for the jet families of the Pati–Salam branch.

1. The holomorphic route is closed (EXACT).  Tree-level exchange of heavy fermions between superpotential vertices gives W_eff = Q u^c H Phi_1 ... Phi_k / M^k
   with holomorphic vev insertions.  The only charged vev is the Pati–Salam condensate Sigma in (4,1,2) (or its conjugate orientation Sigmabar in (4bar,1,2)):
     * one insertion is never SU(4)-invariant (4 x 4bar x 4 has no singlet);
     * Sigma^n alone is invariant only for n = 0 mod 4, through epsilon_{abcd}; every such invariant antisymmetrises at least three Sigma's in SU(4), hence in
       their SU(2)_R doublet indices — impossible for identical bosons (a doublet index takes two values): it vanishes identically;
     * Sigma Sigmabar is invariant, but the D-term vacuum has one orientation only: V_D fixes |Sigma|^2 - |Sigmabar|^2, the sum is a flat direction of V_D, and
       the positive one-loop masses set the conjugate orientation to zero.
   So every holomorphic effective Yukawa of the jet families vanishes.  Together with YK-01 (jet–jet cup products vanish), the tree-level superpotential gives
   the jet families NO masses at any order in the vevs.
2. The non-holomorphic route (STRUCTURAL, exact where stated).  With one orientation condensed the nonabelian D-terms do not vanish: <D^a> ~ x_* (T^a)_{44}
   with x_* = g_U^2 xi / (g_U^2 + C g^2) (C the Casimir factor of the rank-one direction).  The scalar potential then contains g^2 <D^a> (X^dagger T^a H + h.c.)
   for every heavy mode X of the bidoublet block — a bilinear mixing of the light Higgs with heavy bidoublet bosons, theta_X = g^2 <D^a> c_X / M_X^2 with c_X the
   overlap of X with T^a H.  Heavy bosons are not closed forms, so the jet–jet–X couplings are not cup products and do not vanish: the light families get
       Y_ij = sum_X theta_X  int psi_{Q,i} ^ psi_{u,j} ^ X ,
   local at the jets' cusps: diagonal entries set by the heavy modes' values at cusp i, off-diagonal entries by the tails of the harmonic jet representatives
   (small quark mixing), and overall size y ~ g^2 y_X v_R^2 / M_X^2: an order-one top needs the condensate near the first heavy bidoublet scale.
3. The arithmetic skeleton (NUMERICAL, from the frozen period data).  The values of weight-2 forms at the four cusps are their AL eigenvalues times a_1, so the
   Bergman kernel of S_2(Gamma_0(143)) restricted to the cusps is a Klein-group circulant K(c, c') = sum_chi chi(c) chi(c') N_chi with sector weights
   N_chi = a_1^T G_chi^-1 a_1 (G the Petersson Gram from the frozen periods): N_{++} = 2.485, N_{+-} = 4.438, N_{-+} = 4.502, N_{--} = 1.066.  Its eigenvectors
   are the four characters, its eigenvalues 4 N_chi; the three-cusp principal submatrix (any omitted cusp — Atkin–Lehner) has eigenvalues 6.44, 13.16, 17.88
   (0.36 : 0.74 : 1).  This is the evaluation kernel that enters the cusp values of bidoublet modes built from weight-2 forms; it carries no large hierarchy —
   that must come from the mixing angles and the Higgs alignment.
OPEN: the heavy bidoublet spectrum on X0(143) (first C-excited modes), the overlaps c_X, the harmonic jet representatives and their tails (CKM)."""
import itertools
import math
import numpy as np


# ------------------------------------------------------------------ 1. SU(n) invariants by Weyl integration (exact on a fine enough grid)
def su_invariants(n, a, b, M=None):
    """Number of SU(n) invariants in V^{(x)a} (x) Vbar^{(x)b}, by the Weyl integration formula on a uniform torus grid."""
    M = M or 4 * (a + b + n)
    th = 2 * math.pi * np.arange(M) / M
    grids = np.meshgrid(*([th] * (n - 1)), indexing="ij")
    angles = list(grids) + [-sum(grids)]
    z = [np.exp(1j * t) for t in angles]
    chi = sum(z); chib = np.conj(chi)
    delta2 = np.ones_like(chi, dtype=float)
    for j in range(n):
        for k in range(j + 1, n): delta2 = delta2 * np.abs(z[j] - z[k]) ** 2
    val = np.mean(chi ** a * chib ** b * delta2) / math.factorial(n)
    return int(round(val.real))

def holomorphic_orders(max_n=4):
    """SU(4) invariants of Q u^c H Sigma^n (4 x 4bar x 4^n) and of Q u^c H Sigma Sigmabar; the n = 4 epsilon invariant vanishes for identical bosons."""
    out = {f"Sigma^{n}": su_invariants(4, 1 + n, 1) for n in range(0, max_n + 1)}
    out["Sigma Sigmabar"] = su_invariants(4, 2, 2)
    out["epsilon_with_identical_bosons"] = "every order-4 invariant antisymmetrises >= 3 Sigma's in SU(4), hence in SU(2)_R: zero"
    return out


# ------------------------------------------------------------------ 1b/2. the D-term vacuum of the condensate
C4, C2R = 15 / 8, 3 / 4          # sum_a (T^a_{44})^2 for SU(4), sum_a (T^a_{22})^2 for SU(2)_R (Tr T^a T^b = delta/2)

def d_term_vacuum(xi=1.0, gU=0.5, g=0.7, m2=0.05):
    """Rank-one aligned vevs s = |Sigma|^2, sb = |Sigmabar|^2: V = (gU^2/2)(s - sb - xi)^2 + (g^2/2)(C4 + C2R)(s - sb)^2 + m2 (s + sb).  Minimise exactly:
    x = s - sb is fixed by the D-terms, the sum by m2 > 0, so sb = 0.  Returns the vevs and the residual nonabelian D-term."""
    C = C4 + C2R
    x = (gU ** 2 * xi - m2) / (gU ** 2 + C * g ** 2)                # dV/dx = 0 with sb = 0, s = x (the m2 term adds -m2 to the linear coefficient)
    x = max(x, 0.0)
    D_nonabelian = math.sqrt(C) * x
    DU = x - xi
    return {"Sigma_sq": x, "Sigmabar_sq": 0.0, "residual_nonabelian_D": D_nonabelian, "U1_D": DU, "V": 0.5 * gU ** 2 * DU ** 2 + 0.5 * g ** 2 * C * x ** 2 + m2 * x}

def d_term_vacuum_scan(xi=1.0, gU=0.5, g=0.7, m2=0.05, n=201):
    """Brute-force check over (s, sb) >= 0 that the minimum has sb = 0."""
    C = C4 + C2R; best = None
    for s in np.linspace(0, 2 * xi, n):
        for sb in np.linspace(0, 2 * xi, n):
            V = 0.5 * gU ** 2 * (s - sb - xi) ** 2 + 0.5 * g ** 2 * C * (s - sb) ** 2 + m2 * (s + sb)
            if best is None or V < best[0]: best = (V, s, sb)
    return {"V": best[0], "Sigma_sq": best[1], "Sigmabar_sq": best[2]}

def mixing_angle(gsq_D, c_X, M_X_sq):
    """theta_X = g^2 <D^a> c_X / M_X^2 (tree-level Higgs–heavy-bidoublet mixing)."""
    return gsq_D * c_X / M_X_sq


# ------------------------------------------------------------------ 3. the cusp kernel of weight-2 forms
CUSPS = {"inf": (1, 1), "1/13 = W11(inf)": (-1, 1), "1/11 = W13(inf)": (1, -1), "0 = W143(inf)": (-1, -1)}

def _character(chi, cusp):
    e11, e13 = chi; s11, s13 = cusp
    return (e11 if s11 == -1 else 1) * (e13 if s13 == -1 else 1)

def cusp_kernel():
    from mtft.surface.yukawa import load_basis, al_eigenbasis
    from mtft.surface import frozen as FR
    data = load_basis(); forms, labels, dens, vecs = al_eigenbasis(data)
    d = FR.x0143(verify=False); Q = d["period_Q"]; P = Q[:13] + 1j * Q[13:]
    V = np.array([[float(x) for x in v] for v in vecs]); Pe = V @ P
    Icup = np.linalg.inv(d["intersection_cycles"].astype(float)); G = np.real(0.5j * (Pe @ Icup @ Pe.conj().T))
    a1 = np.array([int(f[1]) for f in forms], float)
    sectors = {}
    for i, l in enumerate(labels): sectors.setdefault(tuple(l), []).append(i)
    N = {chi: float(a1[idx] @ np.linalg.solve(G[np.ix_(idx, idx)], a1[idx])) for chi, idx in sectors.items()}
    names = list(CUSPS); K = np.array([[sum(w * _character(chi, CUSPS[ci]) * _character(chi, CUSPS[cj]) for chi, w in N.items()) for cj in names] for ci in names])
    three = {}
    for omit in names:
        idx = [names.index(c) for c in names if c != omit]; three[omit] = np.linalg.eigvalsh(K[np.ix_(idx, idx)])
    return {"sector_weights": N, "sector_sizes": {k: len(v) for k, v in sectors.items()}, "kernel": K, "cusps": names,
            "eigenvalues_four_cusps": np.linalg.eigvalsh(K), "eigenvalues_three_cusps": three}

def report():
    ck = cusp_kernel()
    return {"holomorphic_orders": holomorphic_orders(), "d_term_vacuum": d_term_vacuum(), "d_term_scan": d_term_vacuum_scan(n=81),
            "sector_weights": ck["sector_weights"], "three_cusp_eigenvalues": {k: list(v) for k, v in ck["eigenvalues_three_cusps"].items()}}
