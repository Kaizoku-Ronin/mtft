"""mtft.surface.rrspace — Riemann–Roch spaces on X0(143) with poles at CM points (SM-02 / step 5, v0.30.1).

Exact facts used:
- div(f1 dz) = the 24 Atkin–Lehner fixed points (f1 = the (+,+) differential, i.e. the 143a1 newform): a
  weight-2 form even under an involution vanishes at its fixed points (automorphy factor −1 there), and
  4 + 20 = 24 = deg K.  Verified numerically to 1e-10 with the 400-term expansions.
- f_K = g1^2 w has cuspidal divisor 2A, A = 6(0) + 6(1/11), w the eta quotient with exponents (2,0,−2,0).
Construction: sections of O(A + P1 + P2 + P3) (three W13-fixed CM points) are phi = G / (f1 f_K) with G in
H^0(K^2) (quadratic canonical-ring elements) satisfying ord_q(G) >= 8 at cusps 0 and 1/11 (exact, via AL
signs), G = 0 at the 20 W143-fixed points and at the excluded W13 point P4.  36 − 33 = 3: the sections are 1,
u and a third with simple poles at the two points of the triple sharing a sign of u — the purity theorem's
bookkeeping.  Point conditions are numerical (height-maximising reduction + 400-term q-series, 1e-10).
"""
from __future__ import annotations

import itertools
from typing import Dict

import numpy as np

from .yukawa import load_basis, al_eigenbasis, conv, rank_mod
from . import petersson as PT

W11 = (66, 5, 143, 11)


def _mob(M, t): return (M[0] * t + M[1]) / (M[2] * t + M[3])


def w13_fixed_representatives(amax: int = 40):
    P = []
    for a in range(-amax, amax + 1):
        for c in range(1, amax + 1):
            num = -1 - 13 * a * a
            if num % (11 * c) == 0:
                t = (26 * a + 1j * np.sqrt(52.0)) / (2 * 143 * c); P.append(complex(t.real - np.round(t.real), t.imag))
    return np.array(P)


def w143_fixed_representatives(cmax: int = 14):
    return np.array([complex(a / c, 1 / (c * np.sqrt(143.0))) for c in range(1, cmax + 1) for a in range(c) if (1 + 143 * a * a) % c == 0])


def eigenforms_at(points, prec: int = 400):
    """Values of the 13 AL-eigen weight-2 forms at points of H (weight-2 automorphy included), via reduction."""
    forms, labels, dens, vecs = al_eigenbasis(load_basis(), prec=prec)
    zr, mult, n11, n13 = PT.reduce_points(np.asarray(points, dtype=complex), 1); q = np.exp(2j * np.pi * zr)
    F = []
    for i in range(13):
        v = np.zeros(len(zr), dtype=complex); qp = np.ones(len(zr), dtype=complex)
        for n in range(1, prec + 1): qp = qp * q; v += int(forms[i][n]) * qp
        F.append(mult * (labels[i][0] ** n11) * (labels[i][1] ** n13) * v)
    return np.array(F), labels


def cm_classes(prec: int = 400) -> Dict:
    """The 20 W143-fixed and 4 W13-fixed classes (W11-partners generated explicitly) with eigenform values."""
    Z0 = w143_fixed_representatives(); P0 = w13_fixed_representatives()
    F, labels = eigenforms_at(np.concatenate([Z0, P0]), prec)
    def classes(cols, idxs, tol=1e-6):
        inv = np.array([F[idxs[1], cols] / F[idxs[0], cols], F[idxs[2], cols] / F[idxs[0], cols]]); cl = []
        for j in range(len(cols)):
            for c in cl:
                if np.linalg.norm(inv[:, j] - inv[:, c[0]]) < tol * max(1, np.linalg.norm(inv[:, c[0]])): c.append(j); break
            else: cl.append([j])
        return [cols[c[0]] for c in cl]
    odd143 = [i for i in range(13) if labels[i][0] * labels[i][1] == -1]; odd13 = [i for i in range(13) if labels[i][1] == -1]
    zi = classes(list(range(len(Z0))), odd143); pi = classes(list(range(len(Z0), len(Z0) + len(P0))), odd13)
    allpts = np.concatenate([Z0, P0])
    Z = np.concatenate([allpts[zi], [_mob(W11, allpts[k]) for k in zi]]); P = np.concatenate([allpts[pi], [_mob(W11, allpts[k]) for k in pi]])
    FZ, _ = eigenforms_at(Z, prec); FP, _ = eigenforms_at(P, prec)
    return {"Z": Z, "P": P, "FZ": FZ, "FP": FP, "labels": labels, "f1_max_rel_at_Z": float((np.abs(FZ[0]) / np.abs(FZ).max(axis=0)).max())}


def three_family_sections(exclude: int = 0, prec: int = 400) -> Dict:
    """H^0(S0 (x) O(P1+P2+P3)) as G/(f1 f_K): returns the basis of quadratic products, the nullspace (3 columns),
    the singular values (rank certificate), and |G_k| at the four W13 classes."""
    cm = cm_classes(prec); labels = cm["labels"]
    forms, _, _, _ = al_eigenbasis(load_basis(), prec=130); PREC = 130
    pairs = list(itertools.combinations_with_replacement(range(13), 2)); prods = [conv(forms[i], forms[j]) for i, j in pairs]
    r, piv = rank_mod(prods, PREC); basis = [pairs[c] for c in piv]; B = [prods[c] for c in piv]
    sgn = {"0": [labels[i][0] * labels[i][1] for i in range(13)], "1/11": [labels[i][1] for i in range(13)]}
    rows = [[s[i] * s[j] * int(B[k][n]) for k, (i, j) in enumerate(basis)] for s in sgn.values() for n in range(2, 8)]; ncusp = len(rows)
    FZ, FP = cm["FZ"], cm["FP"]
    for col in range(FZ.shape[1]): rows.append([FZ[i, col] * FZ[j, col] for (i, j) in basis])
    rows.append([FP[i, exclude] * FP[j, exclude] for (i, j) in basis])
    A = np.array(rows, dtype=complex); A /= np.abs(A).max(axis=1, keepdims=True)
    u, s, vh = np.linalg.svd(A); null = vh[-3:].conj().T
    GP = np.abs(np.array([[FP[i, c] * FP[j, c] for (i, j) in basis] for c in range(FP.shape[1])]) @ null); GP /= GP.max()
    return {"basis": basis, "nullspace": null, "singular_values": s, "rank": int(np.sum(s > 1e-8)), "dimension": A.shape[1] - int(np.sum(s > 1e-8)),
            "G_at_W13_classes": GP, "excluded_index": exclude, "f1_max_rel_at_Z": cm["f1_max_rel_at_Z"]}


# ------------------------------------------------ SM-03: derivatives, Higgs target space, up-type Yukawa of M1 (v0.30.2)
import json as _json
from pathlib import Path as _Path
from .petersson import AL, _best_gamma

def ext_inv(dv,cv): return pow(int(dv)%int(cv),-1,int(cv))
def reduce_with_matrix(z, k, K=150):
    z=z.astype(complex); n=len(z); Mm=np.tile(np.eye(2),(n,1,1)); jtot=np.ones(n,dtype=complex); det=np.ones(n); n11=np.zeros(n,int); n13=np.zeros(n,int)
    for it in range(30):
        sh=np.round(z.real); z=z-sh
        T=np.zeros((n,2,2)); T[:,0,0]=1; T[:,1,1]=1; T[:,0,1]=-sh; Mm=T@Mm
        cands=[]
        for Q,(a,b,c,d) in [(1,(1,0,0,1))]+list(AL.items()):
            zq=(a*z+b)/(c*z+d); im,bk,bd,bv=_best_gamma(zq,K); cands.append((Q,(a,b,c,d),zq,im,bk,bd,bv))
        ims=np.array([cd[3] for cd in cands]); choice=np.argmax(ims,axis=0); gain=ims[choice,np.arange(n)]>z.imag*(1+1e-9)
        if not gain.any(): break
        for ci,(Q,(a,b,c,d),zq,im,bk,bd,bv) in enumerate(cands):
            idx=np.where(gain&(choice==ci))[0]
            if len(idx)==0: continue
            if Q!=1:
                j=c*z[idx]+d; jtot[idx]*=j; det[idx]*=Q; z[idx]=zq[idx]; W=np.array([[a,b],[c,d]],float); Mm[idx]=W@Mm[idx]
                if Q in (11,143): n11[idx]+=1
                if Q in (13,143): n13[idx]+=1
            sub=idx[bv[idx]<1]
            if len(sub):
                cc=143*bk[sub]; d0=bd[sub]; a0=np.array([ext_inv(dv,cv) for dv,cv in zip(d0,cc)],dtype=np.int64); b0=(a0*d0-1)//cc
                j=cc*z[sub]+d0; jtot[sub]*=j; z[sub]=(a0*z[sub]+b0)/j
                G=np.zeros((len(sub),2,2)); G[:,0,0]=a0; G[:,0,1]=b0; G[:,1,0]=cc; G[:,1,1]=d0; Mm[sub]=G@Mm[sub]
    return z, Mm, jtot, det, n11, n13
def eigen_vals_and_derivs(points, prec=400):
    forms,labels,dens,vecs=al_eigenbasis(load_basis(),prec=prec)
    w,Mm,jt,det,n11,n13=reduce_with_matrix(np.asarray(points,dtype=complex),1); q=np.exp(2j*np.pi*w)
    cM=Mm[:,1,0]; F=[]; dF=[]
    for i in range(13):
        v=np.zeros(len(w),dtype=complex); dv=np.zeros(len(w),dtype=complex); qp=np.ones(len(w),dtype=complex)
        for nn in range(1,prec+1): qp=qp*q; v+=int(forms[i][nn])*qp; dv+=nn*int(forms[i][nn])*qp
        dv*=2j*np.pi; eps=(labels[i][0]**n11)*(labels[i][1]**n13)
        F.append(det*eps*jt**(-2)*v); dF.append(det*eps*(-2*cM*jt**(-3)*v + jt**(-2)*(det/jt**2)*dv))
    return np.array(F), np.array(dF), labels


def eta_quotient_series(exponents, nterms: int, N: int = 143):
    """q-expansion of prod eta(d tau)^{r_d}: returns (offset, coeffs) with the series = q^offset * sum coeffs[m] q^m."""
    D = [d for d in range(1, N + 1) if N % d == 0]
    off = sum(d * r for d, r in zip(D, exponents)) / 24
    assert abs(off - round(off)) < 1e-12; off = int(round(off))
    s = np.zeros(nterms); s[0] = 1
    for d, r in zip(D, exponents):
        for n_ in range(1, nterms):
            if d * n_ >= nterms: break
            for _ in range(abs(r)):
                new = s.copy()
                if r > 0: new[d * n_:] -= s[:nterms - d * n_]
                else:
                    for m in range(1, nterms // (d * n_) + 1): new[d * n_ * m:] += s[:nterms - d * n_ * m]
                s = new
    return off, s


def _mul(a, b, n):
    out = np.zeros(n, dtype=complex)
    for i in range(min(n, len(a))):
        if a[i] != 0: out[i:] += complex(a[i]) * np.array([complex(x) for x in b[:n - i]])
    return out


def f_K_series(prec: int = 130):
    """f_K = g1^2 w, w = eta(tau)^2/eta(13 tau)^2 (exponents (2,0,-2,0), offset -1): holomorphic weight-2 form with
    cuspidal divisor 2A; coefficients a_0..a_prec (a_0 = 0, a_1 = 1)."""
    g1 = _json.loads((_Path(__file__).parent / "_data" / "x0143_weight1_g1.json").read_text())["coefficients"][: prec + 2]
    off, qw = eta_quotient_series((2, 0, -2, 0), prec + 2); assert off == -1
    h = _mul(g1, g1, prec + 2)
    return _mul(h[1:], qw, prec + 1)


def higgs_target_space(cm: Dict = None, prec: int = 400) -> Dict:
    """H^0(K(2 SigmaP3)) as G'/f1^2, G' in H^0(K^3) (cubic products) vanishing to order 2 at the 20 W143 points and
    at P4: 60 − 42 = 18.  Returns the cubic basis, the 18-column nullspace and the rank certificate."""
    cm = cm or cm_classes(prec); Z, P = cm["Z"], cm["P"]
    F, dF, labels = eigen_vals_and_derivs(np.concatenate([Z, P]), prec); nz = len(Z)
    forms, _, _, _ = al_eigenbasis(load_basis(), prec=130); PREC = 130
    trip = list(itertools.combinations_with_replacement(range(13), 3)); cub = [conv(conv(forms[i], forms[j]), forms[k]) for i, j, k in trip]
    r, piv = rank_mod(cub, PREC); cbasis = [trip[c] for c in piv]
    rows = []
    for c in list(range(nz)) + [nz]:
        rows.append([F[i, c] * F[j, c] * F[k, c] for (i, j, k) in cbasis])
        rows.append([dF[i, c] * F[j, c] * F[k, c] + F[i, c] * dF[j, c] * F[k, c] + F[i, c] * F[j, c] * dF[k, c] for (i, j, k) in cbasis])
    A = np.array(rows, dtype=complex); A /= np.abs(A).max(axis=1, keepdims=True)
    u, s, vh = np.linalg.svd(A)
    return {"cubic_basis": cbasis, "nullspace": vh[-18:].conj().T, "rank": int(np.sum(s > 1e-8)), "dimension": A.shape[1] - int(np.sum(s > 1e-8)), "singular_values": s}


def up_yukawa_M1(prec: int = 400) -> Dict:
    """Up-type Yukawa tensor of model M1: Q_i u^c_j H_k = coordinates of G_i G_j in f_K * (Higgs target basis).
    Returns Y (3 x 3 x 18, symmetric in i,j) and the expansion residual (must be ~1e-12)."""
    cm = cm_classes(prec); sec = three_family_sections(0, prec); tgt = higgs_target_space(cm, prec)
    forms, _, _, _ = al_eigenbasis(load_basis(), prec=130); PREC = 130
    QB = [conv(forms[i], forms[j]) for (i, j) in sec["basis"]]; CB = [conv(conv(forms[i], forms[j]), forms[k]) for (i, j, k) in tgt["cubic_basis"]]
    fK = f_K_series(PREC)
    Gser = [sum(sec["nullspace"][k, i] * np.array([complex(x) for x in QB[k]]) for k in range(len(QB))) for i in range(3)]
    Gp = [sum(tgt["nullspace"][k, j] * np.array([complex(x) for x in CB[k]]) for k in range(len(CB))) for j in range(18)]
    basis8 = np.array([_mul(fK, g, PREC + 1) for g in Gp]).T
    Y = np.zeros((3, 3, 18), dtype=complex); res = 0.0
    for i in range(3):
        for j in range(i, 3):
            rhs = _mul(Gser[i], Gser[j], PREC + 1); coef = np.linalg.lstsq(basis8, rhs, rcond=None)[0]; Y[i, j] = Y[j, i] = coef
            res = max(res, float(np.linalg.norm(basis8 @ coef - rhs) / np.linalg.norm(rhs)))
    return {"Y": Y, "residual": res, "three_family": sec, "target": tgt}


def w13_grading_and_texture(Yres: Dict, prec: int = 130) -> Dict:
    """W13 acts on H^0(S0(SigmaP)) by T(phi) = u (phi o W13) (T^2 = -1/13 on functions; in the slash
    normalisation used here T^2 = -13) and on H^0(K(2 SigmaP)) by the weight-6 slash.  Returns the grades of the
    three sections and the 18 Higgs directions and the selection rule of the Yukawa tensor (SM-04): entries with
    grade_i grade_j grade_k = +1 vanish; a W13-even Higgs direction gives a rank-2 mass matrix."""
    Y = Yres["Y"]; sec = Yres["three_family"]; tgt = Yres["target"]
    forms, labels, _, _ = al_eigenbasis(load_basis(), prec=prec); n = prec + 1; R = n - 7
    s13 = np.array([l[1] for l in labels]); Gnull = sec["nullspace"]; Tnull = tgt["nullspace"]
    Qmat = np.array([np.array([complex(x) for x in conv(forms[i], forms[j])]) for (i, j) in sec["basis"]]).T
    off, us = eta_quotient_series((-1, -1, 1, 1), n + 8)
    def div_u(series):
        a = series[6:]; b = us[:len(a)]; out = np.zeros(len(a), dtype=complex)
        for m in range(len(a)): out[m] = (a[m] - sum(out[k] * b[m - k] for k in range(m))) / b[0]
        return out[:R]
    sq = np.array([s13[a] * s13[b] for (a, b) in sec["basis"]]); T3 = np.zeros((3, 3), dtype=complex); res = 0.0
    for i in range(3):
        TG = div_u(Qmat @ (Gnull[:, i] * sq)); coef = np.linalg.lstsq(Qmat[:R], TG, rcond=None)[0]; res = max(res, float(np.linalg.norm(Qmat[:R] @ coef - TG) / np.linalg.norm(TG)))
        c3 = np.linalg.lstsq(Gnull, coef, rcond=None)[0]; T3[:, i] = c3; res = max(res, float(np.linalg.norm(Gnull @ c3 - coef) / np.linalg.norm(coef)))
    c = (T3 @ T3)[0, 0]; ev3, V3 = np.linalg.eig(T3 / np.sqrt(c))
    sc = np.array([s13[a] * s13[b] * s13[c_] for (a, b, c_) in tgt["cubic_basis"]]); T18 = np.zeros((18, 18), dtype=complex)
    for k in range(18): T18[:, k] = np.linalg.lstsq(Tnull, Tnull[:, k] * sc, rcond=None)[0]
    ev18, V18 = np.linalg.eig(T18)
    Yg = np.einsum("ai,bj,abk,kl->ijl", V3, V3, Y, np.linalg.inv(V18).T); g3 = np.round(ev3.real).astype(int); g18 = np.round(ev18.real).astype(int)
    prod = np.array([[[g3[i] * g3[j] * g18[l] for l in range(18)] for j in range(3)] for i in range(3)]); mag = np.abs(Yg) / np.abs(Yg).max()
    return {"T_squared": complex(c), "closure_residual": res, "grades_sections": g3.tolist(), "grades_higgs": g18.tolist(),
            "max_on_forbidden": float(mag[prod == 1].max()), "max_on_allowed": float(mag[prod == -1].max()), "Y_graded": Yg,
            "rank_even_higgs": int(np.linalg.matrix_rank(np.sum(Yg[:, :, g18 == 1], axis=2), tol=1e-8)), "rank_odd_higgs": int(np.linalg.matrix_rank(np.sum(Yg[:, :, g18 == -1], axis=2), tol=1e-8))}


# ------------------------------------------------ SM-06: the chi_13-twisted (down/lepton) sector (v0.30.4)
def load_chi13():
    """Frozen: bases of S_k(143, chi_13) for k = 2, 4, 6 (PARI mfbasis, coefficients a_0..a_1200, integers) and the
    Atkin–Lehner matrices W_11, W_13, W_143 on the weight-4 and weight-6 spaces (columns = images of basis forms).
    Cusp expansions of a character form at the cusps 1/13, 1/11, 0 are, up to constants, the q-expansions of
    F|W_11, F|W_13, F|W_143 — the route validated by the product test (PARI's mfslashexpansion route was not)."""
    d = np.load(_Path(__file__).parent / "_data" / "x0143_chi13_spaces.npz")
    return {k: d[k] for k in d.files}


def _chi13(d): return 1 if pow(int(d) % 13, 6, 13) == 1 else -1


def reduce_gamma0(z, K_: int = 200):
    """Height-maximising reduction by Gamma_0(143) ONLY (no AL steps: character forms), returning the reduced
    point, the composite integer matrix and the product of automorphy factors."""
    M = [[1, 0], [0, 1]]; jt = 1 + 0j
    for it in range(60):
        n_ = int(np.round(z.real)); z = z - n_; M = [[M[0][0] - n_ * M[1][0], M[0][1] - n_ * M[1][1]], [M[1][0], M[1][1]]]
        best = (1.0, None)
        for kk in range(1, K_):
            c = 143 * kk; d = int(np.round(-c * z.real))
            if d % 13 == 0 or d % 11 == 0: continue
            v = abs(c * z + d)
            if v < best[0]: best = (v, (c, d))
        if best[1] is None: break
        c, d = best[1]; a = pow(d % c, -1, c); b = (a * d - 1) // c; j = c * z + d; jt *= j; z = (a * z + b) / j
        M = [[a * M[0][0] + b * M[1][0], a * M[0][1] + b * M[1][1]], [c * M[0][0] + d * M[1][0], c * M[0][1] + d * M[1][1]]]
    return z, M, jt


def chi13_values(points, k: int, F):
    """Values and derivatives of chi_13-forms of weight k at points of H via Gamma_0-only reduction."""
    vals = np.zeros((len(F), len(points)), dtype=complex); ders = np.zeros_like(vals)
    for pi, p in enumerate(points):
        w, M, jt = reduce_gamma0(complex(p)); q = np.exp(2j * np.pi * w); nn = np.arange(F.shape[1]); qp = q ** nn; ch = _chi13(M[1][1])
        for i in range(len(F)):
            fw = np.sum(F[i] * qp); dfw = 2j * np.pi * np.sum(nn * F[i] * qp)
            vals[i, pi] = ch * jt ** (-k) * fw; ders[i, pi] = ch * (-k * M[1][0] * jt ** (-k - 1) * fw + jt ** (-k - 2) * dfw)
    return vals, ders


def twisted_spaces_chi13(cm: Dict = None) -> Dict:
    """H^0(S0(SigmaP) (x) t3) (dim 3) and H^0(K(2 SigmaP) (x) t3) (dim 18) as chi_13-forms over the same denominators;
    cusp conditions via the AL matrices, point conditions at the 20 W143 points and P4."""
    d = load_chi13(); cm = cm or cm_classes(); Z, P = cm["Z"], cm["P"]; pts = list(Z) + list(P); nz = len(Z)
    out = {}
    for k, target in ((4, 3), (6, 18)):
        F = d[f"F{k}"].astype(float); B = F[:, :20].T; W = {Q: d[f"W{k}_{Q}"] for Q in (11, 13, 143)}
        V, Dv = chi13_values(pts, k, F); rows = []
        if k == 4:
            for Q, nmax in ((143, 7), (13, 7), (11, 1)):
                E = B @ W[Q]
                for n_ in range(1, nmax + 1): rows.append(list(E[n_, :]))
            rows.append(list(B[1, :]))
            for pi in list(range(nz)) + [nz]: rows.append(list(V[:, pi]))
        else:
            for Q in (143, 13, 11):
                E = B @ W[Q]
                for n_ in (1, 2): rows.append(list(E[n_, :]))
            for n_ in (1, 2): rows.append(list(B[n_, :]))
            for pi in list(range(nz)) + [nz]: rows.append(list(V[:, pi])); rows.append(list(Dv[:, pi]))
        A = np.array(rows, dtype=complex); A /= np.abs(A).max(axis=1, keepdims=True); u, s, vh = np.linalg.svd(A)
        out[k] = {"nullspace": vh[-target:].conj().T, "rank": int(np.sum(s > 1e-8)), "conditions": A.shape[0], "dimension": A.shape[1] - int(np.sum(s > 1e-8)), "smallest_sv": float(s[-1])}
    return out


def down_yukawa_M1(prec: int = 400) -> Dict:
    """Down-type (and, by the transposition theorem, charged-lepton) Yukawa tensor of M1 with L_b = L_a (x) t3:
    Q_i (untwisted) x d^c_j (chi_13-twisted) -> f_K x (twisted Higgs target); 3 x 3 x 18, residual ~1e-9."""
    d = load_chi13(); cm = cm_classes(prec); sec = three_family_sections(0, prec); tw = twisted_spaces_chi13(cm)
    forms, _, _, _ = al_eigenbasis(load_basis(), prec=130); PREC = 130; n = PREC + 1
    QB = [conv(forms[i], forms[j]) for (i, j) in sec["basis"]]
    Gser = [sum(sec["nullspace"][k, i] * np.array([complex(x) for x in QB[k]]) for k in range(len(QB))) for i in range(3)]
    Gt = [sum(tw[4]["nullspace"][k, j] * d["F4"][k][:n].astype(complex) for k in range(40)) for j in range(3)]
    Gpt = [sum(tw[6]["nullspace"][k, l] * d["F6"][k][:n].astype(complex) for k in range(68)) for l in range(18)]
    fK = f_K_series(PREC); basis8 = np.array([_mul(fK, g, n) for g in Gpt]).T
    Y = np.zeros((3, 3, 18), dtype=complex); res = 0.0
    for i in range(3):
        for j in range(3):
            rhs = _mul(Gser[i], Gt[j], n); coef = np.linalg.lstsq(basis8, rhs, rcond=None)[0]; Y[i, j] = coef
            res = max(res, float(np.linalg.norm(basis8 @ coef - rhs) / np.linalg.norm(rhs)))
    return {"Y": Y, "residual": res, "twisted": tw, "three_family": sec}


# ------------------------------------------------ SM-07: the cubic-character twist and the lepton tensor (v0.30.5)
_W3 = np.exp(2j * np.pi / 3); _W6 = np.exp(2j * np.pi / 6)


def load_cubic():
    """Frozen: S_4(143, chi) for the cubic Conrey character 133 (coefficients a + b omega, omega = e^{2 pi i/3}) and
    the sextic character 56 = chi_13 * chi_133-bar (a + b zeta_6), a_0..a_1200, with their Atkin–Lehner matrices
    (W_11 same space; W_13, W_143 into the conjugate space, whose basis is the conjugate basis)."""
    d = np.load(_Path(__file__).parent / "_data" / "x0143_cubic_spaces.npz")
    return {133: d["A133"] + d["B133"] * _W3, 100: np.conj(d["A133"] + d["B133"] * _W3), 56: d["A56"] + d["B56"] * _W6, 23: np.conj(d["A56"] + d["B56"] * _W6),
            "W": {(133, 11): d["W133_11"], (133, 13): d["W133_13"], (133, 143): d["W133_143"], (56, 11): d["W56_11"], (56, 13): d["W56_13"], (56, 143): d["W56_143"]}}


def _dlog13(d):
    d = int(d) % 13; x = 1
    for k in range(12):
        if x == d: return k
        x = (x * 2) % 13


def conrey_chi(a: int, d):
    """chi_133(d) = omega^{log_2 d}, chi_100 its conjugate, chi_56 = chi_13 * chi_133-bar, chi_23 = chi_13 * chi_133
    (convention fixed by the direct-vs-reduced evaluation test, residual 1e-14 vs 1.7 for the conjugate)."""
    k = _dlog13(d); c3 = _W3 ** k; s = (-1) ** (k % 2)
    return {133: c3, 100: np.conj(c3), 56: s * np.conj(c3), 23: s * c3}[a]


def cubic_values(a: int, points, k: int = 4):
    F = load_cubic()[a]; V = np.zeros((len(F), len(points)), dtype=complex)
    for pi, p in enumerate(points):
        w, M, jt = reduce_gamma0(complex(p)); q = np.exp(2j * np.pi * w); qp = q ** np.arange(F.shape[1])
        V[:, pi] = np.conj(conrey_chi(a, M[1][1])) * jt ** (-k) * np.array([np.sum(f * qp) for f in F])
    return V


def cubic_family_space(a: int, cm: Dict = None) -> Dict:
    """H^0(S0(SigmaP) (x) t) for the cubic (a = 133) or sextic (a = 56) twist: chi-forms of weight 4 with cusp
    orders (8,8,2,2) via the AL matrices (targets in the conjugate space) and simple vanishing at the 21 CM points."""
    d = load_cubic(); conj_of = {133: 100, 56: 23}; cm = cm or cm_classes(); Z, P = cm["Z"], cm["P"]; pts = list(Z) + list(P); nz = len(Z)
    F = d[a]; B = F[:, :20].T; Bc = d[conj_of[a]][:, :20].T; V = cubic_values(a, pts); rows = []
    E13 = Bc @ d["W"][(a, 13)]; E143 = Bc @ d["W"][(a, 143)]; E11 = B @ d["W"][(a, 11)]
    for n_ in range(1, 8): rows.append(list(E143[n_, :]))
    for n_ in range(1, 8): rows.append(list(E13[n_, :]))
    rows.append(list(E11[1, :])); rows.append(list(B[1, :]))
    for pi in list(range(nz)) + [nz]: rows.append(list(V[:, pi]))
    A = np.array(rows, dtype=complex); A /= np.abs(A).max(axis=1, keepdims=True); u, s, vh = np.linalg.svd(A)
    return {"nullspace": vh[-3:].conj().T, "rank": int(np.sum(s > 1e-8)), "dimension": A.shape[1] - int(np.sum(s > 1e-8)), "smallest_sv": float(s[-1])}


def lepton_yukawa_M1(prec: int = 400) -> Dict:
    """Charged-lepton tensor of M1 with L_b = L_a (x) t3 and L_d = L_c (x) t_cubic: L_i (chi_133-twisted) x
    e^c_j (chi_56-twisted) -> f_K x (chi_13 Higgs space); 3 x 3 x 18, residual ~1e-8."""
    d = load_cubic(); cm = cm_classes(prec); tw = twisted_spaces_chi13(cm); d13 = load_chi13()
    L = cubic_family_space(133, cm); E = cubic_family_space(56, cm); PREC = 130; n = PREC + 1
    Gpt = [sum(tw[6]["nullspace"][k, l] * d13["F6"][k][:n].astype(complex) for k in range(68)) for l in range(18)]
    fK = f_K_series(PREC); basis8 = np.array([_mul(fK, g, n) for g in Gpt]).T
    GL = [sum(L["nullspace"][k, i] * d[133][k][:n] for k in range(40)) for i in range(3)]; GE = [sum(E["nullspace"][k, j] * d[56][k][:n] for k in range(40)) for j in range(3)]
    Y = np.zeros((3, 3, 18), dtype=complex); res = 0.0
    for i in range(3):
        for j in range(3):
            rhs = _mul(GL[i], GE[j], n); coef = np.linalg.lstsq(basis8, rhs, rcond=None)[0]; Y[i, j] = coef; res = max(res, float(np.linalg.norm(basis8 @ coef - rhs) / np.linalg.norm(rhs)))
    return {"Y": Y, "residual": res, "L_space": L, "E_space": E}


def family_equivalence_residual(Ya, Yb, iters: int = 200) -> float:
    """Alternating least-squares fit Ya ~ A (x) B . Yb over GL(3) x GL(3) on the family indices (Higgs index fixed)."""
    A_ = np.eye(3, dtype=complex); B_ = np.eye(3, dtype=complex)
    for it in range(iters):
        T = np.einsum("jb,abk->ajk", B_, Yb); A_ = np.linalg.lstsq(T.reshape(3, -1).T, Ya.reshape(3, -1).T, rcond=None)[0].T
        T = np.einsum("ia,ajk->ijk", A_, Yb); B_ = np.linalg.lstsq(np.transpose(T, (1, 0, 2)).reshape(3, -1).T, np.transpose(Ya, (1, 0, 2)).reshape(3, -1).T, rcond=None)[0].T
    fit = np.einsum("ia,jb,abk->ijk", A_, B_, Yb); return float(np.linalg.norm(fit - Ya) / np.linalg.norm(Ya))


# ------------------------------------------------ SM-13: cubic-class Higgs spaces and mixed twist assignments (model M2)
def load_cubic_higgs():
    """Frozen S_6(143, chi_100) basis to q^1200 (the chi_133 basis is its conjugate) with the Atkin–Lehner matrices of
    both weight-6 spaces and of the weight-4 chi_100 space (W_Q maps chi_100 <-> chi_133 for Q = 13, 143)."""
    import os
    d = np.load(os.path.join(os.path.dirname(__file__), "_data", "x0143_cubic_higgs.npz"))
    return {k: d[k] for k in d.files}

def cubic_family_space_100(cm=None):
    """H^0(S0 (x) O(P1+P2+P3) (x) t_100): chi_100 weight-4 forms with ord >= 7 at cusps 0, 1/11, ord >= 1 at 1/13, infinity,
    vanishing at the 20 W143 points and P4.  Built from the chi_100 space's own AL matrices (conjugating chi_133 sections
    reflects the CM points and gives a different bundle — rejected by the product test, SM-12)."""
    cm = cm or cm_classes(); d = load_cubic(); H = load_cubic_higgs(); F100, F133 = d[100], d[133]
    pts = list(cm["Z"]) + list(cm["P"]); nz = len(cm["Z"]); V = cubic_values(100, pts)
    B = F100[:, :20].T; Bc = F133[:, :20].T; rows = []
    for Q, Bt in ((143, Bc), (13, Bc)):
        E = Bt @ H[f"W{Q}_100_w4"]
        for n_ in range(1, 8): rows.append(list(E[n_, :]))
    rows.append(list((B @ H["W11_100_w4"])[1, :])); rows.append(list(B[1, :]))
    for pi in list(range(nz)) + [nz]: rows.append(list(V[:, pi]))
    A = np.array(rows, dtype=complex); A /= np.abs(A).max(axis=1, keepdims=True); u_, s, vh = np.linalg.svd(A); rank = int(np.sum(s > 1e-8))
    return {"nullspace": vh[-3:].conj().T, "rank": rank, "conditions": A.shape[0], "smallest_sv": float(s[-1]), "F": F100}

def cubic_higgs_space(a, cm=None):
    """H^0(K(2 sum P) (x) t_a) for a in {100, 133}: weight-6 chi_a forms with ord >= 3 at all cusps and double zeros at the
    20 W143 points and P4 (50 conditions on the 68-dimensional space; dim 18)."""
    cm = cm or cm_classes(); H = load_cubic_higgs(); F100 = H["F100_w6"]; F = F100 if a == 100 else np.conj(F100); Fc = np.conj(F)
    pts = list(cm["Z"]) + list(cm["P"]); nz = len(cm["Z"]); dim = F.shape[0]
    V = np.zeros((dim, len(pts)), dtype=complex); D = np.zeros_like(V)
    for pi, p in enumerate(pts):
        w, M, jt = reduce_gamma0(complex(p)); q = np.exp(2j * np.pi * w); nn = np.arange(F.shape[1]); qp = q ** nn; ch = np.conj(conrey_chi(a, M[1][1]))
        fw = F @ qp; dfw = 2j * np.pi * ((nn * F) @ qp); V[:, pi] = ch * jt ** (-6) * fw; D[:, pi] = ch * (-6 * M[1][0] * jt ** (-7) * fw + jt ** (-8) * dfw)
    rows = []; B = F[:, :20].T; Bc = Fc[:, :20].T
    for Q, Bt in ((143, Bc), (13, Bc), (11, B)):
        E = Bt @ H[f"W{Q}_{a}"]
        for n_ in (1, 2): rows.append(list(E[n_, :]))
    for n_ in (1, 2): rows.append(list(B[n_, :]))
    for pi in list(range(nz)) + [nz]: rows.append(list(V[:, pi])); rows.append(list(D[:, pi]))
    A = np.array(rows, dtype=complex); A /= np.abs(A).max(axis=1, keepdims=True); u_, s, vh = np.linalg.svd(A); rank = int(np.sum(s > 1e-8))
    return {"nullspace": vh[-18:].conj().T, "rank": rank, "conditions": A.shape[0], "smallest_sv": float(s[-1]), "F": F}

def _family_series(kind, cm, prec):
    n = prec + 1
    if kind == 0:
        from .yukawa import load_basis, al_eigenbasis, conv
        f, _, _, _ = al_eigenbasis(load_basis(), prec=prec); three = three_family_sections(); qb = three["basis"]; G = three["nullspace"]
        QB = [np.array([complex(x) for x in conv(f[i], f[j])[:n]]) for (i, j) in qb]
        return [sum(G[k, i] * QB[k] for k in range(len(qb))) for i in range(3)]
    if kind == 133:
        sp = cubic_family_space(133, cm); return [sum(sp["nullspace"][k, i] * load_cubic()[133][k][:n] for k in range(40)) for i in range(3)]
    if kind == 100:
        sp = cubic_family_space_100(cm); return [sum(sp["nullspace"][k, i] * sp["F"][k][:n] for k in range(40)) for i in range(3)]
    raise ValueError(kind)

def mixed_yukawa(kind_A, kind_B, higgs_a, cm=None, prec=130):
    """Yukawa tensor of the pair (family class kind_A) x (family class kind_B) -> f_K x (chi_{higgs_a} Higgs space), kinds in
    {0 (untwisted), 133, 100}; character neutrality requires chi_A chi_B = chi_{higgs_a}.  Returns Y (3x3x18) and the
    expansion residual (the product test).  Model M2: up = (133, 0, 133), down = lepton = (133, 133, 100)."""
    cm = cm or cm_classes(); n = prec + 1; GA = _family_series(kind_A, cm, prec); GB = _family_series(kind_B, cm, prec)
    Hs = cubic_higgs_space(higgs_a, cm); Gp = [sum(Hs["nullspace"][k, l] * Hs["F"][k][:n] for k in range(68)) for l in range(18)]; fK = f_K_series(prec)
    def mul(a, b):
        out = np.zeros(n, dtype=complex)
        for i in range(n):
            if a[i] != 0: out[i:] += a[i] * b[:n - i]
        return out
    basis8 = np.array([mul(fK, g) for g in Gp]).T; Y = np.zeros((3, 3, 18), dtype=complex); res = 0.0
    for i in range(3):
        for j in range(3):
            rhs = mul(GA[i], GB[j]); coef = np.linalg.lstsq(basis8, rhs, rcond=None)[0]; Y[i, j] = coef
            res = max(res, float(np.linalg.norm(basis8 @ coef - rhs) / np.linalg.norm(rhs)))
    return {"Y": Y, "residual": res, "higgs": Hs}
