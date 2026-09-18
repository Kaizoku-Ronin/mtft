"""mtft.surface.hym — the compact uniformising metric of X0(143), Green's functions and Hermitian–Yang–Mills
normalisations (SM-08, v0.30.6).

Liouville: with g = e^{2u} g_hyp of curvature −1, u solves  ∇²_hyp u = e^{2u} − 1  on the cusped surface with the
cusp asymptotics u ≈ −2 pi y/w + log(2 pi y/w) + c, imposed as a Neumann flux ∂_n u = −2 pi Y0 + 1 on every horocycle
(Y_c = Y0 · width).  Newton on the SPEC mesh converges in ~5 steps.  GATE (Gauss–Bonnet): ∫ e^{2u} dA_hyp = 48 pi
(151.37 vs 150.80 at h = 0.2, the mesh's own 0.4% area error).  Green's functions of the compact Laplacian,
Δ_c G = δ_p − 1/A_c, for interior points (point loads) and for cusps (unit Neumann flux through the horocycle),
solved with the constant mode fixed by ∫ G e^{2u} = 0.  HYM metric on O(D): h_D = exp(4 pi Σ m_p G_p); kinetic Gram
matrices by quadrature with the compact area element; the Higgs (K(2ΣP)) sector carries y² e^{−2u} for K.
Status: DIAGNOSTIC — no mesh-refinement certification yet; the point-source Green's functions are not corrected
for the unresolved log singularity at the source; the 18-space precision floor (8.7e-5) limits the smallest
normalised couplings.
"""
import numpy as np, scipy.sparse as sps, scipy.sparse.linalg as spla
import mtft.surface as S
from .spectral import face_mesh, assemble_local, R, R2, Evaluator
from . import transport as TR
from .yukawa import load_basis, al_eigenbasis

def assemble(N=143, Y0=2.0, h=0.2, nx=8):
    cx=S.cell_complex(N); width=[len(o) for o in cx.vertices]; faces=[]; key_index={}; nxt=0; rows=[]; cols=[]; Kv=[]; Mv=[]
    bnd={}   # global node -> lumped hyperbolic boundary length (horocycle)
    for f,o in enumerate(cx.faces):
        c_inf,c_0,c_1=cx.vertex_of[o[0]],cx.vertex_of[cx.S(o[0])],cx.vertex_of[cx.S(o[1])]
        Yi,Yz,Yo=Y0*width[c_inf],Y0*width[c_0],Y0*width[c_1]
        nodes,tri,edges=face_mesh(Yi,Yz,Yo,h,nx); Kl,Ml=assemble_local(nodes,tri); gid=-np.ones(len(nodes),dtype=int)
        for i,d in enumerate(o):
            e,sgn=cx.edge_of[d]; idx,s=edges[i]
            for j,sv in zip(idx,s):
                key=(e,round(float(sgn*sv),9))
                if key not in key_index: key_index[key]=nxt; nxt+=1
                gid[j]=key_index[key]
        for j in np.where(gid<0)[0]: gid[j]=nxt; nxt+=1
        Ks,Ms=sps.coo_matrix(Kl),sps.coo_matrix(Ml)
        rows+= [gid[Ks.row],gid[Ms.row]]; cols+=[gid[Ks.col],gid[Ms.col]]; Kv+=[Ks.data,np.zeros(len(Ms.data))]; Mv+=[np.zeros(len(Ks.data)),Ms.data]
        # horocycle boundary nodes of this face: top (y = Yi), cusp-0 circle, cusp-1 circle; lumped lengths ds = dx/y along the top, mapped by R, R2 isometries
        top=np.where(np.abs(nodes.imag-Yi)<1e-7)[0]; top=top[np.argsort(nodes.real[top])]
        for arc_nodes,Y,fmap in ((top,Yi,None),):
            xs=nodes.real[arc_nodes]; L=np.zeros(len(arc_nodes))
            for k in range(len(arc_nodes)-1): seg=(xs[k+1]-xs[k])/Y; L[k]+=seg/2; L[k+1]+=seg/2
            for k,j in enumerate(arc_nodes): bnd[gid[j]]=bnd.get(gid[j],0)+L[k]
        for (Yc,fmap) in ((Yz,R),(Yo,R2)):
            # nodes on the image horocycle: fmap(x + i Yc); recover parameter x by inverse map, sort, same lengths as top
            inv=(lambda z: 1-1/z) if fmap is R else (lambda z: 1/(1-z))       # inverse of R is R2 and vice versa
            cand=np.array([inv(z) for z in nodes]); on=np.where(np.abs(cand.imag-Yc)<1e-6)[0]; on=on[np.argsort(cand.real[on])]
            xs=cand.real[on]; L=np.zeros(len(on))
            for k in range(len(on)-1): seg=(xs[k+1]-xs[k])/Yc; L[k]+=seg/2; L[k+1]+=seg/2
            for k,j in enumerate(on): bnd[gid[j]]=bnd.get(gid[j],0)+L[k]
        faces.append((nodes,tri,gid))
    G=nxt
    K=sps.coo_matrix((np.concatenate(Kv),(np.concatenate(rows),np.concatenate(cols))),shape=(G,G)).tocsc()
    M=sps.coo_matrix((np.concatenate(Mv),(np.concatenate(rows),np.concatenate(cols))),shape=(G,G)).tocsc()
    b=np.zeros(G)
    for j,L in bnd.items(): b[j]=L
    return cx,faces,K,M,b
def solve_liouville(K,M,b,Y0,iters=40):
    G=K.shape[0]; ones=np.ones(G); gn=-2*np.pi*Y0+1.0                       # hyperbolic normal derivative of u on every horocycle (Y_c = Y0 * width)
    bvec=gn*b
    u=np.zeros(G); Mo=M@ones
    for it in range(iters):
        e=np.exp(2*u); F=K@u+M@e-Mo-bvec; J=(K+2*M.multiply(e[None,:])).tocsc()
        du=spla.spsolve(J,-F); lam=1.0
        while lam>1e-4:
            un=u+lam*du; Fn=K@un+M@np.exp(2*un)-Mo-bvec
            if np.linalg.norm(Fn)<np.linalg.norm(F): break
            lam/=2
        u=un
        if np.linalg.norm(Fn)<1e-10*max(1,np.linalg.norm(Mo)): break
    return u, np.linalg.norm(Fn), it
def locate(cx, faces, ev, z):
    """face index, barycentric weights and the three global node ids for a point z of H (via the Evaluator's reduction)."""
    f,w=ev.locate(z); nodes,tri,gid=faces[f]; P=nodes[tri]
    x0,y0=P[:,0].real,P[:,0].imag; x1,y1=P[:,1].real,P[:,1].imag; x2,y2=P[:,2].real,P[:,2].imag
    det=(x1-x0)*(y2-y0)-(x2-x0)*(y1-y0); l1=((w.real-x0)*(y2-y0)-(x2-x0)*(w.imag-y0))/det; l2=((x1-x0)*(w.imag-y0)-(w.real-x0)*(y1-y0))/det; l0=1-l1-l2
    lam=np.column_stack([l0,l1,l2]); ok=np.all(lam>-1e-9,axis=1); j=int(np.where(ok)[0][0]) if ok.any() else int(np.argmax(lam.min(axis=1)))
    return f, np.clip(lam[j],0,1)/np.clip(lam[j],0,1).sum(), gid[tri[j]]


def solve_neumann(K, rhs, Me):
    """K G = rhs (rhs summing to zero on each component): pin one node per connected component of the mesh graph
    (isolated sliver triangles can occur at fine resolution), solve the SPD reduced system, fix int G e^{2u} = 0."""
    import scipy.sparse as _sps, scipy.sparse.linalg as _spla
    from scipy.sparse.csgraph import connected_components
    n=K.shape[0]; nc,lab=connected_components((abs(K)>1e-14).astype(int),directed=False)
    pins=[int(np.where(lab==c)[0][0]) for c in range(nc)]; keep=np.setdiff1d(np.arange(n),pins)
    r=rhs.copy()
    for c in range(nc): idx=lab==c; r[idx]-=r[idx].mean()
    Kp=K[keep][:,keep].tocsc(); G=np.zeros(n); G[keep]=_spla.splu(Kp).solve(r[keep])
    return G-(Me@G)/Me.sum()


def greens_functions(cx, faces, K, M, b, u, points: dict, cusps: dict, Y0: float = 2.0):
    """points: name -> z in H (interior sources); cusps: name -> width of the cusp (source at that cusp).
    Returns dict of nodal Green's functions and the discrete compact area."""
    G = K.shape[0]; e2u = np.exp(2 * u); Me = M @ e2u; Ac = float(np.ones(G) @ Me); ev = Evaluator(cx, faces); out = {}
    for name, p in points.items():
        f, lam, ids = locate(cx, faces, ev, p); rhs = np.zeros(G); rhs[ids] += lam; rhs -= Me / Ac
        out[name] = solve_neumann(K, -rhs, Me)
    if cusps:
        width = [len(o) for o in cx.vertices]; cusp_nodes = {c: set() for c in range(len(width))}
        for fi, o in enumerate(cx.faces):
            nodes, tri, gid = faces[fi]; c_inf, c_0, c_1 = cx.vertex_of[o[0]], cx.vertex_of[cx.S(o[0])], cx.vertex_of[cx.S(o[1])]
            Yi, Yz, Yo = Y0 * width[c_inf], Y0 * width[c_0], Y0 * width[c_1]
            for j in np.where(np.abs(nodes.imag - Yi) < 1e-7)[0]: cusp_nodes[c_inf].add(gid[j])
            for j in range(len(nodes)):
                if abs((1 - 1 / nodes[j]).imag - Yz) < 1e-6: cusp_nodes[c_0].add(gid[j])
                if abs((1 / (1 - nodes[j])).imag - Yo) < 1e-6: cusp_nodes[c_1].add(gid[j])
        by_width = {width[c]: c for c in range(len(width))}
        for name, w in cusps.items():
            bvec = np.zeros(G)
            for j in cusp_nodes[by_width[w]]: bvec[j] = b[j]
            flux = -bvec / bvec.sum(); rhs = -Me / Ac
            out[name] = solve_neumann(K, flux - rhs, Me)
    return out, Ac


QW=np.array([0.109951743655322]*3+[0.223381589678011]*3)/2
QA=np.array([[0.816847572980459,0.091576213509771,0.091576213509771],[0.091576213509771,0.816847572980459,0.091576213509771],[0.091576213509771,0.091576213509771,0.816847572980459],
             [0.108103018168070,0.445948490915965,0.445948490915965],[0.445948490915965,0.108103018168070,0.445948490915965],[0.445948490915965,0.445948490915965,0.108103018168070]])
def quad_with_nodes(cx, faces):
    Z=[];Wt=[];NID=[];LAM=[];Yc=[]
    for f,o in enumerate(cx.faces):
        nodes,tri,gid=faces[f]; a,b,c,d=TR.lift_dart(o[0][0],o[0][1],143); z=nodes[tri]
        A2=np.abs((z[:,1].real-z[:,0].real)*(z[:,2].imag-z[:,0].imag)-(z[:,2].real-z[:,0].real)*(z[:,1].imag-z[:,0].imag))
        for lam,wq in zip(QA,QW):
            w=z@lam; Z.append((a*w+b)/(c*w+d)); Wt.append(wq*A2/w.imag**2); NID.append(gid[tri]); LAM.append(np.tile(lam,(len(tri),1)))
    return np.concatenate(Z), np.concatenate(Wt), np.concatenate(NID), np.concatenate(LAM)


def gauss_bonnet_gate(N: int = 143, Y0: float = 2.0, h: float = 0.2, nx: int = 8) -> dict:
    cx, faces, K, M, b = assemble(N, Y0, h, nx); u, res, it = solve_liouville(K, M, b, Y0)
    ones = np.ones(K.shape[0]); return {"compact_area": float(ones @ (M @ np.exp(2 * u))), "target": 48 * np.pi, "hyperbolic_area": float(ones @ (M @ ones)),
                                       "newton_iterations": it + 1, "residual": float(res), "u": u, "cx": cx, "faces": faces, "K": K, "M": M, "b": b}


def hyp_dist(z, p): return np.arccosh(1 + np.abs(z - p) ** 2 / (2 * z.imag * p.imag))


def up_sector_mass_ratios(h: float = 0.2, nx: int = 8, Y0: float = 2.0, n_dirs: int = 2000, seed: int = 0, log_correction: bool = True) -> dict:
    """SM-08/09 pipeline: compact metric -> Green's functions -> HYM Gram matrices of the up-sector family and Higgs
    spaces -> bilinear-normalised Yukawa tensor -> singular-value ratios over random unit Higgs directions.
    Returns percentiles of m2/m3 and m1/m3, the family Gram spectrum, and the compact area."""
    from . import rrspace as RR, petersson as PT
    from .transport import lift_dart
    cx, faces, K, M, b = assemble(143, Y0, h, nx); u, res, it = solve_liouville(K, M, b, Y0); Gn = K.shape[0]
    cm = RR.cm_classes(); P = cm["P"]
    gr, Ac = greens_functions(cx, faces, K, M, b, u, {"P1": P[1], "P2": P[2], "P3": P[3]}, {"cusp0": 143, "cusp1/11": 13}, Y0)
    Z, Wt, NID, LAM = quad_with_nodes(cx, faces); interp = lambda v: np.sum(v[NID] * LAM, axis=1); uq = interp(u); e2u = np.exp(2 * uq)
    Gq = {k: interp(gr[k]) for k in gr}
    if log_correction:
        ev = Evaluator(cx, faces)
        for name, p in (("P1", P[1]), ("P2", P[2]), ("P3", P[3])):
            f, w = ev.locate(p); nodes, tri, gid = faces[f]; a_, b_, c_, d_ = lift_dart(cx.faces[f][0][0], cx.faces[f][0][1], 143)
            zc = (d_ * Z - b_) / (-c_ * Z + a_); d = hyp_dist(zc, w); logd = np.log(np.maximum(d, 1e-12)); logn = np.log(np.maximum(hyp_dist(nodes, w), 1e-12))
            loc = {gg: j for j, gg in enumerate(gid)}; corr = np.zeros(len(Z))
            for q in np.where((d < 0.6) & np.isfinite(d))[0]:
                ids = NID[q]
                if all(gg in loc for gg in ids): corr[q] = (logd[q] - np.sum(LAM[q] * logn[[loc[gg] for gg in ids]])) / (2 * np.pi)
            Gq[name] = Gq[name] + corr
    forms, labels, dens, vecs = al_eigenbasis(load_basis(), prec=400); F = np.array(PT.eval_forms([(forms[i], labels[i]) for i in range(13)], Z, 1))
    f130, _, _, _ = al_eigenbasis(load_basis(), prec=130); Fq = np.array([[complex(x) for x in f130[i]] for i in range(13)]).T
    cK = np.linalg.lstsq(Fq, RR.f_K_series(130), rcond=None)[0]; f1 = F[0]; fKv = cK @ F
    Yr = RR.up_yukawa_M1(); sec = Yr["three_family"]; tgt = Yr["target"]; Y = Yr["Y"]
    phi = [sum(sec["nullspace"][k, i] * F[a] * F[b2] for k, (a, b2) in enumerate(sec["basis"])) / (f1 * fKv) for i in range(3)]
    psi = [sum(tgt["nullspace"][k, l] * F[a] * F[b2] * F[c2] for k, (a, b2, c2) in enumerate(tgt["cubic_basis"])) / f1 ** 2 for l in range(18)]
    h_fam = np.exp(4 * np.pi * (6 * Gq["cusp0"] + 6 * Gq["cusp1/11"] + Gq["P1"] + Gq["P2"] + Gq["P3"])); h_H = np.exp(8 * np.pi * (Gq["P1"] + Gq["P2"] + Gq["P3"])); wq = Wt * e2u
    Nf = np.array([[np.sum(wq * h_fam * phi[i] * np.conj(phi[j])) for j in range(3)] for i in range(3)]); NH = np.array([[np.sum(Wt * Z.imag ** 2 * h_H * psi[k] * np.conj(psi[l])) for l in range(18)] for k in range(18)])
    Lf = np.linalg.cholesky((Nf + Nf.conj().T) / 2); LH = np.linalg.cholesky((NH + NH.conj().T) / 2); Li = np.linalg.inv(Lf)
    Yt = np.einsum("ijk,kl->ijl", np.einsum("ia,jb,abk->ijk", Li, Li, Y), LH)
    rng = np.random.default_rng(seed); Rr = []
    for _ in range(n_dirs):
        v = rng.standard_normal(18) + 1j * rng.standard_normal(18); v /= np.linalg.norm(v); s = np.sort(np.linalg.svd(np.einsum("ijk,k->ij", Yt, v), compute_uv=False)); Rr.append(s / s[-1])
    Rr = np.array(Rr); ev = np.linalg.eigvalsh((Nf + Nf.conj().T) / 2)
    return {"m2_over_m3_pct": np.percentile(Rr[:, 1], [5, 50, 95]), "m1_over_m3_pct": np.percentile(Rr[:, 0], [5, 50, 95]), "family_gram_spectrum": ev / ev.max(),
            "compact_area": float(np.ones(Gn) @ (M @ np.exp(2 * u))), "Y_normalised": Yt, "class": "DIAGNOSTIC (two-mesh agreement: Gram spectrum 1%, ratio order of magnitude)"}


# ------------------------------------------------ SM-10: twisted sectors (down, lepton) — vectorised character evaluation
from .petersson import _best_gamma

def reduce_gamma0_vec(z, K_=200, iters=60):
    """Vectorised Gamma_0(143)-only height-maximising reduction; returns reduced points, automorphy product j(M,z),
    and the composite matrix modulo 13 (exact), for character evaluation."""
    z=z.astype(complex); n=len(z); jt=np.ones(n,dtype=complex); M=np.tile(np.array([[1,0],[0,1]]),(n,1,1)).astype(np.int64)   # mod 13
    for it in range(iters):
        sh=np.round(z.real).astype(np.int64); z=z-sh
        T=np.tile(np.array([[1,0],[0,1]]),(n,1,1)).astype(np.int64); T[:,0,1]=(-sh)%13; M=np.einsum('nij,njk->nik',T,M)%13
        im,bk,bd,bv=_best_gamma(z,K_); m=bv<1-1e-12
        if not m.any(): break
        idx=np.where(m)[0]; cc=143*bk[idx]; d0=bd[idx].astype(np.int64)
        a0=np.array([pow(int(dv)%int(cv),-1,int(cv)) for dv,cv in zip(d0,cc)],dtype=object); b0=[(int(a)*int(d)-1)//int(c) for a,d,c in zip(a0,d0,cc)]
        j=cc*z[idx]+d0; jt[idx]*=j; z[idx]=(np.array([int(a) for a in a0])*z[idx]+np.array(b0,dtype=float))/j
        Gm=np.zeros((len(idx),2,2),dtype=np.int64); Gm[:,0,0]=[int(a)%13 for a in a0]; Gm[:,0,1]=[int(b)%13 for b in b0]; Gm[:,1,0]=cc%13; Gm[:,1,1]=d0%13
        M[idx]=np.einsum('nij,njk->nik',Gm,M[idx])%13
    return z, jt, M[:,1,1]


def twisted_section_values(Z, cm=None):
    """Values at points Z of the twisted family/Higgs sections: d^c (chi_13), twisted Higgs (chi_13), L (cubic 133),
    e^c (sextic 56), via the vectorised Gamma_0-only reduction (composite matrix tracked mod 13).  Points that Gamma_0
    alone cannot raise above height 0.002 (cusp regions of 1/11 and 1/13) are flagged: their q-series do not converge and
    they carry < 1% of the kinetic integrals (measured with the untwisted sections)."""
    from . import rrspace as RR
    w, jt, dm = reduce_gamma0_vec(np.asarray(Z, dtype=complex)); q = np.exp(2j * np.pi * w); low = w.imag < 0.002
    dlog = np.array([[k for k in range(12) if pow(2, k, 13) == d][0] if d % 13 else 0 for d in dm])
    w3 = np.exp(2j * np.pi / 3); chi13 = np.where(dlog % 2 == 0, 1.0, -1.0); chi133 = w3 ** dlog; chi56 = chi13 * np.conj(chi133)
    cm = cm or RR.cm_classes(); tw = RR.twisted_spaces_chi13(cm); d13 = RR.load_chi13(); dcub = RR.load_cubic()
    Lsp = RR.cubic_family_space(133, cm); Esp = RR.cubic_family_space(56, cm)
    def series_values(coeffs):
        v = np.zeros(len(q), dtype=complex)
        for c in coeffs[::-1]: v = v * q + c
        return v
    out = {"low_height_mask": low}
    for name, ns, F, ch, k in (("dc", tw[4]["nullspace"], d13["F4"].astype(complex), chi13, 4), ("HT", tw[6]["nullspace"], d13["F6"].astype(complex), chi13, 6),
                                ("L", Lsp["nullspace"], dcub[133], chi133, 4), ("E", Esp["nullspace"], dcub[56], chi56, 4)):
        out[name] = np.array([np.conj(ch) * jt ** (-k) * series_values(ns[:, j] @ F) for j in range(ns.shape[1])])
    return out


# ------------------------------------------------ SM-11: flux-divisor scan (family Gram hierarchy for any degree-3 divisor on the 24 AL fixed points)
def family_gram_for_divisor(D_indices, setup=None):
    """HYM Gram spectrum of H^0(S0 (x) O(D)) for D = three of the 24 Atkin–Lehner fixed points (indices 0..19 = W143
    classes, 20..23 = W13 classes P4, P1, P2, P3 in rrspace.cm_classes order).  The sections are G/(f1 f_K) with G in
    H^0(K^2) vanishing at the other 21 points; the metric is e^{4 pi (6 G_0 + 6 G_{1/11} + sum_D G_p)}.  Returns
    (rank of the 33 conditions, spectrum/max, smallest singular value).  `setup` may be reused across divisors
    (one-time metric, eigenform and quadrature evaluation: ~1 min at h = 0.2)."""
    import itertools
    from . import rrspace as RR, petersson as PT
    from .yukawa import conv, rank_mod
    if setup is None:
        Y0, h, nx = 2.0, 0.2, 8
        cx, faces, K, M, b = assemble(143, Y0, h, nx); u, res, it = solve_liouville(K, M, b, Y0); G = K.shape[0]; Me = M @ np.exp(2 * u); Ac = float(np.ones(G) @ Me)
        gr, _ = greens_functions(cx, faces, K, M, b, u, {}, {"cusp0": 143, "cusp1/11": 13}, Y0); ev = Evaluator(cx, faces)
        Z, Wt, NID, LAM = quad_with_nodes(cx, faces); interp = lambda v: np.sum(v[NID] * LAM, axis=1); e2u = np.exp(2 * interp(u)); Gc = interp(gr["cusp0"]) + interp(gr["cusp1/11"])
        forms, labels, _, _ = al_eigenbasis(load_basis(), prec=400); F = np.array(PT.eval_forms([(forms[i], labels[i]) for i in range(13)], Z, 1))
        f130, _, _, _ = al_eigenbasis(load_basis(), prec=130); Fq = np.array([[complex(x) for x in f130[i]] for i in range(13)]).T; cK = np.linalg.lstsq(Fq, RR.f_K_series(130), rcond=None)[0]
        cm = RR.cm_classes(); allpts = list(cm["Z"]) + list(cm["P"]); allvals = np.concatenate([cm["FZ"], cm["FP"]], axis=1)
        pairs = list(itertools.combinations_with_replacement(range(13), 2)); prods = [conv(f130[i], f130[j]) for i, j in pairs]; r, piv = rank_mod(prods, 130); qb = [pairs[c] for c in piv]; QB = [prods[c] for c in piv]
        sgn = {"0": [l[0] * l[1] for l in labels], "1/11": [l[1] for l in labels]}
        cusp_rows = [[sg[i] * sg[j] * int(QB[k][n]) for k, (i, j) in enumerate(qb)] for sg in sgn.values() for n in range(2, 8)]
        setup = dict(cx=cx, faces=faces, K=K, M=M, Me=Me, Ac=Ac, ev=ev, interp=interp, F=F, f1=F[0], fKv=cK @ F, allpts=allpts, allvals=allvals, qb=qb, cusp_rows=cusp_rows, wq=Wt * e2u, Gc=Gc)
    S = setup; rows = list(S["cusp_rows"])
    for p in range(24):
        if p not in D_indices: rows.append([S["allvals"][i, p] * S["allvals"][j, p] for (i, j) in S["qb"]])
    A = np.array(rows, dtype=complex); A /= np.abs(A).max(axis=1, keepdims=True); U, sv, vh = np.linalg.svd(A); null = vh[-3:].conj().T; rank = int(np.sum(sv > 1e-8))
    phi = [sum(null[k, i] * S["F"][a] * S["F"][b_] for k, (a, b_) in enumerate(S["qb"])) / (S["f1"] * S["fKv"]) for i in range(3)]
    Gsum = np.zeros(len(S["wq"]))
    for p in D_indices:
        f_, lam, ids = locate(S["cx"], S["faces"], S["ev"], S["allpts"][p]); rhs = np.zeros(S["K"].shape[0]); rhs[ids] += lam; rhs -= S["Me"] / S["Ac"]; Gsum += S["interp"](solve_neumann(S["K"], -rhs, S["Me"]))
    hf = np.exp(4 * np.pi * (6 * S["Gc"] + Gsum)); Nf = np.array([[np.sum(S["wq"] * hf * phi[i] * np.conj(phi[j])) for j in range(3)] for i in range(3)])
    evs = np.linalg.eigvalsh((Nf + Nf.conj().T) / 2)
    return rank, evs / evs.max(), float(sv[-1]), setup


# ------------------------------------------------ SM-14: Higgs directions, CKM, and frozen M2 tensors (DIAGNOSTIC, h = 0.2)
def load_m2_tensors():
    """Frozen HYM-normalised Yukawa tensors of model M2 (up: cubic Q x untwisted u^c -> chi_133 Higgs; down = lepton:
    cubic x cubic -> chi_100 Higgs), in Gram-orthonormal family bases and the 14 well-conditioned Higgs directions,
    computed on the h = 0.2 mesh (DIAGNOSTIC).  The neutrino Dirac tensor equals the up tensor (same spaces, same H_u).
    Also the best-fit Higgs directions (masses only) and the joint fit (masses + CKM)."""
    import os
    d = np.load(os.path.join(os.path.dirname(__file__), "_data", "x0143_m2_tensors_h02.npz"))
    return {k: d[k] for k in d.files}

def mass_matrix(Yn, v):
    """M = sum_k v_k Y_n[:, :, k] for a unit Higgs direction v."""
    v = np.asarray(v, dtype=complex); v = v / np.linalg.norm(v); return np.einsum("ijk,k->ij", Yn, v)

def mass_ratios(Yn, v):
    s = np.sort(np.linalg.svd(mass_matrix(Yn, v), compute_uv=False)); return s / s[-1]

def ckm(Yu_n, Yd_n, vu, vd):
    """|V_CKM| = |U_u^dag U_d| from the left (shared Q index) singular vectors of the up and down mass matrices, ordered
    light -> heavy."""
    Uu, _, _ = np.linalg.svd(mass_matrix(Yu_n, vu)); Ud, _, _ = np.linalg.svd(mass_matrix(Yd_n, vd))
    return np.abs(Uu[:, ::-1].conj().T @ Ud[:, ::-1])

def fit_higgs_direction(Yn, target_ratios, seeds=20, rng=None):
    """Find a unit Higgs direction whose (m1/m3, m2/m3) match target_ratios in log10 (Nelder–Mead from random seeds).
    Returns (v, ratios, log10-cost).  With the Higgs direction a flat modulus this is a fit, not a prediction."""
    from scipy.optimize import minimize
    rng = rng or np.random.default_rng(1); n = Yn.shape[2]; t = np.log10(np.asarray(target_ratios, dtype=float)); best = None
    def cost(x):
        v = x[:n] + 1j * x[n:]; r = mass_ratios(Yn, v); return float(np.sum((np.log10(r[:2]) - t) ** 2))
    for _ in range(seeds):
        res = minimize(cost, rng.standard_normal(2 * n), method="Nelder-Mead", options={"maxiter": 4000, "xatol": 1e-8, "fatol": 1e-12})
        if best is None or res.fun < best.fun: best = res
    v = best.x[:n] + 1j * best.x[n:]; v = v / np.linalg.norm(v); return v, mass_ratios(Yn, v), float(best.fun)


# ------------------------------------------------ CC-26: canonical HYM normalisation (correct contraction) and its invariance test
def normalise_yukawa(Y, N_A, N_B, N_H, regularise=1e-9):
    """Yukawa tensor in kinetic-orthonormal bases.  With Gram matrices N = L L^dag (Cholesky) the orthonormal sections are
    phi' = A phi with A = L^{-1} (A N A^dag = I), so Y'_{ijl} = sum_{abk} A_{ia} B_{jb} Y_{abk} (L_H)_{kl}.  The v0.30.6–0.31.1
    releases contracted the TRANSPOSE (index pattern 'ai,bj'), which is not a kinetic normalisation and produced spurious
    hierarchies (CC-26, found by Astra's Atlas audit V0311-C01).  The Higgs index is expressed in the Gram-orthonormal
    directions with eigenvalue > regularise * max (numerically null directions dropped); pass regularise=None for the
    Cholesky factor instead."""
    chol = lambda N: np.linalg.cholesky((N + N.conj().T) / 2)
    A = np.linalg.inv(chol(N_A)); B = np.linalg.inv(chol(N_B)); Yn = np.einsum("ia,jb,abk->ijk", A, B, Y)
    if regularise is None: return np.einsum("ijk,kl->ijl", Yn, chol(N_H))
    evH, VH = np.linalg.eigh((N_H + N_H.conj().T) / 2); keep = evH / evH.max() > regularise; Vk = VH[:, keep] / np.sqrt(evH[keep])
    return np.einsum("ijk,kl->ijl", Yn, N_H @ Vk)

def normalisation_is_basis_invariant(Y, N_A, N_B, N_H, seed=0, tol=1e-6):
    """Re-express both family bases by random nonsingular matrices (Y -> G (x) H Y, N -> G N G^dag) and check that the
    singular-value ratios of the normalised mass matrix are unchanged (measured deviation ~1e-8 from the regularised
    ill-conditioned Higgs Gram; tolerance 1e-6 leaves a 100x margin — CI-margin rule, v0.31.1)."""
    rng = np.random.default_rng(seed); G = rng.standard_normal((3, 3)) + 1j * rng.standard_normal((3, 3)); H = rng.standard_normal((3, 3)) + 1j * rng.standard_normal((3, 3))
    Y2 = np.einsum("ia,jb,abk->ijk", G, H, Y); NA2 = G @ N_A @ G.conj().T; NB2 = H @ N_B @ H.conj().T
    Y1n = normalise_yukawa(Y, N_A, N_B, N_H); Y2n = normalise_yukawa(Y2, NA2, NB2, N_H); v = rng.standard_normal(Y1n.shape[2]) + 1j * rng.standard_normal(Y1n.shape[2])
    s1 = np.sort(np.linalg.svd(np.einsum("ijk,k->ij", Y1n, v), compute_uv=False)); s2 = np.sort(np.linalg.svd(np.einsum("ijk,k->ij", Y2n, v), compute_uv=False))
    return float(np.abs(s1 / s1[-1] - s2 / s2[-1]).max()) < tol


# ------------------------------------------------ SM-15: W13-graded normalisation, degeneracy theorem, epsilon admixture
def w13_graded_normalisation(Yres, N_fam, N_H):
    """Normalise the M1 up tensor in the W13 eigenbases (grading preserved: the HYM metric is W13-invariant, so the Grams
    are block-diagonal in the grading up to the mesh error).  Returns Ye (3x3xn_even), Yo (3x3xn_odd) in
    kinetic-orthonormal family and Gram-orthonormal even/odd Higgs bases, and the gradings."""
    from .rrspace import w13_grading_and_texture
    tex = w13_grading_and_texture(Yres); gf = np.array(tex["grades_sections"]); g18 = np.array(tex["grades_higgs"]); Yg = tex["Y_graded"]; V3 = tex["V_sections"]; V18 = tex["V_higgs"]
    Nf_g = V3.T @ N_fam @ np.conj(V3); NH_g = V18.T @ N_H @ np.conj(V18)
    def blockdiag(N, gr):
        B = np.zeros_like(N)
        for s_ in (1, -1): m = gr == s_; B[np.ix_(m, m)] = N[np.ix_(m, m)]
        return B
    chol = lambda N: np.linalg.cholesky((N + N.conj().T) / 2); A = np.linalg.inv(chol(blockdiag(Nf_g, gf))); Yn = np.einsum("ia,jb,abk->ijk", A, A, Yg)
    def orth(mask):
        Nb = NH_g[np.ix_(mask, mask)]; ev, V = np.linalg.eigh((Nb + Nb.conj().T) / 2); keep = ev / ev.max() > 1e-9; return Nb @ (V[:, keep] / np.sqrt(ev[keep]))
    return {"Ye": np.einsum("ijk,kl->ijl", Yn[:, :, g18 == 1], orth(g18 == 1)), "Yo": np.einsum("ijk,kl->ijl", Yn[:, :, g18 == -1], orth(g18 == -1)), "grades_sections": gf, "grades_higgs": g18,
            "gram_offgrade_family": float(np.abs(Nf_g[gf[:, None] != gf[None, :]]).max() / np.abs(Nf_g).max()), "gram_offgrade_higgs": float(np.abs(NH_g[g18[:, None] != g18[None, :]]).max() / np.abs(NH_g).max())}

def w13_epsilon_scan(Ye, Yo, eps_values=(1e-4, 1e-3, 1e-2, 1e-1, 0.3, 1.0), n=300, seed=0):
    """Singular-value ratios of M = M(v_even) + eps M(v_odd) over random unit directions.  Degeneracy theorem: at eps = 0
    the spectrum is (1, 1, 0) up to mesh noise; the admixture gives m1/m3 ~ eps and m2/m3 = 1 - O(eps), so a
    W13-symmetric vacuum cannot separate the two heavy families (m_c = m_t)."""
    rng = np.random.default_rng(seed); rnd = lambda k: (lambda v: v / np.linalg.norm(v))(rng.standard_normal(k) + 1j * rng.standard_normal(k)); out = {}
    for eps in eps_values:
        R = np.array([(lambda s: s / s[-1])(np.sort(np.linalg.svd(np.einsum("ijk,k->ij", Ye, rnd(Ye.shape[2])) + eps * np.einsum("ijk,k->ij", Yo, rnd(Yo.shape[2])), compute_uv=False))) for _ in range(n)])
        out[eps] = {"median_m1_m3": float(np.median(R[:, 0])), "median_m2_m3": float(np.median(R[:, 1])), "band_m1_m3": (float(np.percentile(R[:, 0], 5)), float(np.percentile(R[:, 0], 95)))}
    return out
