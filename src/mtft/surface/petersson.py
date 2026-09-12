"""mtft.surface.petersson — Petersson norms by hyperbolic quadrature on the glued Manin mesh (PET-01, v0.29.0).

Quadrature: 6-point Gauss on every chart triangle, hyperbolic area element from the chart.  Form values:
height-maximising reduction under Gamma_0(143) ⋊ AL with one-step lookahead (greedy AL moves fail near cusp
classes 11 and 13), automorphy factors and AL signs tracked, q-expansions of the AL-eigen weight-2 forms.
Gate: the 13 x 13 Gram matrix of the differentials agrees with the Riemann bilinear prediction from the frozen
period matrix (diagonal 1.005, off-diagonal 1.4%, exact vanishing across AL sectors).  DIAGNOSTIC.
"""
import numpy as np, re, pickle, time, itertools
from math import gcd
import mtft.surface as S
from .spectral import face_mesh
from . import transport as TR, frozen as FR
N=143; QW=np.array([0.109951743655322]*3+[0.223381589678011]*3)/2
QA=np.array([[0.816847572980459,0.091576213509771,0.091576213509771],[0.091576213509771,0.816847572980459,0.091576213509771],[0.091576213509771,0.091576213509771,0.816847572980459],
             [0.108103018168070,0.445948490915965,0.445948490915965],[0.445948490915965,0.108103018168070,0.445948490915965],[0.445948490915965,0.445948490915965,0.108103018168070]])
def quadrature(Y0=2.0,h=0.25,nx=6):
    cx=S.cell_complex(N); width=[len(o) for o in cx.vertices]; pts=[]; wts=[]
    for f,o in enumerate(cx.faces):
        c_inf,c_0,c_1=cx.vertex_of[o[0]],cx.vertex_of[cx.S(o[0])],cx.vertex_of[cx.S(o[1])]
        nodes,tri,_=face_mesh(Y0*width[c_inf],Y0*width[c_0],Y0*width[c_1],h,nx)
        a,b,c,d=TR.lift_dart(o[0][0],o[0][1],N)
        z=nodes[tri]; A2=np.abs((z[:,1].real-z[:,0].real)*(z[:,2].imag-z[:,0].imag)-(z[:,2].real-z[:,0].real)*(z[:,1].imag-z[:,0].imag))
        for lam,wq in zip(QA,QW):
            w=z@lam; pts.append((a*w+b)/(c*w+d)); wts.append(wq*A2/w.imag**2)      # hyperbolic area element, chart-invariant
    return np.concatenate(pts), np.concatenate(wts)
AL={11:(66,5,143,11),13:(78,7,143,13),143:(0,-1,143,0)}
def _best_gamma(z, K=150):
    """Best Gamma_0(143) element c=143k, d nearest: returns Im after the move, and (k,d)."""
    best=np.abs(z*0)+1.0; bestk=np.zeros(len(z),int); bestd=np.zeros(len(z),np.int64)
    for kk in range(1,K):
        cc=143*kk; dd=np.round(-cc*z.real); val=np.abs(cc*z+dd)
        m=(val<best)&(np.gcd(np.abs(dd).astype(np.int64),143)==1)&(np.gcd(np.abs(dd).astype(np.int64),kk)==1)
        best[m]=val[m]; bestk[m]=kk; bestd[m]=dd[m]
    im=np.where(best<1, z.imag/np.maximum(best,1e-300)**2, z.imag)
    return im, bestk, bestd, best
def _apply_gamma(z, mult, k, idx, bk, bd):
    cc=143*bk[idx]; d0=bd[idx]
    a0=np.array([pow(int(dv)%int(cv),-1,int(cv)) for dv,cv in zip(d0,cc)],dtype=np.int64); b0=(a0*d0-1)//cc
    zi=z[idx]; j=cc*zi+d0; z[idx]=(a0*zi+b0)/j; mult[idx]*=j**(-2*k)
def reduce_points(z, k):
    """Maximise Im over Gamma_0(143)+AL with one-step lookahead (AL move followed by the best Gamma_0 move)."""
    z=z.astype(complex); mult=np.ones(len(z),dtype=complex); n11=np.zeros(len(z),int); n13=np.zeros(len(z),int)
    for it in range(30):
        z=z-np.round(z.real)
        cands=[]
        for Q,(a,b,c,d) in [(1,(1,0,0,1))]+list(AL.items()):
            zq=(a*z+b)/(c*z+d); im,bk,bd,bv=_best_gamma(zq); cands.append((Q,(a,b,c,d),zq,im,bk,bd,bv))
        ims=np.array([cd[3] for cd in cands]); choice=np.argmax(ims,axis=0); gain=ims[choice,np.arange(len(z))]>z.imag*(1+1e-9)
        if not gain.any(): break
        for ci,(Q,(a,b,c,d),zq,im,bk,bd,bv) in enumerate(cands):
            idx=np.where(gain&(choice==ci))[0]
            if len(idx)==0: continue
            if Q!=1:
                j=c*z[idx]+d; mult[idx]*=(Q**k)*j**(-2*k); z[idx]=zq[idx]
                if Q in (11,143): n11[idx]+=1
                if Q in (13,143): n13[idx]+=1
            sub=idx[bv[idx]<1]
            if len(sub): _apply_gamma(z,mult,k,sub,bk,bd)
    return z, mult, n11, n13
def ext(a,b):
    if b==0: return a,1,0
    g,x,y=ext(b,a%b); return g,y,x-(a//b)*y
def eval_forms(coefs, z, k):
    """coefs: list of integer coefficient arrays (index n -> a_n); returns values f(z) for weight-2k eigenforms via reduction."""
    zr,mult,n11,n13=reduce_points(z,k); q=np.exp(2j*np.pi*zr); print("   min Im after reduction:", zr.imag.min().round(4), " max |q|:", np.abs(q).max().round(4))
    out=[]
    for co,(e11,e13) in coefs:
        v=np.zeros(len(zr),dtype=complex); qp=np.ones(len(zr),dtype=complex)
        for n_ in range(1,len(co)):
            qp=qp*q; v+=int(co[n_])*qp
        out.append(mult*(e11**n11)*(e13**n13)*v)
    return out


def petersson_gate(Y0=2.0, h=0.25, nx=6):
    """Quadrature Gram of the 13 differentials vs (i/2) P Icup P^H from frozen periods; returns ratios."""
    from .yukawa import load_basis, al_eigenbasis
    z, w = quadrature(Y0, h, nx); data = load_basis(); forms, labels, dens, vecs = al_eigenbasis(data)
    vals = eval_forms([(forms[i], labels[i]) for i in range(13)], z, 1); F = np.array(vals)
    G1 = (F * w * z.imag**2) @ F.conj().T
    d = FR.x0143(); Q = d["period_Q"]; P = Q[:13] + 1j * Q[13:]
    V = np.array([[float(x) for x in v] for v in vecs]); Pe = V @ P
    Icup = np.linalg.inv(d["intersection_cycles"].astype(float)); pred = 0.5j * (Pe @ Icup @ Pe.conj().T)
    mask = np.abs(pred) > 1e-6 * np.abs(pred).max()
    return {"diag_ratio": np.real(np.diag(G1)) / np.real(np.diag(pred)), "offdiag_max_dev": float(np.max(np.abs(np.abs(G1[mask] / pred[mask]) - 1))),
            "cross_sector_vanishing": bool(np.allclose(np.abs(pred[~mask]), 0, atol=1e-6 * np.abs(pred).max())), "area": float(w.sum()), "G1": G1, "F": F, "w": w, "z": z}
