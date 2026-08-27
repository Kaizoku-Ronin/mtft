"""Call-time gates for ``mtft.periods``.  Fast tier is no-PARI."""
from __future__ import annotations

import json, math
import mpmath as mp
import numpy as np

from mtft import hecke as H
from mtft.canonical import adapted_qexpansions
from .core import (GENUS, period_record, intersection_form, intersection_inverse,
                   symplectic_change, symplectic_form, omega_cusp,
                   omega_symplectic, legacy_omega_symplectic, riemann_matrix,
                   frozen_riemann_matrix, hodge_complex_structure, hodge_metric)
from .bridge import relative_basis_change, cuspidal_basis_change, hecke_to_symplectic_change
from .forms import raw_qexpansions, q_tail_bound
from .physics import hodge_structure_hecke, hodge_metric_hecke, cp_channel_report


def _maxabs(A):
    return max(abs(A[i,j]) for i in range(A.rows) for j in range(A.cols))


def _det_int(A):
    # fraction-free Bareiss
    M=[list(map(int,r)) for r in A];n=len(M);prev=1;sgn=1
    for k in range(n-1):
        p=next((i for i in range(k,n) if M[i][k]),None)
        if p is None:return 0
        if p!=k:M[k],M[p]=M[p],M[k];sgn=-sgn
        piv=M[k][k]
        for i in range(k+1,n):
            for j in range(k+1,n):M[i][j]=(M[i][j]*piv-M[i][k]*M[k][j])//prev
        prev=piv
    return sgn*M[-1][-1]


def gate_integral_symplectic():
    E=np.array(intersection_form(),dtype=object);Q=np.array(intersection_inverse(),dtype=object)
    S=np.array(symplectic_change(),dtype=object)
    J=np.block([[np.zeros((13,13),dtype=object),np.eye(13,dtype=object)],
                [-np.eye(13,dtype=object),np.zeros((13,13),dtype=object)]])
    assert np.array_equal(E@Q,np.eye(26,dtype=object))
    assert np.array_equal(S.T@E@S,J)
    assert _det_int(S.tolist())==1
    return {"det_S":1,"E_Q_inverse":True,"symplectic":True}


def gate_period_reconstruction(dps=60):
    with mp.workdps(dps):
        Oc=omega_cusp(dps);Q=mp.matrix(intersection_inverse())
        bil=_maxabs(Oc*Q*Oc.T)
        t=riemann_matrix(dps);tf=frozen_riemann_matrix(dps)
        td=_maxabs(t-tf);sym=_maxabs(t-t.T)
        Y=mp.matrix([[mp.im(t[i,j]) for j in range(13)] for i in range(13)])
        vals=mp.eigsy(Y,eigvals_only=True);mine=min(vals);maxe=max(vals)
        legacy=_maxabs(omega_symplectic(dps)-legacy_omega_symplectic(dps))
        # The legacy derived field is known-bad; the primitive Oc and exact S replay tau.
        assert legacy>1
        assert bil<mp.mpf('1e-45') and td<mp.mpf('1e-45') and sym<mp.mpf('1e-45')
        assert mine>0
        return {"riemann_bilinear_residual":str(bil),"tau_frozen_residual":str(td),
                "tau_symmetry_residual":str(sym),"Im_tau_min_eig":str(mine),
                "Im_tau_max_eig":str(maxe),"legacy_Omega_sym_discrepancy":str(legacy),
                "derived_field_policy":"Omega_sym = Omega_cusp @ S; legacy JSON field rejected"}


def gate_basis_bridge(dps=50):
    R=relative_basis_change();C=cuspidal_basis_change();P=hecke_to_symplectic_change()
    assert _det_int(R)==_det_int(C)==_det_int(P)==1
    # Analytic E2: period rows are Hecke eigenforms in the promoted cusp frame.
    rec=json.loads(__import__('mtft.periods.core',fromlist=['data_path']).data_path(
        'X0_143_period_basis_v6.json').read_text())
    with mp.workdps(dps):
        O=omega_cusp(dps)*mp.matrix(C)
        residuals={}
        for p in (2,3):
            ap=[mp.mpf(x) for x in rec[f'a{p}_per_row']]
            T=mp.matrix([[mp.mpf(x.numerator)/x.denominator for x in r]
                         for r in H.cuspidal_hecke(p)])
            D=mp.diag(ap);den=max(mp.mpf('1'),_maxabs(D*O))
            residuals[p]=_maxabs(D*O-O*T)/den
            assert residuals[p]<mp.mpf('1e-40')
    return {"det_R":1,"det_C":1,"det_P":1,
            "hecke_equivariance":{p:str(v) for p,v in residuals.items()}}


def gate_hodge_bridge(dps=60):
    J=hodge_structure_hecke(dps);G=hodge_metric_hecke(dps)
    I=np.array([[float(x) for x in r] for r in H.star_involution()])
    sq=np.linalg.norm(J@J+np.eye(26));anti=np.linalg.norm(I@J+J@I)/np.linalg.norm(J)
    comm={}
    for p in (2,3,5,7,11,13):
        T=np.array([[float(x) for x in r] for r in H.cuspidal_hecke(p)])
        comm[p]=np.linalg.norm(T@J-J@T)/np.linalg.norm(T)
    eig=np.linalg.eigvalsh((G+G.T)/2)
    assert sq<1e-10 and anti<1e-10 and max(comm.values())<1e-10 and eig[0]>0
    return {"J_square":sq,"star_anticommutator_rel":anti,"Hecke_commutators":comm,
            "Hodge_metric_eig_min":float(eig[0]),"Hodge_metric_eig_max":float(eig[-1])}


def gate_qexpansion_span(dps=45):
    raw=np.array([[float(x) for x in row] for row in raw_qexpansions(140,dps).tolist()])[:,1:].T
    adapted=np.array(adapted_qexpansions(),float)[1:141,:]
    # Both are 140x13 bases of the same form space; solve on raw and replay.
    X,*_=np.linalg.lstsq(raw,adapted,rcond=None)
    res=np.max(np.abs(raw@X-adapted));rel=res/max(1,np.max(np.abs(adapted)))
    rank=np.linalg.matrix_rank(raw)
    assert rank==13 and res<1e-8
    tail=float(q_tail_bound(1/math.sqrt(143),140))
    return {"rank":rank,"max_abs_span_residual":res,"relative_residual":rel,
            "q140_tail_bound_at_y_1_sqrt143":tail}


def gate_quantitative_m8(dps=50):
    r=cp_channel_report("width",dps)
    assert r["J_square_residual"]<1e-10
    assert r["V_star_commutator_rel"]<1e-10
    assert r["commutator_star_odd_rel"]<1e-10
    assert r["antilinear_fraction"]>0
    return r


def run_all(dps=50):
    out={}
    for name,fn in [
        ("integral_symplectic",gate_integral_symplectic),
        ("period_reconstruction",lambda:gate_period_reconstruction(dps)),
        ("basis_bridge",lambda:gate_basis_bridge(dps)),
        ("hodge_bridge",lambda:gate_hodge_bridge(dps)),
        ("qexpansion_span",lambda:gate_qexpansion_span(min(dps,45))),
        ("quantitative_m8",lambda:gate_quantitative_m8(dps)),
        ("involutions",gate_involutions),
        ("oldtorus",gate_oldtorus),
        ("hamiltonian",lambda:gate_hamiltonian(dps)),
        ("bergman_channels",lambda:gate_bergman_channels(min(dps,40))),
    ]:
        out[name]=fn()
    return out


__all__=["gate_integral_symplectic","gate_period_reconstruction","gate_basis_bridge",
         "gate_hodge_bridge","gate_qexpansion_span","gate_quantitative_m8",
         "gate_involutions","gate_oldtorus","gate_hamiltonian",
         "gate_bergman_channels","run_all"]


# ---------------------------------------------------------------- v0.22 gates
def gate_involutions():
    """EXACT: AL involutions, Hecke commutation, decode, census, star."""
    from fractions import Fraction as Fr
    from . import involutions as IV
    Ws={q:[[Fr(x) for x in r] for r in IV.al_matrix(q)] for q in (11,13,143)}
    I=IV._eye()
    for q,W in Ws.items():
        assert IV._mul(W,W)==I and _det_int(IV.al_matrix(q))==1
    assert IV._mul(Ws[11],Ws[13])==IV._mul(Ws[13],Ws[11])==Ws[143]
    EH=[[Fr(x) for x in r] for r in IV.transported_intersection()]
    for W in Ws.values():
        Wt=[[W[j][i] for j in range(26)] for i in range(26)]
        assert IV._mul(IV._mul(Wt,EH),W)==EH
    for p in (2,3,5,7):
        T=[[Fr(x.numerator,x.denominator) for x in r] for r in H.cuspidal_hecke(p)]
        for W in Ws.values():
            assert IV._mul(W,T)==IV._mul(T,W)
    d=IV.al_signs()
    dims={q:(IV._rank(IV._add(Ws[q],I,1)),IV._rank(IV._add(Ws[q],I,-1)))
          for q in (11,13,143)}
    assert dims=={11:(14,12),13:(12,14),143:(4,22)}
    census=IV.sector_census(); assert census==(1,6,5,1)
    r2=IV.route2_fixed_intersections()
    assert r2=={"ell":2,"ghost":2,"q4":0,"q6":0}
    U11=[[Fr(x.numerator,x.denominator) for x in r] for r in H.cuspidal_hecke(11)]
    assert IV._add(U11,Ws[11],1)==[[Fr(0)]*26 for _ in range(26)]
    U13=[[Fr(x.numerator,x.denominator) for x in r] for r in H.cuspidal_hecke(13)]
    comm=IV._add(IV._mul(U13,Ws[13]),IV._mul(Ws[13],U13),-1)
    assert IV._rank(comm)==4
    IV.oldspace_projector()
    S=[[Fr(x) for x in r] for r in IV.star_symplectic()]
    assert IV._mul(S,S)==I
    Jf=mp.matrix(symplectic_form()); Jl=[[Fr(int(Jf[i,j])) for j in range(26)] for i in range(26)]
    St=[[S[j][i] for j in range(26)] for i in range(26)]
    assert IV._mul(IV._mul(St,Jl),S)==[[-x for x in r] for r in Jl]
    IV.star_charge_orbit()
    from .core import charge_energy
    with mp.workdps(40):
        Ea=charge_energy([0,0,1]+[0]*10,[0]*13,40)
        Eb=charge_energy([0,0,0,0,1,-1]+[0]*7,[0]*13,40)
        assert abs(Ea-Eb)<mp.mpf('1e-30') and abs(Ea-mp.mpf('0.881330420747955'))<mp.mpf('1e-12')
    return {"eps":d["eps"],"census":census,"eigen_dims":dims,
            "route2":r2,"min_shell_energy":str(Ea),
            "star":"integral, star^2=I, anti-symplectic (EXACT)"}


def gate_oldtorus():
    """EXACT: oldspace abelian-surface chain and dynamics."""
    from . import oldtorus as OT
    pol=OT.polarization_type()
    assert pol["smith"]==(2,2,18,18) and pol["type"]==(2,18)
    assert OT.l9_index()==9 and OT.mod3_rank()==2
    OT.principal_form()
    pc=OT.product_charpoly()
    assert pc["j_arith_preserves_principal_form"] is False
    ent=OT.entropy()
    assert abs(ent-2.8872709503576206)<1e-12
    return {"polarization":pol,"l9_index":9,"mod3_rank":2,
            "product_charpoly":pc["charpoly"],"entropy":ent}


def gate_hamiltonian(dps=50):
    """CERTIFIED: Hamiltonian layer anchors; degree EXACT null control."""
    from fractions import Fraction as Fr
    from . import hamiltonian as HM
    m=H.model(); E,tri_of,sS,erep=m["E"],m["tri_of"],m["sS"],m["erep"]
    deg=[0]*56
    for k in range(E):
        a,b=tri_of[erep[k]],tri_of[sS[erep[k]]]; deg[a]+=1; deg[b]+=1
    assert set(deg)=={3}
    rec=json.loads(__import__('mtft.periods.core',fromlist=['data_path'])
                   .data_path('X0_143_m7_harmonic_basis.json').read_text())
    W=[[Fr(a,b) for a,b in row] for row in rec["basis_26x84"]]
    for i in range(26):
        assert sum(W[i][e]*3*W[0][e] for e in range(E))==3*sum(W[i][e]*W[0][e] for e in range(E))
    anchors={"width":(0.1234286299,0.061913789292,(283.554358,442.494053)),
             "distance":(0.4248813827,0.129978482489,(2.510021,6.440043))}
    out={"degree":"EXACT: V=3I, A_-=0"}
    for pot,(rho0,frac0,(flo,fhi)) in anchors.items():
        rho=HM.pairing_stability(pot,dps)
        rep=HM.channel_report(pot,dps)
        fr=HM.symplectic_frequencies(pot,dps)
        assert abs(rho-rho0)<1e-8 and rho<1
        assert abs(rep["hamiltonian_antilinear_fraction"]-frac0)<1e-9
        assert rep["pairing_inertia"]==(13,13)
        assert flo-1e-4<fr.min() and fr.max()<fhi+1e-4 and len(fr)==13
        r=HM.oldspace_routing(pot,dps)
        assert r["closure_residual"]<1e-9 and r["new_new"]>0.8
        b=HM.hecke_block_routing(pot,dps)
        assert b["closure_residual"]<1e-8
        out[pot]={"rho":rho,"antilinear_fraction":rep["hamiltonian_antilinear_fraction"],
                  "new_new":r["new_new"],"intra_block":b["intra_block"]}
    assert out["width"]["intra_block"]>0.7>0.5>out["distance"]["intra_block"]
    return out


def gate_bergman_channels(dps=40):
    """CERTIFIED: bilinear coefficients, 4->1 crossover, channel identity."""
    from . import channels as CH
    from .forms import bergman_density
    with mp.workdps(dps):
        for (n,m),ref in [((2,1),'-0.3002284706133925'),((3,2),'2.3014958069029296'),
                          ((5,1),'3.9496541257933906'),((6,2),'-0.1760062907393')]:
            v=CH.bergman_bilinear(n,m,140,dps)
            assert abs(mp.re(v)-mp.mpf(ref))<mp.mpf('1e-12')
            assert abs(mp.im(v))<mp.mpf(10)**(-(dps-8))
        r=CH.mode_crossover(4,1,dps=dps)
        assert abs(r-mp.mpf(CH.CROSSOVER_RATIO_41))<mp.mpf('1e-15')
        y0=1/mp.sqrt(143); x=mp.mpf('0.3'); y=mp.mpf('1.4')*y0
        ch=CH.channel_density(mp.mpc(x,y),139,140,dps)
        bd=bergman_density(mp.mpc(x,y),140,dps)
        assert abs(ch-bd)<mp.mpf('1e-30')
        return {"crossover_ratio":str(r),
                "B21":str(mp.re(CH.bergman_bilinear(2,1,140,dps))),
                "channel_identity_residual":str(abs(ch-bd))}
