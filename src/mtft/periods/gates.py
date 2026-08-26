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
    ]:
        out[name]=fn()
    return out


__all__=["gate_integral_symplectic","gate_period_reconstruction","gate_basis_bridge",
         "gate_hodge_bridge","gate_qexpansion_span","gate_quantitative_m8","run_all"]
