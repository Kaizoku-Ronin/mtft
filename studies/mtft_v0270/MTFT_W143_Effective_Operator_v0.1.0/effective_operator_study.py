#!/usr/bin/env python3
"""W143 finite effective-operator experiment. See frozen PROTOCOL.md."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import sys
import time
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import numpy as np
from scipy.linalg import expm


def norm(a):return float(np.linalg.norm(a))


def complex_parts(m):
    n=len(m)//2;a,b,c,d=m[:n,:n],m[:n,n:],m[n:,:n],m[n:,n:]
    return (a+d+1j*(c-b))/2,(a-d+1j*(c+b))/2


def herm(a):return (a+a.conj().T)/2


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',type=Path,required=True)
    parser.add_argument('--previous',type=Path,required=True)
    args=parser.parse_args();source=args.source.resolve();previous=args.previous.resolve()
    root=Path(__file__).resolve().parent
    sys.path.insert(0,str(source/'src'))
    import mtft
    from mtft.periods.involutions import al_matrix
    started=time.perf_counter()
    path=previous/'coupling_matrices.npz';old=np.load(path)
    r,ri,p=old['R'],old['Ri'],old['P_active_complex']
    wr=r@np.asarray(al_matrix(143),float)@ri
    raw_w,anti=complex_parts(wr)
    j=old['J']
    input_checks={'real_involution':norm(wr@wr-np.eye(26))/np.sqrt(26),
                  'real_selfadjoint':norm(wr-wr.T)/norm(wr),
                  'real_complex_linearity':norm(wr@j-j@wr)/norm(wr),
                  'antilinear_component':norm(anti)/norm(raw_w),
                  'complex_selfadjoint':norm(raw_w-raw_w.conj().T)/norm(raw_w),
                  'complex_involution':norm(raw_w@raw_w-np.eye(13))/np.sqrt(13),
                  'P_idempotence':norm(p@p-p),'P_selfadjoint':norm(p-p.conj().T)}
    assert max(input_checks.values())<1e-9
    # Explicitly remove only the roundoff anti-Hermitian part after gating.
    w=herm(raw_w)
    pe,pu=np.linalg.eigh(herm(p));u=pu[:,pe>.5];f=pu[:,pe<=.5]
    p_reconstructed=u@u.conj().T
    a=herm(u.conj().T@w@u);b=u.conj().T@w@f;d=herm(f.conj().T@w@f)
    left,s,vh=np.linalg.svd(b,full_matrices=True);right=vh.conj().T
    sn=s/norm(w);retained=int(np.sum(sn>=1e-7))
    ambiguous=bool(np.any((sn>1e-9)&(sn<1e-7)))
    assert not ambiguous
    mode_records=[];columns=[]
    for k in range(retained):
        x=u@left[:,k];y=f@right[:,k]
        aj=float(np.vdot(x,w@x).real);dj=float(np.vdot(y,w@y).real)
        bj=float(s[k]);pair=np.column_stack([x,y])
        target=np.array([[aj,bj],[bj,-aj]],complex)
        columns.extend([x,y])
        mode_records.append({'mode':k+1,'a':aj,'b':bj,'d':dj,
            'a_squared_plus_b_squared':aj*aj+bj*bj,'a_plus_d':aj+dj,
            'coupled_pair_invariance_residual':norm(w@pair-pair@target),
            'positive_sign_overlap_rho':(1+aj)/2,
            'principal_angle_degrees':float(np.degrees(np.arccos(np.sqrt((1+aj)/2)))),
            'maximum_toy_complement_norm_fraction':bj*bj})
    ua=u@left[:,retained:];vf=f@right[:,retained:]
    columns.extend(ua.T);columns.extend(vf.T)
    transform=np.column_stack(columns)
    canon=transform.conj().T@w@transform
    predicted=np.zeros((13,13),complex)
    for k,m in enumerate(mode_records):
        predicted[2*k:2*k+2,2*k:2*k+2]=[[m['a'],m['b']],[m['b'],-m['a']]]
    ua_eig=np.linalg.eigvalsh(herm(ua.conj().T@w@ua))
    vf_eig=np.linalg.eigvalsh(herm(vf.conj().T@w@vf))
    # Record uncoupled restrictions before using their measured signs.
    predicted[2*retained:2*retained+ua.shape[1],2*retained:2*retained+ua.shape[1]]=ua.conj().T@w@ua
    predicted[2*retained+ua.shape[1]:,2*retained+ua.shape[1]:]=vf.conj().T@w@vf
    block_checks={'A2_BB_identity':norm(a@a+b@b.conj().T-np.eye(8)),
                  'D2_BstarB_identity':norm(d@d+b.conj().T@b-np.eye(5)),
                  'AB_BD_zero':norm(a@b+b@d),
                  'canonical_transform_unitarity':norm(transform.conj().T@transform-np.eye(13)),
                  'canonical_reconstruction':norm(canon-predicted),
                  'P_reconstruction':norm(p-p_reconstructed),
                  'hermitian_cleanup_size':norm(w-raw_w)}
    i8=np.eye(8);i5=np.eye(5);i13=np.eye(13)

    def response(z):
        full=u.conj().T@np.linalg.solve(z*i13-w,u)
        sigma=b@np.linalg.solve(z*i5-d,b.conj().T)
        denominator=z*i8-a-sigma
        eff=np.linalg.solve(denominator,i8)
        closed=(z*i8+a)/(z*z-1)
        naive=np.linalg.solve(z*i8-a,i8)
        modes=[]
        for k,m in enumerate(mode_records):
            v=left[:,k]
            gf=np.vdot(v,full@v);ge=np.vdot(v,eff@v);gn=np.vdot(v,naive@v)
            modes.append({'full_real':float(gf.real),'full_imag':float(gf.imag),
                          'effective_real':float(ge.real),'effective_imag':float(ge.imag),
                          'naive_real':float(gn.real),'naive_imag':float(gn.imag),
                          'spectral_response_full':float(-gf.imag/np.pi),
                          'spectral_response_naive':float(-gn.imag/np.pi)})
        sv=np.linalg.svd(sigma,compute_uv=False)
        sigma_rank=int(np.sum(sv>=max(sv[0],1)*1e-9))
        row={'z_real':float(z.real),'z_imag':float(z.imag),
             'full_norm':norm(full),'sigma_norm':norm(sigma),
             'sigma_rank_at_relative_floor_1e_9':sigma_rank,
             'complement_condition':float(np.linalg.cond(z*i5-d)),
             'Schur_denominator_condition':float(np.linalg.cond(denominator)),
             'full_condition':float(np.linalg.cond(z*i13-w)),
             'effective_relative_error':norm(eff-full)/norm(full),
             'involution_relative_error':norm(closed-full)/norm(full),
             'naive_relative_error':norm(naive-full)/norm(full),
             'modes':modes}
        return row

    regular=[]
    for eta in [.05,.2,.7]:
        for x in np.linspace(-1.5,1.5,301):regular.append(response(complex(x,eta)))
    for z in [-2+0j,0+0j,2+0j,1j]:regular.append(response(z))
    stress=[];real_poles=[]
    interior=np.linalg.eigvalsh(d)
    interior=interior[np.minimum(abs(interior-1),abs(interior+1))>1e-7]
    for dj in interior:
        k=int(np.argmin([abs(m['d']-dj) for m in mode_records]))
        for power in range(2,13):
            z=complex(dj,10.**(-power));row=response(z)
            row['complement_eigenvalue']=float(dj);row['mode']=k+1
            m=mode_records[k];gcan=(z+m['a'])/(z*z-1)
            row['canonical_channel_response_abs']=float(abs(gcan))
            row['canonical_channel_response_real']=float(gcan.real)
            row['canonical_channel_response_imag']=float(gcan.imag)
            stress.append(row)
        z=complex(dj,0)
        full=u.conj().T@np.linalg.solve(z*i13-w,u)
        closed=(z*i8+a)/(z*z-1)
        v=left[:,k];value=np.vdot(v,full@v)
        real_poles.append({'z':float(dj),'mode':k+1,'full_projected_norm':norm(full),
                           'full_system_condition':float(np.linalg.cond(z*i13-w)),
                           'channel_response_abs':float(abs(value)),
                           'involution_relative_error':norm(closed-full)/norm(full),
                           'Schur_evaluation':'SKIPPED: zI-D singular; use continued projected formula'})
    times=np.linspace(0,2*np.pi,257)
    dynamics={'interpretation':'Chosen dimensionless U(t)=exp(-itW); no physical clock identification',
              'grid_t':times.tolist(),
              'mode_complement_norm_fractions':[(m['b']**2*np.sin(times)**2).tolist() for m in mode_records],
              'matrix_exponential_checks':[]}
    for t in [0,.17,np.pi/4,np.pi/2,np.pi,2*np.pi]:
        ew=expm(-1j*t*w);prop=u.conj().T@ew@u
        trig=np.cos(t)*i8-1j*np.sin(t)*a
        expected=np.sin(t)**2*b@b.conj().T
        loss=i8-prop.conj().T@prop
        naive=expm(-1j*t*a)
        fractions=[]
        for k,m in enumerate(mode_records):
            psi=u@left[:,k]
            fractions.append(norm(f.conj().T@ew@psi)**2)
        dynamics['matrix_exponential_checks'].append({'t':float(t),
            'full_unitarity':norm(ew.conj().T@ew-i13),
            'projected_trig_identity':norm(prop-trig),
            'loss_identity':norm(loss-expected),
            'naive_active_propagator_difference':norm(naive-prop),
            'mode_complement_norm_fractions':fractions,
            'mode_fraction_errors':[abs(x-m['b']**2*np.sin(t)**2) for x,m in zip(fractions,mode_records)]})
    result={'status':'DIAGNOSTIC float64 geometry; exact involution/Schur identities checked separately',
            'mtft_version':mtft.__version__,'input_checks':input_checks,'block_checks':block_checks,
            'coupling_singular_values':s.tolist(),'normalized_singular_values':sn.tolist(),
            'coupling_rank_complex':retained,'rank_ambiguous':ambiguous,'modes':mode_records,
            'uncoupled_active_eigenvalues':ua_eig.tolist(),'uncoupled_fixed_eigenvalues':vf_eig.tolist(),
            'regular_response':regular,'pole_stress':stress,'real_compression_poles':real_poles,
            'dimensionless_unitary_control':dynamics,
            'summary':{'regular_point_count':len(regular),'pole_stress_point_count':len(stress),
                'max_regular_effective_relative_error':max(x['effective_relative_error'] for x in regular),
                'max_regular_involution_relative_error':max(x['involution_relative_error'] for x in regular),
                'max_regular_naive_relative_error':max(x['naive_relative_error'] for x in regular),
                'max_stress_effective_relative_error':max(x['effective_relative_error'] for x in stress),
                'max_stress_involution_relative_error':max(x['involution_relative_error'] for x in stress),
                'max_stress_Schur_condition':max(x['Schur_denominator_condition'] for x in stress)},
            'protocol_sha256':hashlib.sha256((root/'PROTOCOL.md').read_bytes()).hexdigest(),
            'input_sha256':{path.name:hashlib.sha256(path.read_bytes()).hexdigest()},
            'source_sha256':{'src/mtft/periods/involutions.py':hashlib.sha256((source/'src/mtft/periods/involutions.py').read_bytes()).hexdigest()},
            'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'elapsed_seconds':time.perf_counter()-started}
    (root/'effective_operator_results.json').write_text(json.dumps(result,indent=2)+'\n')
    np.savez_compressed(root/'effective_operator_matrices.npz',W=w,W_raw=raw_w,P=p_reconstructed,
        U=u,F=f,A=a,B=b,D=d,coupling_active_vectors=left,coupling_fixed_vectors=right,
        canonical_transform=transform,W_canonical=canon)
    print(json.dumps({'modes':mode_records,'summary':result['summary'],'elapsed_seconds':result['elapsed_seconds']},indent=2))


if __name__=='__main__':main()
