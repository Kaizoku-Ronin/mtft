#!/usr/bin/env python3
"""T2 feedback and memory on the frozen MTFT active/complement split."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import sys
import time
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
import numpy as np
from scipy.linalg import expm


def norm(a):
    return float(np.linalg.norm(a))


def herm(a):
    return (a+a.conj().T)/2


def linear(a):
    n=len(a)//2
    return (a[:n,:n]+a[n:,n:]+1j*(a[n:,:n]-a[:n,n:]))/2


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rank_record(s, scale):
    s=np.asarray(s)/max(float(scale),1e-300)
    return {'normalized_singular_values':s.tolist(),
            'resolved_rank':int(np.sum(s>=1e-7)),
            'ambiguous':bool(np.any((s>1e-9)&(s<1e-7)))}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',type=Path,required=True)
    parser.add_argument('--previous',type=Path,required=True)
    args=parser.parse_args()
    root=Path(__file__).resolve().parent
    source=args.source.resolve();previous=args.previous.resolve()
    sys.path.insert(0,str(source/'src'))
    import mtft
    from mtft.hecke import cuspidal_hecke
    start=time.perf_counter()
    old=np.load(previous/'coupling_matrices.npz')
    r,ri,p=old['R'],old['Ri'],old['P_active_complex']
    exact=cuspidal_hecke(2)
    real=r@np.asarray(exact,float)@ri
    raw=linear(real)
    j=old['J']
    n13=np.eye(13);n8=np.eye(8);n5=np.eye(5)
    raw_checks={
        'frame_inverse':norm(r@ri-np.eye(26))/np.sqrt(26),
        'real_selfadjoint':norm(real-real.T)/norm(real),
        'real_complex_linearity':norm(real@j-j@real)/norm(real),
        'complex_selfadjoint':norm(raw-raw.conj().T)/norm(raw),
        'P_idempotence':norm(p@p-p)/norm(p),
        'P_selfadjoint':norm(p-p.conj().T)/norm(p),
        'fresh_cached_operator_agreement':norm(real-old['operator_T2'])/norm(real),
    }
    if max(raw_checks.values())>1e-9:
        (root/'RAW_GATE_FAILURE.json').write_text(json.dumps(raw_checks,indent=2)+'\n')
        raise RuntimeError('Raw geometry gate failed; no Hermitian experiment performed')
    h=herm(raw)
    ep,up=np.linalg.eigh(herm(p))
    u=up[:,ep>.5];f=up[:,ep<=.5]
    if (u.shape[1],f.shape[1])!=(8,5):
        raise RuntimeError('Frozen sector dimension mismatch')
    a=herm(u.conj().T@h@u)
    b=u.conj().T@h@f
    d=herm(f.conj().T@h@f)
    left,s,vh=np.linalg.svd(b,full_matrices=True)
    right=vh.conj().T
    coupling_rank=rank_record(s,norm(h))
    if coupling_rank['ambiguous']:
        raise RuntimeError('Ambiguous coupling rank; do not choose bright/dark split')
    k=coupling_rank['resolved_rank']
    dark=left[:,k:]
    delayed=b.conj().T@a@dark
    reach=np.column_stack([np.linalg.matrix_power(d,jj)@b.conj().T for jj in range(5)])
    columns=[];pairs=[]
    for jj in range(k):
        x=u@left[:,jj];y=f@right[:,jj]
        pair=np.column_stack([x,y]);block=pair.conj().T@h@pair
        columns.extend([x,y])
        pairs.append({'index':jj+1,'a':float(block[0,0].real),
                      'd':float(block[1,1].real),'b':float(s[jj]),
                      'invariance_relative_residual':norm(h@pair-pair@block)/norm(h)})
    columns.extend((u@dark).T)
    if k<5:
        columns.extend((f@right[:,k:]).T)
    transform=np.column_stack(columns)
    paired=transform.conj().T@h@transform
    target=np.zeros_like(paired)
    for jj in range(k):
        target[2*jj:2*jj+2,2*jj:2*jj+2]=paired[2*jj:2*jj+2,2*jj:2*jj+2]
    target[2*k:2*k+dark.shape[1],2*k:2*k+dark.shape[1]]=paired[2*k:2*k+dark.shape[1],2*k:2*k+dark.shape[1]]
    if k<5:
        tail=2*k+dark.shape[1];target[tail:,tail:]=paired[tail:,tail:]
    bbt=b@b.conj().T;btb=b.conj().T@b
    lam,ev=np.linalg.eigh(h)
    de,dv=np.linalg.eigh(d)
    grouping_tolerance=1e-9*max(1,float(np.max(abs(lam))))
    groups=[]
    for jj,x in enumerate(de):
        if groups and abs(x-de[groups[-1][0]])<=grouping_tolerance:
            groups[-1].append(jj)
        else:
            groups.append([jj])
    residues=[];poles=[]
    for indices in groups:
        vectors=b@dv[:,indices]
        residue=vectors@vectors.conj().T
        residues.append(residue)
        dj=float(np.mean(de[indices]))
        poles.append({'eigenvalue':dj,'multiplicity':len(indices),
                      'residue_trace':float(np.trace(residue).real),
                      'residue_eigenvalues':np.linalg.eigvalsh(herm(residue)).tolist(),
                      'distance_from_full_spectrum':float(np.min(abs(lam-dj)))})
    coefficient=ev.conj().T@u
    regular_z=[complex(x,eta) for eta in [.05,.2,.7] for x in np.linspace(-3,3,601)]
    regular_z.extend([-4+0j,4+0j,2j,.17+.37j])

    def spectral_response(z):
        return coefficient.conj().T@((1/(z-lam))[:,None]*coefficient)

    def response(z):
        full=u.conj().T@np.linalg.solve(z*n13-h,u)
        sigma=b@np.linalg.solve(z*n5-d,b.conj().T)
        schur_matrix=z*n8-a-sigma
        effective=np.linalg.solve(schur_matrix,n8)
        spectral=spectral_response(z)
        naive=np.linalg.solve(z*n8-a,n8)
        row={'z_real':float(z.real),'z_imag':float(z.imag),
             'full_norm':norm(full),'sigma_norm':norm(sigma),
             'full_condition':float(np.linalg.cond(z*n13-h)),
             'complement_condition':float(np.linalg.cond(z*n5-d)),
             'Schur_condition':float(np.linalg.cond(schur_matrix)),
             'effective_relative_error':norm(effective-full)/norm(full),
             'spectral_relative_error':norm(spectral-full)/norm(full),
             'naive_relative_error':norm(naive-full)/norm(full),
             'average_spectral_response_full':float(-np.trace(full).imag/(8*np.pi)),
             'average_spectral_response_naive':float(-np.trace(naive).imag/(8*np.pi))}
        return row,full

    rows=[];regular_full=[];admissible_z=[];skipped=[]
    for z in regular_z:
        # Real endpoints are fixed; do not silently invert a spectral point.
        if z.imag==0 and min(np.min(abs(lam-z.real)),np.min(abs(de-z.real)),
                            np.min(abs(np.linalg.eigvalsh(a)-z.real)))<1e-12:
            skipped.append({'z_real':z.real,'reason':'real endpoint is spectral'})
            continue
        row,full=response(z);rows.append(row);regular_full.append(full);admissible_z.append(z)
    stress=[];real_poles=[];coincident=[]
    for jj,pole in enumerate(poles):
        dj=pole['eigenvalue']
        if pole['distance_from_full_spectrum']<=1e-7*max(1,float(max(abs(lam)))):
            coincident.append(pole);continue
        for power in range(2,13):
            row,_=response(complex(dj,10.**(-power)))
            row['pole_index']=jj;row['compression_pole']=dj
            stress.append(row)
        g=u.conj().T@np.linalg.solve(dj*n13-h,u)
        v=b@dv[:,groups[jj]]
        real_poles.append({'pole_index':jj,'z':dj,
            'projected_norm':norm(g),'full_condition':float(np.linalg.cond(dj*n13-h)),
            'null_relation_relative_residual':norm(g@v)/(norm(g)*norm(v)),
            'projected_singular_values':np.linalg.svd(g,compute_uv=False).tolist(),
            'spectral_relative_error':norm(spectral_response(dj)-g)/norm(g),
            'Schur_evaluation':'SKIPPED: singular complementary inverse'})

    def spectral_evolution(t):
        return (ev*np.exp(-1j*t*lam))@ev.conj().T

    def memory_spectral(t):
        return sum(np.exp(-1j*t*pole['eigenvalue'])*residue
                   for pole,residue in zip(poles,residues))

    times=np.linspace(0,2*np.pi,257)
    dynamics=[]
    for t in times:
        evolution=spectral_evolution(t)
        c=f.conj().T@evolution@u
        cd=c@dark
        pp=u.conj().T@evolution@u
        dynamics.append({'t':float(t),'average_transfer':norm(c)**2/8,
                         'worst_input_transfer':float(np.linalg.norm(c,2)**2),
                         'dark_average_transfer':norm(cd)**2/dark.shape[1] if dark.shape[1] else None,
                         'dark_worst_input_transfer':float(np.linalg.norm(cd,2)**2) if dark.shape[1] else None,
                         'active_unitarity_defect':norm(n8-pp.conj().T@pp),
                         'memory_norm':norm(memory_spectral(t))})
    exp_controls=[]
    for t in [0,.17,np.pi/4,np.pi/2,np.pi,2*np.pi]:
        full=expm(-1j*t*h);spectral=spectral_evolution(t)
        pp=u.conj().T@full@u;c=f.conj().T@full@u
        kernel=b@expm(-1j*t*d)@b.conj().T
        exp_controls.append({'t':float(t),
            'spectral_full_relative_error':norm(full-spectral)/norm(full),
            'full_unitarity':norm(full.conj().T@full-n13),
            'projected_spectral_relative_error':norm(pp-u.conj().T@spectral@u)/norm(pp),
            'norm_balance_residual':norm(n8-pp.conj().T@pp-c.conj().T@c),
            'memory_spectral_relative_error':norm(kernel-memory_spectral(t))/max(norm(kernel),1),
            'isolated_active_relative_error':norm(expm(-1j*t*a)-pp)/norm(pp)})
    short=[]
    ordinary_coefficient=norm(b)**2/8
    dark_coefficient=norm(delayed)**2/(4*dark.shape[1]) if dark.shape[1] else None
    for jj in range(6):
        t=.1/(2**jj);evolution=expm(-1j*t*h)
        c=f.conj().T@evolution@u;cd=c@dark
        short.append({'t':t,
            'ordinary_amplitude_scaled_error':norm(c/t+1j*b.conj().T)/norm(b),
            'ordinary_average_over_t_squared':norm(c)**2/(8*t*t),
            'dark_amplitude_scaled_error':norm(cd/(t*t)+delayed/2)/max(norm(delayed/2),1e-300),
            'dark_average_over_t_fourth':norm(cd)**2/(dark.shape[1]*t**4) if dark.shape[1] else None})

    summary={
        'regular_count':len(rows),'skipped_regular_count':len(skipped),
        'stress_count':len(stress),'compression_poles_tested':len(real_poles),
        'max_regular_effective_relative_error':max(x['effective_relative_error'] for x in rows),
        'max_regular_spectral_relative_error':max(x['spectral_relative_error'] for x in rows),
        'max_regular_naive_relative_error':max(x['naive_relative_error'] for x in rows),
        'max_stress_effective_relative_error':max((x['effective_relative_error'] for x in stress),default=None),
        'max_stress_spectral_relative_error':max((x['spectral_relative_error'] for x in stress),default=None),
        'max_stress_Schur_condition':max((x['Schur_condition'] for x in stress),default=None),
        'max_real_pole_null_residual':max((x['null_relation_relative_residual'] for x in real_poles),default=None),
    }
    for key in ['average_transfer','worst_input_transfer','dark_average_transfer','dark_worst_input_transfer']:
        if dynamics[0][key] is not None:
            peak=max(dynamics,key=lambda x:x[key])
            summary['sampled_max_'+key]={'value':peak[key],'t':peak['t']}
    structure={
        'Hermitian_cleanup_size':norm(h-raw),
        'coupling_singular_values':s.tolist(),'coupling_rank':coupling_rank,
        'full_eigenvalues':lam.tolist(),'complement_eigenvalues':de.tolist(),
        'active_eigenvalues':np.linalg.eigvalsh(a).tolist(),
        'involution_defect_relative':norm(h@h-n13)/np.sqrt(13),
        'SVD_pair_records':pairs,
        'independent_pair_decomposition_relative_residual':norm(paired-target)/norm(h),
        'A_BBstar_commutator_relative':norm(a@bbt-bbt@a)/(norm(a)*norm(bbt)),
        'D_BstarB_commutator_relative':norm(d@btb-btb@d)/(norm(d)*norm(btb)),
        'dark_complex_dimension':dark.shape[1],
        'instantaneous_dark_coupling_relative':norm(b.conj().T@dark)/norm(h),
        'delayed_dark_coupling_norm':norm(delayed),
        'delayed_dark_singular_values':np.linalg.svd(delayed,compute_uv=False).tolist(),
        'delayed_dark_rank':rank_record(np.linalg.svd(delayed,compute_uv=False),max(norm(h)**2,1)),
        'complement_reachable_rank':rank_record(np.linalg.svd(reach,compute_uv=False),norm(reach)),
        'residue_sum_relative_residual':norm(sum(residues)-bbt)/norm(bbt),
        'ordinary_short_time_coefficient':ordinary_coefficient,
        'dark_short_time_coefficient':dark_coefficient}
    result={'status':'DIAGNOSTIC finite-dimensional T2 experiment',
        'mtft_version':mtft.__version__,'raw_checks':raw_checks,'structure':structure,
        'compression_poles':poles,'coincident_full_poles':coincident,
        'regular_response':rows,'skipped_regular':skipped,'pole_stress':stress,
        'real_compression_poles':real_poles,'dimensionless_dynamics':dynamics,
        'exponential_controls':exp_controls,'short_time_controls':short,'summary':summary,
        'interpretation':'T2 is a chosen dimensionless generator; no physical clock, Hamiltonian, mass, QFT, or density-matrix partial trace identified.',
        'protocol_sha256':sha(root/'PROTOCOL.md'),
        'input_sha256':{'coupling_matrices.npz':sha(previous/'coupling_matrices.npz')},
        'source_sha256':{'src/mtft/hecke.py':sha(source/'src/mtft/hecke.py')},
        'script_sha256':sha(Path(__file__)),
        'elapsed_seconds':time.perf_counter()-start}
    (root/'t2_effective_operator_results.json').write_text(json.dumps(result,indent=2)+'\n')
    np.savez_compressed(root/'t2_effective_operator_matrices.npz',H=h,H_raw=raw,H_real=real,
        P=u@u.conj().T,U=u,F=f,A=a,B=b,D=d,left_singular_vectors=left,
        right_singular_vectors=right,dark_basis=dark,delayed_dark=delayed,
        eigenvalues=lam,eigenvectors=ev,complement_eigenvectors=dv,
        residue_matrices=np.asarray(residues),regular_z=np.asarray(admissible_z),
        regular_projected_responses=np.asarray(regular_full))
    print(json.dumps({'structure':structure,'summary':summary,'elapsed_seconds':result['elapsed_seconds']},indent=2))


if __name__=='__main__':
    main()
