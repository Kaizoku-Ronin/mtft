#!/usr/bin/env python3
"""Two-phase, structure-preserving T2 memory approximation study."""
import argparse
import hashlib
import itertools
import json
import os
from pathlib import Path
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import numpy as np

ROOT=Path(__file__).resolve().parent
BUDGETS=[.01,.05,.10]
FAMILIES=['svd','poles','snapshots']


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def herm(a):
    return (a+a.conj().T)/2


def fro(a):
    return float(np.linalg.norm(a))


def input_blocks():
    m=np.load(ROOT/'inputs'/'T2_blocks.npz')
    a,b,d=m['A'],m['B'],m['D']
    h=np.block([[a,b],[b.conj().T,d]])
    frame=np.column_stack([m['U'],m['F']])
    checks={'selfadjoint':fro(h-h.conj().T)/fro(h),
            'frame_orthonormal':fro(frame.conj().T@frame-np.eye(13)),
            'block_reconstruction':fro(frame@h@frame.conj().T-m['H'])/fro(m['H'])}
    if max(checks.values())>1e-9:
        raise RuntimeError('Frozen block geometry gate failed')
    return a,b,d,h,m['F'],checks


def reduced(a,b,d,q):
    br=b@q;dr=herm(q.conj().T@d@q)
    return np.block([[a,br],[br.conj().T,dr]])


def eigen_data(h):
    w,v=np.linalg.eigh(herm(h))
    return w,v[:8,:]


def spectral_response(eig,z):
    w,v=eig
    return np.einsum('ai,ni,bi->nab',v,1/(np.asarray(z)[:,None]-w),v.conj(),optimize=True)


def spectral_propagator(eig,t):
    w,v=eig
    return np.einsum('ai,ni,bi->nab',v,np.exp(-1j*np.asarray(t)[:,None]*w),v.conj(),optimize=True)


def errors(approx,truth):
    delta=approx-truth
    rel=np.linalg.norm(delta,axis=(1,2))/np.linalg.norm(truth,axis=(1,2))
    absolute=np.linalg.svd(delta,compute_uv=False)[:,0]
    return rel,absolute


def training_grid():
    return np.asarray([complex(x,eta) for eta in [.05,.2,.7] for x in np.linspace(-3,3,61)])


def heldout_grid():
    return np.asarray([complex(-2.995+.01*j,eta) for eta in [.05,.2,.7] for j in range(600)])


def times_grid():
    return np.unique(np.concatenate([(np.arange(256)+.5)*2*np.pi/256,[0,.25,1,np.pi,2*np.pi]]))


def train():
    a,b,d,h,f,checks=input_blocks()
    grid=training_grid();truth=spectral_response(eigen_data(h),grid)
    _,sv,vh=np.linalg.svd(b,full_matrices=True)
    svbasis=vh.conj().T
    de,dv=np.linalg.eigh(d)
    gram=np.zeros((5,5),complex)
    for z in grid:
        x=np.linalg.solve(z*np.eye(5)-d,b.conj().T)
        gram+=x@x.conj().T/len(grid)
    mw,mv=np.linalg.eigh(herm(gram));mw=mw[::-1];mv=mv[:,::-1]
    candidates=[];poles_by_rank={};basis_data={'training_z':grid,'snapshot_gram':gram}
    for r in range(6):
        entries=[]
        for subset in itertools.combinations(range(5),r):
            q=dv[:,subset] if subset else np.zeros((5,0),complex)
            response=spectral_response(eigen_data(reduced(a,b,d,q)),grid)
            rel,_=errors(response,truth)
            row={'rank':r,'pole_indices':list(subset),'training_maximum_relative_error':float(np.max(rel))}
            candidates.append(row);entries.append((row,q,rel))
        best=min(x[0]['training_maximum_relative_error'] for x in entries)
        ties=[x for x in entries if x[0]['training_maximum_relative_error']<=best+1e-12]
        selected=min(ties,key=lambda x:tuple(x[0]['pole_indices']))
        poles_by_rank[r]=selected
    models={}
    for family in FAMILIES:
        for r in range(6):
            name=f'{family}_r{r}'
            if family=='svd':
                q=svbasis[:,:r]
                gap=None if r in [0,5] else float((sv[r-1]-sv[r])/sv[0])
                indices=None
            elif family=='snapshots':
                q=mv[:,:r]
                gap=None if r in [0,5] else float((mw[r-1]-mw[r])/max(abs(mw[0]),1e-300))
                indices=None
            else:
                entry,q,_=poles_by_rank[r]
                indices=entry['pole_indices']
                retained=set(indices);discarded=set(range(5))-retained
                cross=[abs(de[i]-de[j])/max(1,float(max(abs(de)))) for i in retained for j in discarded]
                gap=float(min(cross)) if cross else None
            model_h=reduced(a,b,d,q)
            response=spectral_response(eigen_data(model_h),grid)
            rel,_=errors(response,truth)
            models[name]={'family':family,'rank':r,'pole_indices':indices,
                          'training_maximum_relative_error':float(np.max(rel)),
                          'cutoff_relative_gap':gap,'ambiguous_cutoff':bool(gap is not None and gap<=1e-9),
                          'basis_orthonormality':fro(q.conj().T@q-np.eye(r)),
                          'reduced_Hermitian_residual':fro(model_h-model_h.conj().T)}
            basis_data['Q_'+name]=q
            basis_data['training_relative_'+name]=rel
            basis_data['projector_'+name]=f@q@q.conj().T@f.conj().T
    choices={}
    for budget in BUDGETS:
        by_family={}
        for family in FAMILIES:
            eligible=[m for m in models.values() if m['family']==family
                      and m['training_maximum_relative_error']<=budget and not m['ambiguous_cutoff']]
            selected=min(eligible,key=lambda m:m['rank']) if eligible else None
            by_family[family]=None if selected is None else {
                'model_id':f"{family}_r{selected['rank']}",'rank':selected['rank'],
                'training_maximum_relative_error':selected['training_maximum_relative_error']}
        choices[str(budget)]=by_family
    np.savez_compressed(ROOT/'selected_models.npz',**basis_data)
    selection={'status':'FROZEN_BEFORE_HELDOUT_EVALUATION',
        'selection_uses':'183 training spectral points only; no held-out spectral or time errors',
        'raw_checks':checks,'budgets':BUDGETS,'training_count':len(grid),
        'B_singular_values':sv.tolist(),'D_eigenvalues':de.tolist(),
        'snapshot_gram_eigenvalues_descending':mw.tolist(),
        'pole_subset_training_scores':candidates,'models':models,'training_choices':choices,
        'protocol_sha256':sha(ROOT/'PROTOCOL.md'),
        'input_sha256':sha(ROOT/'inputs'/'T2_blocks.npz'),
        'selected_bases_sha256':sha(ROOT/'selected_models.npz'),
        'script_sha256':sha(Path(__file__))}
    (ROOT/'selection_frozen.json').write_text(json.dumps(selection,indent=2)+'\n')
    print(json.dumps({'phase':'training_selection_frozen','training_choices':choices,
                      'snapshot_eigenvalues':mw.tolist()},indent=2))


def test():
    frozen_path=ROOT/'selection_frozen.json'
    frozen=json.loads(frozen_path.read_text())
    if (sha(ROOT/'selected_models.npz')!=frozen['selected_bases_sha256'] or
        sha(ROOT/'PROTOCOL.md')!=frozen['protocol_sha256'] or
        sha(ROOT/'inputs'/'T2_blocks.npz')!=frozen['input_sha256']):
        raise RuntimeError('Selection or inputs changed before holdout')
    a,b,d,h,f,_=input_blocks();basis=np.load(ROOT/'selected_models.npz')
    grid=heldout_grid();times=times_grid()
    truth=spectral_response(eigen_data(h),grid)
    prop=spectral_propagator(eigen_data(h),times)
    full_transfer=1-np.linalg.norm(prop,axis=(1,2))**2/8
    arrays={'heldout_z':grid,'times':times,'full_average_transfer':full_transfer}
    summaries={}
    for name,info in frozen['models'].items():
        q=basis['Q_'+name];r=info['rank']
        hm=reduced(a,b,d,q);eig=eigen_data(hm)
        g=spectral_response(eig,grid);gp=spectral_propagator(eig,times)
        response_rel,response_abs=errors(g,truth)
        time_rel,time_abs=errors(gp,prop)
        transfer=1-np.linalg.norm(gp,axis=(1,2))**2/8
        transfer_error=abs(transfer-full_transfer)
        if r==5:
            qd=np.zeros((5,0),complex)
        elif r==0:
            qd=np.eye(5,dtype=complex)
        else:
            qd=np.linalg.qr(q,mode='complete')[0][:,r:]
        direct=b@qd;mix=q.conj().T@d@qd
        coupling=np.vstack([direct,mix])
        delta=float(np.linalg.norm(coupling,2)) if qd.shape[1] else 0.
        direct_norm=float(np.linalg.norm(direct,2)) if qd.shape[1] else 0.
        mix_norm=float(np.linalg.norm(mix,2)) if mix.size else 0.
        spectral_bound=delta/(grid.imag**2)
        time_bound=np.minimum(2,abs(times)*delta)
        # Floating arithmetic leaves ~1e-14 error even for the exact r=5 control.
        bound_checks={'spectral_maximum_violation':float(max(0,np.max(response_abs-spectral_bound))),
                      'time_maximum_violation':float(max(0,np.max(time_abs-time_bound)))}
        j=int(np.argmax(response_rel));k=int(np.argmax(time_rel))
        per_eta={str(eta):float(max(response_rel[np.isclose(grid.imag,eta)])) for eta in [.05,.2,.7]}
        prefixes={str(horizon):{'maximum_relative_amplitude_error':float(max(time_rel[times<=horizon+1e-14])),
                               'maximum_average_transfer_absolute_error':float(max(transfer_error[times<=horizon+1e-14]))}
                  for horizon in [.25,1,2*np.pi]}
        summary={**info,'spectral_maximum_relative_error':float(max(response_rel)),
            'spectral_worst_point':{'real':float(grid[j].real),'imag':float(grid[j].imag)},
            'time_maximum_relative_error':float(max(time_rel)),'time_worst_point':float(times[k]),
            'average_transfer_maximum_absolute_error':float(max(transfer_error)),
            'spectral_error_by_eta':per_eta,'time_prefixes':prefixes,
            'discarded_direct_coupling_norm':direct_norm,'retained_discarded_D_mixing_norm':mix_norm,
            'embedding_error_operator_norm':delta,'conditional_absolute_bound_checks':bound_checks,
            'sampled_reduced_transfer_range':[float(min(transfer)),float(max(transfer))],
            'heldout_budget_passes':{str(budget):bool(max(response_rel)<=budget and max(time_rel)<=budget)
                                     for budget in BUDGETS}}
        summaries[name]=summary
        for key,values in [('response_relative',response_rel),('response_absolute',response_abs),
                           ('time_relative',time_rel),('time_absolute',time_abs),
                           ('transfer_error',transfer_error),('average_transfer',transfer),
                           ('spectral_absolute_bound',spectral_bound),('time_absolute_bound',time_bound)]:
            arrays[key+'_'+name]=values
    verdicts={}
    for budget,choices in frozen['training_choices'].items():
        verdicts[budget]={}
        for family,choice in choices.items():
            verdicts[budget][family]=None if choice is None else {
                **choice,'heldout_spectral_error':summaries[choice['model_id']]['spectral_maximum_relative_error'],
                'heldout_time_error':summaries[choice['model_id']]['time_maximum_relative_error'],
                'validated':summaries[choice['model_id']]['heldout_budget_passes'][budget]}
    scan={}
    for budget in map(str,BUDGETS):
        scan[budget]={}
        for family in FAMILIES:
            eligible=[m['rank'] for m in summaries.values() if m['family']==family and m['heldout_budget_passes'][budget]]
            scan[budget][family]=min(eligible) if eligible else None
    np.savez_compressed(ROOT/'approximation_pointwise.npz',**arrays)
    result={'status':'DIAGNOSTIC held-out approximation assessment',
        'heldout_count':len(grid),'time_count':len(times),'models':summaries,
        'frozen_choice_verdicts':verdicts,
        'descriptive_smallest_passing_rank_scan':scan,
        'scan_warning':'Ranks chosen after inspecting held-out errors are descriptive scans, not independently validated replacements.',
        'selection_sha256':sha(frozen_path),'selection_basis_sha256':sha(ROOT/'selected_models.npz'),
        'protocol_sha256':sha(ROOT/'PROTOCOL.md'),'input_sha256':sha(ROOT/'inputs'/'T2_blocks.npz'),
        'script_sha256':sha(Path(__file__)),
        'scope':'Fixed numerical spectral/time grids for a dimensionless finite operator. No global relative bound or physical identification.'}
    (ROOT/'approximation_results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'frozen_choice_verdicts':verdicts,'descriptive_scan':scan,
                      'rank4':{f:summaries[f+'_r4'] for f in FAMILIES}},indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--phase',choices=['train','test','all'],default='all')
    args=parser.parse_args()
    if args.phase in ['train','all']:
        train()
    if args.phase in ['test','all']:
        test()
