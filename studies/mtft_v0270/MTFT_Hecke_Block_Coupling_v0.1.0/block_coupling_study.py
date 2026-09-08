#!/usr/bin/env python3
"""Frozen Hecke-block study; see PROTOCOL.md for observables and gates."""
import argparse
from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path
import sys
import time
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import numpy as np

BLOCKS=['ell','old','q4','q6']
NAMES=['T2','T3','T5','T7','W11','W13','W143','STAR','T17','T19']


def norm(a): return float(np.linalg.norm(a))


def realify(a): return np.block([[a.real,-a.imag],[a.imag,a.real]])


def basis_from_projector(p):
    ev,u=np.linalg.eigh((p+p.T)/2)
    return u[:,ev>.5]


def geometry(p,es):
    rows={}
    for name in BLOCKS:
        u=basis_from_projector(es[name])
        ep=u.T@p@u
        ev=np.linalg.eigvalsh((ep+ep.T)/2)
        distance=np.minimum(abs(ev),abs(1-ev))
        rows[name]={'real_dimension':u.shape[1],
                    'active_overlap_real':float(np.trace(p@es[name])),
                    'active_overlap_complex':float(np.trace(p@es[name]))/2,
                    'active_fraction_of_block':float(np.trace(p@es[name]))/u.shape[1],
                    'fixed_overlap_real':u.shape[1]-float(np.trace(p@es[name])),
                    'principal_cosines_squared':ev.tolist(),
                    'active_intersection_real_diagnostic':int(np.sum(abs(1-ev)<=1e-8)),
                    'fixed_intersection_real_diagnostic':int(np.sum(abs(ev)<=1e-8)),
                    'ambiguous_endpoints':bool(np.any((distance>1e-8)&(distance<1e-6))),
                    'J_pairing_max_difference':float(max(abs(ev[::2]-ev[1::2]))),
                    'projector_idempotence':norm(es[name]@es[name]-es[name]),
                    'projector_symmetry':norm(es[name]-es[name].T)}
    pd=sum(es[b]@p@es[b] for b in BLOCKS)
    pblocks=np.array([[norm(es[b]@p@es[c])**2 for c in BLOCKS] for b in BLOCKS])
    return {'blocks':rows,
            'active_overlap_sum_real':sum(v['active_overlap_real'] for v in rows.values()),
            'fixed_overlap_sum_real':sum(v['fixed_overlap_real'] for v in rows.values()),
            'P_block_squared_norms':pblocks.tolist(),
            'P_block_norm_sum':float(pblocks.sum()),
            'block_diagonal_truncation_idempotence_defect':norm(pd@pd-pd),
            'discarded_offdiagonal_P_norm_squared':norm(p-pd)**2,
            'truncation_trace_defect':float(np.trace(pd-p)),
            'orthogonal_projectors_sum_identity':norm(sum(es.values())-np.eye(26))}


def measurements(a,p,es):
    q=np.eye(26)-p
    u=basis_from_projector(p)
    f=basis_from_projector(q)
    s=np.linalg.svd(f.T@a@u,compute_uv=False)
    scale=norm(a)
    sn=s/scale
    c=p@a-a@p
    c2=norm(c)**2
    asplit={b:es[b]@a@es[b] for b in BLOCKS}
    matrix=np.array([[norm(es[b]@c@es[d])**2/(scale*scale) for d in BLOCKS] for b in BLOCKS])
    identity_error=max(norm(es[b]@c@es[d]-(es[b]@p@es[d]@asplit[d]-asplit[b]@es[b]@p@es[d]))/scale for b in BLOCKS for d in BLOCKS)
    ks=[q@asplit[b]@p for b in BLOCKS]
    gram=np.array([[float(np.sum(kb*kc))/(scale*scale) for kc in ks] for kb in ks])
    leak2=norm(q@a@p)**2/(scale*scale)
    diag=float(np.trace(matrix))
    mixing=leak2>=1e-12
    pairs=[]
    for i,b in enumerate(BLOCKS):
        for j,d in enumerate(BLOCKS[i:],i):
            value=matrix[i,j] if i==j else matrix[i,j]+matrix[j,i]
            pairs.append({'blocks':[b,d],'squared_norm_over_operator_squared':float(value),
                          'fraction_of_commutator':float(value/(c2/(scale*scale))) if mixing else None})
    return {'coupling_singular_values':s.tolist(),'singular_values_over_operator_norm':sn.tolist(),
            'coupling_rank_real':int(np.sum(sn>=1e-7)),
            'rank_ambiguous':bool(np.any((sn>1e-9)&(sn<1e-7))),
            'active_to_fixed_norm_fraction':float(np.sqrt(leak2)),
            'commutator_norm_squared_over_operator_squared':c2/(scale*scale),
            'commutator_block_squared_fractions':matrix.tolist(),
            'within_block_commutator_squared_fraction':diag,
            'between_block_commutator_squared_fraction':float(matrix.sum()-diag),
            'between_block_share_of_commutator':float((matrix.sum()-diag)/matrix.sum()) if mixing else None,
            'block_pair_attribution':pairs,
            'commutator_attribution_sum_error':abs(float(matrix.sum())-c2/(scale*scale)),
            'block_commutator_identity_error':identity_error,
            'selfadjoint_cross_norm_identity_error':abs(c2/(scale*scale)-2*leak2),
            'full_operator_offblock_squared_fraction':sum(norm(es[b]@a@es[d])**2/(scale*scale) for b in BLOCKS for d in BLOCKS if b!=d),
            'signed_coupling_Gram_over_operator_squared':gram.tolist(),
            'coupling_Gram_sum':float(gram.sum()),
            'coupling_Gram_sum_error':abs(float(gram.sum())-leak2),
            'Gram_diagonal_sum':float(np.trace(gram)),
            'Gram_offdiagonal_sum':float(gram.sum()-np.trace(gram)),
            'component_coupling_reconstruction_error':norm(sum(ks)-q@a@p)/scale,
            'normalized_attribution_defined':mixing}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',required=True,type=Path)
    parser.add_argument('--prior',required=True,type=Path)
    args=parser.parse_args()
    source=args.source.resolve();prior=args.prior.resolve()
    root=Path(__file__).resolve().parent
    sys.path.insert(0,str(source/'src'))
    import mtft
    from mtft import hecke,liealg
    from mtft.periods.involutions import al_matrix
    started=time.perf_counter()
    exact=json.loads((root/'exact_projectors.json').read_text())
    natives={b:np.array([[float(Fraction(x)) for x in row] for row in exact['projectors'][b]]) for b in BLOCKS}
    ops={f'T{p}':np.array(hecke.cuspidal_hecke(p),float) for p in (2,3,5,7,17,19)}
    ops.update({f'W{p}':np.array(al_matrix(p),float) for p in (11,13,143)})
    ops['STAR']=np.array(hecke.star_involution(),float)
    baseline=np.load(prior/'coupling_matrices.npz')
    result={'status':'DIAGNOSTIC float64 for all period-dependent quantities',
            'release':mtft.__version__,'block_order':BLOCKS,'operator_order':NAMES,
            'protocol_sha256':hashlib.sha256((root/'PROTOCOL.md').read_bytes()).hexdigest(),
            'inputs_sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [prior/'coupling_matrices.npz',root/'exact_projectors.json']},
            'runs':[]}
    base_measured=None
    matrices={}
    for seed in [143,2026,31415]:
        if seed==143:
            r,ri,p=baseline['R'],baseline['Ri'],baseline['P_active_real']
        else:
            r,ri,bs=liealg.x0143_fixed_channels(dps=50,seed=seed)
            cs=[bs[i]@bs[j].conj()-bs[j]@bs[i].conj() for i in range(3) for j in range(i+1,3)]
            _,s,vh=np.linalg.svd(np.vstack(cs),full_matrices=False)
            rank=int(np.sum(s/s[0]>=1e-7))
            assert rank==8
            w=vh.conj().T[:,:rank]
            p=realify(w@w.conj().T)
        es={b:r@natives[b]@ri for b in BLOCKS}
        adapted={name:r@ops[name]@ri for name in NAMES}
        geo=geometry(p,es)
        obs={name:measurements(adapted[name],p,es) for name in NAMES}
        run={'frame_seed':seed,'period_input_dps':50,'geometry':geo,'operators':obs}
        if base_measured is None:
            base_measured=(geo,obs)
            matrices.update(P=p,R=r,Ri=ri)
            matrices.update({f'E_{b}':es[b] for b in BLOCKS})
        else:
            g0,o0=base_measured
            run['baseline_difference']={
                'max_overlap_real':max(abs(geo['blocks'][b]['active_overlap_real']-g0['blocks'][b]['active_overlap_real']) for b in BLOCKS),
                'max_principal_cosine_squared':max(float(np.max(np.abs(np.array(geo['blocks'][b]['principal_cosines_squared'])-g0['blocks'][b]['principal_cosines_squared']))) for b in BLOCKS),
                'max_attribution_fraction':max(float(np.max(np.abs(np.array(obs[name]['commutator_block_squared_fractions'])-o0[name]['commutator_block_squared_fractions']))) for name in NAMES)}
        result['runs'].append(run)
        print(f'frame {seed} completed',flush=True)
    result['elapsed_seconds']=time.perf_counter()-started
    result['script_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    (root/'block_coupling_results.json').write_text(json.dumps(result,indent=2)+'\n')
    np.savez_compressed(root/'block_coupling_matrices.npz',**matrices)
    geo,obs=base_measured
    print(json.dumps({'overlaps_complex':{b:geo['blocks'][b]['active_overlap_complex'] for b in BLOCKS},
                      'operators':{name:{'rank_real':x['coupling_rank_real'],'between_block_share':x['between_block_share_of_commutator'],
                                         'Gram_diagonal_sum':x['Gram_diagonal_sum'],'Gram_offdiagonal_sum':x['Gram_offdiagonal_sum']}
                                   for name,x in obs.items()}},indent=2))


if __name__=='__main__':main()
