#!/usr/bin/env python3
"""Frozen active-module experiment, MTFT 0.26.2; see PROTOCOL.md."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import sys
import time

os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
import numpy as np

SMALL, LARGE = 1e-9, 1e-7
PRIMARY = ['T2', 'T3', 'T5', 'T7', 'W11', 'W13', 'W143', 'STAR']
REPLICATION = ['T17', 'T19']


def norm(a):
    return float(np.linalg.norm(a))


def realify(a):
    return np.block([[a.real, -a.imag], [a.imag, a.real]])


def rank_decision(s):
    t = s / s[0] if len(s) and s[0] else s
    return {'rank': int(np.sum(t >= LARGE)),
            'ambiguous': bool(np.any((t > SMALL) & (t < LARGE))),
            'singular_values': s.tolist(), 'normalized_singular_values': t.tolist()}


def direct_projector(seeds):
    _, s, vh = np.linalg.svd(np.vstack(seeds), full_matrices=False)
    decision = rank_decision(s)
    if decision['ambiguous']:
        raise RuntimeError('Ambiguous primary common-kernel rank')
    w = vh.conj().T[:, :decision['rank']]
    return w @ w.conj().T, w, decision


def measure(a, p):
    q = np.eye(len(p)) - p
    den = norm(a)
    blocks = {'active_to_active': p @ a @ p, 'fixed_to_active': p @ a @ q,
              'active_to_fixed': q @ a @ p, 'fixed_to_fixed': q @ a @ q}
    sizes = {k: norm(v) for k, v in blocks.items()}
    fractions = {k: v*v/(den*den) for k, v in sizes.items()}
    leak = sizes['active_to_fixed']/den
    classification = ('consistent_with_preservation' if leak <= 1e-8 else
                      'resolved_mixing' if leak >= 1e-6 else 'ambiguous')
    return {'operator_norm': den, 'block_norms': sizes,
            'block_squared_fractions': fractions,
            'fraction_sum_error': abs(sum(fractions.values())-1),
            'active_to_fixed_over_operator': leak,
            'fixed_to_active_over_operator': sizes['fixed_to_active']/den,
            'active_to_fixed_over_active_input': sizes['active_to_fixed']/norm(a@p)
                if norm(a@p) > 1e-15*den else None,
            'fixed_to_active_over_fixed_input': sizes['fixed_to_active']/norm(a@q)
                if norm(a@q) > 1e-15*den else None,
            'commutator_over_operator': norm(a@p-p@a)/den,
            'cross_block_adjoint_residual': norm(blocks['active_to_fixed'].T-
                                                blocks['fixed_to_active'])/den,
            'classification': classification}


def invariant_hull(w, operators, j):
    q = w.copy()
    history = [q.shape[1]]
    spectra = []
    normalized = {k: a/norm(a) for k,a in operators.items() if norm(a)}
    for _ in range(27):
        augmented = np.hstack([q] + [a@q for a in normalized.values()])
        u,s,_ = np.linalg.svd(augmented, full_matrices=False)
        decision = rank_decision(s)
        spectra.append(decision)
        if decision['ambiguous']:
            return {'status':'ambiguous','growth':history,'step_spectra':spectra},q
        qnew = u[:, :decision['rank']]
        history.append(qnew.shape[1])
        if qnew.shape[1] == q.shape[1]:
            q = qnew
            break
        q = qnew
    p = q@q.T
    return {'status':'resolved_numerically', 'growth':history,
            'real_dimension':q.shape[1], 'step_spectra':spectra,
            'terminal_leakage':{k:norm((np.eye(26)-p)@a@p)/norm(a)
                                for k,a in operators.items()},
            'J_leakage':norm((np.eye(26)-p)@j@p)/norm(j)}, q


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    args = parser.parse_args()
    source = args.source.resolve()
    sys.path.insert(0,str(source/'src'))
    import mtft
    from mtft import liealg, hecke
    from mtft.periods import physics
    root = Path(__file__).resolve().parent
    started = time.perf_counter()
    native_ops = {f'T{p}':np.array(hecke.cuspidal_hecke(p),float)
                  for p in [2,3,5,7,17,19]}
    native_ops.update(liealg.x0143_symmetry_ops())
    names = PRIMARY+REPLICATION
    n=13
    j = np.block([[np.zeros((n,n)),-np.eye(n)], [np.eye(n),np.zeros((n,n))]])
    result = {'status':'DIAGNOSTIC float64, not exact or interval-certified',
              'package_version':mtft.__version__, 'numpy_version':np.__version__,
              'source':str(source), 'primary':PRIMARY,'replication':REPLICATION,
              'protocol_sha256':hashlib.sha256((root/'PROTOCOL.md').read_bytes()).hexdigest(),
              'rank_bands':{'small_max':SMALL,'large_min':LARGE}, 'runs':[]}
    base = None
    raw = {}
    for dps, seed in [(50,143),(50,2026),(50,31415),(80,143)]:
        r,ri,fixed = liealg.x0143_fixed_channels(dps=dps,seed=seed)
        seeds = [fixed[i]@fixed[k].conj()-fixed[k]@fixed[i].conj()
                 for i in range(3) for k in range(i+1,3)]
        pc,wc,rank = direct_projector(seeds)
        p,w = realify(pc),realify(wc)
        ops = {name:r@native_ops[name]@ri for name in names}
        try:
            closure = liealg.close_lie(seeds)
            rep = liealg.rep_summary(closure['basis'])
            closure_info = {'closure_status':'passed_package_gate',
                            'direct_vs_closure_projector':norm(pc-rep['P_active']),
                            'closure_dimension':len(closure['basis']),
                            'closure_growth':closure['growth']}
        except liealg.LieGateAmbiguous as exc:
            closure_info = {'closure_status':'ambiguous_package_gate',
                            'closure_error':str(exc),
                            'direct_vs_closure_projector':None,
                            'closure_dimension':None,'closure_growth':None}
        g = physics.hodge_metric_hecke(dps)
        checks = {'R_Ri':norm(r@ri-np.eye(26)),
                  'metric_whitening':norm(ri.T@g@ri-np.eye(26))/np.sqrt(26),
                  'P_idempotence':norm(p@p-p), 'P_symmetry':norm(p.T-p),
                  'P_J_commutator':norm(p@j-j@p), **closure_info,
                  'seed_fixed_annihilation':max(norm(a@(np.eye(13)-pc))/norm(a) for a in seeds)}
        measures = {name:measure(ops[name],p) for name in names}
        star = ops['STAR']
        qbasis = np.linalg.svd(np.eye(26)-p)[0][:,:26-w.shape[1]]
        starchecks = {'square_identity':norm(star@star-np.eye(26))/np.sqrt(26),
                      'J_anticommutator':norm(star@j+j@star)/norm(star),
                      'isometry':norm(star.T@star-np.eye(26))/np.sqrt(26),
                      'active_restriction_eigenvalues':np.linalg.eigvalsh(w.T@star@w).tolist(),
                      'fixed_restriction_eigenvalues':np.linalg.eigvalsh(qbasis.T@star@qbasis).tolist()}
        run = {'dps_period_inputs':dps,'frame_seed':seed,
               'complex_active_rank':rank['rank'],'real_active_rank':w.shape[1],
               'direct_seed_spectrum':rank,'checks':checks,'operators':measures,
               'STAR':starchecks}
        if base is None:
            base=(r,ri,p,ops,measures)
            groups={'T2_only':['T2'], 'good_Hecke':['T2','T3','T5','T7'],
                    'Atkin_Lehner':['W11','W13','W143'], 'all_primary':PRIMARY}
            result['hulls']={}
            for group,selected in groups.items():
                hr,hq=invariant_hull(w,{name:ops[name] for name in selected},j)
                hp=hq@hq.T
                hr['replication_leakage']={name:norm((np.eye(26)-hp)@ops[name]@hp)/norm(ops[name])
                                          for name in REPLICATION}
                result['hulls'][group]=hr
                raw['hull_'+group]=hq
            result['controls']={'identity':measure(np.eye(26),p),
                                'full_space_T2':measure(ops['T2'],np.eye(26))}
            for name in names: raw['operator_'+name]=ops[name]
            raw.update(P_active_real=p, P_active_complex=pc, active_basis_real=w,
                       R=r,Ri=ri,J=j,three_fixed_channels=np.stack(fixed),
                       three_Lie_seeds=np.stack(seeds))
        else:
            rb,rbi,pb,ob,mb=base
            trans=rb@ri
            invtrans=r@rbi
            run['baseline_comparison']={
                'transported_P_difference':norm(trans@p@invtrans-pb),
                'transition_isometry':norm(trans.T@trans-np.eye(26)),
                'max_operator_transport_difference':max(norm(trans@ops[name]@invtrans-ob[name])/norm(ob[name]) for name in names),
                'max_block_fraction_difference':max(abs(measures[name]['block_squared_fractions'][k]-mb[name]['block_squared_fractions'][k]) for name in names for k in mb[name]['block_squared_fractions'])}
        result['runs'].append(run)
        print(f'dps={dps} seed={seed}: active={rank["rank"]}, closure={checks["closure_status"]}, direct/closure={checks["direct_vs_closure_projector"]}',flush=True)
        (root/'coupling_results.partial.json').write_text(json.dumps(result,indent=2)+'\n')
    sourcefiles=['src/mtft/liealg.py','src/mtft/hecke.py','src/mtft/periods/physics.py',
                 'src/mtft/periods/hamiltonian.py','src/mtft/periods/involutions.py']
    result['source_sha256']={f:hashlib.sha256((source/f).read_bytes()).hexdigest() for f in sourcefiles}
    result['script_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    result['elapsed_seconds']=time.perf_counter()-started
    (root/'coupling_results.json').write_text(json.dumps(result,indent=2)+'\n')
    (root/'coupling_results.partial.json').unlink(missing_ok=True)
    np.savez_compressed(root/'coupling_matrices.npz',**raw)
    print(json.dumps({'coupling':{k:(v['classification'],v['active_to_fixed_over_operator']) for k,v in result['runs'][0]['operators'].items()},
                      'hulls':{k:v['growth'] for k,v in result['hulls'].items()},
                      'elapsed_seconds':result['elapsed_seconds']},indent=2))


if __name__=='__main__':
    main()
