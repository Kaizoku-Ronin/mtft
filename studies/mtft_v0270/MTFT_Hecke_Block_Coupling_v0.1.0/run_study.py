#!/usr/bin/env python3
"""Replay the exact block controls and numerical coupling study."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys


def main():
    root=Path(__file__).resolve().parent
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source',type=Path,required=True)
    p.add_argument('--prior',type=Path,default=root/'inputs'/'previous')
    args=p.parse_args()
    source=args.source.resolve();prior=args.prior.resolve()
    env=os.environ.copy()
    env['PYTHONPATH']=str(source/'src')+os.pathsep+env.get('PYTHONPATH','')
    env['OPENBLAS_NUM_THREADS']='1'
    commands=[['exact_blocks.py','--source',str(source)],
              ['block_coupling_study.py','--source',str(source),'--prior',str(prior)],
              ['independent_block_geometry.py','--source',str(source),'--previous',str(prior)],
              ['render_block_map.py']]
    with (root/'run_log.txt').open('w') as log:
        for command in commands:
            r=subprocess.run([sys.executable,*command],cwd=root,env=env,capture_output=True,text=True)
            log.write('\n$ python '+' '.join(command)+'\n'+r.stdout+r.stderr)
            log.flush()
            print(command[0]+': '+('completed' if r.returncode==0 else 'FAILED'),flush=True)
            if r.returncode: raise SystemExit(r.returncode)
    exact=json.loads((root/'exact_blocks_results.json').read_text())
    main=json.loads((root/'block_coupling_results.json').read_text())
    independent=json.loads((root/'independent_block_geometry_results.json').read_text())
    runs=main['runs'];base=runs[0]
    operators=base['operators']
    checks={
        'release_is_0_26_2':main['release']=='0.26.2',
        'exact_projector_partition':exact['projectors_sum_to_identity_exact'] and exact['projectors_pairwise_algebraically_orthogonal_exact'] and all(b['idempotent_exact'] for b in exact['blocks'].values()),
        'exact_independent_CRT_and_basis_agreement':all(b['matches_package_nullspace_projector_exact'] for b in exact['blocks'].values()),
        'all_ten_operators_preserve_every_block_exactly':len(exact['operators'])==10 and all(a['all_cross_hecke_blocks_zero_exact'] for a in exact['operators'].values()),
        'exact_W143_positive_projector_rank_four':exact['W143_rank_explanation']['plus_projector_rank_exact']==4,
        'independent_route_checks_passed':independent['all_checks_passed'],
        'independent_geometry_agreement':independent['primary_comparison']['maximum_absolute_difference_all_runs']<1e-9,
        'all_coupling_ranks_decided':all(not a['rank_ambiguous'] for run in runs for a in run['operators'].values()),
        'W143_coupling_rank_saturates_bound':all(run['operators']['W143']['coupling_rank_real']==4 for run in runs),
        'STAR_coupling_rank_zero':all(run['operators']['STAR']['coupling_rank_real']==0 for run in runs),
        'W143_q4_q6_commutator_component_zero_to_tolerance':operators['W143']['commutator_block_squared_fractions'][2][3]<1e-20 and operators['W143']['commutator_block_squared_fractions'][3][2]<1e-20,
        'W11_diagonal_commutator_components_zero_to_tolerance':operators['W11']['within_block_commutator_squared_fraction']<1e-20,
        'commutator_attribution_reconstructs':all(a['commutator_attribution_sum_error']<1e-9 and a['block_commutator_identity_error']<1e-9 for run in runs for a in run['operators'].values()),
        'signed_Gram_reconstructs_in_absolute_error':all(a['coupling_Gram_sum_error']<1e-9 for run in runs for a in run['operators'].values()),
        'active_fixed_overlap_sums':all(abs(run['geometry']['active_overlap_sum_real']-16)<1e-9 and abs(run['geometry']['fixed_overlap_sum_real']-10)<1e-9 for run in runs),
    }
    validation={'checks':checks,'all_checks_passed':all(checks.values()),
                'independent_check_count':len(independent['checks']),
                'scope':'Exact rational packaged-matrix controls plus numerical study consistency checks. Not a full package test suite or an interval certificate.',
                'known_prior_closure_ambiguities':'Not rerun or repaired; no current endpoint depends on numerical Lie closure.'}
    (root/'VALIDATION.json').write_text(json.dumps(validation,indent=2)+'\n')
    print(json.dumps(validation,indent=2))
    if not validation['all_checks_passed']:raise SystemExit(1)


if __name__=='__main__':main()
