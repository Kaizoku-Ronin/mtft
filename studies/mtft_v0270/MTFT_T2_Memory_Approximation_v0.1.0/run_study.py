#!/usr/bin/env python3
"""Replay fixed T2 memory approximations and compare independent routes."""
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess
import sys

ROOT=Path(__file__).resolve().parent
PROTOCOL_SHA='f784e5be3f6620f3bd9d21f6e5ba7501990f1f32035ab66340ec6fce226b5f8e'
INPUT_SHA={
    'T2_blocks.npz':'af01726b204310d959241020af12594b3687814eeb2cf29f85072fbabfd97f8b',
    'T2_independent_geometry.npz':'2990bd3f4ec917f8ae1a264da3688bc235388a716124d9e371ce1026a232237f',
    'mtft-0.26.2.tar.gz':'35cfc00c32877bd6c12b16464ddc8190bcca9e2f054474f4ab0a6c605fe9da38',
}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(name):
    return json.loads((ROOT/name).read_text())


def main():
    os.environ['OPENBLAS_NUM_THREADS']='1'
    import numpy as np
    import scipy
    import sympy
    import mpmath
    import matplotlib
    if sha(ROOT/'PROTOCOL.md')!=PROTOCOL_SHA:
        raise RuntimeError('Frozen protocol hash mismatch')
    for name,expected in INPUT_SHA.items():
        if sha(ROOT/'inputs'/name)!=expected:
            raise RuntimeError(f'Frozen input hash mismatch: {name}')
    commands=[['memory_approximation.py','--phase','train'],
              ['memory_approximation.py','--phase','test'],
              ['independent_approximation.py'],['truncation_identities.py'],
              ['render_memory_approximation.py']]
    freeze_hash=None
    with (ROOT/'run_log.txt').open('w') as log:
        for command in commands:
            result=subprocess.run([sys.executable,*command],cwd=ROOT,env=dict(os.environ),text=True,
                                  stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
            log.write(json.dumps(command)+'\n'+result.stdout+'\n');log.flush()
            if result.returncode:
                raise RuntimeError(f'{command[0]} failed; see run_log.txt')
            if command==commands[0]:
                freeze_hash=sha(ROOT/'selection_frozen.json')
            elif command==commands[1] and sha(ROOT/'selection_frozen.json')!=freeze_hash:
                raise RuntimeError('Training selection changed during holdout')
    selection=read('selection_frozen.json');primary=read('approximation_results.json')
    independent=read('independent_approximation_results.json');isel=read('independent_selection.json')
    symbolic=read('truncation_identity_results.json')
    pp=np.load(ROOT/'approximation_pointwise.npz')
    ip=np.load(ROOT/'independent_approximation_matrices.npz')
    bases=np.load(ROOT/'selected_models.npz')
    names=list(primary['models']);tol=1e-9
    differences={}
    comparisons={}
    for metric in ['response_relative','response_absolute','time_relative','time_absolute','transfer_error']:
        values={name:float(np.max(abs(pp[metric+'_'+name]-ip[metric+'_'+name]))) for name in names}
        differences[metric+'_maximum_difference']=max(values.values())
        comparisons[metric]=values
    differences['ambient_projector_maximum_Frobenius_difference']=max(
        float(np.linalg.norm(bases['projector_'+name]-ip['projector_'+name])) for name in names)
    differences['training_error_maximum_difference']=max(
        float(np.max(abs(bases['training_relative_'+name]-ip['training_relative_'+name]))) for name in names)
    grids_match=bool(np.array_equal(pp['heldout_z'],ip['heldout_z']) and
                     np.array_equal(pp['times'],ip['times']) and
                     np.array_equal(bases['training_z'],ip['training_z']))
    min_separation=float(np.min(abs(bases['training_z'][:,None]-pp['heldout_z'][None,:])))
    choices_match=all(choice['rank']==independent['training_selected_rank_by_budget'][family][budget]
        for budget,byfamily in selection['training_choices'].items() for family,choice in byfamily.items())
    pole_choices_match=all(selection['models'][f'poles_r{r}']['pole_indices']==isel['selected_pole_indices'][str(r)] for r in range(6))
    full_controls={name:{'response':m['spectral_maximum_relative_error'],'time':m['time_maximum_relative_error']}
                   for name,m in primary['models'].items() if m['rank']==5}
    bound_excess=max(max(m['conditional_absolute_bound_checks'].values()) for m in primary['models'].values())
    pass_calculations_consistent=all(m['heldout_budget_passes'][str(budget)]==
        (m['spectral_maximum_relative_error']<=budget and m['time_maximum_relative_error']<=budget)
        for m in primary['models'].values() for budget in selection['budgets'])
    all_values=[pp[key] for key in pp.files]
    checks={
        'parent_study_validation_passed':read('inputs/T2_PARENT_VALIDATION.json')['status']=='PASS',
        'primary_raw_geometry':max(selection['raw_checks'].values())<=tol,
        'training_selection_precedes_and_survives_holdout':freeze_hash==primary['selection_sha256']==sha(ROOT/'selection_frozen.json'),
        'independent_selection_survives_holdout':independent['selection_unchanged_after_holdout']
            and independent['selection_sha256']==sha(ROOT/'independent_selection.json'),
        'basis_freeze_matches':selection['selected_bases_sha256']==sha(ROOT/'selected_models.npz'),
        'protocol_hashes_match':all(x['protocol_sha256']==PROTOCOL_SHA for x in [selection,primary,independent,isel]),
        'prescribed_point_counts':selection['training_count']==183 and primary['heldout_count']==1800 and primary['time_count']==261,
        'training_and_holdout_grids_disjoint':min_separation>1e-12,
        'independent_grids_agree':grids_match,
        'all_32_pole_subsets_scored':len(selection['pole_subset_training_scores'])==32,
        'all_18_family_rank_models_assessed':len(names)==18,
        'independent_training_choices_agree':choices_match,
        'independent_pole_subset_choices_agree':pole_choices_match,
        'independent_errors_and_subspaces_agree':max(differences.values())<=tol,
        'no_ambiguous_selected_cutoffs':not any(m['ambiguous_cutoff'] for m in selection['models'].values()),
        'reduced_bases_orthonormal':max(m['basis_orthonormality'] for m in selection['models'].values())<=tol,
        'reduced_generators_Hermitian':max(m['reduced_Hermitian_residual'] for m in selection['models'].values())<=tol,
        'full_rank_controls_reproduce_full_model':all(max(v.values())<=tol for v in full_controls.values()),
        'independent_identity_and_bound_checks':independent['all_identity_and_bound_checks_passed'],
        'conditional_absolute_bounds_respected_to_roundoff':bound_excess<=tol,
        'budget_labels_match_actual_errors':pass_calculations_consistent,
        'forty_four_symbolic_controls':symbolic['checks_passed']==symbolic['checks_total']==44
            and all(symbolic['checks'].values()),
        'finite_pointwise_arrays':all(np.all(np.isfinite(v)) for v in all_values),
        'raw_transfer_within_unitarity_bounds_to_roundoff':all(
            m['sampled_reduced_transfer_range'][0]>=-tol and m['sampled_reduced_transfer_range'][1]<=1+tol
            for m in primary['models'].values()),
    }
    # NumPy scalar predicates are deliberately converted for JSON serialization.
    checks={key:bool(value) for key,value in checks.items()}
    outcome={str(budget):[name for name,m in primary['models'].items()
                         if m['rank']<5 and m['heldout_budget_passes'][str(budget)]] for budget in selection['budgets']}
    validation={
        'status':'PASS' if all(checks.values()) else 'FAIL',
        'scope':'Execution/identity/independence checks, not a successful-compression verdict.',
        'research_outcome':'No tested reduction below five coordinates passes any of the three joint budgets.' if not any(outcome.values()) else 'See per-budget reduction results.',
        'passing_reduced_models_by_budget':outcome,
        'checks':checks,'checks_passed':sum(checks.values()),'checks_total':len(checks),
        'cross_route_maximum_differences':differences,'cross_route_by_model':comparisons,
        'training_holdout_minimum_spectral_separation':min_separation,
        'full_rank_controls':full_controls,'maximum_absolute_bound_excess':bound_excess,
        'frozen_choices':selection['training_choices'],'heldout_verdicts':primary['frozen_choice_verdicts'],
        'protocol_sha256':PROTOCOL_SHA,'input_sha256':{name:sha(ROOT/'inputs'/name) for name in INPUT_SHA},
        'runtime':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,
                   'sympy':sympy.__version__,'mpmath':mpmath.__version__,'matplotlib':matplotlib.__version__},
        'scope_limitations':'Frozen upstream T2 geometry, finite numerical grids, prescribed candidate families, dimensionless evolution; no global approximation impossibility or physical identification.'}
    (ROOT/'VALIDATION.json').write_text(json.dumps(validation,indent=2)+'\n')
    print(json.dumps({'status':validation['status'],'checks_passed':validation['checks_passed'],
                      'checks_total':validation['checks_total'],'research_outcome':validation['research_outcome'],
                      'cross_route_maximum_differences':differences},indent=2))
    if not all(checks.values()):
        raise SystemExit(1)


if __name__=='__main__':
    main()
