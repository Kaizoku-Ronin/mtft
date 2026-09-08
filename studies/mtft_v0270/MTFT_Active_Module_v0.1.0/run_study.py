#!/usr/bin/env python3
"""Reproduce the active-module study against an extracted MTFT release."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source',type=Path,required=True)
    args=p.parse_args()
    source=args.source.resolve()
    root=Path(__file__).resolve().parent
    env=os.environ.copy()
    env['PYTHONPATH']=str(source/'src')+os.pathsep+env.get('PYTHONPATH','')
    env['OPENBLAS_NUM_THREADS']='1'
    commands=[['frame_check.py'],['kernel_check.py'],
              ['coupling_study.py','--source',str(source)],
              ['independent_coupling_check.py'],['render_coupling.py']]
    with (root/'run_log.txt').open('w') as log:
        for command in commands:
            result=subprocess.run([sys.executable,*command],cwd=root,env=env,
                                  capture_output=True,text=True)
            log.write('\n$ python '+' '.join(command)+'\n')
            log.write(result.stdout)
            log.write(result.stderr)
            log.flush()
            print(command[0]+': '+('completed' if result.returncode==0 else 'FAILED'),flush=True)
            if result.returncode:
                raise SystemExit(result.returncode)
    r=json.loads((root/'coupling_results.json').read_text())
    f=json.loads((root/'frame_check_results.json').read_text())
    k=json.loads((root/'kernel_check_results.json').read_text())
    v=json.loads((root/'independent_coupling_check_results.json').read_text())
    runs=r['runs']
    checks={
        'release_is_0_26_2':r['package_version']=='0.26.2',
        'exact_harmonic_support_rank_four':k['exact_support']['rank']==4 and k['exact_support']['reconstruction_exact'],
        'frame_metric_and_complex_structure_consistent':f['frame']['Ri_transpose_G_Ri_vs_I_rel']<1e-10 and f['frame']['adapted_J_vs_standard_rel']<1e-10,
        'all_direct_projectors_rank_eight_without_ambiguous_singular_values':all(x['complex_active_rank']==8 and not x['direct_seed_spectrum']['ambiguous'] for x in runs),
        'projectors_stable_across_registered_frames':all(x.get('baseline_comparison',{}).get('transported_P_difference',0)<1e-9 for x in runs),
        'block_fractions_stable_across_registered_frames':all(x.get('baseline_comparison',{}).get('max_block_fraction_difference',0)<1e-10 for x in runs),
        'every_registered_non_STAR_operator_resolves_mixing':all(o['classification']=='resolved_mixing' for x in runs for name,o in x['operators'].items() if name!='STAR'),
        'STAR_preserves_split_at_diagnostic_threshold':all(x['operators']['STAR']['classification']=='consistent_with_preservation' for x in runs),
        'independent_localization_agrees':v['checks_passed'] and v['maximum_absolute_difference_from_main']<1e-10,
        'T2_one_augmentation_fills_real_26_with_large_singular_gap':v['T2_one_augmentation']['rank_at_relative_1e_9']==26 and v['T2_one_augmentation']['smallest_normalized_singular_value']>1e-3,
        'all_registered_hulls_fill_real_26':all(h['status']=='resolved_numerically' and h['real_dimension']==26 for h in r['hulls'].values()),
        'identity_and_fullspace_controls':r['controls']['identity']['active_to_fixed_over_operator']<1e-12 and r['controls']['full_space_T2']['active_to_fixed_over_operator']<1e-12,
        'two_exact_fullspace_commutations':all(x['all_entries_exactly_zero'] for x in v['exact_full_commutators'].values())}
    output={'checks':checks,'all_checks_passed':all(checks.values()),
            'scope':'Study replay and consistency gates; not full package tests or interval certification.',
            'unresolved_package_closure_runs':[{'seed':x['frame_seed'],'dps':x['dps_period_inputs'],'error':x['checks'].get('closure_error')}
                for x in runs if x['checks']['closure_status']!='passed_package_gate']}
    (root/'VALIDATION.json').write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps(output,indent=2))
    if not output['all_checks_passed']: raise SystemExit(1)


if __name__=='__main__': main()
