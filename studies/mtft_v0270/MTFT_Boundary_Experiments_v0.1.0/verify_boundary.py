"""Cross-check saved experiment results, browser calculations, and CLI behavior.

Construction and all-subset gates run in boundary_model.py and
run_boundary_experiments.py. This script checks the delivery interfaces against
those computations; it does not repeat the full MTFT package test suite.
"""
import json
import hashlib
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import numpy as np
from boundary_model import BoundaryModel, read, write
from ising_extensions import graph

ROOT=Path(__file__).resolve().parent


def main():
    model=BoundaryModel();results=read('boundary_results.json')
    import mtft
    assert mtft.__version__=='0.26.0'
    _,n,edges=graph(143)
    assert n==56 and [list(e) for e in edges]==model.part['edges']
    tar=ROOT/'mtft-0.26.0.tar.gz'
    if tar.exists():
        assert hashlib.sha256(tar.read_bytes()).hexdigest()=='758f29a7a0c2ca08896a15964b3a607c7143ac1fd335bcc8f4d84a4880e53cdf'
    fixtures=[];maximum_saved_error=0.
    for name in ('ferro','one_reversed_bond'):
        for beta in (0.,.64,1.,3.):
            state=model.evaluate(beta,name);rows=model.all_sensor_masks(state)
            for mask,values in ((0,0),(255,255),(255,85),(1,1),(15,15)):
                obs=model.observation(beta,mask,values,name)
                fixtures.append({'model':name,'beta':beta,'mask':mask,'values':values,
                    'HA':state['H_A'],'floor':state['H_A_given_contacts'],'HB':state['H_B'],
                    'logZ':state['log_Z'],'information':rows[mask]['mutual_information'],
                    'meanEntropy':rows[mask]['H_A_given_observations'],
                    'spinAccuracy':rows[mask]['bayes_spin_accuracy'],
                    'fullAccuracy':rows[mask]['full_configuration_MAP_accuracy'],
                    'event':obs['event_probability'],'entropy':obs['H_A_given_this_observation'],
                    'accuracy':obs['bayes_spin_accuracy_given_this_observation'],
                    'bestFull':obs['best_full_configuration_probability'],
                    'meanPlus':obs['posterior_spin_plus']})
        saved=next(r for r in results['temperature_rows'][name] if r['beta']==.64)
        current=model.evaluate(.64,name)
        for key in ('H_A','H_B','all_boundary_MI','log_Z'):
            error=abs(saved[key]-current[key]);maximum_saved_error=max(maximum_saved_error,error)
            assert error<1e-11,(name,key,error)
    node=shutil.which('node') or os.environ.get('CODEX_PRIMARY_RUNTIME_NODE')
    browser_check={'status':'skipped','reason':'Node.js unavailable'}
    if node:
        with tempfile.TemporaryDirectory() as temp:
            fixture_path=Path(temp)/'fixtures.json'
            fixture_path.write_text(json.dumps(fixtures))
            proc=subprocess.run([node,str(ROOT/'verify_boundary_engine.js'),str(fixture_path)],
                                capture_output=True,text=True,check=True,cwd=ROOT)
            browser_check=json.loads(proc.stdout)
    commands=[['--beta','.64','--observe','6=+','18=-','24=+'],
              ['--beta','.64','--model','one_reversed_bond','--observe','all=+'],
              ['--beta','.64','--best-sensors','2']]
    cli=[]
    for args in commands:
        proc=subprocess.run([sys.executable,str(ROOT/'boundary_cli.py'),*args],
                            capture_output=True,text=True,check=True,cwd=ROOT)
        cli.append(json.loads(proc.stdout))
    assert cli[0]['observation']['mask']==11 and cli[0]['observation']['values']==9
    assert cli[1]['observation']['model']=='one_reversed_bond'
    assert cli[2]['masks_within_1e_minus9_bits']==[{'mask':5,'vertices':[6,22]}]
    # Independently verify spin reversal of outcome-specific reconstructions.
    maximum_flip_error=0.
    for name in ('ferro','one_reversed_bond'):
        for mask,values in ((255,85),(15,7),(1,1)):
            a=model.observation(.64,mask,values,name)
            b=model.observation(.64,mask,values^mask,name)
            error=float(np.max(abs(np.array(a['posterior_spin_plus'])+b['posterior_spin_plus']-1)))
            maximum_flip_error=max(maximum_flip_error,error)
            assert error<1e-12
            assert abs(a['event_probability']-b['event_probability'])<1e-12
    output={'status':'passed','mtft_version':mtft.__version__,'graph_matches_MTFT_source':True,
            'saved_result_max_abs_error':maximum_saved_error,
            'javascript_engine':browser_check,'CLI_examples_passed':len(cli),
            'spin_reversal_max_abs_error':maximum_flip_error,
            'browser_layout_QA':'Not completed: workspace browser URL policy blocked local HTTP and file previews.',
            'scope':'Interface verification; construction and 24,576 sensor cases have separate recorded gates.'}
    write('verification.json',output);print(json.dumps(output,indent=2))


if __name__=='__main__':main()
