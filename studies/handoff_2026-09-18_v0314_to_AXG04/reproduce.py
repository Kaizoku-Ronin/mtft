#!/usr/bin/env python3
"""Run frozen studies on isolated copies with the current Python interpreter."""
from pathlib import Path
import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tarfile
import time

ROOT = Path(__file__).resolve().parents[1]
CATALOG = json.loads((ROOT/'study_catalog.json').read_text())

def source_tree(work):
    archive = ROOT/'baseline/mtft-0.31.4.tar.gz'
    expected = json.loads((ROOT/'provenance.json').read_text())['baseline_sha256']
    if hashlib.sha256(archive.read_bytes()).hexdigest() != expected:
        raise ValueError('Baseline archive digest differs from provenance.')
    destination = work/'baseline'
    if destination.exists():
        return destination/'mtft-0.31.4'
    destination.mkdir()
    with tarfile.open(archive, 'r:gz') as t:
        for entry in t.getmembers():
            p = Path(entry.name)
            if p.is_absolute() or '..' in p.parts or not p.parts or p.parts[0] != 'mtft-0.31.4':
                raise ValueError('Unsafe baseline path: '+entry.name)
            if not (entry.isfile() or entry.isdir()):
                raise ValueError('Unsupported baseline member: '+entry.name)
        t.extractall(destination, filter='data')
    return destination/'mtft-0.31.4'

def command(study, target, source):
    sid = study['id']
    if sid == 'V0314-AUDIT':
        return [sys.executable, str(target/'audit_v0314.py'), '--source', str(source), '--output', str(target/'reproduced_results.json')]
    if sid == 'SC7-01':
        return [sys.executable, str(target/'investigate_spin_circle.py'), '--source', str(source),
                '--archive', str(ROOT/'baseline/mtft-0.31.4.tar.gz'), '--output', str(target/'reproduced_results.json')]
    if sid == 'STUDY-EDITION':
        return [sys.executable, str(target/'verify_hand_calculations.py')]
    if sid == 'HOPF-02':
        return [sys.executable, str(target/'investigate.py')]
    return [sys.executable, str(target/'run_all.py')]

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--list', action='store_true', help='List available studies without executing them.')
    parser.add_argument('--study', choices=[s['id'] for s in CATALOG]+['all'])
    parser.add_argument('--workdir', type=Path, default=ROOT.parent/'MTFT_reproduced')
    args = parser.parse_args()
    if args.list:
        for study in CATALOG:
            count = study['recorded_exact_check_count']
            print(study['id']+': '+(str(count)+' recorded exact checks' if count else 'see original result ledger'))
        return 0
    if not args.study:
        parser.error('Choose --study or --list.')
    work = args.workdir.resolve()
    if work.is_relative_to(ROOT):
        parser.error('--workdir must be outside the preserved handoff folder.')
    selected = CATALOG if args.study == 'all' else [s for s in CATALOG if s['id'] == args.study]
    for study in selected:
        if (work/study['id']).exists():
            parser.error('Refusing to overwrite existing study folder: '+str(work/study['id']))
    work.mkdir(parents=True, exist_ok=True)
    source = source_tree(work) if any(s['id'] in {'V0314-AUDIT','SC7-01'} for s in selected) else None
    status_path = work/'rerun_status.json'
    status = json.loads(status_path.read_text()) if status_path.exists() else {'python':sys.version,'runs':[]}
    env = dict(os.environ)
    env['PYTHONDONTWRITEBYTECODE'] = '1'
    for study in selected:
        target = work/study['id']
        shutil.copytree(ROOT/study['path'], target)
        log = target/'handoff_rerun.log'
        started = time.monotonic()
        print('Running '+study['id']+'; log: '+str(log), flush=True)
        with log.open('w') as out:
            result = subprocess.run(command(study,target,source), cwd=target, stdout=out, stderr=subprocess.STDOUT, env=env)
        status['runs'].append({'study':study['id'], 'returncode':result.returncode,
                              'elapsed_seconds':round(time.monotonic()-started,3), 'log':str(log),
                              'meaning':'Process status; inspect result ledger for mathematical and physical gates.'})
        status_path.write_text(json.dumps(status,indent=2)+'\n')
        print(study['id']+': return code '+str(result.returncode),flush=True)
        if result.returncode:
            return result.returncode
    return 0

if __name__ == '__main__':
    sys.exit(main())
