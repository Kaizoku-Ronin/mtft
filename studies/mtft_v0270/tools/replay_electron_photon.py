"""Replay the finite electron-photon study on an explicitly selected source.

Copies only required original scripts/inputs into a new destination. The frozen
study is never modified; the old local-patch check is replaced by CC-21/CC-22.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys


SCRIPTS = ['common.py', 'charge_study.py', 'local_qed_covariance.py',
           'photon_atomic.py', 'magnetic_and_coupling.py', 'render_results.py']
INPUTS = ['inputs/x0143_certified.json', 'inputs/sector_certificate.json']


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--study', type=Path, required=True)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--render', action='store_true')
    args = parser.parse_args()
    study, source, output = args.study.resolve(), args.source.resolve(), args.output.resolve()
    checker = Path(__file__).with_name('check_release_contracts.py')
    required = [study / name for name in SCRIPTS + INPUTS]
    required += [source / 'src/mtft/__init__.py', checker]
    if any(not path.is_file() for path in required):
        raise SystemExit('Required study inputs, source tree, or adjacent contract checker are missing.')
    if output.exists() or output.is_relative_to(study) or output.is_relative_to(source):
        raise SystemExit('Output must be a new directory outside the frozen study and source tree.')
    output.mkdir(parents=True)
    copied = {}
    for name in SCRIPTS + INPUTS:
        destination = output / name
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(study / name, destination)
        copied[name] = sha(destination)
    source_hashes = {str(path.relative_to(source)): sha(path)
                     for path in sorted((source / 'src/mtft').rglob('*.py'))}
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE='1', OPENBLAS_NUM_THREADS='1')
    env['PYTHONPATH'] = str(source / 'src') + os.pathsep + env.get('PYTHONPATH', '')
    commands = [[sys.executable, str(checker), '--source', str(source),
                 '--output', str(output / 'release_contract_results.json')]]
    commands += [[sys.executable, name] for name in SCRIPTS[1:5]]
    if args.render:
        commands.append([sys.executable, 'render_results.py'])
    record = {'schema_version': 1,
              'scope': 'New-source replay of the finite electron-photon experiment; supplied QED structure and frozen physical references retain their original limitations.',
              'historical_study': str(study), 'source': str(source),
              'copied_file_sha256': copied, 'source_python_sha256': source_hashes,
              'replacement': 'Historical check_corrections.py is superseded by the adjacent CC-21/CC-22 checker; the old source patch is not applied.',
              'commands': [], 'all_commands_pass': False}
    (output / 'REPLAY_README.md').write_text(
        '# New-source electron-photon replay\n\n'
        'This directory contains fresh numerical outputs from the source selected in REPLAY_PROVENANCE.json. '
        'It is not the original frozen study and does not replace its report or manifest. '
        'The original physical inputs, supplied QED structure, and finite experiment scope are unchanged. '
        'The convention check uses CC-21/CC-22 accessors. No historical source patch was applied.\n')
    with (output / 'replay_log.txt').open('w') as log:
        for command in commands:
            print('Running', Path(command[1]).name, flush=True)
            result = subprocess.run(command, cwd=output, env=env, text=True,
                                    stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            row = {'argv': command, 'returncode': result.returncode}
            record['commands'].append(row)
            log.write(json.dumps(row) + '\n' + result.stdout + '\n')
            log.flush()
            (output / 'REPLAY_PROVENANCE.json').write_text(json.dumps(record, indent=2) + '\n')
            if result.returncode:
                raise SystemExit(result.returncode)
    record['all_commands_pass'] = True
    record['outputs_sha256'] = {path.name: sha(path) for path in sorted(output.iterdir())
                               if path.is_file() and path.name not in SCRIPTS + ['REPLAY_PROVENANCE.json']}
    (output / 'REPLAY_PROVENANCE.json').write_text(json.dumps(record, indent=2) + '\n')
    print(json.dumps({'all_commands_pass': True, 'commands': len(commands)}))


if __name__ == '__main__':
    main()
