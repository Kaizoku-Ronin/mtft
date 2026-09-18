#!/usr/bin/env python3
"""Verify handoff bytes and original-archive correspondence; no math rerun."""
from pathlib import Path
import hashlib
import json
import sys
import zipfile

ROOT = Path(__file__).resolve().parents[1]

def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()

def confined(relative):
    p = Path(relative)
    if p.is_absolute() or '..' in p.parts:
        raise ValueError('Unsafe relative path: '+relative)
    target = (ROOT/p).resolve()
    if not target.is_relative_to(ROOT):
        raise ValueError('Path outside handoff: '+relative)
    return target

def main():
    failures = []
    count = 0
    for line in (ROOT/'SHA256SUMS').read_text().splitlines():
        expected, relative = line.split('  ', 1)
        p = confined(relative)
        count += 1
        if not p.is_file() or sha(p) != expected:
            failures.append('Payload hash: '+relative)
    provenance = json.loads((ROOT/'provenance.json').read_text())
    for source in provenance['sources']:
        for relative in source['included_paths']:
            if sha(confined(relative)) != source['sha256']:
                failures.append('Source hash: '+relative)
    grouped = {}
    for record in provenance['archive_members']:
        grouped.setdefault(record['archive'], []).append(record)
    members = 0
    for archive, records in grouped.items():
        with zipfile.ZipFile(confined(archive)) as z:
            bad = z.testzip()
            if bad:
                failures.append('ZIP CRC: '+archive+'/'+bad)
            for record in records:
                payload = z.read(record['member'])
                if payload != confined(record['path']).read_bytes():
                    failures.append('Archive member differs: '+record['path'])
                members += 1
    print(json.dumps({'payload_files_checked':count, 'original_archives_checked':len(grouped),
                      'original_members_compared':members, 'failures':failures,
                      'meaning':'Integrity verification only; mathematics not rerun.'}, indent=2))
    return bool(failures)

if __name__ == '__main__':
    sys.exit(main())
