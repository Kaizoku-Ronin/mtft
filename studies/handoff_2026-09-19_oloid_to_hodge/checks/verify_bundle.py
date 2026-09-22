"""Verify delivery integrity; this does not rerun the research calculations."""
from hashlib import sha256
import json
from pathlib import Path
import sys
from zipfile import ZipFile


def digest(path):
    with path.open('rb') as handle:
        return sha256(handle.read()).hexdigest()


def main():
    root = Path(__file__).resolve().parents[1]
    manifest = json.loads((root / 'MANIFEST.json').read_text(encoding='utf-8'))
    errors = []
    for name, expected in manifest['files'].items():
        path = root / name
        if not path.is_file():
            errors.append('Missing: ' + name)
        elif digest(path) != expected['sha256']:
            errors.append('Hash mismatch: ' + name)
        elif path.stat().st_size != expected['bytes']:
            errors.append('Size mismatch: ' + name)

    historical_count = 0
    for path in sorted((root / 'studies').glob('*/MANIFEST.json')):
        data = json.loads(path.read_text(encoding='utf-8'))
        for name, value in data['files'].items():
            expected = value if isinstance(value, str) else value['sha256']
            target = path.parent / name
            historical_count += 1
            if not target.is_file() or digest(target) != expected:
                errors.append('Original manifest mismatch: ' + str(target.relative_to(root)))

    archive_count = 0
    for archive in sorted((root / 'original_archives').glob('*.zip')):
        archive_count += 1
        with ZipFile(archive) as zipped:
            if zipped.testzip() is not None:
                errors.append('ZIP CRC failure: ' + archive.name)
            for member in zipped.infolist():
                if member.is_dir():
                    continue
                expanded = root / 'studies' / member.filename
                if not expanded.is_file() or expanded.read_bytes() != zipped.read(member):
                    errors.append('Expanded archive mismatch: ' + member.filename)

    source = root / 'sources' / 'mtft-0.32.0.tar.gz'
    expected_source = '46ffc563a0cd81a7a82d975497bf64d82f8471eb424409e3d154532b5d8d6069'
    if not source.is_file() or digest(source) != expected_source:
        errors.append('Reviewed source archive mismatch')

    # The text checksum list also covers MANIFEST.json; it cannot cover itself.
    for line in (root / 'SHA256SUMS.txt').read_text(encoding='utf-8').splitlines():
        expected, name = line.split('  ', 1)
        path = root / name
        if not path.is_file() or digest(path) != expected:
            errors.append('Checksum list mismatch: ' + name)

    report = {
        'status': 'FAIL' if errors else 'PASS',
        'payload_files_checked': len(manifest['files']),
        'historical_payload_hashes_checked': historical_count,
        'original_archives_checked': archive_count,
        'research_experiments_rerun': False,
        'errors': errors,
    }
    print(json.dumps(report, indent=2))
    return bool(errors)


if __name__ == '__main__':
    sys.exit(main())
