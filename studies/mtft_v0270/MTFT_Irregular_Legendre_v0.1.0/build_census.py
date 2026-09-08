"""Build an exact regular/irregular prime census and preserve all witnesses."""
import argparse
import hashlib
import json
from math import isqrt
from pathlib import Path
import platform
import time
import numpy as np
from irregular_arithmetic import prime_sieve, bernoulli_residues, witness_indices

ROOT = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--n-max', type=int, default=100)
    parser.add_argument('--prior-census', type=Path)
    args = parser.parse_args()
    if not 32 <= args.n_max <= 200:
        parser.error('n-max must be between 32 and 200')
    started = time.monotonic()
    limit = (args.n_max + 1) ** 2
    flags = prime_sieve(2 * limit + 3)
    primes = np.flatnonzero(flags[:limit]).astype(np.int64)
    square_n = np.array([isqrt(int(p)) for p in primes], dtype=np.int64)
    rows, chunks, offsets = [], [], [0]
    indices = np.full(len(primes), -1, dtype=np.int16)
    next_progress = time.monotonic() + 10
    for j, pp in enumerate(primes):
        p = int(pp)
        residues = np.array([], dtype=np.uint16) if p == 2 else bernoulli_residues(p).astype(np.uint16)
        ks = witness_indices(residues)
        if p > 2:
            indices[j] = len(ks)
        chunks.append(residues)
        offsets.append(offsets[-1] + len(residues))
        rows.append({'p': p, 'square_n': int(square_n[j]), 'residue_mod30': p % 30,
                     'regularity': 'not_applicable' if p == 2 else ('irregular' if ks else 'regular'),
                     'irregularity_index': None if p == 2 else len(ks), 'witness_indices': ks})
        if time.monotonic() >= next_progress:
            print(json.dumps({'completed_primes': j + 1, 'total_primes': len(primes), 'through_prime': p}), flush=True)
            next_progress = time.monotonic() + 10
    labels = indices > 0
    regular = indices == 0
    n = np.arange(1, args.n_max + 1, dtype=np.int64)
    total = np.bincount(square_n, minlength=args.n_max + 1)[1:]
    irregular_counts = np.bincount(square_n[labels], minlength=args.n_max + 1)[1:]
    regular_counts = np.bincount(square_n[regular], minlength=args.n_max + 1)[1:]
    exceptions = np.bincount(square_n[primes == 2], minlength=args.n_max + 1)[1:]
    assert np.array_equal(total, irregular_counts + regular_counts + exceptions)
    quotients = np.array([-1 if p == 2 else (pow(2, int(p) - 1, int(p) ** 2) - 1) // int(p) for p in primes], dtype=np.int64)
    marks = np.column_stack([flags[primes - 2] | flags[primes + 2], flags[2 * primes + 1],
                             (primes > 2) & flags[(primes - 1) // 2], (primes > 2) & (quotients == 0)])
    np.savez_compressed(ROOT / 'census.npz', p=primes, square_n=square_n, residue_mod30=primes % 30,
                        position_numerator=primes - square_n ** 2, position_denominator=2 * square_n + 1,
                        irregularity_index=indices, irregular=labels, regular=regular,
                        bernoulli_offsets=np.array(offsets, dtype=np.int64),
                        bernoulli_residues=np.concatenate(chunks), n=n, interval_total=total,
                        interval_irregular=irregular_counts, interval_regular=regular_counts,
                        interval_exceptions=exceptions, other_properties=marks,
                        fermat_quotient_base2=quotients)
    (ROOT / 'prime_witnesses.json').write_text(json.dumps(rows, indent=2) + '\n')
    prior = None
    if args.prior_census:
        old = np.load(args.prior_census)
        match = old['p'] < limit
        assert np.array_equal(old['p'][match], primes)
        assert np.array_equal(old['square_n'][match], square_n)
        assert np.array_equal(old['properties'][match], marks)
        assert np.array_equal(old['fermat_quotient_base2'][match], quotients)
        prior = {'filename': args.prior_census.name, 'sha256': hashlib.sha256(args.prior_census.read_bytes()).hexdigest(),
                 'all_overlapping_primes_positions_and_marks_match': True}
    summary = {'n_max': args.n_max, 'prime_upper_bound_exclusive': limit,
               'prime_count': len(primes), 'odd_prime_count': int(np.sum(primes > 2)),
               'irregular_count': int(labels.sum()), 'regular_count': int(regular.sum()),
               'witness_pair_count': int(indices[indices >= 0].sum()),
               'irregularity_index_histogram': {str(k): int(np.sum(indices == k)) for k in np.unique(indices) if k >= 0},
               'irregular_free_intervals': n[irregular_counts == 0].tolist(),
               'regular_free_intervals_n_ge_2': n[(n >= 2) & (regular_counts == 0)].tolist(),
               'prime_free_intervals': n[total == 0].tolist(),
               'wieferich_records': [rows[i] for i in np.flatnonzero(marks[:, 3])],
               'prior_census_reconciliation': prior, 'python_version': platform.python_version(),
               'numpy_version': np.__version__, 'elapsed_seconds': time.monotonic() - started,
               'arithmetic_status': 'exact finite modular arithmetic; not a Lean proof'}
    (ROOT / 'census_summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == '__main__':
    main()
