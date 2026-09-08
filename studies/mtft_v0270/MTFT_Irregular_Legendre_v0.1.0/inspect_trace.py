"""Inspect a prime or square interval in the exact finite census."""
import argparse
import json
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent
MARK_NAMES = ('twin_member', 'sophie_germain', 'safe', 'wieferich_base2')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    choose = parser.add_mutually_exclusive_group(required=True)
    choose.add_argument('--prime', type=int)
    choose.add_argument('--interval', type=int)
    args = parser.parse_args()
    data = np.load(ROOT / 'census.npz')
    records = json.loads((ROOT / 'prime_witnesses.json').read_text())

    def row(j):
        result = dict(records[j])
        result['other_prime_properties'] = dict(zip(MARK_NAMES, data['other_properties'][j].tolist()))
        result['position_in_square_interval'] = [int(data['position_numerator'][j]), int(data['position_denominator'][j])]
        q = int(data['fermat_quotient_base2'][j])
        result['fermat_quotient_base2_mod_p'] = None if q < 0 else q
        return result

    if args.prime is not None:
        found = np.flatnonzero(data['p'] == args.prime)
        if not len(found):
            parser.error('that number is not a prime in this finite census')
        result = row(int(found[0]))
    else:
        n = args.interval
        if not 1 <= n <= int(data['n'][-1]):
            parser.error('interval is outside this finite census')
        result = {'n': n, 'open_interval': [n * n, (n + 1) ** 2],
                  'prime_count': int(data['interval_total'][n - 1]),
                  'irregular_count': int(data['interval_irregular'][n - 1]),
                  'regular_count': int(data['interval_regular'][n - 1]),
                  'exception_count': int(data['interval_exceptions'][n - 1]),
                  'primes': [row(int(j)) for j in np.flatnonzero(data['square_n'] == n)]}
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
