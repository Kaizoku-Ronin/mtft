"""Conditional label randomizations for irregular primes in square intervals."""
import argparse
from fractions import Fraction
import json
from math import comb
from pathlib import Path
import time
import numpy as np

ROOT = Path(__file__).resolve().parent
STAT_NAMES = ('irregular_free_intervals', 'longest_irregular_free_run', 'thin_fraction_difference')


def longest_run(mask):
    """Longest True run per row, for a 1D or 2D Boolean array."""
    mask = np.asarray(mask, dtype=bool)
    one = mask.ndim == 1
    if one:
        mask = mask[None, :]
    running = np.zeros(len(mask), dtype=np.int64)
    maximum = np.zeros(len(mask), dtype=np.int64)
    for column in mask.T:
        running = (running + 1) * column
        maximum = np.maximum(maximum, running)
    return int(maximum[0]) if one else maximum


def make_groups(square_n, residues, model):
    if model == 'global':
        return [np.arange(len(square_n), dtype=np.int64)]
    if model != 'stratified':
        raise ValueError(model)
    keys = ((square_n - 6) // 20) * 30 + residues
    return [np.flatnonzero(keys == key) for key in np.unique(keys)]


def exact_expectations(labels, bin_ids, groups, n_bins):
    """Finite-population expectation, variance and exact zero probabilities."""
    means = [Fraction(0) for _ in range(n_bins)]
    variances = [Fraction(0) for _ in range(n_bins)]
    zeros = [Fraction(1) for _ in range(n_bins)]
    strata = []
    for ids in groups:
        size = len(ids)
        marked = int(labels[ids].sum())
        occupancy = np.bincount(bin_ids[ids], minlength=n_bins)
        q = Fraction(marked, size)
        strata.append({'prime_count': size, 'irregular_count': marked,
                       'interval_occupancy': occupancy.tolist()})
        for j in np.flatnonzero(occupancy):
            m = int(occupancy[j])
            means[j] += m * q
            if size > 1:
                variances[j] += m * q * (1 - q) * Fraction(size - m, size - 1)
            zeros[j] *= Fraction(comb(size - marked, m) if m <= size - marked else 0, comb(size, m))
    return (np.array([float(x) for x in means]), np.array([float(x) for x in variances]),
            zeros, strata)


def statistics(counts, total_counts, thin):
    counts = np.asarray(counts)
    one = counts.ndim == 1
    if one:
        counts = counts[None, :]
    empty = counts == 0
    delta = (counts[:, thin].sum(axis=1) / total_counts[thin].sum()
             - counts[:, ~thin].sum(axis=1) / total_counts[~thin].sum())
    out = np.column_stack([empty.sum(axis=1), longest_run(empty), delta])
    return out[0] if one else out


def holm(pvalues):
    pvalues = np.asarray(pvalues, dtype=float)
    order = np.argsort(pvalues, kind='stable')
    adjusted = np.empty(len(pvalues))
    adjusted[order] = np.minimum(1., np.maximum.accumulate(pvalues[order] * np.arange(len(pvalues), 0, -1)))
    return adjusted


def run_scenario(data, n_min, model, repetitions, seed):
    started = time.monotonic()
    n_max = int(data['n'][-1])
    keep = (data['square_n'] >= n_min) & (data['square_n'] <= n_max)
    ns = np.arange(n_min, n_max + 1)
    if len(ns) < 4:
        raise ValueError('Each analysis range needs at least four intervals to define a bottom quarter.')
    bin_ids = data['square_n'][keep] - n_min
    labels = data['irregular'][keep]
    groups = make_groups(data['square_n'][keep], data['residue_mod30'][keep], model)
    totals = np.bincount(bin_ids, minlength=len(ns))
    observed = np.bincount(bin_ids[labels], minlength=len(ns))
    density_reference = (2 * ns + 1) / np.log(ns * ns + ns + .5)
    relative_density = totals / density_reference
    thin = np.zeros(len(ns), dtype=bool)
    thin[np.lexsort((ns, relative_density))[:len(ns) // 4]] = True
    expected, variance, zero_probs, strata = exact_expectations(labels, bin_ids, groups, len(ns))
    thin_center = (expected[thin].sum() / totals[thin].sum()
                   - expected[~thin].sum() / totals[~thin].sum())
    observed_stats = statistics(observed, totals, thin)
    rng = np.random.default_rng(seed)
    counts = np.empty((repetitions, len(ns)), dtype=np.uint16)
    randomized = np.empty(len(labels), dtype=bool)
    for r in range(repetitions):
        if model == 'global':
            randomized[:] = rng.permutation(labels)
        else:
            for ids in groups:
                randomized[ids] = rng.permutation(labels[ids])
        counts[r] = np.bincount(bin_ids[randomized], minlength=len(ns))
    null_stats = statistics(counts, totals, thin)
    null_quantiles = np.quantile(null_stats, [.025, .5, .975], axis=0, method='inverted_cdf')
    tests = []
    for j, name in enumerate(STAT_NAMES):
        if j < 2:
            tail = int(np.count_nonzero(null_stats[:, j] >= observed_stats[j]))
            alternative = 'excess; upper tail including ties'
            center = None
        else:
            threshold = abs(observed_stats[j] - thin_center)
            tail = int(np.count_nonzero(np.abs(null_stats[:, j] - thin_center) >= threshold - 1e-14))
            alternative = 'two-sided around analytic conditional expectation; including ties'
            center = float(thin_center)
        tests.append({'statistic': name, 'observed': float(observed_stats[j]),
                      'null_mean_mc': float(null_stats[:, j].mean()),
                      'null_central_95_percent': [float(null_quantiles[0, j]), float(null_quantiles[2, j])],
                      'analytic_center_for_two_sided_test': center, 'tail_count': tail,
                      'alternative': alternative, 'p_value_mc': (tail + 1) / (repetitions + 1)})
    interval_quantiles = np.quantile(counts, [.025, .5, .975], axis=0, method='inverted_cdf')
    summary = {'key': f'n{n_min}_{model}', 'n_min': n_min, 'n_max': n_max, 'model': model,
               'prime_count': len(labels), 'irregular_count': int(labels.sum()), 'stratum_count': len(groups),
               'repetitions': repetitions, 'seed': seed, 'rng': 'NumPy Generator PCG64',
               'thin_interval_n': ns[thin].tolist(),
               'thin_pool': {'prime_count': int(totals[thin].sum()), 'irregular_count': int(observed[thin].sum())},
               'other_pool': {'prime_count': int(totals[~thin].sum()), 'irregular_count': int(observed[~thin].sum())},
               'irregular_free_interval_n': ns[observed == 0].tolist(),
               'expected_irregular_free_count_exact_rational': str(sum(zero_probs)),
               'expected_irregular_free_count': float(sum(zero_probs)),
               'zero_count_mc_standard_error': float(null_stats[:, 0].std(ddof=1) / np.sqrt(repetitions)),
               'tests': tests, 'elapsed_seconds': time.monotonic() - started,
               'intervals': [{'n': int(n), 'total': int(totals[j]), 'irregular': int(observed[j]),
                              'expected_irregular': float(expected[j]), 'variance_irregular': float(variance[j]),
                              'zero_probability': float(zero_probs[j]),
                              'zero_probability_exact': str(zero_probs[j]),
                              'null_pointwise_quantiles_025_50_975': interval_quantiles[:, j].tolist(),
                              'relative_prime_density': float(relative_density[j]), 'thin': bool(thin[j])}
                             for j, n in enumerate(ns)], 'strata': strata}
    return summary, null_stats


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repetitions', type=int, default=49999)
    parser.add_argument('--seed', type=int, default=20260906)
    args = parser.parse_args()
    if args.repetitions < 999:
        parser.error('use at least 999 randomizations')
    data = np.load(ROOT / 'census.npz')
    scenarios, draws = [], {}
    for index, (n_min, model) in enumerate([(6, 'global'), (6, 'stratified'), (32, 'global'), (32, 'stratified')]):
        scenario, stats = run_scenario(data, n_min, model, args.repetitions, args.seed + index)
        scenarios.append(scenario)
        draws[scenario['key']] = stats
        print(json.dumps({k: scenario[k] for k in ('key', 'expected_irregular_free_count', 'irregular_free_interval_n', 'tests', 'elapsed_seconds')}), flush=True)
    all_tests = [test for scenario in scenarios for test in scenario['tests']]
    adjusted = holm([test['p_value_mc'] for test in all_tests])
    for test, adj in zip(all_tests, adjusted):
        test['p_value_holm_12'] = float(adj)
    result = {'scenarios': scenarios, 'holm_family_size': len(all_tests),
              'earlier_explored_bound_exclusive': 1024,
              'inference_scope': 'Conditional organization of labels on this fixed finite prime census; not a test or proof of Legendre.'}
    (ROOT / 'analysis_results.json').write_text(json.dumps(result, indent=2) + '\n')
    np.savez_compressed(ROOT / 'permutation_statistics.npz', **draws)


if __name__ == '__main__':
    main()
