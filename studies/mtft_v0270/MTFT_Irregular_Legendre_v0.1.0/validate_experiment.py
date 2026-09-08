"""Arithmetic and null-model checks independent of classification decisions."""
from fractions import Fraction
from itertools import combinations, product
import json
from math import isqrt
from pathlib import Path
import time
import numpy as np
from sympy import bernoulli
from irregular_arithmetic import verify_series_product, power_sum_residue, witness_indices
from analyze_traces import exact_expectations, longest_run, holm

ROOT = Path(__file__).resolve().parent


def check_null_formulas():
    bins = np.array([0, 0, 1, 1, 2, 2])
    passed = 0
    for groups, labels in [([np.arange(6)], np.array([1, 1, 0, 0, 0, 0], dtype=bool)),
                           ([np.array([0, 2, 4]), np.array([1, 3, 5])], np.array([1, 1, 0, 0, 0, 0], dtype=bool))]:
        expected, variance, zero_probs, _ = exact_expectations(labels, bins, groups, 3)
        choices = [list(combinations(ids.tolist(), int(labels[ids].sum()))) for ids in groups]
        outcomes = []
        for assignment in product(*choices):
            marked = np.array([j for part in assignment for j in part], dtype=int)
            outcomes.append(np.bincount(bins[marked], minlength=3))
        outcomes = np.array(outcomes)
        assert np.allclose(outcomes.mean(axis=0), expected, rtol=0, atol=1e-14)
        assert np.allclose(outcomes.var(axis=0), variance, rtol=0, atol=1e-14)
        for j in range(3):
            assert Fraction(int(np.sum(outcomes[:, j] == 0)), len(outcomes)) == zero_probs[j]
        passed += len(outcomes)
    assert longest_run([False, True, True, False, True]) == 2
    assert longest_run([True, True, True]) == 3
    assert longest_run([False, False]) == 0
    assert np.allclose(holm([.001, .04, .03]), [.003, .06, .06])
    return passed


def main():
    started = time.monotonic()
    data = np.load(ROOT / 'census.npz')
    ps, offsets, all_residues = data['p'], data['bernoulli_offsets'], data['bernoulli_residues']
    prime_set = set(ps.tolist())
    limit = (int(data['n'][-1]) + 1) ** 2
    for n in range(limit):
        trial = n >= 2 and all(n % d for d in range(2, isqrt(n) + 1))
        assert trial == (n in prime_set), n
    rational_cache = {k: bernoulli(k) for k in range(2, 1022, 2)}
    small_entries = products_checked = positive_checks = 0
    next_progress = time.monotonic() + 10
    for j, pp in enumerate(ps):
        p = int(pp)
        if p == 2:
            assert data['irregularity_index'][j] == -1
            continue
        residues = all_residues[offsets[j]:offsets[j + 1]]
        assert verify_series_product(p, residues), ('series_product', p)
        products_checked += 1
        ks = witness_indices(residues)
        assert len(ks) == int(data['irregularity_index'][j])
        if p < 1024:
            for index, actual in enumerate(residues):
                k = 2 + 2 * index
                exact = rational_cache[k]
                expected = int(exact.p) * pow(int(exact.q), -1, p) % p
                assert int(actual) == expected, (p, k)
                small_entries += 1
        for k in ks:
            assert power_sum_residue(p, k) == 0, ('positive_power_sum', p, k)
            positive_checks += 1
        if time.monotonic() >= next_progress:
            print(json.dumps({'series_products_checked': products_checked, 'through_prime': p, 'positive_checks': positive_checks}), flush=True)
            next_progress = time.monotonic() + 10
    rng = np.random.default_rng(2026090624)
    candidate_ids = np.flatnonzero(all_residues != 0)
    negative_ids = rng.choice(candidate_ids, size=min(128, len(candidate_ids)), replace=False)
    for flat in negative_ids:
        j = int(np.searchsorted(offsets, flat, side='right') - 1)
        p = int(ps[j]); k = int(2 + 2 * (flat - offsets[j]))
        assert power_sum_residue(p, k) == int(all_residues[flat]), ('negative_power_sum', p, k)
    anchors = {37: [32], 59: [44], 67: [58], 83: [], 89: [], 97: [], 101: [68], 103: [24], 691: [12, 200]}
    for p, ks in anchors.items():
        j = int(np.searchsorted(ps, p))
        assert witness_indices(all_residues[offsets[j]:offsets[j + 1]]) == ks
    assert np.array_equal(data['interval_total'], data['interval_regular'] + data['interval_irregular'] + data['interval_exceptions'])
    null_assignments = check_null_formulas()
    analysis_checks = []
    analysis_path = ROOT / 'analysis_results.json'
    if analysis_path.exists():
        analysis = json.loads(analysis_path.read_text())
        draws = np.load(ROOT / 'permutation_statistics.npz')
        for scenario in analysis['scenarios']:
            stats = draws[scenario['key']]
            se = scenario['zero_count_mc_standard_error']
            difference = abs(float(stats[:, 0].mean()) - scenario['expected_irregular_free_count'])
            assert difference <= 6 * se + 1e-12, ('MC versus exact zero expectation', scenario['key'])
            assert len(stats) == scenario['repetitions']
            for test in scenario['tests']:
                assert 0 < test['p_value_mc'] <= test['p_value_holm_12'] <= 1
                assert test['p_value_mc'] == (test['tail_count'] + 1) / (scenario['repetitions'] + 1)
            analysis_checks.append(scenario['key'])
    result = {'status': 'passed', 'integers_checked_by_trial_division': limit,
              'full_modular_series_product_certificates': products_checked,
              'bernoulli_residues_checked_against_exact_rationals_below_1024': small_entries,
              'all_positive_witnesses_checked_by_independent_power_sums': positive_checks,
              'negative_entries_checked_by_independent_power_sums': len(negative_ids),
              'known_prime_witness_anchors': len(anchors), 'count_partition_all_intervals': True,
              'exhaustively_enumerated_toy_null_assignments': null_assignments,
              'monte_carlo_expectation_and_tail_checks': analysis_checks,
              'elapsed_seconds': time.monotonic() - started,
              'limit': 'Finite exact computation and independent checks; not a Lean kernel proof.'}
    (ROOT / 'validation_results.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2), flush=True)


if __name__ == '__main__':
    main()
