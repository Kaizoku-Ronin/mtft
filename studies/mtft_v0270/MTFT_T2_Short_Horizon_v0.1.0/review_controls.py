"""Independent mathematical and numerical controls for short-horizon bounds."""
from pathlib import Path
import argparse
import json
import math
import sys

import mpmath as mp
import numpy as np
import sympy as sp
from scipy.linalg import expm

ROOT = Path(__file__).resolve().parent
mp.mp.dps = 80


def mp_exact_float(x):
    n, d = float(x).as_integer_ratio()
    return mp.mpf(n) / d


def mp_matrix(x):
    return mp.matrix([[mp.mpc(mp_exact_float(v.real), mp_exact_float(v.imag))
                       for v in row] for row in x])


def fnorm(x):
    return mp.sqrt(mp.fsum(abs(x[i, j]) ** 2
                          for i in range(x.rows) for j in range(x.cols)))


def run(study):
    data = np.load(study / 'inputs/T2_blocks.npz')
    q = np.load(study / 'inputs/selected_models.npz')['Q_svd_r4']
    a, b, d = data['A'], data['B'], data['D']
    br = b @ q
    dr_raw = q.conj().T @ d @ q
    dr = (dr_raw + dr_raw.conj().T) / 2
    hf = np.block([[a, b], [b.conj().T, d]])
    hr = np.block([[a, br], [br.conj().T, dr]])
    assert np.array_equal(hf, hf.conj().T)
    assert np.array_equal(hr, hr.conj().T)
    mf, mr = mp_matrix(hf), mp_matrix(hr)
    pf, pr = mp.eye(13), mp.eye(12)
    coeff, tail = [], []
    for k in range(14):
        coeff.append(fnorm(pf[:8, :8] - pr[:8, :8]))
        tail.append(fnorm(pf[:, :8]) + fnorm(pr[:, :8]))
        pf, pr = pf * mf, pr * mr
    orders = [4, 6, 8, 10, 12]
    horizons = [mp.mpf(1)/8, mp.mpf(1)/4, mp.mpf(1)/2, mp.mpf(3)/4, mp.mpf(1)]

    def bound(t, m):
        return (mp.fsum(t**k * coeff[k] / mp.factorial(k) for k in range(m+1))
                + t**(m+1) * tail[m+1] / mp.factorial(m+1)) / mp.sqrt(3)

    table = [{'T': str(t), 'order': m, 'relative_bound': mp.nstr(bound(t,m), 65),
              'passes_one_percent_numerically': bool(bound(t,m) <= mp.mpf('0.01'))}
             for t in horizons for m in orders]
    checks = []
    def exact(name, condition):
        checks.append({'name': name, 'passed': bool(condition), 'kind': 'exact symbolic'})

    ai = sp.Matrix([[2, 1-sp.I], [1+sp.I, -1]])
    bi = sp.Matrix([1, sp.I])
    di = sp.Matrix([[3]])
    hi = ai.row_join(bi).col_join(bi.conjugate().T.row_join(di))
    for k in [0,1]:
        exact(f'discarded complex toy moment {k}', (hi**k)[:2,:2] == ai**k)
    exact('discarded complex toy second moment', (hi**2)[:2,:2]-ai**2 == bi*bi.conjugate().T)
    third = (hi**3)[:2,:2]-ai**3-(ai*bi*bi.conjugate().T + bi*bi.conjugate().T*ai + bi*di*bi.conjugate().T)
    exact('discarded complex toy third moment', third.applyfunc(sp.simplify) == sp.zeros(2))
    chain = sp.Matrix([[0,1,0],[1,0,1],[0,1,0]])
    reduced = sp.Matrix([[0,1],[1,0]])
    for k in range(4):
        exact(f'indirect coupling toy moment {k}', (chain**k)[0,0] == (reduced**k)[0,0])
    exact('indirect coupling toy fourth moment', (chain**4)[0,0] - (reduced**4)[0,0] == 1)
    t = sp.Symbol('t', real=True)
    series = sp.series((1+sp.cos(sp.sqrt(2)*t))/2-sp.cos(t), t, 0, 6).removeO()
    exact('indirect coupling fourth-order propagator', sp.expand(series) == t**4/24)
    block_diagonal = sp.diag(ai, 3)
    exact('zero coupling all moments through 13', all((block_diagonal**k)[:2,:2] == ai**k for k in range(14)))
    exact('dimension denominator toy', sp.simplify((sp.cos(t)**2+1)-1) == sp.cos(t)**2)
    assert all(c['passed'] for c in checks)

    sampled = []
    denominator_min = float('inf')
    unit_singular_residual = 0.0
    for time in np.linspace(0,1,161):
        uf = expm(-1j*time*hf)[:8,:8]
        ur = expm(-1j*time*hr)[:8,:8]
        denom = np.linalg.norm(uf)
        denominator_min = min(denominator_min, float(denom))
        sv = np.linalg.svd(uf, compute_uv=False)
        unit_singular_residual = max(unit_singular_residual, float(np.max(abs(sv[:3]-1))))
        observed = float(np.linalg.norm(uf-ur)/denom)
        sample_bounds = [float(bound(mp_exact_float(time), m)) for m in orders]
        sampled.append({'t': float(time), 'relative_error': observed,
                        'order_12_bound': sample_bounds[-1]})
        assert all(observed <= bb + 2e-14 for bb in sample_bounds)
        assert denom >= math.sqrt(3)-2e-14

    out = {
        'scope': 'Independent 80-digit diagnostics and exact toy controls; continuous bound proof is in the review.',
        'precision_decimal_digits': mp.mp.dps,
        'full_exact_float_hermitian': bool(np.array_equal(hf, hf.conj().T)),
        'reduced_exact_float_hermitian': bool(np.array_equal(hr, hr.conj().T)),
        'active_moment_difference_norms': [mp.nstr(v,65) for v in coeff],
        'sum_active_column_power_norms': [mp.nstr(v,65) for v in tail],
        'bound_table': table,
        'exact_controls': checks,
        'exact_controls_passed': sum(c['passed'] for c in checks),
        'exact_controls_total': len(checks),
        'sampled_times': len(sampled),
        'sampled_bounds_passed': True,
        'minimum_sampled_full_active_frobenius': denominator_min,
        'maximum_three_unit_singular_values_residual': unit_singular_residual,
        'sampled_evolution': sampled,
    }
    (ROOT / 'review_control_results.json').write_text(json.dumps(out, indent=2) + '\n')
    print(json.dumps({'exact_controls_passed': out['exact_controls_passed'],
                      'exact_controls_total': out['exact_controls_total'],
                      'sampled_times': len(sampled), 'sampled_bounds_passed': True,
                      'half_horizon_order12_bound': str(bound(mp.mpf('.5'),12))}, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--study', type=Path, default=ROOT,
                        help='Study root; defaults to the directory containing this script')
    run(parser.parse_args().study.resolve())
