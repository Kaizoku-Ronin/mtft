"""Independent rational verification using SymPy's Gaussian rational domain.

The primary implementation uses pairs of integer matrices and common dyadic
scales. This verifier computes in QQ_I and checks the certificate's inequalities
by rational squaring. It does not call primary matrix or norm computations.
"""
from pathlib import Path
import argparse
from fractions import Fraction as F
from math import factorial, isqrt
import copy
import importlib.util
import json

import sympy as sp
from sympy.polys.domains import QQ_I
from sympy.polys.matrices import DomainMatrix

ROOT = Path(__file__).resolve().parent


def fraction(x):
    return F(int(x.numerator), int(x.denominator))


def rat(x):
    return F(int(x['numerator']), int(x['denominator']))


def rec(x):
    return {'numerator': str(x.numerator), 'denominator': str(x.denominator)}


def matrix(encoded):
    scale = 2 ** encoded['exponent']
    m = sp.Matrix([[sp.Rational(int(re), scale) + sp.I*sp.Rational(int(im), scale)
                    for re, im in zip(rr, ii)]
                   for rr, ii in zip(encoded['real'], encoded['imag'])])
    assert m == m.conjugate().T
    return DomainMatrix.from_Matrix(m).convert_to(QQ_I)


def square_norm(entries):
    return sum((fraction(z.x)**2 + fraction(z.y)**2 for row in entries for z in row), F())


def independent_sqrt_upper(x):
    # Different rational grid from the primary verifier; no floating arithmetic.
    scale = 2**128
    n, d = x.numerator, x.denominator
    radicand = n*d*scale*scale
    root = isqrt(radicand)
    if root*root != radicand:
        root += 1
    upper = F(root, d*scale)
    assert upper*upper >= x
    return upper


def run(study):
    models = json.loads((study/'FROZEN_MODELS.json').read_text())
    cert = json.loads((study/'BOUND_CERTIFICATE.json').read_text())
    h, hr = matrix(models['H']), matrix(models['Hr'])
    ident_h = DomainMatrix.eye((13,13), QQ_I)
    ident_r = DomainMatrix.eye((12,12), QQ_I)
    ph, pr = ident_h, ident_r
    diff_squared, column_squared = [], []
    for k in range(14):
        lh, lr = ph.to_list(), pr.to_list()
        diff = [[lh[i][j]-lr[i][j] for j in range(8)]for i in range(8)]
        diff_squared.append(square_norm(diff))
        column_squared.append((square_norm([row[:8] for row in lh]),
                               square_norm([row[:8] for row in lr])))
        if k < 13:
            ph, pr = ph*h, pr*hr
    coefficient_checks = []
    for k, squared in enumerate(diff_squared):
        supplied = rat(cert['difference_norm_uppers'][k])
        coefficient_checks.append(supplied >= 0 and supplied*supplied >= squared)
    assert all(coefficient_checks)
    row_checks = []
    for row in cert['rows']:
        horizon, m = rat(row['horizon']), row['order']
        finite = sum(horizon**k*rat(cert['difference_norm_uppers'][k])/factorial(k)
                     for k in range(m+1))
        assert finite == rat(row['finite_sum_upper'])
        a,b = column_squared[m+1]
        upper_sum = rat(row['column_tail_upper'])*factorial(m+1)/horizon**(m+1)
        # S >= sqrt(a)+sqrt(b) iff S >= 0, S^2-a-b >= 0, and
        # (S^2-a-b)^2 >= 4ab. All quantities here are rational.
        residual = upper_sum**2-a-b
        tail_ok = upper_sum >= 0 and residual >= 0 and residual**2 >= 4*a*b
        denom = rat(cert['denominator_lower'])
        denominator_ok = denom > 0 and denom*denom <= 3
        expected_relative = (finite+rat(row['column_tail_upper']))/denom
        relative_ok = expected_relative == rat(row['relative_bound_columns'])
        row_checks.append(tail_ok and denominator_ok and relative_ok)
    assert all(row_checks)
    t, m = F(1,2), 12
    independent_finite = sum(t**k*independent_sqrt_upper(diff_squared[k])/factorial(k)
                             for k in range(m+1))
    independent_tail = t**(m+1)*sum(independent_sqrt_upper(v) for v in column_squared[m+1])/factorial(m+1)
    numerator = independent_finite+independent_tail
    assert numerator*numerator <= F(3,10000)

    # Negative-input controls target the primary API, independently of the proof.
    spec = importlib.util.spec_from_file_location('primary_dyadic_bounds', study/'dyadic_bounds.py')
    primary = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(primary)
    rejected = []
    def reject(name, call):
        try:
            call()
        except (ValueError, TypeError):
            rejected.append({'name': name, 'passed': True})
        else:
            rejected.append({'name': name, 'passed': False})
    bad = copy.deepcopy(models['H']); bad['imag'][0][0] = '1'
    reject('imaginary diagonal', lambda: primary.decode(bad))
    bad2 = copy.deepcopy(models['H']); bad2['real'][0][1] = str(int(bad2['real'][0][1])+1)
    reject('non-Hermitian off-diagonal', lambda: primary.decode(bad2))
    bad3 = copy.deepcopy(models['H']); bad3['exponent'] = -1
    reject('negative dyadic exponent', lambda: primary.decode(bad3))
    bad4 = copy.deepcopy(models['H']); bad4['imag'].pop()
    reject('shape mismatch', lambda: primary.decode(bad4))
    reject('nonpositive dimension denominator', lambda: primary.certify(models['H'],models['Hr'],6,[F(1,2)],[4],F(1,100)))
    reject('negative horizon', lambda: primary.certify(models['H'],models['Hr'],8,[F(-1,2)],[4],F(1,100)))
    reject('negative Taylor order', lambda: primary.certify(models['H'],models['Hr'],8,[F(1,2)],[-1],F(1,100)))
    reject('negative square root', lambda: primary.sqrt_upper(F(-1)))
    assert all(item['passed'] for item in rejected)
    result = {
        'implementation_route': 'SymPy DomainMatrix over QQ_I; rational squared-inequality checks',
        'coefficient_enclosures_checked': len(coefficient_checks),
        'coefficient_enclosures_passed': all(coefficient_checks),
        'column_tail_and_relative_bounds_checked': len(row_checks),
        'column_tail_and_relative_bounds_passed': all(row_checks),
        'independent_order12_half_horizon_absolute_upper': rec(numerator),
        'independent_exact_one_percent_comparison': {'left': rec(numerator*numerator), 'right': rec(F(3,10000)), 'passed': True},
        'negative_input_controls': rejected,
        'negative_input_controls_passed': sum(item['passed'] for item in rejected),
    }
    (ROOT/'exact_verification_results.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['independent_order12_half_horizon_absolute_upper','independent_exact_one_percent_comparison','negative_input_controls']},indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--study', type=Path, default=ROOT,
                        help='Study root; defaults to the directory containing this script')
    run(parser.parse_args().study.resolve())
