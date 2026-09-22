#!/usr/bin/env python3
"""Superstable principal period-doubling cascade of x -> x*x + c.

Only Python's standard library is used. Decimal bisection uses floating-point
function evaluations: this is reproducible high-precision numerical evidence,
NOT interval arithmetic and NOT a rigorous certification of the roots/limits.

c_n denotes the center with critical orbit of primitive period 2**n. These are
superstable centers, not the parameter values where stable cycles bifurcate.
"""

from decimal import Decimal, localcontext
from pathlib import Path
import csv
import json
import time

HERE = Path(__file__).resolve().parent


def iterate(c, steps):
    x = Decimal(0)
    for _ in range(steps):
        x = x*x + c
    return x


def solve_root(n, previous, preceding, precision):
    """Find the first sign-changing root away from the preceding center.

The search interval is a fixed fraction of the preceding parameter gap, not
chosen using a tabulated Feigenbaum constant. Scanning toward decreasing c
avoids returning the lower-period center c_(n-1), itself a root of this iterate.
The resulting primitive period is checked by evaluating the halfway iterate.
"""
    gap = preceding - previous
    steps = 2**n
    # Principal roots lie in this modest interval for all orders computed here.
    fractions = [Decimal(k)/Decimal(400) for k in range(20, 161)]
    x_old = previous - gap*fractions[0]
    f_old = iterate(x_old, steps)
    bracket = None
    for fraction in fractions[1:]:
        x_new = previous - gap*fraction
        f_new = iterate(x_new, steps)
        if f_old*f_new < 0:
            bracket = (x_new, x_old, f_new, f_old)
            break
        x_old, f_old = x_new, f_new
    if bracket is None:
        raise ArithmeticError(f'No sign-changing root found at n={n}')
    lo, hi, f_lo, f_hi = bracket
    original_bracket = [str(lo), str(hi)]
    tolerance = Decimal(10)**(-(precision-12))
    iterations = 0
    while hi-lo > tolerance:
        mid = (lo+hi)/2
        f_mid = iterate(mid, steps)
        if f_mid == 0:
            # Decimal rounding may evaluate exactly zero; retain a numerical
            # precision-sized bracket rather than claim an exact algebraic root.
            lo, hi = mid-tolerance/2, mid+tolerance/2
            break
        if f_lo*f_mid < 0:
            hi, f_hi = mid, f_mid
        else:
            lo, f_lo = mid, f_mid
        iterations += 1
    root = (lo+hi)/2
    halfway = iterate(root, steps//2)
    residual = iterate(root, steps)
    if abs(halfway) < Decimal('1e-20'):
        raise ArithmeticError(f'Lower-period root suspected at n={n}')
    return root, halfway, residual, original_bracket, str(hi-lo), iterations


def run(precision, max_n=10):
    with localcontext() as context:
        context.prec = precision
        rows = [dict(n=0, period=1, c=Decimal(0), d=None,
                     residual=Decimal(0), root_status='exact'),
                dict(n=1, period=2, c=Decimal(-1), d=Decimal(-1),
                     residual=Decimal(0), root_status='exact')]
        for n in range(2, max_n+1):
            root, halfway, residual, bracket, width, count = solve_root(
                n, rows[-1]['c'], rows[-2]['c'], precision)
            rows.append(dict(n=n, period=2**n, c=root, d=halfway,
                             residual=residual, initial_bracket=bracket,
                             numerical_bracket_width=width,
                             bisection_steps=count,
                             root_status='numerical_not_interval_certified'))
        for n, row in enumerate(rows):
            c = row['c']
            # c=r(2-r)/4, branch r=1+sqrt(1-4c): x^2+c is affinely
            # conjugate to the logistic family, via z=r*(1/2-x).
            row['logistic_r'] = 1+(1-4*c).sqrt()
            row['delta_ratio'] = ((rows[n-1]['c']-rows[n-2]['c']) /
                                  (c-rows[n-1]['c'])) if n >= 2 else None
            row['alpha_magnitude_ratio'] = (
                abs(rows[n-1]['d']/row['d'])) if n >= 2 else None
            row['alpha_signed_ratio'] = (
                rows[n-1]['d']/row['d']) if n >= 2 else None
            row['arithmetic_precision_digits'] = precision
        return rows


def strings(value):
    if isinstance(value, Decimal):
        return str(value)
    if isinstance(value, dict):
        return {k: strings(v) for k, v in value.items()}
    if isinstance(value, list):
        return [strings(v) for v in value]
    return value


def main():
    started = time.time()
    low, high = run(60), run(90)
    with localcontext() as context:
        context.prec = 100
        consistency = []
        for left, right in zip(low, high):
            diffs = {field: abs(left[field]-right[field])
                     for field in ['c', 'd', 'delta_ratio',
                                   'alpha_magnitude_ratio', 'logistic_r']
                     if left.get(field) is not None}
            consistency.append(dict(n=right['n'], period=right['period'],
                                    absolute_precision_rerun_differences=diffs))
    result = dict(
        claim_status='high_precision_numerical_evidence_not_interval_certified',
        family='f_c(z)=z^2+c',
        selection='first sign-changing primitive-period root moving left from prior principal center within 0.05 to 0.4 of prior gap',
        parameter_kind='superstable centers, not bifurcation boundaries',
        indexing='c_0=0 has period 1; c_1=-1 has period 2; c_n has period 2^n',
        delta_definition='(c_(n-1)-c_(n-2))/(c_n-c_(n-1)) for n>=2',
        spatial_definition='d_n=f_(c_n)^(2^(n-1))(0), n>=1; alpha_n=abs(d_(n-1)/d_n), n>=2',
        interpretation='Finite ratios approach the universal constants; the last row is not an estimate with a rigorous error bound.',
        primitive_period_check='Every proper divisor of 2^n divides 2^(n-1); the numerically nonzero halfway return excludes these at the computed precision.',
        precision_rerun_digits=[60, 90],
        consistency=consistency, results=high,
        elapsed_seconds=round(time.time()-started, 3))
    (HERE/'scaling_results.json').write_text(json.dumps(strings(result), indent=2)+'\n')
    columns=['n', 'period', 'c', 'logistic_r', 'd', 'delta_ratio',
             'alpha_magnitude_ratio', 'alpha_signed_ratio', 'residual',
             'root_status', 'arithmetic_precision_digits']
    with (HERE/'scaling_results.csv').open('w', newline='') as handle:
        writer=csv.DictWriter(handle, fieldnames=columns, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(strings(high))
    print(json.dumps(dict(elapsed_seconds=result['elapsed_seconds'],
                          last_row=strings(high[-1]),
                          final_precision_differences=strings(consistency[-1])), indent=2))


if __name__ == '__main__':
    main()
