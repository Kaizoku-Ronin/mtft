"""Exact bounds for active blocks of finite Hermitian dyadic matrices.

No floating arithmetic is used after matrix entries are encoded. The JSON
matrix format is dependency-free: real and imaginary integer arrays / 2**exponent.
This is research support code, not a claim about exact MTFT period geometry.
"""
from fractions import Fraction as F
from math import isqrt, factorial
import re

PRECISION = 10**30


def rational(value):
    return F(int(value['numerator']), int(value['denominator']))


def record(value):
    return {'numerator': str(value.numerator), 'denominator': str(value.denominator)}


def sqrt_lower(x):
    if x < 0:
        raise ValueError('Square root requires nonnegative input')
    return F(isqrt(x.numerator * PRECISION**2 // x.denominator), PRECISION)


def sqrt_upper(x):
    lo = sqrt_lower(x)
    return lo if lo*lo == x else lo+F(1, PRECISION)


def encode(a):
    """Encode a numpy-like complex array exactly, including signed dyadics."""
    entries = [[[float(z.real).as_integer_ratio(), float(z.imag).as_integer_ratio()]
                for z in row] for row in a]
    exponent = max(d.bit_length()-1 for row in entries for z in row for _, d in z)
    scale = 1 << exponent
    return {'exponent': exponent,
            'real': [[str(z[0][0]*(scale//z[0][1])) for z in row] for row in entries],
            'imag': [[str(z[1][0]*(scale//z[1][1])) for z in row] for row in entries]}


def decode(encoded):
    def integer(x):
        if type(x) is int:
            return x
        if isinstance(x, str) and re.fullmatch(r'-?[0-9]+', x):
            return int(x)
        raise ValueError('Encoded entries must be decimal integer strings or integers')
    real = [[integer(x) for x in row] for row in encoded['real']]
    im = [[integer(x) for x in row] for row in encoded['imag']]
    n = len(real)
    if n == 0 or len(im) != n or any(len(row) != n for row in real+im):
        raise ValueError('Expected nonempty square complex matrix')
    if type(encoded['exponent']) is not int or encoded['exponent'] < 0:
        raise ValueError('Invalid dyadic exponent')
    if any(real[i][j] != real[j][i] or im[i][j] != -im[j][i]
           for i in range(n) for j in range(n)):
        raise ValueError('Matrix must be exactly Hermitian')
    return real, im


def multiply(a, b):
    ar, ai = a; br, bi = b; n = len(ar)
    rr = [[0]*n for _ in range(n)]; ii = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            rr[i][j] = sum(ar[i][k]*br[k][j]-ai[i][k]*bi[k][j] for k in range(n))
            ii[i][j] = sum(ar[i][k]*bi[k][j]+ai[i][k]*br[k][j] for k in range(n))
    return rr, ii


def norm_squared(a, denominator, rows=None, columns=None):
    re, im = a; n = len(re)
    return F(sum(re[i][j]**2+im[i][j]**2
                 for i in range(n if rows is None else rows)
                 for j in range(n if columns is None else columns)), denominator**2)


def maximum_row_bound(a, exponent):
    re, im = a; scale = 1 << exponent
    return max(sum(sqrt_upper(F(x*x+y*y, scale*scale)) for x, y in zip(r, s))
               for r, s in zip(re, im))


def moment_data(full, reduced, active, maximum_order):
    if not isinstance(active, int) or not 0 < active <= min(len(full['real']), len(reduced['real'])):
        raise ValueError('Invalid active dimension')
    if not isinstance(maximum_order, int) or maximum_order < 0:
        raise ValueError('Invalid Taylor order')
    matrices = [decode(full), decode(reduced)]
    exponents = [full['exponent'], reduced['exponent']]
    sizes = [len(a[0]) for a in matrices]
    powers = [([[int(i == j) for j in range(n)] for i in range(n)],
               [[0]*n for _ in range(n)]) for n in sizes]
    difference_norms = []; column_norms = [[], []]
    for k in range(maximum_order+2):
        exponent = max(exponents)*k
        factors = [1 << (exponent-e*k) for e in exponents]
        diff = tuple([[powers[0][part][i][j]*factors[0]-powers[1][part][i][j]*factors[1]
                       for j in range(active)] for i in range(active)] for part in [0, 1])
        difference_norms.append(sqrt_upper(norm_squared(diff, 1 << exponent)))
        for j in range(2):
            column_norms[j].append(sqrt_upper(norm_squared(powers[j], 1 << (exponents[j]*k), columns=active)))
        if k < maximum_order+1:
            powers = [multiply(powers[j], matrices[j]) for j in range(2)]
    rows = [maximum_row_bound(matrices[j], exponents[j]) for j in range(2)]
    return {'differences': difference_norms, 'columns': column_norms, 'rows': rows}


def certify(full, reduced, active, horizons, orders, budget):
    if any(t < 0 for t in horizons) or not orders or min(orders) < 0 or budget < 0:
        raise ValueError('Invalid horizon, order, or budget')
    deficit = 2*active-len(full['real'])
    if deficit <= 0:
        raise ValueError('Global denominator requires active > complement')
    data = moment_data(full, reduced, active, max(orders))
    denominator = sqrt_lower(F(deficit))
    rows = []
    for t in horizons:
        for m in orders:
            finite = sum(t**k*data['differences'][k]/factorial(k) for k in range(m+1))
            tail_columns = t**(m+1)*sum(c[m+1] for c in data['columns'])/factorial(m+1)
            tail_rows = sqrt_upper(F(active))*t**(m+1)*sum(r**(m+1) for r in data['rows'])/factorial(m+1)
            bc = (finite+tail_columns)/denominator
            br = (finite+tail_rows)/denominator
            rows.append({'horizon': record(t), 'order': m,
                         'finite_sum_upper': record(finite), 'column_tail_upper': record(tail_columns),
                         'row_tail_upper': record(tail_rows), 'relative_bound_columns': record(bc),
                         'relative_bound_rows': record(br), 'relative_bound': record(min(bc, br)),
                         'certifies_budget': min(bc, br) <= budget})
    eligible = [rational(r['horizon']) for r in rows if r['certifies_budget']]
    return {'arithmetic': 'exact Gaussian integers and rational outward square-root enclosures',
            'active_dimension': active, 'full_dimension': len(full['real']),
            'reduced_dimension': len(reduced['real']), 'sqrt_enclosure_denominator': str(PRECISION),
            'denominator_lower': record(denominator), 'budget': record(budget),
            'difference_norm_uppers': [record(x) for x in data['differences']],
            'row_norm_uppers': [record(x) for x in data['rows']],
            'rows': rows, 'largest_certified_listed_horizon': record(max(eligible)) if eligible else None}
