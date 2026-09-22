#!/usr/bin/env python3
"""Exact finite checks and high-precision diagnostics for the research note.

Python standard library only. Global identities are proved in the report;
finite arithmetic checks below supplement those proofs. Decimal stability
is a numerical diagnostic, not an interval certificate.
"""
from decimal import Decimal as D, localcontext
from fractions import Fraction as F
from pathlib import Path
import csv
import json
import math

OUT = Path(__file__).resolve().parent
PI = D('3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117067982148086513282306647')


def trim(p):
    p = list(p)
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return p


def mul(p, q):
    r = [0] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        for j, b in enumerate(q):
            r[i+j] += a*b
    return trim(r)


def iterate_poly(m, c):
    p = [0, 1]
    for _ in range(m):
        p = mul(p, p)
        p[0] += c
    return trim(p)


def bareiss(a):
    a = [row[:] for row in a]
    n = len(a)
    previous, sign = 1, 1
    for k in range(n-1):
        if a[k][k] == 0:
            swap = next((i for i in range(k+1, n) if a[i][k]), None)
            if swap is None:
                return 0
            a[k], a[swap] = a[swap], a[k]
            sign = -sign
        pivot = a[k][k]
        for i in range(k+1, n):
            for j in range(k+1, n):
                numerator = pivot*a[i][j] - a[i][k]*a[k][j]
                assert numerator % previous == 0
                a[i][j] = numerator // previous
        for i in range(k+1, n):
            a[i][k] = 0
        previous = pivot
    return sign*a[-1][-1]


def resultant(p, q):
    p, q = trim(p), trim(q)
    m, n = len(p)-1, len(q)-1
    rows = []
    for i in range(n):
        rows.append([0]*i + p[::-1] + [0]*(n-1-i))
    for i in range(m):
        rows.append([0]*i + q[::-1] + [0]*(m-1-i))
    return bareiss(rows)


def discriminant(p):
    degree = len(p)-1
    assert p[-1] == 1
    dp = [i*p[i] for i in range(1, len(p))]
    return (-1)**(degree*(degree-1)//2)*resultant(p, dp)


def critical_values(m, c):
    p, result = 0, []
    for _ in range(m):
        p = p*p+c
        result.append(p)
    return result


def predicted_disc(m, c):
    degree = 2**m
    value = (-1)**(degree*(degree-1)//2)*2**(m*degree)
    for j, p in enumerate(critical_values(m, c), 1):
        value *= p**(2**(m-j))
    return value


def agm(a, b, precision):
    for _ in range(100):
        if abs(a-b) < D(10)**(-precision+8):
            return (a+b)/2
        a, b = (a+b)/2, (a*b).sqrt()
    raise ArithmeticError('AGM failed to converge')


def period_table(precision):
    rows = []
    with localcontext() as context:
        context.prec = precision
        previous = None
        target = 1/(2*D(2).sqrt())
        for exponent in [2, 4, 8, 12, 20, 30]:
            eps = D(10)**(-exponent)
            s = (1+eps).sqrt()
            a = (s*(s+1)).sqrt()
            # Stable identity s-1=eps/(s+1) avoids cancellation.
            b = (s*eps/(s+1)).sqrt()
            modulus_complement = b/a
            elliptic_k = PI/(2*agm(D(1), modulus_complement, precision))
            integral = elliptic_k/a
            log_inverse_eps = -eps.ln()
            asymptote = log_inverse_eps/(2*D(2).sqrt()) + D(8).ln()/D(2).sqrt()
            slope = None if previous is None else (integral-previous[1])/(log_inverse_eps-previous[0])
            rows.append(dict(epsilon=str(eps), log_inverse_epsilon=str(log_inverse_eps),
                             period_integral=str(integral), leading_asymptote=str(asymptote),
                             asymptote_error=str(integral-asymptote),
                             secant_slope=None if slope is None else str(slope),
                             predicted_slope=str(target)))
            previous = log_inverse_eps, integral
    return rows


def main():
    disc_rows = []
    for m in range(1, 5):
        for c in [-3, -2, -1, 0, 1, 2, 3]:
            measured = discriminant(iterate_poly(m, c))
            expected = predicted_disc(m, c)
            assert measured == expected, (m, c, measured, expected)
            disc_rows.append(dict(iteration_depth=m, parameter=c, discriminant=str(measured)))

    # p_m(c) is monic and p_m'(c)=1 modulo 2: squarefree in characteristic 0.
    p = [0]
    for m in range(1, 9):
        p = mul(p, p)
        if len(p) < 2:
            p += [0]*(2-len(p))
        p[1] += 1
        derivative_mod2 = trim([(i*p[i]) % 2 for i in range(1, len(p))])
        assert derivative_mod2 == [1]

    hodge_rows = []
    for m in range(2, 8):
        genus, old_genus = 2**(m-1)-1, 2**(m-2)-1
        plus = [j for j in range(genus) if (-1)**(j+1) == 1]
        minus = [j for j in range(genus) if (-1)**(j+1) == -1]
        assert len(plus) == old_genus
        assert len(minus) == old_genus+1
        # Pullbacks 2*x*(x^2+c)^k have distinct odd leading powers.
        assert [2*k+1 for k in range(old_genus)] == plus
        hodge_rows.append(dict(depth=m, genus=genus, inherited_h10=len(plus), new_h10=len(minus),
                               h1_total=2*genus))

    modular_rows = []
    for c in map(F, [-4, -3, -2, 1, 2, 3]):
        t = -64*(c+1)
        j = 64*(4*c+3)**3/(c+1)
        j2 = 64*(3-c)**3/(c+1)**2
        assert j == (t+16)**3/t
        assert j2 == (t+256)**3/t**2
        fricke = -c/(c+1)
        assert -fricke/(fricke+1) == c
        assert 64*(4*fricke+3)**3/(fricke+1) == j2
        # Classical level-2 modular polynomial, evaluated with exact fractions.
        phi2 = (j**3+j2**3-j*j*j2*j2+1488*j*j2*(j+j2)
                -162000*(j*j+j2*j2)+40773375*j*j2
                +8748000000*(j+j2)-157464000000000)
        assert phi2 == 0
        # Quotient (u,v)=(x^2,x*y) satisfies v^2=u*(u^2+2*c*u+c^2+c).
        for x in map(F, [-3, -1, 0, 2]):
            y_squared = (x*x+c)**2+c
            u = x*x
            assert x*x*y_squared == u*(u*u+2*c*u+c*c+c)
        modular_rows.append(dict(c=str(c), t=str(t), j=str(j), j_2isogenous=str(j2), fricke_c=str(fricke)))

    low, high = period_table(60), period_table(90)
    max_period_difference = max(abs(D(a['period_integral'])-D(b['period_integral'])) for a,b in zip(low,high))
    assert max_period_difference < D('1e-50')
    result = dict(exact_finite_checks=dict(discriminants=disc_rows,
                                          critical_polynomial_squarefree_mod2_depths=list(range(1,9)),
                                          hodge_decompositions=hodge_rows, modular_2isogenies=modular_rows),
                  periods=dict(classification='NUMERICAL DIAGNOSTIC, not interval-certified',
                               precisions=[60,90], max_precision_rerun_difference=str(max_period_difference),
                               normalization='I(c)=integral_b^a dx/sqrt((a^2-x^2)(x^2-b^2)); half a closed period up to phase',
                               rows=high))
    (OUT/'geometry_results.json').write_text(json.dumps(result,indent=2)+'\n')
    with (OUT/'period_results.csv').open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=high[0].keys())
        writer.writeheader()
        writer.writerows(high)
    print(json.dumps(dict(exact_discriminant_comparisons=len(disc_rows),
                          hodge_depths=len(hodge_rows), exact_modular_examples=len(modular_rows),
                          period_precision_rerun_difference=str(max_period_difference),
                          final_period_slope=high[-1]['secant_slope'],
                          predicted_period_slope=high[-1]['predicted_slope']),indent=2))


if __name__ == '__main__':
    main()
