"""Exact rational specialization checks accompanying DYNAMICS_AND_FRICKE.md.

The note gives symbolic proofs. These independent numerical specializations use
Fraction arithmetic and do not purport to prove identities by finite sampling.
"""
from fractions import Fraction as Q
from pathlib import Path
import json


def mm(a, b):
    return [[sum(a[i][k] * b[k][j] for k in range(2))
             for j in range(2)] for i in range(2)]


def add(a, b):
    return [[a[i][j] + b[i][j] for j in range(2)] for i in range(2)]


def scale(s, a):
    return [[s * x for x in row] for row in a]


def w(c):
    return -c / (c + 1)


def A(c):
    return [[-1/(4*c), 1/(4*c*(c+1))], [Q(-1, 4), 1/(4*c)]]


def D(c):
    return [[Q(1), Q(0)], [c/(c+1), 1/(c+1)]]


def dD(c):
    return [[Q(0), Q(0)], [1/(c+1)**2, -1/(c+1)**2]]


def main():
    values = [Q(-5), Q(-3), Q(-2), Q(-3, 2), Q(-1, 2),
              Q(1, 3), Q(1), Q(2), Q(4), Q(7)]
    connection_checks = isogeny_checks = roundtrip_checks = 0
    identity = [[Q(1), Q(0)], [Q(0), Q(1)]]
    for c in values:
        cp = w(c)
        assert w(cp) == c
        assert mm(D(cp), D(c)) == identity
        roundtrip_checks += 1
        lhs = add(add(dD(c), scale(1/(2*(c+1)), D(c))), mm(D(c), A(c)))
        rhs = scale(-1/(c+1)**2, mm(A(cp), D(c)))
        explicit = [[(c-1)/(4*c*(c+1)), 1/(4*c*(c+1))],
                    [1/(2*(c+1)**2), 1/(4*c*(c+1)**2)]]
        assert lhs == rhs == explicit
        connection_checks += 1
        s2 = 2*(c+1)
        for x in [Q(-4), Q(-1), Q(1, 2), Q(1), Q(3), Q(5)]:
            y2 = x**3 - 4*c*x**2 - 4*c*x
            xp = (x-4*c-4*c/x)/s2
            yp2 = y2*(1+4*c/x**2)**2/s2**3
            assert yp2 == xp**3-4*cp*xp**2-4*cp*xp
            isogeny_checks += 1
    t = [[Q(1), Q(1)], [Q(0), Q(1)]]
    u = [[Q(1), Q(0)], [Q(-2), Q(1)]]
    uinv = [[Q(1), Q(0)], [Q(2), Q(1)]]
    tu = mm(t, u)
    assert mm(tu, tu) == scale(Q(-1), identity)
    h = mm(t, uinv)
    assert h == [[Q(3), Q(1)], [Q(2), Q(1)]]
    assert h[0][0]*h[1][1]-h[0][1]*h[1][0] == 1
    result = {
        "method": "exact Fraction specializations; symbolic proofs in note",
        "fricke_involution_and_normalized_roundtrip_checks": roundtrip_checks,
        "connection_intertwining_checks": connection_checks,
        "isogeny_equation_checks": isogeny_checks,
        "modular_word_checks": "(TU)^2=-I; TU^-1=[[3,1],[2,1]]; determinant=1",
        "all_passed": True,
    }
    dest = Path(__file__).with_name("fricke_checks.json")
    dest.write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
