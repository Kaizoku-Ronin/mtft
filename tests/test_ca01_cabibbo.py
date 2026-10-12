"""CA-01: the Cabibbo angle and two measurements of the curve's size."""
from mtft.matter import cabibbo as C


def test_two_measurements():
    lo, hi = C.r_for_cabibbo(1); assert 90 < lo < 92 and 890 < hi < 900
    v = C.r_from_vcb(0.1781, 0.5629); assert 18 < v[0] < 19.5 and 27 < v[1] < 29 and C.gap(v, (lo, hi)) > 3.2


def test_fricke_and_pairs():
    f = C.r_from_vcb(0.107, 0.0); assert 78 < f[0] < 81 and 103 < f[1] < 106 and f[1] > C.r_for_cabibbo(1)[0]
    lo_c, mid, data = C.c_bracket(); assert 0.19 < lo_c < 0.21 and lo_c < data < mid
    n, m = C.pairs_and_scale(1); assert n == 5 and 2.5e3 < m < 4.5e3
