"""TH-01: thresholds and the scale tension."""
from fractions import Fraction as F
from mtft.matter import thresholds as T


def test_left_right_pairs_cancel():
    assert T.pair_b("421") == (F(4, 3), F(8, 3), F(8, 9)) and T.pair_b("412") == (F(4, 3), F(0), F(32, 9))
    s1, s2 = T.slopes("421"), T.slopes("412"); assert s1 == (F(4, 3), F(-4, 3)) and tuple(a + b for a, b in zip(s1, s2)) == (0, 0)


def test_requirements():
    N, D = T.N_D(); assert abs(N - 56.51) < 0.05 and abs(D - 15.36) < 0.05
    assert abs(T.required_nl(7 / 3) - 29.2) < 0.2 and abs(T.required_nl(1.0) - 96.9) < 0.3 and T.adjoint_requirement(7 / 3) > 190


def test_counts_and_dkl():
    assert [T.vector_like_pairs(c) for c in (0, 1, 2)] == [1, 5, 17] and abs(T.dkl_factor() - 0.3527) < 1e-3
