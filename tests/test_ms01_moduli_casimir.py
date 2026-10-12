"""MS-01: the Casimir energy and the curve's size."""
from mtft.matter import moduli_casimir as M


def test_lattice_sums_match_matter05():
    S0, S1 = M.odd_sums("1/2"); assert abs(S1 - 2.107219) < 1e-5 and abs(S1 / S0 - 1.3261) < 1e-3
    P0, P1 = M.odd_sums("tau/2"); assert abs(P1 - 3.238510) < 1e-5 and abs(P1 / P0 - 1.0810) < 1e-3


def test_minimum_condition_and_sign():
    assert abs(M.minimum_ratio(-4.0, "1/2") - 11.11) < 0.02 and abs(M.minimum_ratio(-4.0, "tau/2") - 9.06) < 0.02
    a = M.a1_counts(); assert a["bosons"] == -8.0 and a["fermions"] == 16.0 and a["vector_multiplet"] == 8.0
    assert M.minimum_ratio(a["vector_multiplet"] / 16) is None


def test_full_twisted_spectrum():
    for d in (-6, 0, 11, 13, 24):
        a0, a1, s = M.dolbeault_a1(-24, d); assert abs(s - 4.0) < 1e-12 and abs(a0 - (-4 + d / 2)) < 1e-12
    a1, N = M.block_a1(-24); assert a1 == 32.0 and N == 16 and a1 / N == 2.0
    assert not M.minimum_possible(-24) and M.minimum_possible(2)
    assert M.total_curvature(-24, 1.0, (0.225, 0.55)) < M.total_curvature(-24, 1.0, ())
