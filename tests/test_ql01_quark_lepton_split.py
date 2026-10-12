"""QL-01: the quark–lepton split."""
from mtft.matter import quark_lepton_split as Q


def test_b_tau_bound():
    assert 3.09 < Q.z_e() < 3.10 and 0.32 < Q.dressing_bound() < 0.33
    assert 1.39 < Q.tau_over_b(8.9e11) < 1.41 and 1.52 < Q.tau_over_b(2.1e14) < 1.54
    assert 4.3 < Q.shortfall(8.9e11) < 4.8


def test_targets_and_insertions():
    assert 0.29 < Q.first_over_third() < 0.31
    f1, f2, f3 = Q.factors(8.9e11); assert 0.41 < f1 < 0.43 and 4.3 < f2 < 4.4 and 1.39 < f3 < 1.41
    k1, k2, k3 = Q.insertions(8.9e11); assert k1 < 0 < k3 < 1 < k2
    lo, hi = Q.wilson_line_limit(); assert 2.6 < lo < 3.1 < hi < 3.4
