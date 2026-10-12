"""DP-01: dressing strengths and the quark–lepton conflict."""
from mtft.matter import dressing_profiles as D


def test_weights_and_strengths():
    assert abs(D.weight(0.0) - 1.0) < 1e-12 and abs(D.weight(3.0 ** 0.5) - 0.5) < 1e-12
    m2, m3 = D.masses_for_z(0.25, 1.0); assert abs(m2 - 3.873) < 1e-3 and m3 == 0.0
    m2, m3 = D.masses_for_z(0.25, 0.8); assert abs(m2 - 4.899) < 1e-3 and abs(m3 - 0.75) < 1e-9


def test_conflict():
    assert 3.05 < D.z_leptons() < 3.12 and 12.0 < D.conflict() < 12.7 and 0.058 < D.tail_ratio_if_shared(D.z_leptons()) < 0.063
