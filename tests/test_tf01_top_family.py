"""TF-01: the minimal-breaking theorem, the two-fold anisotropy, the two-plus-one design."""
from mtft.matter import top_family as T


def test_minimal_breaking_units():
    assert [T.min_breaking_units(e, e) for e in (1, 2, 3)] == [2, 3, 4] and T.min_breaking_units(1, 3) == 2 and T.min_breaking_units(3, 1) == 4


def test_two_fold_anisotropy_and_shortfall():
    tf = T.two_fold_anisotropy()
    assert abs(tf["E2hat"] - 0.106) < 0.003 and 0.018 < tf["delta2_edge"] < 0.020 and 0.07 < tf["delta2_across"] < 0.08
    assert T.report()["shortfall"] > 8


def test_design():
    d = T.flux_family_design()
    assert d["Q"][2] == 1 and d["u_c"][2] == 1 and d["higgs_torus_flux"] == 0 and d["jet_breaking_units"] == 3
