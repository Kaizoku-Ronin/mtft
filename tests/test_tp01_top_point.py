"""TP-01: flux families come in 24s; the constant-jet top; the two-point design; the unique Higgs."""
from mtft.matter import top_point as TP


def test_flux_families_come_in_twenty_fours():
    assert TP.flux_family_chirality(11, 1) == 24 and TP.flux_family_chirality(13, -1) == -24 and TP.flux_family_chirality(0, 1) == 24


def test_two_point_design():
    f = TP.two_point_families(); assert f["P1"] == 2 and f["P2"] == 1 and f["higgs_exceptional"] == (4, 2)
    assert abs(TP.global_cost(1 / 3) - 0.424) < 0.002 and abs(TP.global_cost(0.55) - 0.700) < 0.002 and abs(TP.global_cost(1 / 3, False) - 0.106) < 0.002
    assert TP.d_flat_relation(1, 3.0, 1.0, 1.0) == 0.0


def test_unique_higgs():
    assert TP.higgs_count(24) == 1 and TP.higgs_count(0) == 13 and TP.higgs_count(30) == 0 and TP.higgs_count(12) == 0
