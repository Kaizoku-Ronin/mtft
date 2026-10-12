"""RT-01: the ratio-3 background — leading ratio, first correction, jet-sector neutrinos."""
from mtft.matter import ratio_three as R


def test_leading_and_corrected_ratio():
    assert R.leading_ratio() == 3
    w = R.corrected_ratio_window()
    for v in w.values():
        assert 2.8 < v["ratio_at_r_u"] < 3.0 and 2.8 < v["ratio_at_r_d"] < 3.0
        assert v["ratio_at_r_u"] > 1.2 * 2.3049                                   # the 25 per cent excess over the data
    rs = sorted(v["r_from_eps_u"] for v in w.values()); assert 4000 < rs[0] < 5500 and 40000 < rs[1] < 60000


def test_jet_sector_neutrinos_are_type_one_at_vR():
    ms = [R.jet_neutrino_estimate(t) for t in (0.1, 0.3, 1.0, 3.0)]
    assert max(ms) / min(ms) < 1.01                                               # independent of the shared D-term suppression
    assert 0.01 < ms[0] < 0.05                                                    # GeV: tens of MeV, ~10^9 above 0.05 eV
