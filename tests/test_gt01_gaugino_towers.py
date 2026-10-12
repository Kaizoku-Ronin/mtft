"""GT-01: the two gaugino towers."""
from mtft.matter import gaugino_towers as G


def test_where_the_wr_tower_lives():
    d = G._setup()
    assert G.share(G.tau_at(0.17), 0, d) == 1.0
    assert G.share(G.tau_at(0.17), 2, d) < 0.003 and 0.05 < G.share(G.tau_at(0.39), 2, d) < 0.06
    assert 0.70 < G.share(complex(0.18363, 0.03346), 2, d) < 0.74 and 0.19 < G.share(complex(0, 143 ** -0.5), 2, d) < 0.21


def test_strength_bound():
    Z3u, m2 = G.reversal_requirements(); assert 0.32 < Z3u < 0.33 and 12 < m2 < 13
    assert 3.9 < G.mass_ratio_needed(1.0) < 4.0 and 4.5 < G.mass_ratio_needed(0.72) < 4.8 and G.lepton_bound(0.8) == 1.25
