"""TL-01: the gamma tilt is a rotation."""
from mtft.matter import tilt as T


def test_block_tilts_cancel():
    vus, vub = T.block_tilts(); assert vus < 1e-12 and vub < 1e-12


def test_c_is_the_down_share_of_vcb():
    for args in ((0.25, 0.03, 0.03, 0.02), (0.25, 0.03, 0.016, 0.025), (0.2, -0.04, 0.02, 0.03, 2.0)):
        c, share = T.down_share_check(*args); assert abs(c / share - 1) < 0.05


def test_programme_numbers():
    lo, _ = T.c_range(T.Z_DP01); assert 0.56 < lo < 0.57 and lo > T.WINDOW[1]
    assert 0.83 < T.c_range(1.0)[0] < 0.85 and 0.155 < T.z_for_window() < 0.159
