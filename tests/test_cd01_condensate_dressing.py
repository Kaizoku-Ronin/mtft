"""CD-01: the condensate's u^c-but-not-d^c dressing."""
from mtft.matter import condensate_dressing as C


def test_gauge_selection():
    assert C.dressed() == ["e^c", "nu^c", "u^c"] and C.undressed() == ["d^c"]


def test_required_dressing():
    z1, _ = C.z_required(0.381); assert abs(z1 - 1.0) < 0.01
    zq, (q, d, u) = C.z_required(0.117); assert 0.24 < zq < 0.26 and q < u < d and d / q < 1.5
