"""VS-01: the condensate's feedback on the (4,2,1) scalars."""
from mtft.matter import vacuum_stability as V


def test_feedback_masses():
    m = V.scalar_masses_4L(0.0, 1.0)
    assert abs(m["squark"] + 7 / 13) < 1e-12 and abs(m["slepton"] + 11 / 13) < 1e-12 and not V.stable(0.0, 1.0)
    assert V.stable(V.stability_bound(), 1.0) and not V.stable(V.stability_bound() - 1e-9, 1.0)


def test_ordering_and_asymmetry():
    assert not V.ordering(11, 13) and V.ordering(13, 11) and abs(V.scalar_masses_4L(11.0, 13.0)["slepton"]) < 1e-12
    assert abs(V.radiative_asymmetry(0.4472) - 0.2) < 1e-3
