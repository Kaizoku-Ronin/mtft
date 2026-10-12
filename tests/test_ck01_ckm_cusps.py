"""CK-01: inter-cusp distances, the 1/143 mixing bound, the excluded split assignments, the jet-charge texture."""
import math
import numpy as np
from mtft.matter import ckm_cusps as CK


def test_distances_and_bounds():
    d = CK.cusp_distances()
    assert abs(d[("inf", "0")] - math.log(143)) < 1e-12 and abs(d[("inf", "1/11")] - math.log(143 ** 2 / 13)) < 1e-12 and abs(d[("inf", "1/13")] - math.log(143 ** 2 / 11)) < 1e-12
    assert d[("1/11", "1/13")] == d[("inf", "0")] and d[("0", "1/13")] == d[("inf", "1/11")]       # Atkin–Lehner isometries
    b = CK.mixing_bounds(); assert abs(b[("inf", "0")] - 1 / 143) < 1e-15 and abs(b[("inf", "1/11")] - 13 / 143 ** 2) < 1e-15


def test_only_the_single_point_survives_quark_mixing():
    s = CK.split_assignments()
    assert [k for k, v in s.items() if v["allowed"]] == ["3 (one point)"]
    assert s["1+1+1 (one per cusp)"]["exceed 1/143"] == ["Vus", "Vcb"] and s["2+1 (third family alone)"]["exceed 1/143"] == ["Vcb"]
    assert CK.DATA["Vub"] < 1 / 143 < CK.DATA["Vcb"] < CK.DATA["Vus"]


def test_texture_scaling_and_scan():
    Y = CK.texture(0.1, np.ones((3, 3)))
    assert np.allclose(np.abs(Y), [[1, 0.1, 0.01], [0.1, 0.01, 1e-3], [0.01, 1e-3, 1e-4]])
    eu, ed = CK.fit_eps(CK.DATA["mc_over_mt"]), CK.fit_eps(CK.DATA["ms_over_mb"])
    assert abs(eu - 0.0602) < 1e-3 and abs(ed - 0.1387) < 1e-3
    mc = CK.monte_carlo(eu, ed, n=1500, seed=3)
    assert 1e-5 < mc["mu_over_mt"]["median"] < 5e-5 and 0.1 < mc["Vus"]["median"] < 0.25
    assert mc["Vcb"]["median"] > 2 * CK.DATA["Vcb"] and mc["Vub"]["median"] > 2 * CK.DATA["Vub"]        # the symmetric-texture overprediction, recorded
