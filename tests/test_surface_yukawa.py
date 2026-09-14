"""YUK-01 / COND-01 gates, GP-free from the frozen weight-2 basis."""
import numpy as np

from mtft.surface import yukawa as YK, condensation as CD


def test_yukawa_tensor_gates_and_condensation():
    Y = YK.yukawa_tensor()
    g = Y["gates"]
    assert g["cubic_rank"] == 60 and g["quartic_rank"] == 84
    assert g["cubic_sectors"] == [12, 18, 17, 13] and g["quartic_sectors"] == [24, 18, 19, 23]
    assert g["selection_rule_violations"] == 0 and g["nonzero"] == 7543
    rng = np.random.default_rng(0)
    c = CD.extension_cohomology(Y, rng.standard_normal(84))
    assert (c["rank_delta"], c["h0_V"], c["h1_V"], c["index"]) == (13, 48, 0, 48)
    e = np.array([rng.standard_normal() if tuple(s) == (1, 1) else 0.0 for s in Y["quartic_sectors"]])
    c2 = CD.extension_cohomology(Y, e)
    assert [c2["sectors"][s]["rank"] for s in ((1, 1), (1, -1), (-1, 1), (-1, -1))] == [1, 6, 5, 1]
    assert CD.tachyon_mass2(-72) == -3.0 and abs(CD.condensation_energy(-72)["split"] - 216 * np.pi) < 1e-9
