"""mtft.surface.bimodule + frozen U/block extension (v0.26.0). No GP needed."""
import numpy as np

from mtft.surface import bimodule as BM, frozen as FR


def test_frozen_extension_gates_and_block_invariants():
    d = FR.x0143()
    assert all(d["gates"].values()) and len(d["gates"]) >= 19
    inv = FR.block_invariants()
    assert inv["ell"]["polarization_elementary_divisors"] == [4, 4]
    assert inv["ghost"]["polarization_elementary_divisors"] == [2, 2, 18, 18]
    assert inv["q4"]["polarization_elementary_divisors"] == [2] * 6 + [18, 18]
    assert inv["q6"]["polarization_elementary_divisors"] == [1, 1] + [2] * 8 + [4, 4]
    for k, v in inv.items():                       # covolume^2 = polarization det exactly
        assert abs(v["hodge_covolume"] ** 2 - v["polarization_det"]) < 1e-8 * v["polarization_det"]


def test_af09_census_exact_pattern():
    c = BM.x0143_census()
    status = {r["twist"]: r["order_zero_status"] for r in c["rows"]}
    assert status == {"untwisted": "FAIL", "W11": "FAIL", "W13": "PASS", "W143": "PASS"}
    assert all(r["max_one_form_size"] < 1e-10 and r["one_form_dimension"] == 0 for r in c["rows"])
    assert c["sectors_V"] == {"(+,+)": 2, "(+,-)": 12, "(-,+)": 10, "(-,-)": 2}
    assert c["sectors_VxV"] == {"(+,+)": 252, "(+,-)": 88, "(-,+)": 88, "(-,-)": 248}
    assert c["AL_adjoint_identity_U13"] < 1e-10 and c["AL_adjoint_identity_U11"] < 1e-10


def test_random_alphabet_control():
    d = FR.x0143()
    Jint = d["intersection_cycles"].astype(float)
    G = Jint @ d["J_true"]; G = (G + G.T) / 2
    G = -G if np.linalg.eigvalsh(G)[0] < 0 else G
    rows = BM.census(BM.random_alphabet(26, 2, True, 1), G, {"untwisted": None})
    assert rows[0]["order_zero_status"] == "FAIL"          # non-commutative alphabet
    assert rows[0]["one_form_dimension"] == 0
