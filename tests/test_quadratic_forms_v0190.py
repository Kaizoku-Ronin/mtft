from mtft import quadratic_forms as QF


def test_gate_all_layers():
    out = QF.gate()
    assert out["layers"] == "L1-L7 PASS"
    assert out["eureka_to"] == 20000
