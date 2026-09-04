"""CC-19 (alpha^-1 target -> CODATA 2022) and CC-20 (LN_MONSTER from exact factorisation)."""
import math

import mtft.constants as C
import mtft.falsify as F


def test_cc20_ln_monster_exact():
    order = 1
    for p, e in C.MONSTER_ORDER_FACTORIZATION.items():
        order *= p ** e
    assert order == C.MONSTER_ORDER == 808017424794512875886459904961710757005754368000000000
    assert abs(C.LN_MONSTER - 124.12642336632464) < 1e-13
    assert abs(C.LN_MONSTER - C.LN_MONSTER_TRUNCATED_RETIRED) < 1e-9      # retired value stays exposed
    assert abs(C.LN_MONSTER - math.log(order)) < 1e-12


def test_cc19_alpha_target_is_codata_2022():
    assert C.ALPHA_INV_CODATA2022 == 137.035999177 and C.ALPHA_INV_CODATA2022_ERR == 2.1e-8
    assert C.ALPHA_INV_CODATA2018_RETIRED == 137.035999084
    preds = [p for p in F.predictions() if "α⁻¹" in p.relation] if hasattr(F, "predictions") else []
    for p in preds:
        assert p.observed == C.ALPHA_INV_CODATA2022


def test_cc19_status_bookkeeping():
    g = C._MTFTGauge()
    s3 = abs(g.alpha_inv - C.ALPHA_INV_CODATA2022) / C.ALPHA_INV_CODATA2022_ERR
    s4 = abs(g.alpha_inv_4term - C.ALPHA_INV_CODATA2022) / C.ALPHA_INV_CODATA2022_ERR
    assert 7.0 < s3 < 8.5                       # 3-term: FAIL under both targets (12.1 -> 7.7 sigma)
    assert s4 < 1.0                             # 4-term: 3.9 sigma (2018) -> 0.55 sigma (2022); NOT registered
    sm = abs(g.alpha_inv_monster - C.ALPHA_INV_CODATA2022) / C.ALPHA_INV_CODATA2022_ERR
    assert sm > 2e4                             # Monster expression rejected (Wave M1)
