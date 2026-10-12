"""CI-01: the condensate insertion."""
import math
from mtft.matter import condensate_insertion as C


def test_entry_pattern_and_tail_theorem():
    assert 1.75 < C.r_B() < 1.77 and 0.95 < C.r_A() < 0.97 and 0.91 < C.product_relation() < 0.93
    assert C.tail_shortfall(8.9e11) > 2.45 and C.tail_shortfall(2.1e14) > 2.6 and 1 / C.r_B() < 0.57


def test_lever_arm_bound():
    assert 7.6 < C.lever_needed(1.401) < 7.8 and 5.8 < C.lever_needed(1.527) < 5.95
    assert C.lever_needed(1.0) == math.inf and 1.04 < C.threshold_f3() < 1.05
    assert 2.3 < C.lever_arm(0.95) < 2.4 and 2.95 < C.lever_arm(1.21) < 3.05
    s, lam = 0.4, 3.0                                   # brute-force check of the AM–GM optimum
    best = max(s * (math.sin(t) ** 2 + lam / (2 * math.sqrt(2)) * math.sin(2 * t)) for t in [k * math.pi / 2000 for k in range(2001)])
    assert abs(best - C.max_excess(s, lam)) < 1e-4
