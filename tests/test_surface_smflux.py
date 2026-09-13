"""SM-01: four-stack no-go, five-stack solutions, ledger of the representative model, purity of its chiral sectors."""
from fractions import Fraction as Fr

from mtft.surface import smflux as SF

YS = [Fr(k, 6) for k in range(-6, 7)]


def test_four_stack_no_go_and_five_stack_solutions():
    assert SF.search(["c", "L", "a", "b"], 8, {"c": [Fr(1, 6)], "L": [Fr(0), Fr(1, 3)], "a": YS, "b": YS}) == []
    s5 = SF.search(["c", "L", "a", "b", "d"], 6, {"c": [Fr(1, 6)], "L": [Fr(0), Fr(1, 3)], "a": YS, "b": YS, "d": YS})
    assert len(s5) == 30 and all(s["degrees"]["L"] == -3 for s in s5) and min(s["vector_like_pairs"] for s in s5) == 6
    assert all(s["y"]["a"] == s["y"]["b"] or s["y"]["a"] == s["y"]["d"] or s["y"]["b"] == s["y"]["d"] for s in s5)


def test_representative_ledger_and_purity():
    M1 = {"degrees": {"c": 0, "L": -3, "a": 3, "b": 3, "d": 0}, "y": {"c": Fr(1, 6), "L": Fr(0), "a": Fr(-1, 2), "b": Fr(1, 2), "d": Fr(-1, 2)}}
    led = SF.anomaly_ledger(M1)
    assert led["ledger"]["Y"] == {"SU3^2": "0", "SU2^2": "0", "grav": "0"} and led["Y_in_anomaly_free_span"] and len(led["anomaly_free_combinations"]) == 3
    rows = SF.rep_table(M1); mult = {(r["rep"][0], r["rep"][1], r["rep"][2]): r["multiplicity"] for r in rows}
    assert mult[("3", "2", Fr(1, 6))] == 3 and mult[("3bar", "1", Fr(-2, 3))] == 3 and mult[("1", "1", Fr(1))] == 3
    assert SF.purity_with_S0(3, "cm")["pure"] and SF.purity_with_S0(-6, "cm2")["pure"] and not SF.purity_with_S0(0, "trivial")["pure"]
