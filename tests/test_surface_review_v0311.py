import numpy as np
from fractions import Fraction
from mtft.surface import smflux as SF, arithspin as AS

def test_exact_anomaly_extension_and_majorana_obstruction():
    S = SF.M1_STACKS; N, m = S["N"], S["m"]; phase = {s: 1 for s in N}
    for t in ((1, 0, 0), (0, 1, 0), (0, 0, 1), (2, -3, 5), (1, 1, 1)):
        q = {s: t[0] * Fraction(phase[s]) + t[1] * Fraction(S["Y"][s]) + t[2] * Fraction(S["B-L"][s]) for s in N}
        assert SF.abelian_anomaly_polynomial(N, m, q) == 0
        assert all(v == 0 for v in SF.mixed_nonabelian_anomalies(N, m, q).values())
    assert SF.abelian_anomaly_polynomial(N, m, {"c": 0, "L": 0, "a": 1, "b": 0, "d": 0}) != 0     # a generic direction is anomalous (the c-only direction happens to cancel)
    r = SF.majorana_obstruction(); assert r["nu_c_B_minus_L"] == 1 and not r["bare_majorana_allowed"]

def test_atkin_lehner_lift_on_S0_is_D8():
    L = AS.al_lift_on_S0(); I = np.eye(2)
    assert np.allclose(L["A2"], I) and np.allclose(L["B2"], -I) and np.allclose(L["commutator"], -I) and np.allclose(L["AB2"], I)
    gens = [L["A"], L["B"]]; seen = [I]; frontier = [I]
    while frontier:
        nxt = []
        for g in frontier:
            for h in gens:
                k = g @ h
                if not any(np.allclose(k, s) for s in seen): seen.append(k); nxt.append(k)
        frontier = nxt
    assert len(seen) == 8 and sum(np.allclose(np.linalg.matrix_power(g, 4), I) and not np.allclose(np.linalg.matrix_power(g, 2), I) for g in seen) == 2
