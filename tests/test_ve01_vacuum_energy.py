"""VE-01 (2026-10-17): the one-loop vacuum energy of the twisted compactification is minus a square, |Tr U^w|^2; it vanishes exactly at the
Bose–Fermi-balanced Z/2 point, which is the global maximum; the curve corrections to Matter-05's closed forms are per-cent (CC-70)."""
import math

import numpy as np

from mtft.matter import ledger as L
from mtft.matter import puzzle as Z
from mtft.matter import vacuum_energy as V


def test_perfect_square_and_its_zero():
    rng = np.random.default_rng(7)
    for _ in range(3):
        th = rng.random((8, 2))
        assert abs(V.potential(th) - V.block_sum(th)) < 1e-9 and V.potential(th) <= 0
    S = V.lattice_sums()
    assert abs(S["S0"] - 1.589032) < 2e-6 and abs(S["S1"] - 2.107219) < 2e-6                   # PG-6's and Matter-05's sums
    jet = V.z2_theta({"4": 0, "L": 0, "R": 0}); bal = V.balanced_point()
    assert abs(V.potential(jet) + 64 * S["S0"]) < 1e-9 and abs(V.potential(bal)) < 1e-12
    sad = V.z2_theta({"4": 1, "L": 0, "R": 0}, eps=(1, 0))                                    # a half-period that is not the twist point
    assert abs(V.potential(sad) + 32 * S["S0"]) < 1e-9


def test_z2_table():
    rows = V.z2_table()
    assert len(rows) == 8 and all(r["identity"] for r in rows)
    assert sorted(r["TrW_squared"] for r in rows) == [0, 0, 16, 16, 16, 16, 64, 64]
    for r in rows:
        if r["TrW_squared"] == 64: assert r["stability"].startswith("minimum")
        if r["TrW_squared"] == 0: assert r["stability"].startswith("maximum")
        if r["TrW_squared"] == 16: assert r["stability"].startswith("flat")
    assert next(r for r in rows if r["wilson"] == {"4": 1, "L": 0, "R": 0})["bosons_untwisted"] == 32


def test_balanced_point_is_the_top():
    bal = V.balanced_point(); jet = V.z2_theta({"4": 0, "L": 0, "R": 0})
    sb = V.signature(V.hessian(bal)); sj = V.signature(V.hessian(jet))
    assert (sb["negative"], sb["zero"], sb["positive"]) == (2, 14, 0) and (sj["negative"], sj["zero"], sj["positive"]) == (0, 2, 14)
    d = V.t_direction().ravel()
    assert abs(d @ V.hessian(jet) @ d + d @ V.hessian(bal) @ d) < 1e-8                         # PG-4: the two Hessians are negatives
    ev = sorted(set(np.round(sj["eigenvalues"][2:], 3)))
    assert abs(ev[1] / ev[0] - 166.935 / 28.893) < 2e-3                                        # PG-6's theta-coordinate Hessian ratio
    q = V.quartic_profile(V.su4_cartan_direction())
    assert q["sign"] == "maximum" and 3.8 < q["fitted_power"] < 4.1
    assert all(abs(x) < 1e-20 for x in V.quartic_profile(np.array([[0, 1]] * 8, float))["values"])   # the overall U(1) is exactly flat


def test_rigidity_obstructions():
    a = V.t_anomalies()
    assert a["anomalous"] and a["SU2L^2"] == 12 and a["SU2R^2"] == -12 and a["SU4^2"] == 0 and a["grav^2"] == 0
    assert V.dparity_projection()["tachyon_removed"] is False and V.dparity_projection(2.0)["projection_on_T"] == -8.0
    s = V.stueckelberg_criterion()
    assert 1.8 < s["M_St_over_M_KK_light_direction"][0] < s["M_St_over_M_KK_light_direction"][1] < 2.3 and not s["lifts_nonabelian_directions"]


def test_curve_corrections_are_per_cent():
    for r in (230.0, 280.0):
        c = V.curve_corrections(r)
        assert 0.06 < c["expansion_parameter_<t>/R_C^2"] < 0.09 and c["56pi_over_r"] > 0.6
        assert 0.0095 < c["a1_relative_V"] < 0.0125 and all(0.01 < x < 0.018 for x in c["a1_relative_mass2"])
        assert 0 < c["a2_relative_V_estimate"] < 1e-3 and 0 < c["geodesic_relative_V_per_class"] < 1e-5
        assert abs(c["systole"] - 2 * math.acosh(2)) < 1e-12
    assert abs(V.curve_corrections(250.0)["a1_relative_V"] * 250 - 2.78) < 0.01


def test_kk_rung_and_record():
    for r in (230.0, 280.0):
        k = V.curve_kk_rung(r)
        assert 0.06 < k["continuum_threshold_1/(2R_C)"] < 0.071 and 7 < k["M_SS_over_threshold"] < 8
        assert k["first_maass_sqrt(lambda_1)/R_C"] > k["continuum_threshold_1/(2R_C)"]
    for row in V.energy_ledger():
        assert -4.6 < row["V_loop_jet_in_M_KK4"] < -3.6 and row["curvature_over_loop"] > 1e5 and 94 < row["log10_loop_over_rho_Lambda"] < 106
    e = {x["id"]: x for x in L.E}
    for cc, cls in (("VE-1", "EXACT"), ("VE-2", "EXACT"), ("VE-3", "EXACT"), ("VE-4", "EXACT"), ("VE-5", "EXACT"), ("VE-6", "STRUCTURAL"), ("CC-70", "EXACT")):
        assert e[cc]["class"] == cls and e[cc]["test"] == "test_ve01_vacuum_energy.py" and e[cc]["wave"] == ("correction" if cc == "CC-70" else "VE-01")
    assert e["SC-4"]["status"] == "caveat" and "CC-70" in e["SC-4"]["note"]
    sub = Z.subs()["I.6"][1]
    assert sub.status == "open" and "VE-2" in sub.claims and "VE01_VACUUM_ENERGY.md" in sub.docs and "Tr W" in sub.settle
    assert Z.WAVE_DATES["VE-01"] == "2026-10-17" and Z.WAVE_REPORT["VE-01"] == "VE01_VACUUM_ENERGY.md"
    assert all(Z.check().values()) and all(L.check().values())
