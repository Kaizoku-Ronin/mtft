"""PG-01: the weights as the divisor-gas energy at the Hagedorn wall, the two energies, the Glaisher identities, the clock as zeta(-1); the
critical-point atlas of the vacuum, Hubble friction, the toy run and the simulation gate."""
import contextlib
import io
import json
import math
from fractions import Fraction as Fr

import numpy as np
import pytest

from mtft import divisor_gas as D
from mtft.matter import ledger as L
from mtft.matter import puzzle as Z
from mtft.matter import vacuum_atlas as VA
from mtft.matter.__main__ import main as matter_main
from mtft.legend import REGISTRY


# ------------------------------------------------------------------------------------------------------------------------------ the arithmetic clock
def test_weight_is_the_divisor_gas_energy():
    from mtft.marked_gas import weights_sieve
    w = weights_sieve(300)
    for n in (2, 12, 143, 144, 255, 300):
        assert abs(D.divisor_energy(n) - w[n]) < 1e-12
    r = D.weight_is_energy(120)
    assert r["agree"] and r["max_abs_difference"] < 1e-9
    g = D.divisor_gas(143)
    assert g["divisors"] == 4 and abs(g["energy_numerator"] - w[143]) < 1e-12 and g["is_weight"] and 0 < g["mean_energy"] < g["energy_numerator"]
    assert abs(D.divisor_gas(144)["energy_numerator"] - 2.809505) < 1e-5
    assert abs(D.divisor_partition(6, 1.0) - (1 + 1 / 2 + 1 / 3 + 1 / 6)) < 1e-15


def test_two_energies_and_the_wall():
    e, u = D.dirichlet_ensemble(2.0), D.unmarked_gas(2.0)
    assert abs(e["Z"] - 0.3259046) < 1e-6 and abs(e["mean_log_n"] - 1.780032) < 1e-5 and e["variance_log_n"] > 0 and e["heat_capacity"] > 0
    assert abs(u["mean_energy"] - 0.569961) < 1e-5 and u["rh_sensitive"] and e["rh_sensitive"]       # CC-65: both energies see the zeros
    vm = D.von_mangoldt_check(2.0, 5000)
    assert abs(vm["difference"]) < 1e-3 and vm["difference"] > 0                                  # the prime sum approaches -zeta'/zeta from below
    w = D.hagedorn_wall()
    assert w["beta_wall"] == 1.0 and w["residue_equals_minus_zeta_prime_2"]
    for b in (1.05, 1.5, 3.0):
        assert D.dirichlet_ensemble(b)["mean_log_n"] > D.unmarked_gas(b)["mean_energy"]           # the marking raises the energy
        dc = D.marked_energy_decomposition(b)
        assert dc["agree"] and abs(dc["U_D"] - D.dirichlet_ensemble(b)["mean_log_n"]) < 1e-9 and abs(dc["U"] - D.unmarked_gas(b)["mean_energy"]) < 1e-9
    rz = D.residue_at_first_zero()                                                                 # the marked energy keeps the pole at rho_1, residue -1
    assert rz["is_minus_one"] and abs(rz["rho_1"].imag - 14.134725) < 1e-5


def test_glaisher_identities():
    m = D.mean_weight_glaisher()
    assert m["agree"] and abs(m["minus_zeta_prime_2"] - 0.9375482543) < 1e-9
    d1, d2 = D.partial_sum_density(2000), D.partial_sum_density(20000)
    assert abs(d2["density"] - d2["limit"]) < 3e-3 < abs(d1["density"] - d1["limit"])                 # the error shrinks like log^2 x / x
    for d in (d1, d2):
        assert d["error"] < 0 and abs(d["error"] - d["predicted_error"]) < 0.05 * abs(d["predicted_error"])  # -log^2 x/(4x) - (ln 2pi) log x/(2x)
    nv = D.alphabet_negative_values(7)
    assert all(nv["checks"].values()) and all(r["agree"] for r in nv["rows"])
    assert abs(nv["rows"][0]["W"] + math.log(2 * math.pi) / 24) < 1e-15 and nv["rows"][1]["W"] == 0 and nv["rows"][3]["W"] == 0
    assert abs(float(D.bendersky_constant_log(0)) - 0.5 * math.log(2 * math.pi)) < 1e-15                # A_0 = sqrt(2 pi)
    import mpmath as mp
    assert abs(float(D.bendersky_constant_log(1)) - float(mp.log(D.GLAISHER_A))) < 1e-15               # A_1 = A


def test_euler_characteristic_is_zeta_minus_one_and_bost_connes():
    c = D.euler_characteristic_is_zeta_minus_one(143)
    assert c["equal"] and c["chi_SL2Z"] == Fr(-1, 12) and c["chi_PSL2Z"] == Fr(-1, 6) and c["equals_minus_zeta_minus_1_over_2"]
    assert c["matches"] and c["g_minus_1_from_zeta"] == 12 and c["index"] == 168 and c["cusps"] == 4 and "coincidence" in c["abelianisation_match"]
    assert D._elliptic_points(143) == (0, 0) and D._elliptic_points(13) == (2, 2) and D._elliptic_points(65) == (4, 0)
    for N in (57, 65, 91, 145):
        cN = D.euler_characteristic_is_zeta_minus_one(N); assert cN["matches"], N
    bc = D.bost_connes()
    assert bc["wall"] == 1.0 and "unbroken" in bc["hot_phase"] and "broken" in bc["cold_phase"] and "no map" in bc["physical_temperature"]
    assert "not a KMS state" in bc["marked_gas"] and "nothing follows" in bc["weights_at_the_wall"]
    for name in ("divisor_gas_energy", "two_primon_energies", "mean_weight_glaisher", "alphabet_negative_values", "clock_is_zeta_minus_one"):
        e = REGISTRY[name]; assert e.exactness == "EXACT" and all(u in REGISTRY for u in e.upstream)
    assert "keeps every pole" in REGISTRY["two_primon_energies"].nature and "-1/zeta(-1)" in REGISTRY["mean_weight_glaisher"].nature   # CC-65 wording
    assert "coincidence of size" in REGISTRY["clock_is_zeta_minus_one"].nature and "euler_characteristic" in REGISTRY["clock_is_zeta_minus_one"].example


# ------------------------------------------------------------------------------------------------------------------------------ the matter clock
def test_classifier_covers_the_eight_types():
    assert VA.classify([-1, -2, -3])[0] == 2 and VA.classify([1, 2, 3])[0] == 1
    assert VA.classify([-1, 2, 3])[0] == 3 and VA.classify([-1, -2, 3])[0] == 4
    assert VA.classify([complex(-0.3, 1), complex(-0.3, -1), -2])[0] == 5 and VA.classify([complex(0.3, 1), complex(0.3, -1), 2])[0] == 6
    assert VA.classify([complex(0.3, 1), complex(0.3, -1), -2])[0] == 7 and VA.classify([complex(-0.3, 1), complex(-0.3, -1), 2])[0] == 8
    assert VA.classify([0.0, 1, 2])[0] is None and VA.classify([-1.0])[0] == 2 and VA.classify([-1.0, 2.0])[0] is None


def test_shape_well_and_volume():
    sw = VA.shape_well()
    assert sw["cosh_identity"] and sw["curvature_at_floor_phi"] == 1 and (sw["A"], sw["B"]) == (19, 11) and sw["type"][0] == 2
    assert float(sw["r_star"]) == pytest.approx(math.sqrt(19 / 11)) and float(sw["hessian_r"]) > 0
    vr = VA.volume_runaway()
    assert vr["dV_dVol_negative"] and not vr["finite_critical_point"] and "runaway" in vr["pv01"]


def test_condensate_and_wilson_line():
    cd = VA.condensate_dterm(m2=Fr(1, 20))
    assert cd["gradient_zero"] and cd["degenerate"] and cd["type_unlifted"][0] is None and cd["hessian_eigenvalues"][0] == pytest.approx(0.0, abs=1e-12)
    assert cd["vbar_eigenvalue_is_4m2"] and cd["type_lifted"][0] == 2 and cd["v_lifted"] ** 2 == pytest.approx(8 / 13 * (1 - 1 / 20))
    for eps in ((0, 1), (1, 0), (1, 1)):
        wl = VA.wilson_line_hessians(eps)
        assert wl["type_jet"][0] == 2 and wl["type_untwisted"][0] == 1 and wl["minus_each_other"] and min(wl["eigs_jet"]) > 0 and max(wl["eigs_untwisted"]) < 0
    from mtft.matter.path_a import hessian_sign
    wl = VA.wilson_line_hessians((0, 1))
    H = np.array(wl["hessian_jet"]); d = np.array([0.0, 0.5])                                     # along the straight path to theta_eta = eps/2
    assert float(d @ H @ d) == pytest.approx(hessian_sign((0, 1))["at_jet_vacuum"], rel=1e-9)       # PA-01's second derivative along the straight path


def test_atlas_friction_toy_and_gate():
    A = VA.atlas()
    assert len(A["rows"]) == 7 and A["saddles"] == 0 and A["spirals_without_friction"] == 0
    assert A["rows"][0][2].startswith("sink") and "runaway" in A["rows"][1][2] and A["rows"][4][2].startswith("source")
    f = VA.friction_type(1.0, 0.2); assert f["type"][0] == 5 and "spiral" in f["kind"] and f["threshold_H"] == pytest.approx(2 / 3)
    assert VA.friction_type(1.0, 1.0)["type"][0] == 2 and "node" in VA.friction_type(1.0, 1.0)["kind"]
    assert VA.friction_type(1.0, 2 / 3)["kind"].startswith("critically")
    r = VA.toy_run(0.2, t_end=50.0); assert r["captured"] and r["floor_crossings"] == 15 and r["class"] == "DIAGNOSTIC"
    r60 = VA.toy_run(0.2, t_end=60.0); assert r60["floor_crossings"] == 18 > r["floor_crossings"]     # the count grows with the run (CC-65)
    for run in (r, r60):
        assert run["energy"]["non_increasing"] and run["energy"]["fraction_remaining"] < 1e-6         # dE/dt = -3H phi'^2 <= 0: the capture argument
    r = VA.toy_run(1.0, t_end=50.0); assert r["captured"] and r["floor_crossings"] == 0 and r["energy"]["non_increasing"]
    r = VA.toy_run(0.4, mode="radiation", t_end=60.0); assert r["captured"] and r["floor_crossings"] >= 10 and r["energy"]["non_increasing"]
    v = VA.volume_companion(t_end=40.0); assert v["monotone_increasing"] and v["final_log_vol"] > 2
    v = VA.volume_companion(t_end=1000.0, n=4000); assert 0.5 < v["late_time_slope_dlogVol_dlogt"] < 2 / 3 and v["asymptotic_slope"] == 2 / 3   # -> (2/3) log t
    g = VA.simulation_gate(); assert [k for k, _ in g["needs"]][:3] == ["I.1", "I.6", "GS-1"] and g["status"].startswith("open")


def test_command_lines_and_record():
    def run(*argv):
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            code = matter_main(["atlas", *argv])
        return code, buf.getvalue()
    code, out = run(); assert code == 0 and "cosh well" in out and "simulation gate" in out
    code, out = run("--json"); d = json.loads(out); assert code == 0 and d["atlas"]["saddles"] == 0
    assert "divisor gas of 143" in D.text() and "residue of U_D at the first zero" in D.text()
    pg = [e for e in L.E if e["wave"] == "PG-01"]
    assert [e["id"] for e in pg] == ["PG-1", "PG-2", "PG-3", "PG-4", "PG-5"] and all(e["module"] == "mtft.matter.vacuum_atlas" for e in pg)
    assert [e["status"] for e in pg] == ["active", "superseded", "superseded", "active", "active"]                 # CC-65, append-only
    assert all("CC-65" in e["note"] for e in pg if e["status"] == "superseded")
    cc = next(e for e in L.E if e["id"] == "CC-65")
    assert cc["status"] == "active" and "PG-4" in cc["statement"] and "PG-5" in cc["statement"] and "keeps every pole" in cc["statement"]
    a = Z.assignment(); assert all(a[f"PG-{i}"] == "IV" for i in (1, 2, 3, 4, 5)) and a["CC-65"] == "IV"
    subs = Z.subs()
    assert subs["IV.7"][1].status == "answered" and set(subs["IV.7"][1].claims) == {"PG-1", "PG-4", "PG-5", "CC-65"}
    assert subs["I.7"][1].status == "open" and "I.6" in subs["I.7"][1].settle and subs["I.7"][1].claims == ("PG-5",)
    assert any(n["id"] == "N22" and "PG-4" in n["claims"] and "PG-2" not in n["claims"] for n in Z.NEGATIVES)
