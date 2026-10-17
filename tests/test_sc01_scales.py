"""SC-01 (2026-10-17): the scale ladder — every scale is M_KK times a number or times sqrt(lambda); on the measured background every one-loop
modulus is heavier than the KK scale, and the Wilson line decays linearly into messenger pairs."""
import math

from mtft.matter import ledger as L
from mtft.matter import puzzle as Z
from mtft.matter import scales as S


def test_ladder_exact_ratios():
    l = S.ladder(); t = S.tau()
    assert abs(l["M_SS"] - 0.5 / math.sqrt(t.imag)) < 1e-12 and abs(l["pair_threshold"] - 2 * l["M_SS"]) < 1e-12
    assert abs(l["pa01_messenger_M_sqrtA"] - 3.1057) < 2e-4                                   # PA-01's messenger scale, second route
    assert abs(l["saddle_gap"] - 0.562936) < 1e-6
    assert abs(l["omega_hat_over_sqrt_lambda"][0] - 0.21747) < 1e-5 and abs(l["omega_hat_over_sqrt_lambda"][1] - 0.40389) < 1e-5
    assert abs(l["m_phi_hat_over_sqrt_lambda"] - 0.13242) < 1e-5
    assert abs(l["omega_hat_over_sqrt_lambda"][1] / l["omega_hat_over_sqrt_lambda"][0] - math.sqrt(128.985 / 37.394)) < 1e-3   # PG-6's metric curvatures


def test_measured_background_puts_every_modulus_above_M_KK():
    r = S.report()
    assert 70 < r["lambda_range"][0] < 75 and 100 < r["lambda_range"][1] < 106
    assert 1.8 < r["omega_hat_range"][0][0] and r["omega_hat_range"][0][1] < 2.3 and 3.4 < r["omega_hat_range"][1][0] and r["omega_hat_range"][1][1] < 4.2
    assert 1.1 < r["m_phi_hat_range"][0] and r["m_phi_hat_range"][1] < 1.4 and r["all_moduli_above_M_KK"]
    for b in r["background"]:
        assert 0.6 < b["expansion_parameter_56pi_over_r"] < 0.8 and 4 < b["loop_parameter_lambda_over_16"] < 7


def test_dissipation_regime():
    for lam in (72.0, 103.0):
        d = S.dissipation(lam)
        assert d["both_modes_above_pair_threshold"] and not d["broad_resonance_possible"] and d["regime"].startswith("linear decay")
        assert 0.2 < d["width_over_mass_light_mode"][0] < 0.3 and 2.5 < d["width_over_mass_light_mode"][1] < 3.3
    assert S.dissipation(1.0)["regime"] == "broad parametric resonance" and S.dissipation(10.0)["regime"] == "below the messenger-pair threshold"


def test_record():
    e = {x["id"]: x for x in L.E}
    for cc, cls in (("SC-1", "EXACT"), ("SC-2", "NUMERICAL"), ("SC-3", "EXACT"), ("SC-4", "STRUCTURAL"), ("SC-5", "STRUCTURAL")):
        assert e[cc]["class"] == cls and e[cc]["module"] == "mtft.matter.scales" and e[cc]["test"] == "test_sc01_scales.py"
    sub = Z.subs()["III.5"][1]
    assert sub.status == "open" and "SC-2" in sub.claims and "heavier than the KK scale" in sub.answer + sub.settle
    assert all(Z.check().values()) and all(L.check().values())
