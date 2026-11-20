"""TE-01: the twist landscape — the relative heat trace, and the energy and pin of every twist (records FF-14, FF-15, FF-20)."""
import numpy as np
import pytest

from mtft.lines import twist_energy as TE

# FF-15's energies (ledger units) of its named twists at r = 230 and 280 (record ff15_numbers.json, "selection table").
FF15_ENERGY = {
    6665056: (-5.965621297963878e-05, -3.306473860490074e-05),      # φ*η
    14025081: (-0.0001516849310759383, -8.407194262518517e-05),     # 143a1's pullback (0, 1)
    11776537: (-0.00015168493107593786, -8.407194262518539e-05),    # 143a1's pullback (1, 1)
    52741932: (-0.002105175658998131, -0.0011621601731760963),      # X0(11)'s pullbacks
    64041154: (-0.002105175658998131, -0.0011621601731760963),
    16119790: (-0.003934899299766736, -0.0021716501167429575),
    45796554: (-0.004091338357100325, -0.002258355343722602),       # FF-14's hardest pin: the energy minimum
    394080: (-3.8122674971724724e-05, -2.1129673179971348e-05),     # FF-14's weakest pin
    39774856: (-0.00399480163944004, -0.002204849621026865),        # the cuspidal twists φA·φB, φ1·φB, φ1·φA
    33196637: (-0.0002733448310348769, -0.00015150088718321724),
    61104341: (-0.0040094439751875, -0.0022129666801225943),
}


@pytest.fixture(scope="module")
def curves():
    return {z: TE.relative_heat_trace(z) for z in (394080, 45796554)}


def test_fwht_against_the_direct_transform():
    a = np.random.default_rng(143).normal(size=256)
    F = TE.fwht(a.copy())
    direct = [sum(a[y] * (-1) ** bin(x & y).count("1") for y in range(256)) for x in range(256)]
    assert np.allclose(F, direct, rtol=0, atol=1e-10)
    with pytest.raises(ValueError):
        TE.fwht(np.zeros(6))


def test_classes_against_ff14_and_the_frozen_trace_weights():
    L0, Y, _ = TE.classes(14.0)
    assert len(L0) == 93720 and int((Y == 0).sum()) == 248 and float(L0[Y == 0].min()) == pytest.approx(9.897, abs=5e-4)
    chk = TE.per_trace_check(14.0)
    assert (chk["traces with classes"], chk["empty traces (both routes)"], chk["empty in one route only"]) == (321, 773, 0)
    assert chk["traces found here but absent from the frozen table"] == 0 and chk["max |ratio - 1|"] < 1e-12
    assert TE.switch_check(230.0) < 1e-12 and TE.switch_check(280.0) < 1e-12


def test_the_relative_heat_trace(curves):
    lo, ce = TE.relative_heat_trace(0)
    assert not lo.any() and not ce.any()
    assert TE.heat_terms([2.633915793849633, 9.9, 14.0]).min() >= 0
    lo, ce = curves[394080]
    assert lo.min() >= 0 and (ce >= lo).all() and ce[-1] == pytest.approx(1.0, abs=1e-12)
    # FF-14's J = ∫ t̃⁻² D for its weakest and hardest twists and for the mean twist (record ff14_numbers.json, T4-T5)
    assert (TE.J(ce), TE.J(lo)) == (pytest.approx(0.11183864192916144, rel=1e-9), pytest.approx(0.0369746989156233, rel=1e-9))
    lo, ce = curves[45796554]
    assert (TE.J(ce), TE.J(lo)) == (pytest.approx(1.1489525772212512, rel=1e-9), pytest.approx(1.074088634207713, rel=1e-9))
    lo, ce = TE.mean_relative_heat_trace()
    assert (TE.J(ce), TE.J(lo)) == (pytest.approx(0.6776921189970775, rel=1e-9), pytest.approx(0.6028281759835392, rel=1e-9))


def test_energies_of_ff15s_named_twists():
    for z, (v230, v280) in FF15_ENERGY.items():
        assert TE.energy(z, 230.0) == pytest.approx(v230, rel=1e-9)
        assert TE.energy(z, 280.0) == pytest.approx(v280, rel=1e-9)
    assert TE.energy(0) == 0.0 and TE.pin(0) == 0.0
    w, t = TE.energy_weights(230.0)
    assert (w < 0).all() and t < 0                         # every weight lowers the energy, and so does the tail
    assert TE.energy(45796554) == pytest.approx(TE.energy(45796554, central=False) + t, rel=1e-12)


def test_two_routes_inside_the_tool(curves):
    for z, (lo, ce) in curves.items():
        pins = TE.pin_matrix(ce, 230.0)["per-cycle pins"]
        assert TE.pin(z, 230.0, "1") == pytest.approx(pins[0], rel=1e-9)
        assert TE.pin(z, 230.0, "tau") == pytest.approx(pins[1], rel=1e-9)
        assert TE.pin(z, 230.0, "1", central=False) == pytest.approx(TE.pin_matrix(lo, 230.0)["per-cycle pins"][0], rel=1e-9)
        assert TE.energy(z, 230.0) == pytest.approx(TE.hold(ce, 230.0)["V (w = 0 included)"], rel=1e-9)


def test_pins_and_the_hold_against_ff14(curves):
    _, ce = curves[45796554]
    pm = TE.pin_matrix(ce, 230.0)
    assert pm["per-cycle pins"] == [pytest.approx(0.3141647944463996, rel=1e-9), pytest.approx(0.2390672652887669, rel=1e-9)]
    assert pm["twelve: 64π² eig(M_A)"] == [pytest.approx(0.07566178788347801, rel=1e-9), pytest.approx(0.20095424198410528, rel=1e-9)]
    h = TE.hold(ce, 230.0)
    assert h["V (without w = 0, FF-14's print)"] == pytest.approx(-0.0033578210220447215, rel=1e-9)
    assert h["|grad V|_hyp (analytic)"] == pytest.approx(0.00011608940307752983, rel=1e-8)
    assert h["|grad V|_hyp (differences)"] == pytest.approx(h["|grad V|_hyp (analytic)"], rel=1e-6)
    assert h["sinh d"][0.31] == pytest.approx(0.00037128891200956444, rel=1e-8)
    _, ce = curves[394080]
    h = TE.hold(ce, 230.0)
    assert h["V (without w = 0, FF-14's print)"] == pytest.approx(-3.74159114751194e-05, rel=1e-9)
    assert h["V (w = 0 included)"] == pytest.approx(-3.8122674971724724e-05, rel=1e-9)      # FF-15's convention


def test_the_landscape_is_invariant_under_the_involutions():
    from mtft.lines import twists as T
    for z in [394080, 45796554] + np.random.default_rng(20).integers(1, 1 << 26, 5).tolist():
        for g in ("W11", "W13", "W143", "mirror"):
            assert TE.energy(T.act(g, z)) == pytest.approx(TE.energy(z), rel=1e-12)
            assert TE.pin(T.act(g, z), 280.0, "tau") == pytest.approx(TE.pin(z, 280.0, "tau"), rel=1e-12)
    v = np.full(16, 2.0)
    v[[3, 9]] = 1.0
    v[0] = np.nan
    s = TE.summary(v)
    assert s["twists at the min"] == [3, 9] and s["argmin"] == 3 and s["argmax"] == 1


def test_bands_on_a_small_space():
    z = np.arange(256)
    k = (z & 1) + ((z >> 1) & 1)
    V = -(k + 0.1 * z / 256) - 0.05
    V[0] = np.nan
    b = TE.bands(V, systole_classes=(1, 2), chunk=64)
    assert [b["bands"][j]["twists"] for j in range(3)] == [63, 128, 64] and b["overlap"] is False
    assert b["bands"][2]["depth max"] == 1.0
    with pytest.raises(ValueError):
        TE.landscape(quantity="mass")
