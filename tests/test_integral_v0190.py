"""v0.19.0 integral-model arc: certificates v1-v9 as live gates."""
import pytest

from mtft.canonical import integral_gates as G


def test_frame():
    assert G.gate_frame()["ideal_basis_frame"] == "s2"


def test_integral_model_gate():
    assert G.gate_integral_model()["saturation_steps"][2] == 25


def test_counts_mod2_model_triple():
    assert G.gate_counts_mod2() == {"saturated": 4, "packaged_s2": 7,
                                    "adapted_mixed": 3}


def test_ci_a_both_frames():
    out = G.gate_ci_a()
    assert out["a_codifferent"] == -637


def test_cusp_bijection_mod2():
    assert G.gate_cusps(2)["bijection"]


@pytest.mark.slow
def test_counts_mod3():
    assert G.gate_counts_mod3()["saturated"] == 4


@pytest.mark.slow
def test_cusp_bijection_mod3():
    assert G.gate_cusps(3)["bijection"]


@pytest.mark.slow
def test_product_chain_13_49_637():
    assert G.gate_product_chain()["index"] == 49


@pytest.mark.slow
def test_al_splitting():
    q = G.gate_al_splitting()
    assert q["W13_in_W11minus"] == [1, 1, 1, 1, 1, 26]
