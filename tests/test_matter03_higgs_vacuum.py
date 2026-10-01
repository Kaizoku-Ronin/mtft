"""Independent exact checks of the paired branch and its declared limits."""
from fractions import Fraction as F
import json
import pytest

from mtft.matter.higgs_vacuum import (
    HD_VECTOR, HU_VECTOR, POWERS, Q, paired_higgs_witness,
    single_arrow_gate, classical_gauge_certificate, report,
)


def test_frozen_sections_reconstructed_from_archived_basis_independently():
    from mtft.canonical import s2_qexpansions
    from mtft.matter.collinear import yukawa_witness
    matrix = s2_qexpansions()
    down = tuple(sum(F(x)*v for x, v in zip(row, HD_VECTOR)) for row in matrix)
    up = tuple(sum(F(x)*v for x, v in zip(row, HU_VECTOR)) for row in matrix)
    pair = paired_higgs_witness()
    assert down == pair.Hd_q_coefficients == yukawa_witness().canonical_q_coefficients
    assert up == pair.Hu_q_coefficients
    assert down[1:6] == (0, 0, 0, 1, 0)
    assert up[1:5] == (0, 0, 1, 0)


def test_two_coordinate_charts_give_same_odd_exceptional_product_order():
    import sympy as sp
    pair = paired_higgs_witness()
    u, v = sp.symbols("u v")
    # Compute actual Jacobians and substitutions in both charts, rather
    # than assuming the blow-up Jacobian order supplied by the theorem.
    for z, w in ((u, u*v), (u*v, u)):
        jacobian = sp.Matrix([[sp.diff(z, u), sp.diff(z, v)],
                              [sp.diff(w, u), sp.diff(w, v)]]).det()
        product = sp.expand(z**3*z**2*jacobian**2)
        product_order = min(monomial[0] for monomial, _ in sp.Poly(product, u, v).terms())
        assert product_order == pair.product_exceptional_valuation == 7
        assert product_order % 2 == 1
    assert pair.product_curve_order == 5


@pytest.mark.parametrize("t", [F(1), F(2, 3), F(-7, 5)])
def test_nonzero_common_scaling_changes_spectral_coefficient_not_parity(t):
    pair = paired_higgs_witness()
    product = {}
    for i, a in enumerate(pair.Hd_q_coefficients[:10]):
        for j, b in enumerate(pair.Hu_q_coefficients[:10]):
            if a and b:
                product[i+j] = product.get(i+j, F(0))+t*a*t*b
    assert min(n for n, value in product.items() if value) == 7
    assert product[7] == t*t
    assert min(n for n, value in product.items() if value)-2+2 == 7


@pytest.mark.parametrize("eps", [F(1, 2), F(1), F(3, 2)])
def test_spectator_constraints_force_wall_but_single_arrow_is_still_excluded(eps):
    from mtft.surface.blowup import slope
    gate = single_arrow_gate(3, 3, eps)
    assert gate["A1_slope"] == slope((2, 0, -6), 3, 3, eps)
    assert gate["A2_slope"] == slope((-3, 0, 9), 3, 3, eps)
    assert gate["spectator_necessary_gate"] == (eps == 1)
    assert not gate["nonzero_single_Hd_D_flat"]
    assert gate["on_wall"] == (eps == 1)


@pytest.mark.parametrize("args", [(3, 3, 0), (3, 3, 3), (0, 3, 1)])
def test_nonample_polarizations_rejected(args):
    with pytest.raises(ValueError):
        single_arrow_gate(*args)


@pytest.mark.parametrize("bad", [True, 3.0, "3", None])
def test_exact_polarization_types_are_explicit(bad):
    with pytest.raises(TypeError):
        single_arrow_gate(bad, 3, 1)


def test_gauge_dimension_and_Q_with_direct_matrix_multiplication():
    def unit(i, j):
        return [[int(a == i and b == j) for b in range(8)] for a in range(8)]
    def mul(a, b):
        return [[sum(a[i][k]*b[k][j] for k in range(8)) for j in range(8)] for i in range(8)]
    qmat = [[Q[i] if i == j else F(0) for j in range(8)] for i in range(8)]
    for root in (unit(3, 7), unit(7, 3)):
        assert mul(qmat, root) == mul(root, qmat)
    assert sum(POWERS) == 0
    gauge = classical_gauge_certificate()
    assert gauge["dimension"] == 13
    assert gauge["dimension_SU8"] == 12
    assert gauge["Q_commutes_with_both_arrows"]
    assert gauge["axionic_photon_kernel"] == "OPEN"


def test_stated_existence_is_conditional_without_metric_or_mass_selection():
    output = report()
    json.dumps(output)
    assert output["paired_higgs"]["active_rank2_Higgs_pair_stable"]
    assert output["paired_higgs"]["full_rank8_polystable_on_wall"]
    assert output["paired_higgs"]["Hermitian_metric_exists_by_surface_HK"]
    assert not output["paired_higgs"]["Hermitian_metric_computed"]
    assert not output["paired_higgs"]["physical_mass_values_derived"]
    assert not output["branch"]["selected_scale"]
    assert not output["branch"]["Higgs_modulus_stabilized"]
    assert not output["classical_gauge"]["unique_photon"]
    assert "https://arxiv.org/pdf/1312.2673" in [r["url"] for r in output["primary_references"]]
