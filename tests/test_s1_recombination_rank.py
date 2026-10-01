"""Regression gates for exceptional S1 torus-pairing ranks."""
import mpmath as mp
import pytest
import sympy as sp
from mtft.research import product_surface as PS, theta_torus as TT


def test_s1_gluing_generic_and_exceptional_counts():
    assert [PS.s1_colour_recombination_dims(k)["dim_Ext1_E_k_L_L"] for k in (1,2,3)] == [17,2,2]
    assert PS.s1_colour_recombination_dims(1, pairing_rank=0)["dim_Ext1_E_k_L_L"] == 32
    assert PS.s1_colour_recombination_dims(2, pairing_rank=1)["dim_Ext1_E_k_L_L"] == 17
    assert PS.s1_colour_recombination_dims(1)["possible_gluing_dimensions"] == (32,17)
    assert PS.s1_colour_recombination_dims(2)["possible_gluing_dimensions"] == (17,2)
    assert PS.s1_colour_recombination_dims(3)["possible_gluing_dimensions"] == (2,)
    assert PS.s1_colour_recombination_dims(1)["rank_assumption"].startswith("generic")


def test_s1_explicit_kernel_and_independent_rank_one_plane():
    # Formula theta_product_coefficient is invariant under shift -> -shift.
    # Rows j=1,2 have shifts +/-1/3 (column 0), +/-1/6 (column 1).
    B = TT.torus_factor_closed_form(1,2,dps=60)["B"]
    xi = mp.matrix([0,1,-1]); e0 = mp.matrix([1,0,0])
    assert sp.Matrix.hstack(sp.Matrix([0,1,-1]), sp.Matrix([1,0,0])).rank() == 2
    assert mp.norm(B.T * xi) < mp.mpf("1e-50")
    assert mp.norm(B.T * e0) > mp.mpf("0.1")
    # Thus span(xi) has pairing rank0, while span(xi,e0) has rank1.
    assert PS.s1_colour_recombination_dims(1,0)["dim_Ext1_E_k_L_L"] == 32
    assert PS.s1_colour_recombination_dims(2,1)["dim_Ext1_E_k_L_L"] == 17


@pytest.mark.parametrize("k,rank", [(1,2),(2,0),(3,1),(4,None),(-1,None),(1.5,None)])
def test_s1_rejects_ranks_incompatible_with_independence(k,rank):
    with pytest.raises(ValueError):
        PS.s1_colour_recombination_dims(k,rank)
