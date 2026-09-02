"""mtft.surface gates (v0.25.0).  GP-dependent tests skip without PARI/GP."""
import numpy as np
import pytest

import mtft.surface as S
from mtft.surface import gauge as G, hodge as H

GP = __import__("mtft.gprun", fromlist=["find_gp"]).find_gp()


def test_invariants_143_and_crt_projective_line():
    inv = S.invariants(143)
    assert (inv.index, inv.genus, inv.cusps, inv.nu2, inv.nu3) == (168, 13, 4, 0, 0)
    assert inv.cusp_widths == (1, 11, 13, 143)
    assert len(S.manin.projective_line(143)) == 168
    assert S.invariants(6).genus == 0 and S.invariants(6).nu3 == 0          # CC-MSL-01
    assert S.invariants(105).genus == 13 and S.invariants(105).cusps == 8


def test_elliptic_levels_refused():                                        # CC-MSL-02
    for N in (5, 7, 13):
        with pytest.raises(ValueError, match="elliptic"):
            S.cell_complex(N)


@pytest.mark.parametrize("N", [6, 11, 35, 105, 143])
def test_exact_layer_gates(N):
    cx = S.cell_complex(N)
    S.assert_gates(cx.gates())
    cb = S.tree_cotree(cx)
    S.assert_gates(cb.gates)
    assert cb.rank == 2 * cx.inv.genus and abs(cb.unimodular_determinant) == 1


def test_certified_hodge_and_whitney_family_143():
    cx = S.cell_complex(143); cb = S.tree_cotree(cx)
    uh = H.unweighted_hodge(cx, cb)
    assert uh["nullity"] == 26 and abs(uh["intersection_det"]) == 1 and uh["intersection_rounding"] < 1e-10
    fam = H.refinement_family(cx, cb, uh["intersection_cycles"], 2)
    assert [lv["euler"] for lv in fam] == [-24, -24, -24]
    assert [lv["cells"] for lv in fam] == [(4, 84, 56), (88, 336, 224), (424, 1344, 896)]
    ratios = [lv["branch"]["bracket_ratio"] for lv in fam]
    assert ratios[0] > ratios[1] > ratios[2] > 2.0                           # non-vacuity: not converged
    ex = [lv["branch"]["exact"][0] for lv in fam]
    assert ex[0] >= ex[1] >= ex[2]                                           # conforming branch non-increasing
    # r=0 identity: H^T M1 H = I/sqrt(3) on the l2-harmonic subspace -> same star
    d01 = H.siegel_distance(fam[0]["G"], fam[1]["G"])[0]
    assert d01 > 0.05


def test_gauge_exact_formulas_and_ast_census():
    inv = S.invariants(143)
    assert abs(G.closed_area(inv) - 48 * np.pi) < 1e-12
    assert abs(G.flux_action(inv, 1) - np.pi / 24) < 1e-14
    assert G.flat_connection_torus_dimension(inv) == 26 and G.flat_connection_torus_dimension(inv, True) == 29
    assert G.line_bundle_index(inv, 12) == 0
    z = G.ym_partition_function(inv, "SU(2)", 1.0)
    assert z["converged"] and abs(z["Z"] - 1.0) < 1e-6                       # (2j+1)^-24 suppression
    ctl = G.same_genus_control(105, 143, 0.01)
    assert ctl["closed_equal"] and ctl["cusps"] == (8, 4) and ctl["cusp_relative_difference"] > 1e-12
    c = G.line_operator_census(35)
    assert all(g["status"] == "PASS" for g in c["gates"]) and c["theta_orbit_sizes"] == [1, 5, 7, 35]
    assert G.manin_theta_orbits(35) == [1, 5, 7, 35]
    c12 = G.line_operator_census(12)
    assert c12["lagrangian_count"] == 28 and S.invariants(12).index == 24     # 4 non-cyclic Lagrangians


@pytest.mark.skipif(GP is None, reason="PARI/GP not found")
def test_transport_periods_and_j_invariant_143():
    from mtft.surface import hodge_structure as HS, transport as TR
    rep = S.report(143, max_refinement=2)
    assert S.all_pass(rep)
    tr = rep["transport"]
    assert str(tr.charpoly_route_A[2].as_expr().factor()).count("**2") >= 2
    assert TR.sector_census(tr.atkin_lehner, 11, 13) == {"(+,+)": 2, "(+,-)": 12, "(-,+)": 10, "(-,-)": 2}
    d = [x["siegel_distance"] for x in rep["family_distances"]]
    assert d[0] > d[1] > d[2] > 0.1
    assert abs(d[0] - 0.3764) < 2e-3 and abs(d[2] - 0.1755) < 2e-3          # lab v0.4.0 regression
    r = HS.elliptic_block_check(rep["hodge_structure"], tr.hecke[2], 0, [0, -1, 1, -1, -2])
    assert r["status"] == "PASS" and r["closest_rel_error"] < 1e-8           # j(143a1) from periods
    import mtft.periods as P
    assert P.sector_census() == (1, 6, 5, 1)                                  # third route agrees
