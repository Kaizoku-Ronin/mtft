from fractions import Fraction
from mtft.research import unified_parent as UP, compactification as CP

def test_parent_admissibility_and_family_blocks():
    for g in ("E6", "E7", "E8", "G2", "F4"): a = UP.parent_admissible(g); assert a["no_cubic_casimir"] and a["gauge_anomaly_factorizes"]
    assert not UP.parent_admissible("SU(3)")["no_cubic_casimir"] and not UP.parent_admissible("SO(10)")["gauge_anomaly_factorizes"]
    e7 = UP.e7_three_27_families(); assert e7["net_27_families"] == 3 and e7["pure"] and sum(e7["decomposition"].values()) == 133
    e6 = UP.e6_three_16_families(); assert e6["net_16_families"] == 3 and sum(e6["decomposition"].values()) == 78
    assert UP.e8_adjoint_net_families((3, 3, 3, -9))["net_families"] == 0 and UP.e8_adjoint_net_families((3, 0, -3))["net_families"] == 0
    for g, n in (("E6", 78), ("E7", 133), ("E8", 248)):                       # CC-33: dim(adj)/28 is never an integer for E6, E7, E8
        r = UP.gaugino_p2_integrality(g); assert not r["integral"] and r["net_tensors_needed"] == Fraction(n, 28)
    assert UP.susy_gravitational_condition(133, 1)["H"] == 377 and UP.susy_gravitational_condition(248, 1)["H"] == 492

def test_radion_is_brans_dicke():
    r = CP.radion_brans_dicke(2); assert r["omega_BD"] == Fraction(-1, 2) and not r["passes_cassini_massless"]
