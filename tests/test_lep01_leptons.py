"""LEP-01: the charged-lepton texture."""
from mtft.matter import leptons as LP
from mtft.matter.flux_scan import scan, StackFlux


def test_georgi_jarlskog_at_matching():
    g = LP.gj_at_matching()
    assert 2.9 < g["mmu_over_ms"] < 3.2 and 0.27 < g["me_over_md"] < 0.32
    assert g["eps_ratio_second_family"] > 1.6 and g["eps_ratio_first_family"] < 0.8          # no common rescaling


def test_su4_identity_and_vertex_charges():
    assert LP.su4_identity_holds(scan(3, 1))
    v = LP.vertex_charges(StackFlux((-1, -1, 1, 0, 2), (0, 0, -1, 2, 0)))
    assert all(sum(t) == 0 for t in v.values()) and v["down"] == v["lepton"] and [abs(x) for x in v["up"]] == [1, 0, 1] and [abs(x) for x in v["down"]] == [1, 2, 3]


def test_alignment_fits_koide_and_numerology():
    for m in (LP.UP_MZ, LP.DOWN_MZ, LP.LEPTON_POLE):
        s, res = LP.fit_alignment(m); assert res < 1e-10
    assert abs(LP.fit_alignment(LP.DOWN_MZ)[0][1]) < 0.3
    assert all(v < 2.1 for v in LP.coherent_light_mass_spread().values())
    assert abs(LP.koide(LP.LEPTON_POLE) - 2 / 3) < 1e-5 and LP.koide(LP.DOWN_MZ) > 0.74
    n = LP.numerology_example(); assert -0.06 < n["me_over_mtau_pred"] / n["me_over_mtau_data"] - 1 < -0.05


def test_alignments_are_not_unique():                                                              # LEP-6, CC-61
    sets = {k: LP.all_alignments(tuple(m)) for k, m in (("up", LP.UP_MZ), ("down", LP.DOWN_MZ), ("leptons", LP.LEPTON_POLE))}
    assert len(sets["up"]) >= 3 and len(sets["down"]) >= 2 and len(sets["leptons"]) >= 4            # non-vacuity: several exact solutions
    for k, m in (("up", LP.UP_MZ), ("down", LP.DOWN_MZ), ("leptons", LP.LEPTON_POLE)):
        for s in sets[k]:
            r2, r1, _ = LP.mass_ratios(s)
            assert abs(r2 / (m[1] / m[2]) - 1) < 1e-8 and abs(r1 / (m[0] / m[2]) - 1) < 1e-8             # every member is an exact fit
        canon = LP.fit_alignment(m)[0]; assert canon in sets[k] and max(abs(canon[1]), abs(canon[2])) < 1     # charge-0 dominated
    assert any(s[1] > 3 for s in sets["leptons"]) and any(s[1] < 0.4 for s in sets["leptons"])           # leptons fit both ways
    assert abs(LP.fit_alignment(LP.DOWN_MZ)[0][1] - 0.2207) < 1e-3 and abs(LP.fit_alignment(LP.LEPTON_POLE)[0][1] - 0.3387) < 1e-3


def test_fit_alignment_raises_without_a_solution():                                                 # CC-11: every raise path is triggered
    import pytest
    with pytest.raises(ValueError):
        LP.fit_alignment((3.0, 2.0, 1.0))                                                          # ratios above one cannot be fitted
