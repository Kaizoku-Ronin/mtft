"""Independent root-algebra and cohomology checks for the charged complex."""
from itertools import product

from mtft.matter.charged_modes import (split_charged_census,
    paired_charged_report, certified_paired_charged_report)
from mtft.matter.collinear import cohomology


NAMES = {"c1": 0, "c2": 1, "c3": 2, "L+": 3, "L-": 4,
         "A1": 5, "A2": 6, "B": 7}


def root(row):
    return NAMES[row["target"]], NAMES[row["source"]]


def trace_bracket(a, b, c):
    i, j = a
    k, l = b
    m, n = c
    return int(j == k and l == m and i == n) - int(l == i and j == m and k == n)


def test_complete_split_census_and_qzero_gauge_separation():
    census = split_charged_census()
    assert len(census["rows"]["Qminus"]) == len(census["rows"]["Qplus"]) == 6
    for sign, rows in census["rows"].items():
        for row in rows:
            cc = cohomology(row["exponent"])
            assert (row["a_dim"], row["phi_dim"], row["gauge_dim"]) == (cc.h1, cc.phi_dim, cc.h0)
        assert census["totals"][sign]["chiral_dim"] == 279
        assert census["totals"][sign]["gauge_dim"] == 1
    assert census["totals"]["Qminus"]["a_dim"] == 216
    assert census["totals"]["Qplus"]["a_dim"] == 217
    assert cohomology(-3).holo_chiral_dim == 54


def test_six_blocks_follow_actual_u8_trace_algebra():
    census = split_charged_census()
    actual = []
    for minus, plus in product(census["rows"]["Qminus"], census["rows"]["Qplus"]):
        for label, phi in (("down", (3, 7)), ("up", (7, 3))):
            if trace_bracket(root(minus), root(plus), phi):
                actual.append((minus["target"], label, minus["exponent"],
                               plus["exponent"], (minus["a_dim"], plus["a_dim"])))
    result = paired_charged_report()
    reported = [(b["spectator"], b["higgs_entry"], b["negative_a_exponent"],
                 b["positive_a_exponent"], b["shape"]) for b in result["six_a_mass_blocks"]]
    assert sorted(actual) == sorted(reported)
    assert len(actual) == 6
    assert sum(min(shape) for *_, shape in actual) == 118


def test_weak_gauge_generators_have_nonzero_phi_images_in_pair():
    # Pair Phi=E37+E73. Matrix-unit products give the two charged images.
    def bracket_entries(a, b):
        i, j = a
        k, l = b
        out = {}
        if j == k:
            out[(i, l)] = out.get((i, l), 0) + 1
        if l == i:
            out[(k, j)] = out.get((k, j), 0) - 1
        return out
    for gauge, expected in (((4, 3), {(4, 7): -1}),
                            ((3, 4), {(7, 4): 1})):
        entries = {}
        for phi in ((3, 7), (7, 3)):
            for entry, value in bracket_entries(phi, gauge).items():
                entries[entry] = entries.get(entry, 0) + value
        assert {entry: value for entry, value in entries.items() if value} == expected


def test_bound_applies_before_exact_profiles_are_selected():
    result = paired_charged_report()
    assert result["a_map_rank_upper_bound"] == 118
    assert result["paired_lower_bounds"]["Qminus"] == {
        "a_kernel_lower_bound": 98, "phi_cokernel_dimension": 62,
        "hyperH1_lower_bound": 160}
    assert result["paired_lower_bounds"]["Qplus"] == {
        "a_kernel_lower_bound": 99, "phi_cokernel_dimension": 61,
        "hyperH1_lower_bound": 160}
    assert "a_map_rank" not in result
    assert not result["physical_mass_spectrum"]


def test_exact_profile_dimensions_from_deformation_sequence():
    result = certified_paired_charged_report()
    ranks = [block["certified_rank"] for block in result["six_a_mass_blocks"]]
    assert ranks == [2, 2, 5, 3, 17, 12]
    assert result["a_map_rank"] == 41
    assert result["canonical_vanishing_orders"] == {"down": 3, "up": 2}
    assert result["paired_hyperH1"]["Qminus"] == {
        "a_kernel_dimension": 175, "phi_cokernel_dimension": 62,
        "hyperH1_dimension": 237}
    assert result["paired_hyperH1"]["Qplus"] == {
        "a_kernel_dimension": 176, "phi_cokernel_dimension": 61,
        "hyperH1_dimension": 237}
    assert not result["physical_mass_spectrum"]


def test_reports_return_fresh_mutable_data():
    result = certified_paired_charged_report()
    result["six_a_mass_blocks"][0]["certified_rank"] = 99
    result["split_census"]["rows"]["Qminus"].clear()
    fresh = certified_paired_charged_report()
    assert fresh["six_a_mass_blocks"][0]["certified_rank"] == 2
    assert len(fresh["split_census"]["rows"]["Qminus"]) == 6
