"""mtft.maass.spectrum: the frozen Maass spectrum of Gamma_0(143), 41 untwisted and 42 eta-twisted forms."""
import numpy as np
import pytest

from mtft.maass import hejhal as H
from mtft.maass import spectrum as SP


def test_data_files_match_their_provenance():
    assert SP.verify_provenance() == {"maass_143.json": True, "maass_143_coefficients.npz": True}
    meta = SP.metadata()
    assert meta["level"] == 143 and meta["coefficients_kept"] == 400
    assert meta["setup"]["M"] == 632 and meta["setup"]["samples"] == 788
    assert meta["complete_below_R"] == {"untwisted": 1.98, "twisted": 1.8}


def test_counts_by_sector():
    un, tw = SP.forms(twist=False), SP.forms(twist=True)
    assert len(un) == 41 and len(tw) == 42 and len(SP.forms()) == 83
    assert len({f["id"] for f in SP.forms()}) == 83
    assert SP.sectors(False) == {("even", 1, 1): 4, ("even", 1, -1): 4, ("even", -1, 1): 4, ("even", -1, -1): 7,
                                 ("odd", 1, 1): 6, ("odd", 1, -1): 6, ("odd", -1, 1): 5, ("odd", -1, -1): 5}
    assert sum(SP.sectors(True).values()) == 42 and len(SP.sectors(True)) == 8
    assert max(f["R"] for f in un) < 2.0 and max(f["R"] for f in tw) < 1.8
    assert all(f["c"].shape == (400,) and f["c"][0] == 1.0 for f in SP.forms())
    assert [f["R"] for f in un] == sorted(f["R"] for f in un)
    with pytest.raises(KeyError):
        SP.form("u.o.mm.99")


def test_every_form_was_found_again_by_the_package_solver():
    """The table's own record of its build: each R against the lab's, the two horocycles, the Hecke relations."""
    fs = SP.forms()
    assert max(abs(f["R"] - f["R_lab"]) for f in fs) < 1e-10      # the build: 5.0e-12
    assert max(f["agreement"] for f in fs) < 1e-9                 # 1.2e-11
    assert max(f["hecke"] for f in fs) < 1e-9                     # 2.3e-11
    assert min(f["agreement"] for f in fs) > 0                    # computed, not filled in
    assert all(abs(f["lambda"] - (0.25 + f["R"] ** 2)) < 1e-14 for f in fs)


def test_hecke_relations_hold_on_the_frozen_coefficients():
    """Recomputed here from the coefficient file; the relations were never imposed."""
    worst = max(H.hecke_residual(f["c"], f["twist"]) for f in SP.forms())
    assert 0 < worst < 1e-9
    f = SP.first()
    broken = f["c"].copy()
    broken[5] += 1e-6                                             # c_6 spoiled: the check must see it
    assert H.hecke_residual(broken) > 5e-7


def test_newforms_have_c_p_equal_minus_eps_p_over_sqrt_p():
    new = [f for f in SP.forms() if f["newform"]]
    old = [f for f in SP.forms() if not f["newform"]]
    assert len(new) >= 60 and len(old) >= 4
    for f in new:
        assert abs(f["c11_sqrt11"] + f["eps11"]) < 1e-8 and abs(f["c13_sqrt13"] + f["eps13"]) < 1e-8, f["id"]
    for f in old:                                                 # an oldform misses the rule at one prime at least
        assert abs(f["c11_sqrt11"] + f["eps11"]) > 1e-3 or abs(f["c13_sqrt13"] + f["eps13"]) > 1e-3, f["id"]


def test_the_gap():
    g = SP.gap()
    assert g["id"] == "u.o.mm.1" and g["parity"] == "odd" and g["eps"] == (-1, -1)
    assert abs(g["R_1"] - 0.3738518157) < 1e-9 and abs(g["lambda_1"] - 0.3897651801) < 1e-9
    assert all(f["lambda"] > 0.25 for f in SP.forms())            # nothing below the continuum's edge
    assert abs(SP.first(True)["R"] - 0.4627415219) < 1e-9


def test_oldforms_carry_the_eigenvalue_of_their_lower_level():
    """The level-13 form (LMFDB 13.1), computed at level 13, reappears at level 143 in two sectors."""
    r13 = H.refine(0.97081541, 1, (-1,), False, H.setup(13), half_width=1e-6)["R"]
    twins = [f for f in SP.forms(twist=False, parity=1) if abs(f["R"] - r13) < 1e-6]
    assert sorted(f["eps"] for f in twins) == [(-1, -1), (1, -1)]             # eps_13 = -1 kept, both signs at 11
    assert all(abs(f["R"] - r13) < 1e-9 and not f["newform"] for f in twins)


def test_series_below_the_domain_agrees_with_the_value_carried_into_it():
    """On a horocycle under the Ford domain the Fourier series must equal chi times the value at the reduced point."""
    xs = np.array([0.013, 0.087, 0.2041, 0.3377, 0.4512, -0.1234])
    for fid in ("u.o.mm.1", "u.e.mm.1", "u.o.pp.1", "t.o.mm.1", "t.e.pp.1"):
        f = SP.form(fid)
        direct, carried = SP.series(f, xs, 0.0093), SP.evaluate(f, xs, 0.0093)
        assert np.max(np.abs(direct - carried)) < 1e-8 * max(1.0, np.max(np.abs(direct))), fid
        assert np.max(np.abs(direct)) > 1e-3                      # not zero against zero
    f, g = SP.form("u.o.mm.1"), dict(SP.form("u.o.mm.1"))
    g["eps11"] = 1                                                # a wrong sign at 11: the two must now differ
    assert np.max(np.abs(SP.series(f, xs, 0.0093) - SP.evaluate(g, xs, 0.0093))) > 1e-3
    assert abs(SP.evaluate(f, 0.0, 1.0)[()]) < 1e-15              # an odd form vanishes on the imaginary axis


def test_rederive_one_form_at_call_time():
    r = SP.rederive("u.e.pp.1")
    assert abs(r["difference"]) < 1e-10 and r["agreement"] < 1e-9 and r["hecke"] < 1e-9 and r["id"] == "u.e.pp.1"
