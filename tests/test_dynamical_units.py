"""tests/test_dynamical_units.py — gates for the dynamical-units wave
(studies/du01_two_clock_ledger.py, studies/du02_cycle_space_map.py),
v0.11.4 integration wave.

Both study mains run end-to-end in ~2 s each and assert internally at
every stage; this tier re-runs them with the ledger output redirected
to a tmp dir and then re-asserts the headline EXACT certificates from
the outside:

  du01  A1 tr L = 166 (one self-loop)          [Add. BT leg 1]
        A3 integer spectrum {0,1,2,4,5} simple  [Add. BT leg 1]
        A4 Kirchhoff tau = 3518081582959364640  [Add. BT leg 1]
        A5 b1 = 29 = 2g + (cusps - 1)           [Add. BT leg 1]
        B4 uniform-weight variances 0, 35/16, 10/3 (exact rationals)
                                              [Add. BT leg 2]
        C3 anchor count = 2 (standing H4 negative carried)
  du02  A2 det D = -1 (Lefschetz duality, shipped conventions)
                                              [Add. BT leg 3]
        B1 width census {1:143, 11:13, 13:11, 143:1}
        B4 self-loop = link of the width-1 (divisor-143) cusp
        C2/C5 transported Hecke integral, same 13 lines
        C4 Eisenstein block charpoly (x - 3)^3
        F1 signature basis counts Eis:3 old:4 f1:2 f2:8 f3:12
"""
import importlib.util
import os
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
STUDIES = os.path.join(HERE, "..", "studies")
if STUDIES not in sys.path:
    sys.path.insert(0, STUDIES)  # studies scripts import each other by bare name


def _load(name):
    spec = importlib.util.spec_from_file_location(
        name, os.path.join(STUDIES, name + ".py"))
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def du01_ledger(tmp_path_factory):
    mod = _load("du01_two_clock_ledger")
    mod._HERE = str(tmp_path_factory.mktemp("du01"))
    mod.main()
    return mod.LEDGER


@pytest.fixture(scope="module")
def du02_ledger(tmp_path_factory):
    mod = _load("du02_cycle_space_map")
    mod._HERE = str(tmp_path_factory.mktemp("du02"))
    mod.main()
    return mod.LEDGER


# ---------------------------------------------------------------- du01

def test_du01_graph_clock_exact(du01_ledger):
    L = du01_ledger
    assert L["A1 self_loops (EXACT)"] == 1
    assert L["A1 trace_L (EXACT)"] == 166
    assert [tuple(x) for x in L["A3 integer_eigenvalues (EXACT)"]] == [
        (0, 1), (1, 1), (2, 1), (4, 1), (5, 1)]
    assert L["A4 spanning_trees (EXACT)"] == 3518081582959364640
    assert L["A4 spanning_trees (EXACT)"] == 2**5 * 3 * 5 * 17 * 941 * 101921 * 4495339
    assert L["A5 b1_spine (EXACT)"] == 29


def test_du01_hecke_clock(du01_ledger):
    L = du01_ledger
    assert L["B2 oracle_cross_cert_max_err (Cert)"] < 1e-9
    assert 0.95 < L["B3 ramanujan_margin_max (Cert)"] < 1.0
    assert L["B4 uniform_var (EXACT)"] == {"f1": "0", "f2": "35/16",
                                           "f3": "10/3"}
    assert L["B4 var_E2_err (Cert)"] < 1e-9
    lines = L["B1 T2_lines (Cert, x2 each)"]
    assert len(lines) == 13
    assert abs(sum(lines) + 1.0) < 1e-6   # old -2 x2 + new (0+3+0) = -1; tr T2 = 2*sum = -2


def test_du01_anchor_protocol(du01_ledger):
    assert du01_ledger["C3 anchor_count"] == 2
    assert "H4" in du01_ledger["C3 standing_negative"]


# ---------------------------------------------------------------- du02

def test_du02_lefschetz_duality(du02_ledger):
    assert du02_ledger["A1 relations_killed (EXACT)"] is True
    assert du02_ledger["A2 det_D (EXACT)"] == -1


def test_du02_cusp_links(du02_ledger):
    L = du02_ledger
    assert L["B1 width_census (EXACT)"] == {"1": 143, "11": 13,
                                            "13": 11, "143": 1}
    assert L["B4 selfloop_is_width1_link (EXACT)"] == "143"
    assert L["B5 link_supports (EXACT)"] == {"1": 23, "11": 13,
                                             "13": 11, "143": 1}


def test_du02_transported_hecke(du02_ledger):
    L = du02_ledger
    assert L["C2 T2star_integral (EXACT)"] is True
    assert L["C4 eisenstein_charpoly"] == "(x - 3)**3"
    assert L["C5 charpoly_match (EXACT)"] is True
    assert L["C5 cuspidal_charpoly"] == (
        "x**2*(x + 2)**4*(x**4 - 3*x**3 - x**2 + 5*x + 1)**2*"
        "(x**6 - 10*x**4 + 2*x**3 + 24*x**2 - 7*x - 12)**2")


def test_du02_interface(du02_ledger):
    L = du02_ledger
    assert L["F1 signature_counts (Cert)"] == {"Eis": 3, "old": 4,
                                               "f1": 2, "f2": 8,
                                               "f3": 12}
    assert L["F1 basis_cond (Cert)"] < 1e3
    gram = sorted(L["D4 link_gram_levels (Cert)"])
    assert len(gram) == 3 and abs(gram[0] - 1.30564596) < 1e-6
