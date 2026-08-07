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


# ----------------------------------------------------------------------
# du03 — the dispersion (v0.11.5 candidate wave)
# ----------------------------------------------------------------------

@pytest.fixture(scope="module")
def du03_ledger(tmp_path_factory):
    mod = _load("du03_dispersion")
    mod._HERE = str(tmp_path_factory.mktemp("du03"))
    mod.main()
    return mod.LEDGER


def test_du03_free_channels_closed(du03_ledger):
    L = du03_ledger
    assert float(L["A1 cusp_well_closed (Cert)"]) < 1e-10
    assert L["A2 commutator_nnz (EXACT)"] == 249          # noncommuting


def test_du03_minimal_coupling(du03_ledger):
    L = du03_ledger
    assert abs(L["B2 trace_Vh (Cert)"]) < 1e-9            # traceless
    # mechanism regression: NOT equidistribution (profile-specific)
    assert float(L["C2c equidistribution_dev (Cert)"]) > 1e-3


def test_du03_parity_selection_rule(du03_ledger):
    P = du03_ledger["C5 parity"]
    assert P["anti_rel"] < 1e-6                           # V eta-ODD
    assert abs(P["comm_rel"] - 2.0) < 1e-3
    disp = du03_ledger["C2 dispersion"]
    assert len(disp) == 12
    assert all(abs(d["omega"]) < 1e-10 for d in disp)     # means vanish
    assert {d["m"] for d in disp} == {2, 4}
    f1 = next(d for d in disp if abs(d["a"]) < 1e-9)
    assert abs(f1["delta"] - 0.1141) < 5e-4               # form factor


def test_du03_second_order_dispersion(du03_ledger):
    rows = du03_ledger["C7 second_order_dispersion"]
    assert all(r["parity_pure"] for r in rows)            # eta resolved
    assert all(r["w2_mean"] < 0 for r in rows)            # neg. semidef.
    f1 = next(r for r in rows if abs(r["a"]) < 1e-9)
    assert abs(f1["w2_mean"] + 0.00693959) < 1e-6
    assert abs(f1["w2_plus"][0] + 0.013515) < 1e-5        # 37x asymmetry
    assert abs(f1["w2_minus"][0] + 0.000364) < 1e-5


def test_du03_systole(du03_ledger):
    E = du03_ledger["E systole (EXACT)"]
    assert E["trace"] == 4 and E["trace3_excluded"] is True
    assert abs(E["ell"] - 2.633915793849633) < 1e-11


def test_du03_census_and_pheno(du03_ledger):
    F = du03_ledger["F1 census"]
    assert F["hits"] == 2 and F["p"] > 0.05               # seeded rng
    G = du03_ledger["G"]
    assert abs(G["chi_g_s"] / 3.054870e-11 - 1) < 1e-5
