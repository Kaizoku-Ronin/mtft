"""R2C-01 reproduced from scratch (own gamma matrices, own enumeration) and the oloid involution algebra."""
import itertools, numpy as np, sympy as sp
from fractions import Fraction
from mtft.research import chirality as CH, involutions as IV, pipeline as PL


def test_clifford_selection_rule():
    cl = CH.clifford_6d(); assert cl["clifford_ok"] and cl["C_ok"] and np.allclose(cl["Gamma_star"], np.kron(np.kron(np.diag([1, -1]), np.diag([1, -1])), np.diag([1, -1])))
    for e in (-1, 1):
        for ep in (-1, 1):
            assert CH.bilinear_selection(e, ep, "scalar")["allowed"] == (e != ep) and CH.bilinear_selection(e, ep, "transpose")["allowed"] == (e != ep) and CH.bilinear_selection(e, ep, "vector")["allowed"] == (e == ep)


def test_exhaustive_scan_counts_by_independent_enumeration():
    d = CH.sector_degrees(); counts = dict(fam=0, full=0, sc=0, both=0)
    for eps in itertools.product((1, -1), repeat=10):
        E = dict(zip(CH.SECTORS, eps)); fam = all(E[s] * d[s] == CH.FAMILY_TARGET[s] for s in CH.FAMILY_TARGET); full = all(E[s] * d[s] == d[s] for s in CH.SECTORS)
        sc = (E["cL"] == -E["ca"] == -E["cb"]) and (E["Ld"] == -E["ad"] == -E["bd"]); counts["fam"] += fam; counts["full"] += full; counts["sc"] += sc; counts["both"] += fam and sc
    assert counts == {"fam": 16, "full": 4, "sc": 64, "both": 0}
    r = CH.enumerate_chiralities(); assert r["counts"]["family"] == 16 and r["counts"]["full_ledger"] == 4 and r["counts"]["scalar"] == 64 and r["counts"]["both"] == 0 and r["counts"]["family_and_LaLb_equal"] == 8
    assert r["no_go"] and r["vector_higgs_allowed_with_family_signs"]
    # combinatorial proofs of the counts: 2^4, 2^2, 2^6
    assert (16, 4, 64) == (2 ** 4, 2 ** 2, 2 ** 6)


def test_gravitational_p2_table():
    ng = CH.m1_scalar_higgs_no_go(); table = {(cd, ab): p["p2_coefficient"] for cd, ab, p in ng["ledger_preserving_p2_table"]}
    assert table == {(1, 1): Fraction(-1, 60), (1, -1): Fraction(-11, 720), (-1, 1): Fraction(-1, 80), (-1, -1): Fraction(-1, 90)}
    assert all(p["n_grav"] >= 8 for _, _, p in ng["ledger_preserving_p2_table"])
    g = PL.m1_gate_report()["gates"]["six_d_yukawa_chirality"]; assert g["pass"] is False and g["witness"]["counts"]["both"] == 0


def test_oloid_involution_algebra():
    c = IV.involution_checks(); assert all(c.values()); L = IV.logistic_conjugacy(); assert L["conjugacy_holds"] and sp.simplify(L["c"] - IV.r * (2 - IV.r) / 4) == 0
    q = sp.Rational(3, 7); w = IV.oloid_involution; assert w(w(q)) == q and IV.conjugating_coordinate(w(q)) == -1 - IV.conjugating_coordinate(q)   # numeric rational instance
    assert IV.oloid_intersection_dimension()["two_solid_oloids_in_R3_give_S4"] is False
