"""PS-01: Pati–Salam breaking, anomalies, the gamma direction and the b–tau check."""
from fractions import Fraction
import sympy as sp
from mtft.matter import ps_breaking as P


def test_branching_matches_pati_salam_and_pv02():
    b = P.branching()
    assert b["pv02_consistent"]
    assert b["(4,2,1)"] == {"Q (c<-L)": Fraction(1, 6), "l (A1<-L)": Fraction(-1, 2)}
    assert set(b["(4bar,1,2)"].values()) == {Fraction(1, 3), Fraction(-2, 3), Fraction(1), Fraction(0)}
    assert set(b["(1,2,2)"].values()) == {Fraction(1, 2), Fraction(-1, 2)}


def test_anomalies_nonabelian_cancel_and_only_the_decoupled_trace_u1_is_free():
    A = P.anomalies()
    assert A["SU4^3"] == 0 and A["Witten_SU2L_doublets"] % 2 == 0 and A["Witten_SU2R_doublets"] % 2 == 0 and A["SU4^2 Q4"] == 0
    assert (A["SU2L^2 QL"], A["SU2R^2 QR"], A["SU2L^2 Q4"]) == (-12, 12, 6)
    x, y, z = sp.symbols("x y z")
    eqs = [sum(sp.Rational(str(A[f"{m} {U}"])) * s for U, s in (("Q4", x), ("QL", y), ("QR", z))) for m in ("SU4^2", "SU2L^2", "SU2R^2", "grav^2")]
    assert sp.linsolve(eqs, [x, y, z]) == sp.FiniteSet((z, z, z))                  # only Q4 + QL + QR, which no field carries


def test_wilson_lines_keep_rank_and_the_higgs_pair_needs_NL_equal_NB():
    w = P.wilson_lines({"c": "O", "A1": "N1", "L": "N2", "B": "N2", "A2": "N3"})
    assert w["unbroken"] == ["SU(3)", "U(1)_{B-L}", "SU(2)_L", "U(1)_{3R}"] and w["higgs_pair_survives"] and not w["reaches_SM"]
    assert not P.wilson_lines({"c": "O", "A1": "N1", "L": "N2", "B": "N3", "A2": "N4"})["higgs_pair_survives"]


def test_fi_vacua_delta_gives_the_bidoublet_and_kappa_gives_the_ps_higgs():
    d = P.single_block_vacua(P.fi_terms(delta=0.1))
    assert [v["rep"] for v in d["vacua"]] == ["(1,2,2)"] and abs(d["group_FI"]["4"]) < 1e-12
    k = P.single_block_vacua(P.fi_terms(kappa=-0.01, gamma_charge=P.PS_GAMMA_FLUX))
    assert [v["rep"] for v in k["vacua"]] == ["(4,1,2)/(4bar,1,2)"] and abs(k["group_FI"]["L"]) < 1e-12
    assert sum(P.RANK[s] * P.PS_GAMMA_FLUX[s] for s in P.RANK) == 0                 # traceless gamma flux


def test_gamma_direction_facts():
    g = P.gamma_facts()
    assert g["gamma_sq"] == -8 and g["K_dot_gamma"] == 0 and g["chi_formula_holds"]
    assert g["chi_table"][(1, 0)] == -4 and g["chi_table"][(2, 0)] == -16          # 4c^2 localised modes
    assert g["E_Q_torsion_order_divides"] == 1 and g["cusps_fill_fibre_over_O"]
    from mtft.surface.blowup import net_chirality
    assert all(net_chirality((a, 0, c, e)) == net_chirality((a, 0, 0, e)) for a in (-1, 0, 1) for c in (-2, 1, 2) for e in (-3, 3))


def test_b_tau_against_data():
    bt = P.b_tau_running()
    assert 9e5 < bt["mu_equal"] < 1.4e6 and abs(bt["R_mZ"] - 1.6377) < 1e-3
    assert 5e5 < P.b_tau_running(m_b=2.79)["mu_equal"] < 7e5 and 2e6 < P.b_tau_running(m_b=2.93)["mu_equal"] < 2.6e6
    assert 0.30 < bt["ms_over_mmu_at_mu"] < 0.36 and 3.1 < bt["md_over_me_at_mu"] < 3.7         # the Georgi–Jarlskog pattern, not equality
    assert 0.02 < P.seesaw_estimate() < 0.04                                                   # GeV: far above the 0.45 eV bound
