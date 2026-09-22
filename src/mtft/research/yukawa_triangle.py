"""R2C-06 (2026-09-19): the telescoping theorem for gauge-vertex Yukawas, and the light-Higgs no-go of the class on a curve.

A Yukawa from the 6D gauge vertex psi-bar Gamma^zbar A_zbar psi couples two family blocks (charges q1, q2) to the Higgs block A with
q_H = -(q1 + q2).  Charges are linear in the line-bundle degrees, so d_H = -(d1 + d2): the flux that makes the families chiral fixes the
Higgs block's degree, and by R2C-05 its tachyon mass -|d_H|/24.  Two left-handed families of index 3 give d_H = -6: m_H^2 = -1/4 in
curvature units, and the Higgs block carries index -6 (M1's six extra doublet pairs).  A flux-neutral Higgs (d_H = 0) would need
d1 = -d2 — two 'families' of opposite chirality, never a Standard-Model pair.  Because slope = degree on a curve, chirality and Higgs
mass are the same number: no vector-Higgs parent on a curve has a light Higgs at leading order; with R2C-01 (no elementary-scalar
Higgs for M1's families) and H-19 (a hyper-scalar Higgs in a charged block is massive at +|d|/24), the class has no leading-order
electroweak hierarchy.  Also: E7 -> E6 x U(1) has blocks of charge 0, +-1 only, so two 27_{+1} families have no gauge-vertex Higgs
(charge -2 absent): M3 fails the Yukawa gate.  Escapes: an internal space of complex dimension >= 2, where slope and chiral index are
independent (line bundles with mu(L) = 0 and c1(L)^2 != 0 — the fibred constructions), or non-gauge Yukawas in a non-supersymmetric
charged parent, subject to an R2C-01-type chirality scan."""
import sympy as sp

def triangle_closure(d1, d2): return {"d_H": -(d1 + d2), "higgs_tachyon_m2": -sp.Rational(abs(d1 + d2), 24), "higgs_block_index": -(d1 + d2), "flux_neutral_possible": d1 + d2 == 0}

def forced_higgs_hypercharge(Y1, Y2): return -(sp.Rational(Y1) + sp.Rational(Y2))

def m1_triangles():
    m = {"c": 0, "L": -3, "a": 3, "b": 3, "d": 0}; d = lambda i, j: m[i] - m[j]
    tri = {"up": (d("c", "L"), d("a", "c"), d("L", "a")), "down": (d("c", "L"), d("b", "c"), d("L", "b")), "lepton": (d("d", "L"), d("b", "d"), d("L", "b")), "neutrino": (d("d", "L"), d("a", "d"), d("L", "a"))}
    Y = {"Q": sp.Rational(1, 6), "uc": sp.Rational(-2, 3), "dc": sp.Rational(1, 3), "L": sp.Rational(-1, 2), "ec": 1, "nc": 0}
    return {"degrees": tri, "closed": all(sum(v) == 0 for v in tri.values()), "higgs_degree": -6, "higgs_hypercharges": {"H_u": forced_higgs_hypercharge(Y["Q"], Y["uc"]), "H_d": forced_higgs_hypercharge(Y["Q"], Y["dc"]), "H_d(lep)": forced_higgs_hypercharge(Y["L"], Y["ec"]), "H_u(nu)": forced_higgs_hypercharge(Y["L"], Y["nc"])},
            "extra_doublet_pairs": 6}

def e7_triangle_test(block_charges=(0, 1, -1), family_charge=1):
    need = -2 * family_charge; return {"needed_higgs_charge": need, "present": need in block_charges, "yukawa_gate": "FAIL" if need not in block_charges else "PASS"}

def light_higgs_no_go():
    return {"theorem": "d_H = -(d_1 + d_2) for every gauge-vertex Yukawa triangle; on a curve slope = degree, so chiral families of index 3 force m_H^2 = -1/4",
            "elementary_scalar_higgs": "excluded for M1's families (R2C-01)", "hyper_scalar_higgs": "massive at +|d|/24 in a charged block (H-19)", "vector_higgs": "tachyonic at -|d|/24 (R2C-05)",
            "conclusion": "no leading-order electroweak hierarchy in the class (line-bundle flux on a curve, families as chiral zero modes)",
            "escapes": ["internal complex dimension >= 2: slope and chiral index independent (fibred constructions over X0(143))", "non-gauge Yukawas in a non-supersymmetric charged parent (R2C-01-type scan required)"]}
