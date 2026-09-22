"""R2C-03 (2026-09-19): tensor integrality for the p2 term, and the colour-cubic (2-form x 6-form) obstruction of M1's declared content.

p2 coefficients of I8 per field: complex Weyl (signed dimension) -1/1440; self-dual tensor +1/360; gravitino -49/288.  Cancelling the
irreducible p2 term of a fermion spectrum with n_grav signed dimensions needs n_T - n_T' = n_grav/4 net self-dual tensors: an integer
only if n_grav = 0 mod 4.  For M1's four ledger-preserving assignments: n_grav = 24 (6 tensors) and 16 (4) survive; 22 and 18 are
excluded; a gravitino cannot be balanced by tensors alone (269/4).
Colour-cubic ledger: dI8/dD3 = sum over colour triplet sectors of eps * A(rep) * q * mult/2 (F32: D3 q, F3: D3 q/2; q = f_i - f_j).  The f_L
coefficient is sourced only by the Q sector (the only coloured field with U(1)_L charge among Hom-type bifundamentals (i, j-bar)); it
equals -eps_cL = -1 for family preservation, and it is a 2-form x 6-form term, hence not cancellable by products of 4-forms.  Escape
routes: non-Hom sectors (N_i, N_j) with eps = +1 cancel a U(1)_x term at the cost of three vector-like coloured pairs; or a parent
whose simple factor has no cubic Casimir (R2C-04).  CC-29 (2026-09-19): the earlier "escape route (b)" — an L-stack without its own U(1)
and the flux (3, 0, 0, 0, -3) — was WRONG: that vector flips u^c, d^c, L (d_ca = +3, needs -3); the family-preserving flux without U(1)_L is
(3, 0, 6, 6, 3) = M1 + 3, all sector degrees (Higgs -6 included) unchanged, and removing U(1)_L cures only the f_L term while the u^c and d^c
sectors leave -f_a/2 and -f_b/2.  General statement: EVERY U(1)_x paired with colour in a family sector carries an irreducible
U(1)_x–SU(3)^3 term; since hypercharge needs U(1)_{a,b,d}, no Hom-type unitary-stack flux model with SM bifundamental families is a
consistent 6D gauge theory."""
from fractions import Fraction
import sympy as sp

P2_WEYL = Fraction(-1, 1440); P2_SELF_DUAL_TENSOR = Fraction(1, 360); P2_GRAVITINO = Fraction(-49, 288)

def tensor_integrality(n_grav, gravitino=False):
    p2 = n_grav * P2_WEYL + (P2_GRAVITINO if gravitino else 0); need = -p2 / P2_SELF_DUAL_TENSOR
    return {"n_grav": n_grav, "p2_fermions": n_grav * P2_WEYL, "net_self_dual_tensors_needed": need, "integral": need.denominator == 1}

def m1_tensor_survivors():
    out = []
    for e_cd in (1, -1):
        for e_ab in (1, -1):
            n = 16 + 3 * e_cd + 2 + 2 + e_ab; t = tensor_integrality(n); out.append({"eps_cd": e_cd, "eps_ab": e_ab, **t})
    return {"table": out, "survivors": [(r["eps_cd"], r["eps_ab"], r["net_self_dual_tensors_needed"]) for r in out if r["integral"]], "gravitino_possible": tensor_integrality(24, gravitino=True)["integral"]}

F = dict(zip("cLabd", sp.symbols("f_c f_L f_a f_b f_d")))

def colour_cubic_ledger(sectors):
    """sectors: list of (charge_expr in the f's, eps, weight) with weight = A(rep) * (SU(2) multiplicity)/2: Q-type (3,2): 1; (3,1): 1/2."""
    return sp.expand(sum(sp.nsimplify(w) * e * q for q, e, w in sectors))

def m1_colour_cubic(eps=None, extra=()):
    eps = eps or {"cL": 1, "ca": 1, "cb": 1, "cd": 1}
    sectors = [(F["c"] - F["L"], eps["cL"], 1), (F["c"] - F["a"], eps["ca"], sp.Rational(1, 2)), (F["c"] - F["b"], eps["cb"], sp.Rational(1, 2)), (F["c"] - F["d"], eps["cd"], sp.Rational(1, 2))] + list(extra)
    d = colour_cubic_ledger(sectors); return {"dI8_dD3": d, "coefficients": {k: d.coeff(F[k]) for k in "cLabd"}, "f_L_irreducible": d.coeff(F["L"]) != 0, "term_type": "2-form x 6-form: not a product of 4-forms (no Green–Schwarz cancellation)"}

def hom_type_obstruction():
    base = m1_colour_cubic(); withcd = m1_colour_cubic(extra=[(F["c"] - F["d"], -1, sp.Rational(3, 2))])
    return {"f_L_coefficient": base["coefficients"]["L"], "unchanged_by_extra_cd": withcd["coefficients"]["L"] == base["coefficients"]["L"],
            "statement": "within Hom-type bifundamentals (i, j-bar) the f_L D3 term is sourced by Q alone and is irreducible; M1's declared content has no anomaly-free 6D completion preserving the families",
            "escape_routes": ["non-Hom sectors (N_i, N_j) with eps = +1 (three vector-like coloured pairs per cancellation)", "a parent whose simple factor has no cubic Casimir, with the families as one complex block of U(1)-flux degree 3 (R2C-04)"]}

def family_preserving_flux(m_L=0):
    """Solve the six R2C-01 family targets d_cL = 3, d_ca = -3, d_cb = -3, d_Ld = -3, d_ad = 3, d_bd = 3 for the flux vector:
    (m_c, m_L, m_a, m_b, m_d) = (m_L + 3, m_L, m_L + 6, m_L + 6, m_L + 3) — M1 up to a common shift; the Higgs degrees d_La = d_Lb = -6 are
    shift-invariant.  CC-29: the vector (3, 0, 0, 0, -3) once quoted as an escape route is not family-preserving."""
    m = {"c": m_L + 3, "L": m_L, "a": m_L + 6, "b": m_L + 6, "d": m_L + 3}
    d = {"cL": m["c"] - m["L"], "ca": m["c"] - m["a"], "cb": m["c"] - m["b"], "Ld": m["L"] - m["d"], "ad": m["a"] - m["d"], "bd": m["b"] - m["d"], "La": m["L"] - m["a"], "Lb": m["L"] - m["b"]}
    return {"m": m, "sector_degrees": d, "family_preserved": [d[k] for k in ("cL", "ca", "cb", "Ld", "ad", "bd")] == [3, -3, -3, -3, 3, 3], "higgs_degrees": (d["La"], d["Lb"])}

def colour_cubic_without_U1L():
    """With SU(2)_L (no U(1)_L) the Q sector carries only f_c: the f_L term disappears but -f_a/2 (u^c) and -f_b/2 (d^c) remain."""
    d = colour_cubic_ledger([(F["c"], 1, 1), (F["c"] - F["a"], 1, sp.Rational(1, 2)), (F["c"] - F["b"], 1, sp.Rational(1, 2)), (F["c"] - F["d"], 1, sp.Rational(1, 2))])
    return {"dI8_dD3": d, "irreducible_terms": {k: d.coeff(F[k]) for k in "ab" if d.coeff(F[k]) != 0}, "escape_b_works": False}
