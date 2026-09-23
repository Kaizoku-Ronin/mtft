"""R2C-03 (2026-09-19): tensor integrality for the p2 term, and the colour-cubic (2-form x 6-form) obstruction of M1's declared content.

p2 coefficients of I8 per field (CC-33, 2026-09-23; v0.33.0): complex Weyl (signed dimension) -1/1440; self-dual tensor -7/360
(the -L8/8 genus: L8 = (7 p2 - p1^2)/45 — the +1/360 quoted in R2C-03 is its p1^2 coefficient, not its p2 coefficient); gravitino
-49/288.  In Weyl units the three are 1 : 28 : 245, which is the 273 = 245 + 28 and 29 = 28 + 1 of the (1,0) count H - V + 29 T = 273
and the 1/48 = (28 + 2)/1440 of the (2,0) tensor multiplet.  Cancelling the p2 term of a fermion spectrum with n_grav signed
dimensions needs n_grav/28 net anti-self-dual tensors: an integer only if n_grav = 0 mod 28.  For M1's four ledger-preserving
assignments n_grav = 24, 22, 18, 16 NONE survives (R2C-03's survivors 24 -> 6 and 16 -> 4 are withdrawn); a gravitino cannot be
balanced by tensors alone (269/28).  `tensor_coefficient_witnesses` re-derives the constant on three independent routes.
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

P2_WEYL = Fraction(-1, 1440); P2_SELF_DUAL_TENSOR = Fraction(-7, 360); P2_GRAVITINO = Fraction(-49, 288)      # CC-33: tensor was +1/360
P1SQ_SELF_DUAL_TENSOR = Fraction(1, 360)                        # the coefficient that was mistaken for the p2 one (R2C-03)
P2_WEYL_UNITS = {"weyl": 1, "self_dual_tensor": 28, "gravitino": 245}

def tensor_coefficient_witnesses():
    """Three independent derivations of the self-dual-tensor p2 coefficient (CC-33), each an exact rational identity:
    (i) the genera: A-roof_8 = (7 p1^2 - 4 p2)/5760 gives the Weyl -1/1440, -L8/8 with L8 = (7 p2 - p1^2)/45 gives the tensor -7/360;
    (ii) the (1,0) supergravity count H - V + 29 T = 273: gravity multiplet 273 = 245 (gravitino) + 28 (self-dual tensor), tensor multiplet 29 = 28 + 1;
    (iii) the (2,0) tensor multiplet (one self-dual tensor, two Weyl): p2 magnitude 1/48 = (28 + 2)/1440."""
    p1sq, p2 = sp.symbols("p1sq p2")
    a_roof8 = (7 * p1sq - 4 * p2) / 5760; L8 = (7 * p2 - p1sq) / 45
    weyl = Fraction(str(sp.Rational(a_roof8.coeff(p2)))); tensor = Fraction(str(sp.Rational((-L8 / 8).coeff(p2)))); tensor_p1sq = Fraction(str(sp.Rational((-L8 / 8).coeff(p1sq))))
    units = {"weyl": 1, "self_dual_tensor": tensor / weyl, "gravitino": P2_GRAVITINO / weyl}
    return {"weyl": weyl, "self_dual_tensor": tensor, "tensor_p1sq_coefficient": tensor_p1sq, "weyl_units": units,
            "route_i_genera": weyl == P2_WEYL and tensor == P2_SELF_DUAL_TENSOR and tensor_p1sq == P1SQ_SELF_DUAL_TENSOR,
            "route_ii_273": units["gravitino"] + units["self_dual_tensor"] == 273 and units["self_dual_tensor"] + units["weyl"] == 29,
            "route_iii_2_0_multiplet": (units["self_dual_tensor"] + 2) * (-P2_WEYL) == Fraction(1, 48)}

def tensor_integrality(n_grav, gravitino=False):
    """Net self-dual tensors needed to cancel the p2 term of n_grav signed Weyl dimensions (negative = anti-self-dual): -n_grav/28.
    Integral iff n_grav = 0 mod 28 (CC-33; the earlier n_grav/4 used the p1^2 coefficient of the tensor genus)."""
    p2 = n_grav * P2_WEYL + (P2_GRAVITINO if gravitino else 0); need = -p2 / P2_SELF_DUAL_TENSOR
    return {"n_grav": n_grav, "p2_fermions": n_grav * P2_WEYL, "net_self_dual_tensors_needed": need, "anti_self_dual_tensors_needed": -need,
            "integral": need.denominator == 1, "rule": "n_grav = 0 mod 28 (CC-33)"}

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
