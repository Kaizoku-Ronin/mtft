#!/usr/bin/env python3
"""
mtft.ledger_peel — v0.12.0 addendum entries (peel wave), Entry-schema.
Merge into mtft.ledger at the next consolidation; values with closed forms
are generated at import (self-verifying, rule 4).
"""
from __future__ import annotations
import mpmath as mp
from mpmath import mpf, log, exp, pi
from .ledger import Entry

mp.mp.dps = 50
_z = mp.zeta
def _s(x): return mp.nstr(x, 30)

LN2, LN3 = log(mpf(2)), log(mpf(3))
_minf = log(27*LN2/(8*LN3))

ENTRIES = {
 "peel.bulk_pole_s": Entry("3", "EXACT",
    "studies/peel_2026aug/weights_index.py", "this drop",
    form="rightmost pole of F(s) = -zeta(s-2) zeta'(s-1)",
    supersedes=(("2", "dictionary F(s) = -zeta(s-1) zeta'(s)"),),
    comment="index correction; residue at s=3 is T_inf = -zeta'(2)"),
 "peel.skeleton_C0": Entry(_s(-mp.diff(_z, mpf(-1))/_z(mpf(-1))), "EXACT",
    "studies/peel_2026aug/skeleton_zeros.py", "this drop",
    tol=1e-14, form="-zeta'(-1)/zeta(-1) = 1 - 12 ln A (Glaisher)"),
 "peel.bulk_y_coeff": Entry(_s(-_z(mpf(3))/(240*pi)), "EXACT",
    "studies/peel_2026aug/mu_expansion.py", "this drop",
    form="-zeta(3)/(240 pi): first odd term of mu_bulk"),
 "peel.bulk_y3_coeff": Entry(_s(-_z(mpf(5))/(252*pi)), "EXACT",
    "studies/peel_2026aug/mu_expansion.py", "this drop",
    form="-zeta(5)/(252 pi)"),
 "marked.m_inf_identity": Entry(_s(_minf), "EXACT",
    "studies/peel_2026aug/rung4_coset.py", "this drop",
    form="ln(27 ln2/(8 ln3)) = eps_3(3) - eps_2(3), marked levels at beta_c=3",
    comment="unifies with T_inf: residue and spacing of one spectral point; "
            "sigma=3 site-measure derivation flagged for PR-5.1 audit"),
 "marked.sigma_star": Entry(_s(log(LN3/LN2)/log(mpf(3)/2)), "DERIVED",
    "studies/peel_2026aug/rung4_coset.py", "this drop",
    form="ln(log2 3)/ln(3/2): gap-closing level crossing"),
 "marked.R_inf": Entry(_s(_minf/(1-exp(-_minf))), "CERTIFIED",
    "studies/peel_2026aug/rung4_coset.py", "this drop", tol=1e-11,
    form="m_inf/(1-e^-m_inf); matches ledger EP value 1.42507723746"),
 "gl2.conductor": Entry("143", "EXACT",
    "studies/peel_2026aug/gl2_setup.py", "this drop",
    form="Delta = -1859 = -11*13^2, c4 = 64 coprime to both"),
 "gl2.root_number": Entry("-1", "CERTIFIED",
    "studies/peel_2026aug/gl2_setup.py", "this drop",
    comment="two-route s=4 test: 10.5 digits (eps=-1) vs 0.7 (eps=+1)"),
 "gl2.Lprime1": Entry("0.9456964112", "CERTIFIED",
    "studies/peel_2026aug/gl2_setup.py", "this drop", tol=1e-9,
    form="L'(143a1, 1); nonzero -> analytic rank exactly 1"),
 "gl2.rank_read": Entry("1", "CERTIFIED",
    "studies/peel_2026aug/gl2_peel.py", "this drop", tol=7e-7,
    form="-X (S_f - smooth - zero osc) at five depths",
    comment="third independent route to rank(143a1) = 1"),
 "gl2.first_zero": Entry("3.2930459", "CERTIFIED",
    "studies/peel_2026aug/gl2_peel.py", "this drop", tol=1e-6,
    form="lowest gamma of L(143a1, s) on Re s = 1"),
}

def verify(tol_digits=12):
    """Recompute closed-form entries; return list of (key, ok)."""
    out = []
    checks = {
     "peel.skeleton_C0": lambda: -mp.diff(_z, mpf(-1))/_z(mpf(-1)),
     "peel.bulk_y_coeff": lambda: -_z(mpf(3))/(240*pi),
     "peel.bulk_y3_coeff": lambda: -_z(mpf(5))/(252*pi),
     "marked.m_inf_identity": lambda: 3*log(mpf(3)/2) - log(LN3/LN2),
     "marked.R_inf": lambda: _minf/(1-exp(-_minf)),
    }
    for k, f in checks.items():
        ok = abs(f() - ENTRIES[k].mpf()) < mpf(10)**(-tol_digits)
        out.append((k, bool(ok)))
    return out
