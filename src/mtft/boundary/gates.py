"""Gates for `mtft.boundary`.

Recomputed from the frozen data at call time, not asserted.  The N = 143
results are re-derived in pure Python from the shipped S_4 basis, cusp
functionals and operator matrices; the census is checked for internal
consistency against nu_2 and nu_3 computed from Kronecker symbols.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import combinations_with_replacement

from . import (
    CENSUS, HYPOTHESES, RELATIONS_143, boundary_dim, cusp_functionals,
    nu2, nu3, operator, s4_qexpansions,
)

__all__ = ["verify", "gate_census_consistency", "gate_h0_2k",
           "gate_cusp_functionals", "gate_leakage", "gate_t2_from_qexpansions",
           "gate_hypotheses"]

NCOEF = 140
CUSPS = (1, 11, 13, 143)


def _rank(rows):
    rows = [[Fraction(v) for v in r] for r in rows]
    n = len(rows)
    m = len(rows[0]) if n else 0
    r = 0
    for c in range(m):
        piv = next((i for i in range(r, n) if rows[i][c]), None)
        if piv is None:
            continue
        rows[r], rows[piv] = rows[piv], rows[r]
        pv = rows[r][c]
        rows[r] = [v / pv for v in rows[r]]
        for i in range(n):
            if i != r and rows[i][c]:
                f = rows[i][c]
                rows[i] = [a - f * b for a, b in zip(rows[i], rows[r])]
        r += 1
        if r == n:
            break
    return r


def _solve_in_s4(vec, S4cols):
    """Coordinates of a q-expansion in the S_4 basis (exact, overdetermined)."""
    rows = [list(r) + [v] for r, v in zip(S4cols, vec)]
    n, m = len(rows), len(rows[0]) - 1
    rows = [[Fraction(v) for v in r] for r in rows]
    piv, r = [], 0
    for c in range(m):
        k = next((i for i in range(r, n) if rows[i][c]), None)
        if k is None:
            continue
        rows[r], rows[k] = rows[k], rows[r]
        pv = rows[r][c]
        rows[r] = [v / pv for v in rows[r]]
        for i in range(n):
            if i != r and rows[i][c]:
                f = rows[i][c]
                rows[i] = [a - f * b for a, b in zip(rows[i], rows[r])]
        piv.append(c)
        r += 1
    for i in range(r, n):
        if rows[i][m] != 0:
            raise AssertionError("q-expansion is not in the span of the S_4 basis")
    out = [Fraction(0)] * m
    for i, c in enumerate(piv):
        out[c] = rows[i][m]
    return out


_H0_CACHE = None


def _h0_2k_coords():
    """H^0(2K) as 91 coordinate vectors in the S_4 basis (cached)."""
    global _H0_CACHE
    if _H0_CACHE is not None:
        return _H0_CACHE
    from ..canonical import MONOMIALS, adapted_qexpansions
    E = adapted_qexpansions()
    S4 = s4_qexpansions()
    out = []
    for i, j in MONOMIALS:
        prod = [sum(E[a][i] * E[n - a][j] for a in range(n + 1))
                for n in range(NCOEF + 1)]
        out.append(_solve_in_s4(prod, S4))
    _H0_CACHE = out
    return out


def gate_census_consistency():
    """b_N = 4 + nu_2 + nu_3 for every recorded level, from Kronecker symbols."""
    bad = []
    for N, (p, q, g, bN, n2, n3, good, up, uq) in CENSUS.items():
        if (nu2(N), nu3(N)) != (n2, n3) or boundary_dim(N) != bN:
            bad.append(N)
        if good != (4, n2 + n3):
            bad.append(N)
    return {"levels": len(CENSUS), "mismatches": bad, "ok": not bad}


def gate_h0_2k():
    """dim H^0(2K) = 36 and dim S_4 = 40, from shipped q-expansions."""
    H = _h0_2k_coords()
    d = _rank([list(v) for v in H])
    return {"dim_H0_2K": d, "dim_S4": len(s4_qexpansions()[0]), "ok": d == 36}


def gate_cusp_functionals():
    """All four c_d vanish identically on H^0(2K)."""
    H, C = _h0_2k_coords(), cusp_functionals()
    worst = max(abs(sum(c[k] * v[k] for k in range(len(v))))
                for c in C for v in H)
    return {"max_value_on_H0": str(worst), "ok": worst == 0}


def _leak_rows(op_name):
    H, C = _h0_2k_coords(), cusp_functionals()
    A = operator(op_name)
    n = len(A)
    AH = [[sum(A[r][k] * v[k] for k in range(n)) for r in range(n)] for v in H]
    return [[sum(c[k] * w[k] for k in range(n)) for w in AH] for c in C]


def gate_leakage():
    """Ranks 4, 2, 1 and the exact vanishing/extra relations at N = 143."""
    out, ok = {}, True
    for name, key in (("T2", "T_l (l = 2, 3, 5, 7)"), ("U11", "U_11"), ("U13", "U_13")):
        rows = _leak_rows(name)
        spec = RELATIONS_143[key]
        r = _rank(rows)
        zero = tuple(CUSPS[i] for i in range(4) if all(v == 0 for v in rows[i]))
        extra = None
        if all(a == b for a, b in zip(rows[0], rows[1])):
            extra = "c_1 = c_11"
        good = (r == spec["rank"] and zero == spec["zero"]
                and extra == spec["extra"] and any(v for row in rows for v in row))
        ok &= good
        out[name] = {"rank": r, "zero_cusps": zero, "extra": extra, "ok": good}
    out["ok"] = ok
    return out


def gate_t2_from_qexpansions():
    """Cross-check: rebuild T_2 from the shipped S_4 q-expansions.

    (T_2 F)_n = a_{2n} + 2^3 a_{n/2} for weight 4, needing coefficients to
    2 * 56 = 112 <= 140.  Compares against the shipped operator matrix, so
    a corrupted matrix cannot pass.
    """
    S4 = s4_qexpansions()
    d = len(S4[0])
    sturm = 56
    A = operator("T2")
    worst = 0
    for k in range(d):
        col = [S4[n][k] for n in range(NCOEF + 1)]
        img = [col[2 * n] + (8 * col[n // 2] if n % 2 == 0 else 0)
               for n in range(sturm + 1)]
        rebuilt = [sum(A[r][k] * S4[n][r] for r in range(d))
                   for n in range(sturm + 1)]
        worst = max(worst, max(abs(a - b) for a, b in zip(img, rebuilt)))
    return {"max_discrepancy": str(worst), "sturm": sturm, "ok": worst == 0}


def gate_hypotheses():
    """The recorded verdicts match what the census table actually says."""
    q1 = all(v[6] == (4, v[4] + v[5]) for v in CENSUS.values())
    q2 = all(v[8][0] == 1 for v in CENSUS.values())
    q3 = all((v[7][1] > 0 and v[8][1] == 0)
             for v in CENSUS.values() if v[4] + v[5] > 0)
    q4 = all(v[7][1] <= 1 and v[8][1] <= 1 for v in CENSUS.values())
    got = {"Q1": q1, "Q2": q2, "Q3": q3, "Q4": q4}
    want = {"Q1": True, "Q2": False, "Q3": True, "Q4": False}
    return {"computed": got, "recorded": want, "ok": got == want,
            "note": "Q2 and Q4 must come out FALSE -- they are falsified, "
                    "and a gate that let them pass would mean the census "
                    "table had been edited"}


def verify(full=True):
    out = {
        "census_consistency": gate_census_consistency(),
        "hypotheses": gate_hypotheses(),
    }
    if full:
        out["h0_2k"] = gate_h0_2k()
        out["cusp_functionals"] = gate_cusp_functionals()
        out["leakage"] = gate_leakage()
        out["t2_crosscheck"] = gate_t2_from_qexpansions()
    out["ok"] = all(v["ok"] for v in out.values() if isinstance(v, dict))
    return out
