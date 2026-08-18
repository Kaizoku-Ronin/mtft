"""Gates for `mtft.canonical`.

Every number the canonical-ideal arc claimed is recomputed here from the
frozen q-expansions, in pure Python, at call time.  Nothing is asserted
from memory: if the shipped data were altered, these fail.

Exact-rank policy follows the arc: integer elimination mod two distinct
large primes.  rank_Q >= rank_p always, and each target is capped above by
theory, so agreement at both primes pins the rational rank.  Residual
checks are exact integer arithmetic against zero, not a tolerance.
"""

from __future__ import annotations

from itertools import combinations_with_replacement

from . import (
    COORDINATE_LABELS,
    DESCENT,
    GENUS,
    MONOMIALS,
    PREDICTIONS,
    SECTOR_ORDER,
    adapted_qexpansions,
    ci_a_quadric,
    ideal_basis,
    ideal_by_sector,
    monomial_sector,
    s2_qexpansions,
)

__all__ = [
    "PRIMES", "NCOEF", "ROUTE2_PRIMES_EXPECTED", "verify",
    "gate_petri", "gate_generation", "gate_sector_grading", "gate_bundles",
    "gate_projection", "gate_descent", "gate_route2", "gate_ci_a",
]

#: Two primes for the rank certificates.
PRIMES = (2147483647, 1000003)
#: q-expansions are shipped to q^140; the weight-4 Sturm bound is 56.
NCOEF = 140

_SECTOR_COORDS = {
    s: [k for k, (_, sec) in enumerate(COORDINATE_LABELS) if sec == s]
    for s in SECTOR_ORDER
}


# ---------------------------------------------------------------- helpers

def _product(E, i, j):
    """Coefficient vector of e_i * e_j, q^0..q^NCOEF."""
    return [sum(E[a][i] * E[n - a][j] for a in range(n + 1))
            for n in range(NCOEF + 1)]


def _rank_mod_p(cols, p):
    """Rank over F_p of the matrix whose columns are `cols`."""
    rows = [list(r) for r in zip(*cols)] if cols else []
    rows = [[v % p for v in r] for r in rows]
    n, m = len(rows), (len(rows[0]) if rows else 0)
    rank = 0
    for c in range(m):
        piv = next((i for i in range(rank, n) if rows[i][c]), None)
        if piv is None:
            continue
        rows[rank], rows[piv] = rows[piv], rows[rank]
        inv = pow(rows[rank][c], p - 2, p)
        rows[rank] = [v * inv % p for v in rows[rank]]
        for i in range(n):
            if i != rank and rows[i][c]:
                f = rows[i][c]
                pr = rows[rank]
                rows[i] = [(a - f * b) % p for a, b in zip(rows[i], pr)]
        rank += 1
        if rank == n:
            break
    return rank


def _rank(cols):
    """Rank agreed by both primes, else raise."""
    rs = {_rank_mod_p(cols, p) for p in PRIMES}
    if len(rs) != 1:
        raise AssertionError(f"rank disagreed across primes: {rs}")
    return rs.pop()


def _legendre(a, p):
    a %= p
    if a == 0:
        return 0
    return 1 if pow(a, (p - 1) // 2, p) == 1 else -1


def _eta_11a1(n_terms=NCOEF):
    """g = eta(t)^2 eta(11t)^2 = the 11a1 newform, exact integer series."""
    N = n_terms + 1
    prod = [0] * N
    prod[0] = 1
    for shift in (1, 11):
        for _ in range(2):                       # squared
            for d in range(shift, N, shift):     # multiply by (1 - q^d)
                for k in range(N - 1, d - 1, -1):
                    prod[k] -= prod[k - d]
    return [0] + prod[:N - 1]                    # overall factor q^1


def _i2_in_y():
    """I_2 as 55 dicts {monomial index: coeff}, adapted coordinates."""
    cols = []
    for s in SECTOR_ORDER:
        cols.extend(ideal_by_sector()[s][1])
    return cols


# ------------------------------------------------------------------ gates

def gate_petri():
    """P1-P4: dim S_2, Sym^2, h^0(2K), dim I_2; and residual of all 55."""
    F = s2_qexpansions()
    out = {"P1_dim_S2": len(F[0]), "P2_dim_Sym2": len(MONOMIALS)}
    prods = [_product(F, i, j) for i, j in MONOMIALS]
    out["P3_h0_2K"] = _rank(prods)
    out["P4_dim_I2"] = len(MONOMIALS) - out["P3_h0_2K"]

    I2 = ideal_basis()
    worst = 0
    for c in range(55):
        for n in range(NCOEF + 1):
            worst = max(worst, abs(sum(I2[m][c] * prods[m][n]
                                       for m in range(91))))
    out["max_residual"] = worst
    out["ok"] = (out["P1_dim_S2"] == PREDICTIONS["P1_dim_S2"]
                 and out["P2_dim_Sym2"] == PREDICTIONS["P2_dim_Sym2"]
                 and out["P3_h0_2K"] == PREDICTIONS["P3_h0_2K"]
                 and out["P4_dim_I2"] == PREDICTIONS["P4_dim_I2"]
                 and worst == 0)
    return out


def gate_generation():
    """P5-P8: the 55 quadrics generate I_3 (Petri; excludes trigonal)."""
    triples = list(combinations_with_replacement(range(GENUS), 3))
    tix = {t: i for i, t in enumerate(triples)}
    I2 = ideal_basis()
    cubics = []
    for k in range(GENUS):
        for c in range(55):
            row = [0] * len(triples)
            for m, (i, j) in enumerate(MONOMIALS):
                v = I2[m][c]
                if v:
                    row[tix[tuple(sorted((i, j, k)))]] += v
            cubics.append(row)
    r = _rank(cubics)
    out = {"P5_dim_Sym3": len(triples), "P8_rank_V_I2": r,
           "P6_h0_3K": len(triples) - r, "P7_dim_I3": r,
           "n_cubics": len(cubics)}
    out["ok"] = (out["P5_dim_Sym3"] == PREDICTIONS["P5_dim_Sym3"]
                 and r == PREDICTIONS["P8_rank_V_I2"]
                 and out["P6_h0_3K"] == PREDICTIONS["P6_h0_3K"])
    return out


def gate_sector_grading():
    """26/5/4/20, support confinement, and residual 0 in adapted coords."""
    E = adapted_qexpansions()
    by = ideal_by_sector()
    dims, mons, confined = {}, {}, True
    for s in SECTOR_ORDER:
        idx, cols = by[s]
        dims[s] = len(cols)
        mons[s] = len(idx)
        if any(monomial_sector(m) != s for m in idx):
            confined = False
    worst = 0
    cache = {}
    for s in SECTOR_ORDER:
        for col in by[s][1]:
            for m in col:
                if m not in cache:
                    cache[m] = _product(E, *MONOMIALS[m])
            for n in range(NCOEF + 1):
                worst = max(worst, abs(sum(v * cache[m][n]
                                           for m, v in col.items())))
    out = {"dims": dims, "monomials": mons, "support_confined": confined,
           "max_residual": worst, "total": sum(dims.values())}
    out["ok"] = (dims == DESCENT["ideal_grading"]
                 and mons == DESCENT["monomials_per_sector"]
                 and confined and worst == 0 and out["total"] == 55)
    return out


def gate_bundles():
    """Ten product-bundle rank tests against deg L = (0, 6, 5, 1).

    On a genus-1 base h^0(L) = deg L for positive degree and 1 for the
    trivial bundle, so each rank is min(h^0 of the product, #monomials).
    """
    E = adapted_qexpansions()
    deg = DESCENT["bundle_degrees"]
    rows, ok = [], True
    for a in range(4):
        for b in range(a, 4):
            sa, sb = SECTOR_ORDER[a], SECTOR_ORDER[b]
            ca, cb = _SECTOR_COORDS[sa], _SECTOR_COORDS[sb]
            pairs = (list(combinations_with_replacement(ca, 2)) if a == b
                     else [(i, j) for i in ca for j in cb])
            d = deg[sa] + deg[sb]
            h0 = d if d > 0 else 1
            pred = min(h0, len(pairs))
            got = _rank([_product(E, min(i, j), max(i, j)) for i, j in pairs])
            ok &= (got == pred)
            rows.append({"product": f"{sa} x {sb}", "monomials": len(pairs),
                         "deg": d, "predicted": pred, "computed": got})
    return {"tests": rows, "n": len(rows),
            "passed": sum(r["predicted"] == r["computed"] for r in rows),
            "ok": ok}


def _dim_intersect(cols, keep):
    """dim(I_2 ∩ Sym^2(keep)) for a set of adapted coordinates."""
    keep = set(keep)
    outside = [m for m in range(91)
               if MONOMIALS[m][0] not in keep or MONOMIALS[m][1] not in keep]
    proj = [[c.get(m, 0) for m in outside] for c in cols]
    return 55 - _rank([list(r) for r in zip(*proj)] if outside else [])


def gate_projection():
    """The nine-row projection table, including the three honest excesses."""
    cols = _i2_in_y()
    blk = {"f1": [0], "f3": list(range(1, 7)), "old": [7, 12],
           "f2": list(range(8, 12))}
    expected = {"f3": 9, "f2": 0, "old": 0, "f1+f3": 10, "f1+f2": 1,
                "f2+f3": 32, "newspace": 33, "f2+f3+old": 44, "all13": 55}
    subsets = {
        "f3": blk["f3"], "f2": blk["f2"], "old": blk["old"],
        "f1+f3": blk["f1"] + blk["f3"], "f1+f2": blk["f1"] + blk["f2"],
        "f2+f3": blk["f2"] + blk["f3"],
        "newspace": blk["f1"] + blk["f2"] + blk["f3"],
        "f2+f3+old": blk["f2"] + blk["f3"] + blk["old"],
        "all13": list(range(GENUS)),
    }
    got = {k: _dim_intersect(cols, v) for k, v in subsets.items()}
    return {"table": got, "expected": expected, "ok": got == expected}


def gate_descent():
    """Eta identity for the ghost lines, plus the ramification bookkeeping."""
    E = adapted_qexpansions()
    g = _eta_11a1()
    g13 = [g[n // 13] if n % 13 == 0 else 0 for n in range(NCOEF + 1)]
    e8 = [E[n][7] for n in range(NCOEF + 1)]
    e13 = [E[n][12] for n in range(NCOEF + 1)]
    plus = all(e8[n] == g[n] + 13 * g13[n] for n in range(NCOEF + 1))
    minus = all(e13[n] == g[n] - 13 * g13[n] for n in range(NCOEF + 1))

    d = DESCENT["branch_degrees"]
    deg = DESCENT["bundle_degrees"]
    bidouble = {
        "(+,-)": (d["W13"] + d["W143"]) == 2 * deg["(+,-)"],
        "(-,+)": (d["W11"] + d["W143"]) == 2 * deg["(-,+)"],
        "(-,-)": (d["W11"] + d["W13"]) == 2 * deg["(-,-)"],
    }
    rh = {k: 24 - 2 * (2 * DESCENT["quotient_genera"][k] - 2) ==
          DESCENT["fixed_points"][k] for k in ("W11", "W13", "W143")}
    return {
        "e8_is_g_plus_13g13": plus,
        "e13_is_g_minus_13g13": minus,
        "bidouble_relations": bidouble,
        "riemann_hurwitz": rh,
        "degrees_sum_to_12": sum(deg.values()) == 12,
        "fixed_points_sum_to_degK": sum(DESCENT["fixed_points"].values()) == 24,
        "ok": plus and minus and all(bidouble.values()) and all(rh.values())
        and sum(deg.values()) == 12,
    }


#: The wave's Route 2 ran 32 primes in [3, 149]. In-package the shipped
#: q-expansions stop at q^140, so a_p(143a1) is unavailable for p = 149 and
#: the gate reproduces 31 of the 32. The 32nd needed external curve data and
#: is recorded in studies/ci_2026aug/X0_143_CI_D_REPORT.md, not re-derived here.
ROUTE2_PRIMES_EXPECTED = 31


def gate_route2(pmax=NCOEF + 1):
    """Jac(X0(143)/W143) ~ 143a1 x 11a1, prime by prime.

    a_p(143a1) is read off the shipped e1 (= 72 f_143a1) and a_p(11a1) off
    the eta product, so the test needs no external curve data — at the cost
    of stopping where the shipped q-expansions do. See ROUTE2_PRIMES_EXPECTED.
    """
    E = adapted_qexpansions()
    if any(E[n][0] % 72 for n in range(NCOEF + 1)):
        raise AssertionError("e1 is not 72 * (integral newform)")
    a143 = [E[n][0] // 72 for n in range(NCOEF + 1)]
    a11 = _eta_11a1()

    def count(p):
        tot = 2                                   # the point at infinity
        for x in range(p):
            rhs = (x * x * x - x * x - x - 2) % p
            for y in range(p):
                if (y * y + y) % p != rhs:
                    continue
                if (x - 4) % p == 0 and (y + 7) % p == 0:
                    tot += 1 + _legendre(52, p)
                elif (x - 2) % p == 0 and (y + 1) % p == 0:
                    tot += 1 + _legendre(2, p)
                else:
                    l1 = (y - x + 3) % p
                    if l1 == 0:
                        tot += 1
                    else:
                        tot += 1 + _legendre(l1 * ((y + 3 * x - 5) % p), p)
        return tot

    rows, good = [], 0
    for p in range(3, pmax):
        if p in (11, 13) or any(p % q == 0 for q in range(2, p)):
            continue
        ap = p + 1 - count(p)
        hit = (ap == a143[p] + a11[p])
        good += hit
        rows.append({"p": p, "a_C2": ap, "a_143a1": a143[p],
                     "a_11a1": a11[p], "match": hit})
    return {"primes": len(rows), "matches": good, "rows": rows,
            "coefficient_depth": NCOEF,
            "ok": good == len(rows) and len(rows) >= ROUTE2_PRIMES_EXPECTED}


def gate_ci_a():
    """Q*: residual 0, confined to (+,+), leading coefficient factored."""
    E = adapted_qexpansions()
    q = ci_a_quadric()
    worst = 0
    prods = {m: _product(E, *MONOMIALS[m]) for m in q}
    for n in range(NCOEF + 1):
        worst = max(worst, abs(sum(v * prods[m][n] for m, v in q.items())))
    a = q[MONOMIALS.index((0, 0))]
    sectors = {monomial_sector(m) for m in q}
    return {"max_residual": worst, "a": a,
            "a_equals_minus_7sq_13_1957sq": a == -(7 ** 2) * 13 * 1957 ** 2,
            "sectors": sorted(sectors), "n_terms": len(q),
            "ok": worst == 0 and sectors == {"(+,+)"}
            and a == -(7 ** 2) * 13 * 1957 ** 2}


def verify(full=True):
    """Run the gates.  `full=False` skips the two slowest."""
    out = {
        "petri": gate_petri(),
        "sector_grading": gate_sector_grading(),
        "bundles": gate_bundles(),
        "projection": gate_projection(),
        "descent": gate_descent(),
        "ci_a": gate_ci_a(),
    }
    if full:
        out["generation"] = gate_generation()
        out["route2"] = gate_route2()
    out["ok"] = all(v["ok"] for v in out.values() if isinstance(v, dict))
    return out
