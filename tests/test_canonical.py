"""Gates for the canonical-ideal arc (v0.17.0).

These are not stored assertions: every number is recomputed from the frozen
q-expansions in `mtft/canonical/_data/`.  Altering the data fails the suite.
"""

import hashlib

import pytest

from mtft import canonical
from mtft.canonical import gates


# ------------------------------------------------------------- provenance

def test_frozen_data_matches_provenance_manifest():
    """The shipped data must be byte-identical to the audited wave files."""
    manifest = canonical.data_path("PROVENANCE.txt").read_text().splitlines()
    entries = [ln.split() for ln in manifest
               if ln.strip() and not ln.startswith("#")]
    assert len(entries) == 6
    for sha, name in entries:
        blob = canonical.data_path(name).read_bytes()
        assert hashlib.sha1(blob).hexdigest() == sha, name


def test_data_shapes():
    assert len(canonical.s2_qexpansions()) == 141
    assert len(canonical.adapted_qexpansions()[0]) == canonical.GENUS
    assert len(canonical.adapted_basis()) == canonical.GENUS
    assert len(canonical.ideal_basis()) == 91
    assert len(canonical.ideal_basis()[0]) == 55
    assert len(canonical.MONOMIALS) == 91


# ------------------------------------------------- sector-ordering convention

def test_sector_ordering_is_pinned_and_convertible():
    """Two orderings circulate in the corpus; the map between them is exact."""
    assert canonical.SECTOR_ORDER == ("(+,+)", "(+,-)", "(-,+)", "(-,-)")
    dims = tuple(canonical.DESCENT["sector_dims_S2"][s]
                 for s in canonical.SECTOR_ORDER)
    assert dims == (1, 6, 5, 1)
    # the v0.16.0 CHANGELOG quotes the same data as (1, 5, 6, 1)
    assert canonical.reorder_sectors(dims) == (1, 5, 6, 1)
    assert canonical.reorder_sectors(
        canonical.reorder_sectors(dims),
        frm=canonical.SECTOR_ORDER_CHANGELOG,
        to=canonical.SECTOR_ORDER) == dims


def test_monomial_sector_partition():
    counts = {s: 0 for s in canonical.SECTOR_ORDER}
    for m in range(91):
        counts[canonical.monomial_sector(m)] += 1
    assert counts == canonical.DESCENT["monomials_per_sector"]
    assert sum(counts.values()) == 91


# --------------------------------------------------------------- P1 - P4

def test_gate_petri():
    """P1-P4 and exact residual 0 for all 55 quadrics."""
    r = gates.gate_petri()
    assert r["P1_dim_S2"] == 13
    assert r["P2_dim_Sym2"] == 91
    assert r["P3_h0_2K"] == 36           # 3g - 3
    assert r["P4_dim_I2"] == 55          # (g-2)(g-3)/2, Petri
    assert r["max_residual"] == 0        # exact integers, not a tolerance
    assert r["ok"]


def test_not_hyperelliptic():
    """Hyperelliptic would force dim I_2 = C(12,2) = 66, not 55."""
    assert gates.gate_petri()["P4_dim_I2"] != 66


# ------------------------------------------------------------- the grading

def test_gate_sector_grading():
    r = gates.gate_sector_grading()
    assert r["dims"] == {"(+,+)": 26, "(+,-)": 5, "(-,+)": 4, "(-,-)": 20}
    assert r["total"] == 55
    assert r["support_confined"]
    assert r["max_residual"] == 0
    assert r["ok"]


def test_gate_bundles():
    """deg L = (0, 6, 5, 1) on E = X0(143)*, ten product-rank tests."""
    r = gates.gate_bundles()
    assert r["n"] == 10
    assert r["passed"] == 10
    assert r["ok"]


def test_gate_projection():
    """Nine projections, six matching Riemann-Roch on the quotients."""
    r = gates.gate_projection()
    assert r["table"] == r["expected"]
    # the three honest excesses are part of the record, not bugs
    assert r["table"]["newspace"] == 33      # excess 3: the ghost-only piece
    assert r["table"]["f2+f3+old"] == 44     # excess 2
    assert r["table"]["f1+f2"] == 1          # excess 1
    assert r["ok"]


def test_newspace_deficiency_is_three():
    """A 3-dim piece of H^0(2K) is reachable only through the level-11 ghost."""
    r = gates.gate_projection()
    assert 66 - r["table"]["newspace"] == 33     # rank, not 36
    assert 36 - 33 == 3


# ---------------------------------------------------------------- descent

def test_gate_descent():
    """Ghost sections as eta products, plus the ramification bookkeeping."""
    r = gates.gate_descent()
    assert r["e8_is_g_plus_13g13"]
    assert r["e13_is_g_minus_13g13"]
    assert all(r["bidouble_relations"].values())
    assert all(r["riemann_hurwitz"].values())
    assert r["degrees_sum_to_12"]
    assert r["ok"]


def test_W11_acts_freely():
    d = canonical.DESCENT
    assert d["fixed_points"]["W11"] == 0
    assert d["branch_degrees"]["W11"] == 0
    assert sum(d["fixed_points"].values()) == 24      # = deg K_X


def test_L_minusminus_is_a_generator_point():
    d = canonical.DESCENT["L_minusminus"]
    assert d["point"] == (4, -7)
    a = canonical.DESCENT["weierstrass"]
    x, y = d["point"]
    lhs = y * y + a["a1"] * x * y + a["a3"] * y
    rhs = x ** 3 + a["a2"] * x * x + a["a4"] * x + a["a6"]
    assert lhs == rhs, "P must lie on 143a1"


def test_gate_ci_a():
    """Q* is confined to (+,+) and vanishes exactly."""
    r = gates.gate_ci_a()
    assert r["max_residual"] == 0
    assert r["sectors"] == ["(+,+)"]
    assert r["a"] == -(7 ** 2) * 13 * 1957 ** 2
    assert r["ok"]


# ------------------------------------------------------------------- slow

@pytest.mark.slow
def test_gate_generation():
    """P5-P8: quadrics generate I_3; excludes trigonal and plane quintic."""
    r = gates.gate_generation()
    assert r["P5_dim_Sym3"] == 455
    assert r["P6_h0_3K"] == 60           # 5g - 5
    assert r["P8_rank_V_I2"] == 395      # Enriques-Babbage
    assert r["ok"]


@pytest.mark.slow
def test_gate_route2():
    """Jac(X0(143)/W143) ~ 143a1 x 11a1, prime by prime."""
    r = gates.gate_route2()
    assert r["primes"] >= gates.ROUTE2_PRIMES_EXPECTED
    assert r["matches"] == r["primes"]
    assert r["ok"]
