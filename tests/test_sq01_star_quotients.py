"""SQ-01: star quotients X0(pq)* for supersingular p < q — two routes, the census, what singles out 143, and the record."""
import contextlib
import io
import json

import pytest

from mtft.matter import ledger as L
from mtft.matter import odd_modes as OM
from mtft.matter import puzzle as Z
from mtft.matter import star_quotients as S
from mtft.matter.__main__ import main as matter_main

TORI = [57, 58, 65, 77, 82, 91, 118, 123, 141, 142, 143, 145, 155]
SPHERES = [6, 10, 14, 15, 21, 22, 26, 33, 34, 35, 38, 39, 46, 51, 55, 62, 69, 87, 94, 95, 119]


def test_class_numbers_symbols_and_fixed_points():
    assert [S.class_number(D) for D in (-3, -4, -7, -8, -12, -23, -56, -84, -143, -572)] == [1, 1, 1, 1, 1, 3, 4, 4, 10, 10]
    assert [S.conductor(D) for D in (-3, -12, -27, -28, -52, -572)] == [1, 2, 3, 2, 1, 2]
    assert S.kronecker(-8, 5) == -1 and S.kronecker(-4, 5) == 1 and S.kronecker(-7, 2) == 1 and S.kronecker(-3, 2) == -1
    assert S.eichler(-28, 2) == 1 and S.kronecker(-28, 2) == 0                    # l divides the conductor: Eichler symbol 1
    fp = {Q: S.fixed_points(Q, [l for l in (11, 13) if l not in Q]) for Q in ((11,), (13,), (11, 13))}
    assert fp == {(11,): 0, (13,): 4, (11, 13): 20}
    # X0(10) has genus 0, so tr W_2 = 0 and W_2 has two fixed points; both come from the Z[i] term (the sqrt(-2) term vanishes at 5)
    assert S.genus((2, 5)) == 0 and S.fixed_points((2,), (5,)) == 2 and S.class_number(-8) * (1 + S.kronecker(-8, 5)) == 0


def test_two_routes_agree_wherever_both_run():
    rows = S.census()
    assert len(rows) == 105 and all(isinstance(r["gstar"], int) for r in rows)        # route 2: integral at every level
    both = [r for r in rows if r["routes"] == 2]
    assert len(both) == 99 and all(r["gstar"] == r["gstar_route1"] for r in both)
    assert sorted(r["N"] for r in rows if r["routes"] == 1) == [2201, 2419, 2773, 2911, 3337, 4189]
    for r in rows:                                                                     # route 1's orbits fill the new space the genus formula predicts
        if r["orbits"] is not None:
            assert sum(r["orbits"]) == r["new_dim"] == S.genus((r["p"], r["q"])) - 2 * S.genus((r["p"],)) - 2 * S.genus((r["q"],))
    d = S.dense_census()
    assert len(d) == 75 and all(r["gstar"] == r["gstar_route1"] for r in d)


def test_census():
    assert S.spheres() == SPHERES and S.tori() == TORI
    assert sum(1 for r in S.census() if r["gstar"] >= 2) == 71
    assert S.ranks() == {N: 1 for N in TORI}
    assert {N: S.root_number(N) for N in TORI} == {N: -1 for N in TORI}              # odd rank forced by the signs
    assert S.star_genus((11, 13)) == 1 and S.genus((11, 13)) == 13 and S.index((11, 13)) == 168


def test_signs_are_constant_on_orbits_and_the_sectors_of_143():
    assert S.ROUTE1[143] == ((1, 1, 1), (4, -1, 1), (6, 1, -1))                       # one sign pair per Galois orbit
    sectors = {(a, b): 0 for a in (1, -1) for b in (1, -1)}
    for d, a, b in S.ROUTE1[143]:
        sectors[(a, b)] += d
    sectors[(-1, 1)] += 1; sectors[(-1, -1)] += 1                                      # old forms from 11a: w_11 = -1 (Ogg), W_13 = +1 and -1
    assert sectors == {(1, 1): 1, (1, -1): 6, (-1, 1): 5, (-1, -1): 1} and sum(sectors.values()) == S.genus((11, 13))


def test_what_singles_out_143():
    assert S.twin_ladder() == [(15, 0), (35, 0), (143, 1), (323, 5), (899, 12)]
    assert [N for N, g in S.consecutive_ladder() if g == 1] == [77, 143]
    assert S.new_dim_levels("smaller") == [26, 38, 39, 51, 57, 58, 65, 69, 85, 91, 143]
    larger = S.new_dim_levels("larger")
    assert larger == [221, 247, 299, 377, 403, 533, 611, 767, 923]
    consecutive = {r["N"] for r in S.census() if r["consecutive"]}
    assert [N for N in larger if N in consecutive] == [221]                             # 'the larger prime', read literally, picks 221
    assert S.distinctions(143) == {"only_twin_torus": True, "only_new_dim_smaller_and_consecutive": True, "only_orbit_pattern": True,
                                   "only_consecutive_torus": False, "only_index": False, "only_torus": False}
    assert S.siblings(143) == [123] and S.siblings(123) == [143]
    assert not any(S.distinctions(123).values())
    r123 = next(r for r in S.census() if r["N"] == 123)
    assert (r123["index"], r123["genus"], r123["torsion_free"], r123["new_dim"], r123["orbits"]) == (168, 13, True, 7, (1, 1, 2, 3))


def test_moonshine_cross_check():
    m = S.monster_check()
    assert len(S.MONSTER_ELEMENT_ORDERS) == m["n_orders"] == 73 and max(S.MONSTER_ELEMENT_ORDERS) == 119
    assert m["spheres_are_orders"] and len(m["spheres"]) == 28 and m["orders_not_spheres"] == [57, 93] and not m["143_is_an_order"]


def test_denser_factorizations():
    d = S.dense_census()
    assert [r["N"] for r in d if r["gstar"] == 1] == [102, 114, 130, 138, 174, 182, 190, 195, 210, 231, 238]
    assert [r["N"] for r in d if r["gstar"] == 0] == [30, 42, 66, 70, 78, 105, 110]
    assert {r["cusps"] for r in d} == {8, 16}


def test_flavour_obstruction_does_not_use_the_level():
    cost = S.flavour_transfer()
    assert sorted(cost) == TORI and min(cost.values()) == 24 and max(cost.values()) == 96
    assert cost[143] == OM.family_cost(3) == 72                                         # OM-3's 72 at genus 13


def test_raise_paths():
    for bad in (0, 5, -5):
        with pytest.raises(ValueError):
            S.class_number(bad)
    for bad in ((2, 2), (4,), ()):
        with pytest.raises(ValueError):
            S.genus(bad)
    with pytest.raises(ValueError):
        S.fixed_points((3,), (3,))
    with pytest.raises(ValueError):
        S.star_genus_route1(4189)
    with pytest.raises(ValueError):
        S.root_number(15)
    with pytest.raises(ValueError):
        S.new_dim_levels("middle")
    with pytest.raises(ValueError):
        S.distinctions(144)
    with pytest.raises(ValueError):
        S.text(144)


def test_command_line():
    def run(*argv):
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            code = matter_main(["levels", *argv])
        return code, buf.getvalue()
    code, out = run(); assert code == 0 and "143 -> 1" in out and "siblings of 143: 123" in out
    code, out = run("--json"); d = json.loads(out); assert code == 0 and d["tori"] == TORI and d["routes_agree"] is True
    code, out = run("--level", "143"); assert code == 0 and "only_twin_torus" in out and "(1, 4, 6)" in out
    code, out = run("--level", "144"); assert code == 2 and "not a product" in out


def test_the_record():
    E = {e["id"]: e for e in L.E}
    sq = [e for e in L.E if e["wave"] == "SQ-01"]
    assert [e["id"] for e in sq] == [f"SQ-{i}" for i in range(1, 9)] and all(e["module"] == "mtft.matter.star_quotients" for e in sq)
    assert all(str(N) in E["SQ-1"]["statement"] for N in TORI) and "221" not in E["SQ-1"]["statement"]
    a = Z.assignment()
    assert [a[f"SQ-{i}"] for i in range(1, 9)] == ["II"] * 6 + ["IX", "IV"]
    subs = Z.subs()
    assert all(subs[s][1].status == "answered" for s in ("II.7", "II.8", "II.9")) and subs["IV.6"][1].status == "open"
    assert "SQ-8" in subs["IV.6"][1].claims and "123" in subs["IV.6"][1].settle
    assert any(n["id"] == "N20" and n["piece"] == "II" for n in Z.NEGATIVES)


def test_pari_route1_reproduces_the_frozen_table():
    cypari2 = pytest.importorskip("cypari2")
    pari = cypari2.Pari()
    for N in TORI + [15, 221]:
        p, q = [x for x in S.SUPERSINGULAR if N % x == 0]
        pari(f"mf = mfinit([{N},2],0); S = mfsplit(mf)")
        degs = [int(x) for x in pari("apply(poldegree, S[2])")]
        wp = [int(v[0]) for v in pari(f"mfatkineigenvalues(mf,{p})")]
        wq = [int(v[0]) for v in pari(f"mfatkineigenvalues(mf,{q})")]
        assert tuple(sorted(zip(degs, wp, wq), key=lambda o: (o[0], -o[1], -o[2]))) == S.ROUTE1[N], N
