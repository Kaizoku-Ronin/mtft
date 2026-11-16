"""mtft.lines.cutting, crossings and lorenz: every crossing of the closed geodesics, exactly, by cutting along the
fundamental domain (MP13, MP14) and, for the totals, by Ghys's Lorenz template, which shares no step with it."""
import math
from fractions import Fraction as Fr
from itertools import combinations

import pytest

from mtft.lines import closed as C
from mtft.lines import crossings as X
from mtft.lines import cutting as K
from mtft.lines import darts as D
from mtft.lines import involutions as IV
from mtft.lines import lorenz as LZ
from mtft.lines import periods as P
from mtft.lines import rademacher as RM
from mtft.lines import table as TB

SIX = {**C.MOTIONS, **C.MORE_MOTIONS}


def _point(m, p):
    """The point (x, y^2) of the fundamental domain, carried into the plane by the frame m."""
    z = complex(float(p[0]), math.sqrt(float(p[1])))
    return (m[0][0] * z + m[0][1]) / (m[1][0] * z + m[1][1])


def test_exact_cutting_sequences():
    lab = TB.record("crossings_24")["walks"]
    for name, A in SIX.items():
        w = K.walk(A)
        assert abs(w.length - 2 * math.acosh(abs(D.trace(A)) / 2)) < 1e-12          # the arcs add up to the geodesic's length
        for k, ar in enumerate(w.arcs):                                             # each arc is the axis seen from its frame
            assert ar.N == K.conj(ar.m, K.pos(A)) and K.form(ar.N) == ar.form
            assert all(isinstance(v, Fr) for p in (ar.entry, ar.exit) for v in p)
            for x, y2 in (ar.entry, ar.exit):                                       # the ends lie on the arc and on the boundary of F
                a, b, c = ar.form
                assert a * (x * x + y2) + b * x + c == 0 and K.kind((x, y2)) != "int"
            nxt = w.arcs[(k + 1) % len(w.arcs)]                                     # the arcs join up in the plane, the last to the first
            here = _point(ar.m, ar.exit)                                            # through one period of A
            there = _point(nxt.m if k + 1 < len(w.arcs) else D.mm(K.pos(A), nxt.m), nxt.entry)
            assert abs(here - there) < 1e-9 * max(1.0, abs(here))
        if name in lab:
            assert (len(w.arcs), w.exits, w.through_i, w.through_rho) == (lab[name]["arcs"], lab[name]["exits"], lab[name]["through_i"], lab[name]["through_rho"])
    assert (K.walk(SIX["cat"]).exits, K.walk(SIX["s10"]).exits) == ("LiRi", "pLL")
    assert K.walk(SIX["pell"]).through_i == 0 and len(K.walk(SIX["pell"]).arcs) == 25
    assert [len(K.walk(SIX[k]).arcs) for k in ("uv3", "u2v2", "uv5")] == [6, 6, 8]
    assert all(K.mirror_relation(A)[0] == -1 for A in SIX.values())                 # each orbit is its own mirror image, reversed
    assert K.walk(D.neg(SIX["cat"])) is K.walk(SIX["cat"])
    with pytest.raises(ValueError):
        K.walk(D.mpow(SIX["cat"], 2))                                               # a proper power
    with pytest.raises(ValueError):
        K.walk(D.T)
    p = K.cross((-1, -1, 1), (-1, 1, 1))                                            # the cat's two arcs through i
    assert p == (Fr(0), Fr(1)) and K.kind(p) == "i"
    assert K.cross((1, 0, -2), (1, -1, -1)) is None                                 # they meet at x = 1, outside the domain
    assert K.cross((1, 0, -2), (2, 0, -9)) is None                                  # concentric: no crossing at all


def test_the_24_geodesics_cross_as_the_lab_found():
    Z, lab = X.crossings_143(), TB.record("crossings_24")
    assert len(Z.lines) == 28 and len(Z.geodesics) == 24
    assert [(G["motion"], G["ell"], [n + 1 for n in G["lines"]]) for G in Z.geodesics] == [(G["motion"], G["ell"], G["lines"]) for G in lab["geodesics"]]
    assert Z.I == lab["I"] and Z.iota == lab["iota"] and Z.self_crossings() == lab["self"]
    assert Z.check_conventions()                                                    # the mirrored half-open convention: the same
    n = 24
    assert all(Z.I[g][h] == Z.I[h][g] and Z.iota[g][h] == -Z.iota[h][g] for g in range(n) for h in range(n))
    assert all(abs(Z.iota[g][h]) <= Z.I[g][h] and (Z.I[g][h] - Z.iota[g][h]) % 2 == 0 for g in range(n) for h in range(n) if g != h)
    s = Z.summary()
    assert (s["crossings_between"], s["self_crossings"], s["unlinked_pairs"], s["pairs"], s["signed_total"], s["simple"]) == (3214, 2918, 103, 276, 202, 15)
    assert s["places"] == {"i": 84, "int": 3696, "left": 2184, "rho": 168}
    assert Z.count("mirrored")["places"] == {"i": 84, "int": 3696, "right": 2184, "rho1": 168}
    most = max(Z.I[g][h] for g in range(n) for h in range(n) if g != h)
    assert Z.self_crossings()[20:] == [0, 10, 24, 2794] and most == Z.I[22][23] == 506      # the Pell lines of lengths 13 and 143


def test_linking_numbers_of_the_cat_unions():
    Z = X.crossings_143()
    cat = Z.family("cat")
    Lk = Z.linking()
    assert [Z.label(g) for g in cat] == ["cat7(1-4)", "cat7(2-3)", "cat35(5-6)", "cat35(7-8)"]
    assert [[Lk[a][b] for b in cat] for a in cat] == [[None, 0, -6, -8], [0, None, -8, -6], [-6, -8, None, -24], [-8, -6, -24, None]]
    assert all(v is None or v <= 0 for row in Lk for v in row)                      # never positive
    # the cat geodesics cross only over the order-2 point: 84 crossings, and for each P: sum I(P, Q) + 2 self(P) = 2 ell
    assert all(sum(Z.I[p][q] for q in cat if q != p) + 2 * Z.I[p][p] == 2 * Z.geodesics[p]["ell"] for p in cat)
    assert sum(Z.I[p][q] for p, q in combinations(cat, 2)) + sum(Z.I[p][p] for p in cat) == 84
    sys4 = Z.family("s10")[:4]                                                      # the four systoles: simple and pairwise disjoint
    assert all(Z.geodesics[g]["ell"] == 1 for g in sys4) and all(Z.I[a][b] == 0 for a in sys4 for b in sys4)
    assert (Z.total("cat", "s10"), Z.total("cat", "pell"), Z.total("s10", "pell")) == (336, 336, 1008)
    assert {k: -sum(Z.I[p][d] for p in cat for d in Z.family(k)) // 2 for k in ("s10", "pell")} == TB.record("crossings_24")["totals"]["Lk_cat_B"] == {"s10": -168, "pell": -168}


def test_the_cat_block_by_a_formula_on_the_darts():
    """The lab's second route for the cat geodesics (MP13, P3), with no cutting: they cross only over the order-2 point, once
    at each of its 84 lifts, and at the lift {y, yS} the two branches are the cat geodesics through the darts y and y T^-1."""
    Z = X.crossings_143()
    dd, cat = Z.darts, Z.family("cat")
    for which in (0, 1):                                                            # either dart of the edge gives the same pair
        M = {(a, b): 0 for a in cat for b in cat}
        for o in D.graph(143).edges:
            y = sorted(o)[which]
            g = Z.geodesic_at("cat", dd.points[y])[0]
            h = Z.geodesic_at("cat", dd.points[dd.act(y, D.TI)])[0]
            M[(g, h)] += 1
            if g != h:
                M[(h, g)] += 1
        assert [[M[(a, b)] for b in cat] for a in cat] == [[Z.I[a][b] for b in cat] for a in cat] == [[0, 0, 6, 8], [0, 0, 8, 6], [6, 8, 16, 24], [8, 6, 24, 16]]


def test_the_mirror_and_the_pell_zeros():
    Z, lab = X.crossings_143(), TB.record("crossings_24")
    mir = Z.mirror()
    assert [m[0] for m in mir] == lab["mirror"]["image"] and [m[1] for m in mir] == lab["mirror"]["eps"]
    assert Z.mirror_check()                                                         # counts kept, signs reversed
    rev = Z.reversed_by_mirror()
    assert len(rev) == 10 and set(Z.family("pell")) <= set(rev)
    assert all(Z.iota[g][h] == 0 for g in rev for h in rev)                         # two reversed geodesics: signed count zero
    pell = Z.family("pell")
    assert all(v == 0 for g in pell for v in Z.iota[g])                             # the Pell rows vanish on these 24 ...
    assert all(Z.I[pell[3]][h] >= 6 for h in range(24) if h != pell[3])             # ... though the long one crosses every other
    other = [h for h in range(24) if h not in rev]
    assert len(other) == 14 and all(Z.I[pell[0]][h] == 0 for h in other)            # the Pell line of length 1 never meets the 14
    # classes: the mirror's sign on the class agrees with the geometry
    for g in range(24):
        sign = IV.mirror_sign(Z.lines[Z.geodesics[g]["rep"]].cls)
        assert sign == (mir[g][1] if mir[g][0] == g else 0)


def test_signed_counts_are_a_form_on_the_periods():
    Z = X.crossings_143()
    vec = [P.flat(P.coordinates(Z.lines[G["rep"]].cls)) for G in Z.geodesics]
    rel = P.relations(vec)
    assert P.rank(vec) == 14 and len(rel) == 10
    for v in rel:                                                                   # every relation among the periods kills the signed rows
        assert all(sum(v[g] * Z.iota[g][h] for g in range(24)) == 0 for h in range(24))
    assert P.rank(Z.iota) == 12


def test_the_38_geodesics_of_the_third_record():
    E, lab = X.crossings_143_extended(), TB.record("crossings_38")
    assert len(E.lines) == 44 and len(E.geodesics) == 38
    mp, eps = [], []
    for G in lab["geodesics"]:                                                      # matched by motion and dart
        g, d = E.geodesic_at(G["motion"], tuple(G["dart"]))
        mp.append(g)
        eps.append(d)
        assert E.geodesics[g]["ell"] == G["ell"] and len(E.geodesics[g]["lines"]) == G["lifts"]
    assert sorted(mp) == list(range(38))
    assert all(E.I[mp[a]][mp[b]] == lab["I"][a][b] for a in range(38) for b in range(38))
    assert all(eps[a] * eps[b] * E.iota[mp[a]][mp[b]] == lab["iota"][a][b] for a in range(38) for b in range(38))
    assert E.check_conventions() and E.mirror_check()
    mir = E.mirror()
    kinds = ["reversed" if mir[g] == (g, -1) else ("kept" if mir[g] == (g, 1) else "swapped") for g in range(38)]
    assert all(kinds[mp[a]] == G["mirror"] for a, G in enumerate(lab["geodesics"]))
    rev = E.reversed_by_mirror()
    assert len(rev) == 20 and all(E.iota[g][h] == 0 for g in rev for h in rev)
    # the Pell class is not isotropic: UV^5's short lines cross the Pell line of length 1 once each, with signs
    p1 = E.geodesic_at("pell", (0, 1))[0]
    short = [g for g in E.family("uv5") if E.geodesics[g]["ell"] == 7]
    assert sorted(E.iota[p1][g] for g in short) == [-1, -1, 1, 1] and all(E.I[p1][g] == 1 for g in short)
    # the old 24 sit inside the 38 unchanged
    Z = X.crossings_143()
    old = [E.geodesic_at(Z.geodesics[g]["motion"], Z.lines[Z.geodesics[g]["rep"]].dart) for g in range(24)]
    assert all(E.I[old[a][0]][old[b][0]] == Z.I[a][b] for a in range(24) for b in range(24))
    assert all(old[a][1] * old[b][1] * E.iota[old[a][0]][old[b][0]] == Z.iota[a][b] for a in range(24) for b in range(24))


def test_lorenz_words_and_the_template():
    want = {"cat": ("xy",), "s10": ("xxy", "xyy"), "pell": ("x" * 22 + "y", "x" + "y" * 22),
            "uv3": ("xxxy", "xyyy"), "u2v2": ("xxyy",), "uv5": ("xxxxxy", "xyyyyy")}
    for name, A in SIX.items():
        assert LZ.words(A) == want[name]
        w = LZ.lorenz_word(A)
        assert w.count("x") - w.count("y") == RM.Psi(A)                             # Ghys: #x - #y is Rademacher's number
        assert abs(D.trace(LZ.word_matrix(w))) == abs(D.trace(A))
    assert LZ.is_reciprocal(SIX["cat"]) and LZ.is_reciprocal(SIX["u2v2"]) and not LZ.is_reciprocal(SIX["s10"])
    y22, x22 = "x" + "y" * 22, "x" * 22 + "y"
    vals = [LZ.lorenz_linking(a, b) for a, b in (("xy", "xyy"), ("xy", "xxy"), ("xy", y22), ("xy", x22), ("xyy", y22), ("xxy", x22), ("xyy", x22), ("xxy", y22))]
    assert vals == [1, 1, 1, 1, 2, 2, 1, 1]                                         # the record's pre-registered values
    assert LZ.lorenz_linking("x", "y") == 0 and LZ.lorenz_linking("xyy", "xy") == LZ.lorenz_linking("xy", "xyy")
    assert (LZ.transfer_sum(SIX["cat"], SIX["s10"]), LZ.transfer_sum(SIX["cat"], SIX["pell"]), LZ.transfer_sum(SIX["s10"], SIX["pell"])) == (2, 2, 6)
    with pytest.raises(ValueError):
        LZ.lorenz_linking("xy", "yx")
    with pytest.raises(ValueError):
        LZ.transfer_sum(SIX["cat"], D.inv(SIX["cat"]))


def test_the_transfer_by_two_routes_at_three_levels():
    """Cutting geodesics at level N against psi(N) times the template's sums at level 1: fifteen pairs, three levels."""
    sums = {(a, b): LZ.transfer_sum(SIX[a], SIX[b]) for a, b in combinations(SIX, 2)}
    assert len(sums) == 15 and sorted(set(sums.values())) == [2, 4, 6, 8, 12]
    for N, index in ((143, 168), (11, 12), (23, 24)):
        Z = X.crossings(N, SIX)
        assert Z.darts.size == index and Z.check_conventions()
        for (a, b), s in sums.items():
            assert Z.total(a, b) == index * s, (N, a, b)
    assert X.crossings(143, SIX) is X.crossings_143_extended()


def test_elliptic_levels_are_refused():
    for N in (13, 1, 7):
        with pytest.raises(ValueError):
            X.Crossings(N)
    with pytest.raises(ValueError):
        X.Crossings(143, {"square": D.mpow(SIX["cat"], 2)})
