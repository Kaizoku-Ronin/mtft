"""mtft.lines.periods and mtft.lines.involutions: the 26 periods by Hecke block by two routes, the lattice of 143a1,
the Atkin-Lehner involutions and the mirror on the classes, and the lab's missed conjecture (MP12)."""
import math
import random

import pytest

from mtft.integral_lattice import smith_invariants
from mtft.lines import closed as C
from mtft.lines import darts as D
from mtft.lines import involutions as IV
from mtft.lines import periods as P
from mtft.lines import rademacher as RM
from mtft.lines import table as TB

BLOCKS = ("ell", "old", "q4", "q6")


def test_route_1_gives_the_lab_table_of_periods():
    L, lab = C.lines_143(), TB.lines()
    for l, r in zip(L, lab):
        assert {k: list(v) for k, v in P.coordinates(l.cls).items()} == r["blocks"]
        assert P.coordinates_of(l.gamma) == P.coordinates(l.cls)
    zero = {name: [n + 1 for n, l in enumerate(L) if not any(P.coordinates(l.cls)[name])] for name in BLOCKS}
    assert zero == {"ell": [9, 10, 11, 12, 22, 24, 25, 26, 27, 28], "old": [25, 26, 27, 28], "q4": [1, 2, 3, 4], "q6": []}
    v = P.thirty(C.MOTIONS["pell"])
    assert len(v["vector"]) == 30 and v["vector"][:4] == (-21, 0, -15, -1) and v["periods"]["ell"] == (0, 0)
    with pytest.raises(ValueError):
        P.coordinates((1, 2, 3))


def test_route_2_through_the_package_manin_symbols():
    """mtft.hecke's own model and blocks: no file.  What no choice of basis changes must agree with route 1."""
    L = C.lines_143()
    one = [P.coordinates(l.cls) for l in L]
    two = [P.components(l.gamma) for l in L]
    gen = P.generator_components()
    for name in BLOCKS:
        a, b = [c[name] for c in one], [c[name] for c in two]
        assert [not any(v) for v in a] == [not any(v) for v in b]                   # the same lines have no component
        assert P.rank(a) == P.rank(b) > 0
        assert P.relations(a) == P.relations(b)                                     # the same linear relations among the 28
        for l, comp in zip(L, b):                                                   # the Manin class is linear in cls
            assert comp == tuple(sum(l.cls[e] * gen[name][e][i] for e in range(29)) for i in range(len(comp)))
    assert {name: P.rank([c[name] for c in two]) for name in BLOCKS} == {"ell": 2, "old": 4, "q4": 5, "q6": 8}     # found, not predicted
    whole_1, whole_2 = [P.flat(c) for c in one], [P.flat(c) for c in two]
    assert P.rank(whole_1) == P.rank(whole_2) == 14 and P.relations(whole_1) == P.relations(whole_2)
    assert len(P.relations(whole_1)) == 14
    for name in C.MOTIONS:                                                          # the transfer, read in mtft.hecke's model alone
        fam = [P.flat(c) for l, c in zip(L, two) if l.motion == name]
        assert all(sum(v[i] for v in fam) == 0 for i in range(26))
    assert all(not any(two[n][name]) for n in (24, 25, 26, 27) for name in ("ell", "old"))      # the Pell lines
    with pytest.raises(ValueError):
        P.components(((2, 1), (1, 1)))


def test_the_two_bases_span_one_lattice():
    cb = P.compare_bases()
    for name, dim in zip(BLOCKS, (2, 4, 8, 12)):
        assert cb[name] == {"rank": dim, "consistent": True, "invertible": True, "smith": [1] * dim, "primitive": True}
    co = TB.coord()
    cusp = [[v for name in BLOCKS for v in co[name][e]] for e in range(29)]
    found = TB.record("found")
    assert math.prod(int(x) for x in smith_invariants(cusp)) == found["index_cuspidal"] == 1327104 == 2 ** 14 * 3 ** 4
    phi = TB.record("phi_on_generators")
    big = [list(phi[e]) + cusp[e] for e in range(29)]
    assert math.prod(int(x) for x in smith_invariants(big)) == found["index_product_lattice"] == 5573836800
    G = D.graph(143)                                                                # the cusp loops have no periods
    assert all(not any(P.flat(P.coordinates(G.face_cycle(f)))) for f in range(4))


def test_the_period_lattice_of_143a1_by_two_routes():
    lat, lab = P.lattice_143a1(), TB.record("periods")
    assert abs(lat["Omega1"] - complex(*lab["Omega1"])) < 1e-12 and abs(lat["Omega2"] - complex(*lab["Omega2"])) < 1e-12
    assert abs(lat["Omega1"] - complex(-0.984962, 2.015772)) < 1e-6 and abs(lat["Omega2"] + 1.969923) < 1e-6
    assert 0 < lat["residual"] < 1e-8 and lat["j_error"] < 1e-6                     # CERTIFIED(1e-8); the record: 6e-15
    assert abs(lat["j_curve"] + 141.01344809037116) < 1e-10 and abs(lat["Omega1"].real - lat["Omega2"].real / 2) < 1e-12
    assert P.small_elements() == [tuple(tuple(r) for r in g) for g in lab["small"]]
    assert [list(c) for c in lat["coordinates"]] == lab["ell_small"] and len({c for c in lat["coordinates"]}) > 2
    assert abs(P.period_143a1(C.MOTIONS["pell"])) < 1e-10                           # W_143 reverses the Pell unit
    assert abs(P.period_of(C.lines_143()[0].cls) + lat["Omega2"]) < 1e-12           # a cat line of length 7: (0, -1), real
    assert P.ap_143a1(2) == 0 and [P.ap_143a1(p) for p in (3, 5, 7)] == [-1, -1, -2] and P.an_143a1(200)[:6] == (1, 0, -1, -2, -1, 0)
    # second route: the package's stored periods (PARI/GP, mtft.periods) lie in the lattice and generate it
    s = P.stored_periods_in_lattice()
    assert s["count"] == 26 and s["worst"] < 1e-10 and s["index"] == 1 and len(set(s["integers"])) > 4


def test_atkin_lehner_involutions_and_the_mirror_on_the_classes():
    assert IV.group_law() == {"involutions": True, "commute": True, "product": True, "mirror_involution": True}
    assert IV.block_signs() == {"ell": (1, 1, 1), "old": (-1, 0, 0), "q4": (-1, 1, -1), "q6": (1, -1, -1)}
    assert {k: list(v) for k, v in IV.block_signs().items()} == TB.record("block_signs")
    G = D.graph(143)
    pell = C.MOTIONS["pell"]
    assert IV.al_conjugate(pell, 143) == D.inv(pell) and IV.apply(G.cls(pell), IV.al_matrix(143)) == tuple(-v for v in G.cls(pell))
    assert any(IV.al_matrix(11)[e] != IV.al_matrix(13)[e] for e in range(29))
    rng = random.Random(13)
    sample = [G.random_element(rng, rng.randint(1, 6)) for _ in range(40)] + [l.gamma for l in C.lines_143()]
    moved = {11: 0, 13: 0, 143: 0}
    for g in sample:                                                                # the matrices are the conjugations
        for Q in (11, 13, 143):
            h = IV.al_conjugate(g, Q)
            assert h[1][0] % 143 == 0 and D.trace(h) == D.trace(g) and G.cls(h) == IV.apply(G.cls(g), IV.al_matrix(Q))
            # the fibre number is not kept: it shifts by the winding of (eta(Q tau)/eta(tau))^24
            assert RM.Psi(h) - RM.Psi(g) == IV.psi_shift(g, Q) == 24 * RM.winding(g, {Q: 1, 1: -1})
            assert IV.psi_shift(h, Q) == -IV.psi_shift(g, Q)                        # and the shift is odd under the involution
            moved[Q] += RM.Psi(h) != RM.Psi(g)
    assert min(moved.values()) > 30                                                 # in most cases it does move
    assert RM.Psi(pell) + IV.psi_shift(pell, 143) == 21 == RM.Psi(D.inv(pell))
    assert IV.mirror_laws(sample) == {"psi_reversed": True, "windings_reversed": True, "linear": True, "ell_conjugated": True}
    L, lab = C.lines_143(), TB.lines()
    assert [IV.mirror_sign(l.cls) for l in L] == [r["mirror"] for r in lab]
    assert [IV.mirror_sign(l.cls) for l in L] == [1] * 4 + [-1] * 4 + [0] * 12 + [-1] * 8


def test_the_missed_conjecture_is_a_miss():
    """MP12, P13: 'every line with no 143a1 component is negated by one of the three involutions'.  False, and recorded."""
    cj = IV.conjecture_p13()
    assert cj["no_143a1_component"] == [9, 10, 11, 12, 22, 24, 25, 26, 27, 28]
    assert cj["negated_by_W143"] == cj["negated_by_some"] == [9, 10, 11, 12, 25, 26, 27, 28]
    assert cj["negated_by_none"] == [22, 24] and cj["holds"] is False
    assert cj["converse_holds"] is True                                             # the part that is a theorem: 143a1 has sign + under all three
    L, lab = C.lines_143(), TB.lines()
    assert [[IV.line_signs(l.cls)[Q] for Q in (11, 13, 143)] for l in L] == [r["al"] for r in lab]
    assert (L[21].motion, L[21].ell, L[23].ell) == ("s10", 30, 30)                  # the two zeros without a mechanism


def test_where_an_involution_sends_a_line():
    L = C.lines_143()
    count = 0
    for i, mo, ell, Q, hits in TB.record("al_images"):
        img = IV.al_image(L[i], Q)
        assert (L[i].motion, L[i].ell) == (mo, ell)
        if hits:
            assert sorted(img["hits"]) == sorted((s, j) for s, j in hits)
            count += len(hits)
    assert count >= 20
    first = IV.al_image(L[0], 11)                                                   # W_11 exchanges the cat lines 1 and 2
    assert (1, 1) in first["hits"] and first["psi"] == L[0].psi
    assert IV.al_image(L[24], 143)["hits"] == [(-1, 24)]                            # W_143 returns the Pell line of length 1, reversed
