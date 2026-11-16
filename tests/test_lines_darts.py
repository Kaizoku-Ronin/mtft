"""mtft.lines.darts and mtft.lines.rademacher: the graph on the darts of P^1(Z/143), classes, Dedekind sums, windings."""
import random
from fractions import Fraction as Fr

import numpy as np
import pytest

from mtft.lines import darts as D
from mtft.lines import rademacher as RM
from mtft.lines import table as TB


def _add(u, v):
    return tuple(a + b for a, b in zip(u, v))


def test_graph_of_level_143_and_its_generators():
    G = D.graph(143)
    assert G.summary() == {"N": 143, "darts": 168, "V": 56, "E": 84, "F": 4, "loops": 1, "faces": [1, 11, 13, 143], "b1": 29, "genus": 13}
    assert G.summary() == {**TB.record("graph"), "N": 143, "darts": 168}
    dd = G.darts
    assert dd.cycle_type(D.T) == [1, 11, 13, 143] and dd.cycle_type(D.S) == [2] * 84 and dd.cycle_type(D.R) == [3] * 56
    assert dd.cycle_type(((2, 1), (1, 1))) == [7] * 4 + [35] * 4                    # the cat map
    assert dd.points[dd.base] == (0, 1) and dd.elliptic() == (0, 0)
    assert G.generators == TB.generators()                                          # the lab's matrices, in its order
    for e, g in enumerate(G.generators):
        assert g[0][0] * g[1][1] - g[0][1] * g[1][0] == 1 and g[1][0] % 143 == 0
        assert G.cls(g) == tuple(int(j == e) for j in range(29)) and G.is_cycle(G.crossings(g))
    assert sorted(G.cusp_of_face(f) for f in range(4)) == [1, 11, 13, 143]
    assert {len(G.faces[f]): G.cusp_of_face(f) for f in range(4)} == {1: 143, 143: 1, 13: 11, 11: 13}      # width = N / denominator


def test_class_is_additive_and_conjugation_invariant():
    G, rng = D.graph(143), random.Random(1431)
    seen = set()
    for _ in range(200):
        g, h = G.random_element(rng, rng.randint(1, 6)), G.random_element(rng, rng.randint(1, 6))
        assert G.cls(D.mm(g, h)) == _add(G.cls(g), G.cls(h))
        assert G.cls(D.mm(D.mm(h, g), D.inv(h))) == G.cls(g)
        assert G.cls(D.inv(g)) == tuple(-v for v in G.cls(g))
        seen.add(G.cls(g))
    assert len(seen) > 150                                                          # the sample is not a few classes repeated
    c = tuple(rng.randint(-3, 3) for _ in range(29))
    assert G.cls(G.element(c)) == c
    with pytest.raises(ValueError):
        G.cls(D.S)                                                                  # S is not in the group: the path does not close


def test_cusp_loops_have_the_face_cycles_as_classes():
    G = D.graph(143)
    lab = TB.record("face_cycles")
    total = (0,) * 29
    for f in range(4):
        loop = G.cusp_loop(f)
        assert loop[1][0] % 143 == 0 and abs(D.trace(loop)) == 2                    # parabolic
        assert G.cls(loop) == G.face_cycle(f) == tuple(lab[str(len(G.faces[f]))])
        total = _add(total, G.face_cycle(f))
    assert total == (0,) * 29                                                       # the four cusp loops bound together


def test_words_and_the_character_to_z6():
    rng = random.Random(7)
    gens = (D.S, D.T, D.TI)
    for _ in range(200):
        g = D.I2
        for _ in range(rng.randint(1, 14)):
            g = D.mm(g, rng.choice(gens))
        back = D.unword(D.word(g))
        assert back in (g, D.neg(g))
        h = D.I2
        for _ in range(rng.randint(1, 9)):
            h = D.mm(h, rng.choice(gens))
        assert D.chi6(D.mm(g, h)) == (D.chi6(g) + D.chi6(h)) % 6                    # a homomorphism of PSL(2, Z)
    assert D.chi6(D.T) == 1 and D.chi6(D.S) == 3 and D.chi6(D.R) == 4
    assert D.mpow(D.R, 3) == D.I2 and D.mm(D.S, D.R) == D.T and D.mpow(D.T, -5) == ((1, -5), (0, 1))
    assert D.graph(143).chi_form() == TB.record("chi_generators")
    with pytest.raises(ValueError):
        D.as_matrix(((2, 0), (0, 2)))


def test_dedekind_sums_and_rademacher():
    for k in range(1, 40):
        for h in range(1, k + 1):
            if np.gcd(h, k) == 1:
                assert RM.dedekind(h, k) == RM.dedekind_by_definition(h, k)
    assert RM.dedekind(1, 3) == Fr(1, 18) and RM.dedekind(1, 2) == 0
    motions = {"cat": ((2, 1), (1, 1)), "s10": ((84, -47), (143, -80)), "pell": ((12, 1), (143, 12))}
    assert {k: RM.Psi(A) for k, A in motions.items()} == {"cat": 0, "s10": -1, "pell": -21} == TB.record("psi_motions")
    assert {k: list(RM.positive_word(A)) for k, A in motions.items()} == TB.record("letters")
    rng = random.Random(11)
    for A in motions.values():
        nu, nv = RM.positive_word(A)
        assert RM.Psi(A) == nu - nv and RM.Psi(D.inv(A)) == -RM.Psi(A)
        for _ in range(20):                                                         # a class function of SL(2, Z)
            g = D.I2
            for _ in range(rng.randint(1, 8)):
                g = D.mm(g, rng.choice((D.S, D.T, D.TI)))
            assert RM.Psi(D.mm(D.mm(g, A), D.inv(g))) == RM.Psi(A)
    assert RM.Phi(D.T) == 1 and RM.Phi(D.mpow(D.T, -7)) == -7 and RM.Psi(D.T) == 1
    with pytest.raises(ValueError):
        RM.dedekind(2, 4)


def test_windings_are_integer_linear_forms_of_the_class():
    G, rng = D.graph(143), random.Random(143)
    phi_gen = [RM.windings(g) for g in G.generators]
    assert [list(v) for v in phi_gen] == TB.record("phi_on_generators")
    nonzero = 0
    for _ in range(300):
        g = G.random_element(rng, rng.randint(1, 8))
        c = G.cls(g)
        w = RM.windings(g)
        assert w == tuple(sum(c[j] * phi_gen[j][i] for j in range(29)) for i in range(3))
        nonzero += any(w)
    assert nonzero > 250
    with pytest.raises(ValueError):
        RM.windings(((1, 0), (11, 1)))                                              # in the group of level 11 only


def test_windings_on_the_cusp_loops_are_the_orders_of_the_units():
    """Two routes: Dedekind sums along the loop round a cusp, and Ligozat's formula for the order at the cusp."""
    G = D.graph(143)
    want = {1: (5, -5, -1), 143: (5, 65, 11), 13: (-5, -65, 1), 11: (-5, 5, -11)}       # by the cusp's width
    assert {int(k): tuple(v) for k, v in TB.record("cusp_loops").items()} == want
    for f in range(4):
        width, den = len(G.faces[f]), G.cusp_of_face(f)
        w = RM.windings(G.cusp_loop(f))
        assert w == want[width]
        ligozat = tuple(RM.cusp_order(RM.UNITS[k], den, 143) for k in RM.UNIT_ORDER)
        assert ligozat == tuple(Fr(v) for v in w)
    for k in RM.UNIT_ORDER:                                                         # a unit: its orders add up to zero
        assert sum(RM.cusp_order(RM.UNITS[k], c, 143) for c in (1, 11, 13, 143)) == 0


def _winding_by_integral(g, r, nt=6000):
    """(1/24) times the integral of sum delta r_delta E2(delta z) dz from tau to g tau, by Gauss-Legendre."""
    sig1 = np.zeros(nt + 1)
    for dd in range(1, nt + 1):
        sig1[dd::dd] += dd

    def e2(z):
        q = np.exp(2j * np.pi * z)
        n = min(nt, int(40 / (2 * np.pi * z.imag.min())) + 50)
        qn = np.cumprod(np.repeat(q[:, None], n, axis=1), axis=1)
        return 1 - 24 * np.sum(sig1[1:n + 1][None, :] * qn, axis=1)

    (a, _), (c, d) = g
    z0, z1 = complex(-d, 1) / c, complex(a, 1) / c
    panels = max(8, int(abs(z1 - z0) * abs(c) * 2))
    xg, wg = np.polynomial.legendre.leggauss(24)
    tot = 0j
    for p in range(panels):
        lo, hi = z0 + (z1 - z0) * p / panels, z0 + (z1 - z0) * (p + 1) / panels
        z = lo + (hi - lo) * (xg + 1) / 2
        tot += np.sum(wg * sum(e * dl * e2(dl * z) for dl, e in r.items())) * (hi - lo) / 2
    return tot / 24


def test_windings_by_integrating_the_logarithmic_derivative():
    small = [((1, 0), (143, 1)), ((12, 1), (143, 12)), ((84, -47), (143, -80)), ((1, 0), (286, 1)), ((72, 1), (143, 2)), ((191, 2), (286, 3))]
    worst, values = 0.0, set()
    for g in small:
        for k in RM.UNIT_ORDER:
            exact = RM.winding(g, RM.UNITS[k])
            worst = max(worst, abs(_winding_by_integral(g, RM.UNITS[k]) - float(exact)))
            values.add(exact)
    assert worst < 1e-8 and len(values) > 5                                         # the record: 1e-13; and not all zero


def test_other_levels():
    g11 = D.graph(11)
    assert g11.summary() == {"N": 11, "darts": 12, "V": 4, "E": 6, "F": 2, "loops": 1, "faces": [1, 11], "b1": 3, "genus": 1}
    # one loop at every level: T fixes the dart (0 : 1), so that dart and its partner through S meet at one vertex
    assert all(len(D.graph(N).loops) == 1 for N in (11, 23, 47, 143))
    assert D.graph(23).genus == 2 and D.graph(23).b1 == 5 and D.graph(47).genus == 4
    for g in g11.generators:
        assert g[1][0] % 11 == 0 and sum(abs(v) for v in g11.cls(g)) == 1
    assert D.darts(13).size == 14 and D.darts(13).elliptic() == (2, 2) and D.darts(1).size == 1
    assert D.darts(13).cycle_type(D.T) == [1, 13]
    for N in (1, 13, 2, 7):
        with pytest.raises(ValueError):
            D.Graph(N)
    with pytest.raises(ValueError):
        D.Darts(0)
