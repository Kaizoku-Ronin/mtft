"""mtft.lines.closed: the 28 closed lines of Gamma_0(143)\\PSL(2, R) over three motions, against the lab's table (MP12)."""
import math
import random

import pytest

from mtft.integral_lattice import smith_invariants
from mtft.lines import closed as C
from mtft.lines import darts as D
from mtft.lines import periods as P
from mtft.lines import rademacher as RM
from mtft.lines import table as TB
from mtft.matter import one_loop_adjoints as OL


def test_the_28_lines_equal_the_lab_table():
    assert TB.verify_provenance() == {"lines_143.json": True}
    L, lab = C.lines_143(), TB.lines()
    assert len(L) == len(lab) == 28
    for l, r in zip(L, lab):
        assert (l.motion, l.ell, list(l.dart), l.psi, list(l.windings), list(l.cls)) == (r["motion"], r["ell"], r["dart"], r["psi"], r["phi"], r["cls"])
        assert l.gamma[1][0] % 143 == 0 and l.gamma[0][0] * l.gamma[1][1] - l.gamma[0][1] * l.gamma[1][0] == 1
        assert D.graph(143).cls(l.gamma) == l.cls                              # tracing the cycle = decomposing the large matrix
        assert l.cycle[0] == min(l.cycle) and len(set(l.cycle)) == l.ell
    types = {k: [l.ell for l in L if l.motion == k] for k in C.MOTIONS}
    assert types == {"cat": [7] * 4 + [35] * 4, "s10": [1] * 4 + [5] * 4 + [6] * 4 + [30] * 4, "pell": [1, 11, 13, 143]}
    assert max(len(str(abs(v))) for l in L for row in l.gamma for v in row) > 150   # the Pell line of length 143: about 199 digits
    assert not any(all(v == 0 for v in l.cls) for l in L)                           # no line is null-homologous


def test_transfer_congruence_and_reciprocity():
    G = D.graph(143)
    chi = G.chi_form()
    L = C.lines_143()
    for l in L:
        assert l.psi == l.ell * RM.Psi(l.A)
        assert l.psi % 6 == l.chi == sum(c * x for c, x in zip(l.cls, chi)) % 6     # Psi = chi(cls) modulo 6
        assert abs(l.length - l.ell * 2 * math.acosh(abs(D.trace(l.A)) / 2)) < 1e-12
    for name, A in C.MOTIONS.items():
        t = C.transfer(A)
        assert t["sum_cls"] == (0,) * 29 and t["sum_windings"] == (0, 0, 0)         # the lifts of a motion add up to zero ...
        assert t["sum_psi"] == 168 * RM.Psi(A) and t["darts"] == 168                # ... but for the fibre: 168 Psi(A)
        co = [P.flat(P.coordinates(l.cls)) for l in L if l.motion == name]
        assert all(sum(c[i] for c in co) == 0 for i in range(26))                   # and so do their periods, block by block
    assert [C.transfer(A)["sum_psi"] for A in C.MOTIONS.values()] == [0, -168, -3528]
    cat = [l for l in L if l.motion == "cat"]
    for l in cat:                                                                   # each cat lift and its partner through S
        mate = [m for m in cat if tuple(sorted(m.cycle)) == C.partner(l)]
        assert len(mate) == 1 and mate[0] is not l and mate[0].cls == tuple(-v for v in l.cls)


def test_classes_of_closed_lines_fill_a_lattice_of_index_6():
    G, rng = D.graph(143), random.Random(30)
    chi, vecs = G.chi_form(), []
    while len(vecs) < 60:
        g = G.random_element(rng, rng.randint(2, 9))
        if D.is_hyperbolic(g):
            v = [RM.Psi(g)] + list(G.cls(g))
            assert v[0] % 6 == D.chi6(g) == sum(c * x for c, x in zip(v[1:], chi)) % 6
            vecs.append(v)
    inv = [int(x) for x in smith_invariants(vecs)]
    assert len(inv) == 30 and math.prod(inv) == 6                                   # rank 30, index 6 in Z^30
    assert TB.record("found")["lattice"] == {"rank": 30, "index": 6}


def test_ranks_of_the_28_classes():
    L, found = C.lines_143(), TB.record("found")
    assert P.rank([l.cls for l in L]) == found["rank28"] == 16
    assert P.rank([(l.psi,) + l.cls for l in L]) == found["rank28_with_psi"] == 17
    assert {k: P.rank([l.cls for l in L if l.motion == k]) for k in C.MOTIONS} == found["rank_by_motion"] == {"cat": 4, "s10": 13, "pell": 3}


def test_elements_given_directly_and_the_systole():
    G, rng = D.graph(143), random.Random(5)
    L = C.lines_143()
    for l in (L[0], L[8], L[20], L[24]):
        h = G.random_element(rng, 5)
        v = C.invariants(D.mm(D.mm(h, l.gamma), D.inv(h)))                          # a conjugate, given as a bare matrix
        assert (v["psi"], v["cls"], v["windings"]) == (l.psi, l.cls, l.windings)
        assert abs(v["length"] - l.length) < 1e-9 * max(1.0, l.length)
    s10 = C.MOTIONS["s10"]
    v = C.invariants(s10)
    assert v["trace"] == 4 and abs(v["length"] - 2 * math.acosh(2)) < 1e-15 and abs(v["length"] - 2.633915793849633) < 1e-14
    assert (v["psi"], v["cls"]) == (L[8].psi, L[8].cls) and L[8].gamma == s10       # it fixes the dart (0 : 1): line 9
    # Matter-05's systole (ledger M05-3) is one of the four: its example has the class of line 12
    old = OL.systole_X0_143()
    w = C.invariants(old["example"])
    assert old["trace"] == 4 and w["trace"] == 4 and [n + 1 for n, l in enumerate(L) if (l.psi, l.cls) == (w["psi"], w["cls"])] == [12]
    assert all(l.trace == 4 for l in L[8:12]) and min(l.trace for l in L) == 4      # the four systoles; nothing among the 28 is shorter
    with pytest.raises(ValueError):
        C.invariants(((2, 1), (1, 1)))                                              # the cat map is not in the group
    with pytest.raises(ValueError):
        C.lifts(D.T)                                                                # not hyperbolic


def test_lines_at_other_levels():
    for N, index in ((11, 12), (23, 24), (47, 48)):
        for A in C.MOTIONS.values():
            t = C.transfer(A, N)
            assert sum(t["cycle_type"]) == index == t["darts"]
            assert not any(t["sum_cls"]) and t["sum_psi"] == index * RM.Psi(A)
    assert [l.ell for l in C.lifts(C.MOTIONS["pell"], 11)] == [1, 11] and C.lifts(C.MOTIONS["cat"], 11)[0].windings is None
    with pytest.raises(ValueError):
        C.lifts(C.MOTIONS["cat"], 13)                                               # elliptic points: no graph
