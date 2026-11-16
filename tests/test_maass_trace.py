"""mtft.maass.trace: the Selberg trace formula by sector, its frozen weights, and what it certifies about the table.

Gates of the v1.1.0 plan (Tool 1, second half).  Nothing here uses a Maass form at level 1; level 11 uses three
forms the solver finds there; level 143 uses the frozen table.
"""
import hashlib
import math
from pathlib import Path

import mpmath as mp
import numpy as np
import pytest

from mtft.maass import trace as T

DATA = Path(T.__file__).resolve().parent / "_data"


def test_frozen_weights_and_their_provenance():
    want = [ln.split()[0] for ln in (DATA / "PROVENANCE.txt").read_text(encoding="utf-8").splitlines() if ln.endswith("trace_weights.npz")]
    assert want == [hashlib.sha1((DATA / "trace_weights.npz").read_bytes()).hexdigest()]
    assert T.length_bound() == 17.0 and T.frozen_levels() == (1, 11, 143)
    sizes = {(N, Q, kind): T.classes(N, Q, kind)[0].size for N, Q in ((143, 1), (143, 11), (143, 13), (143, 143)) for kind in ("hyp", "gli")}
    assert sizes == {(143, 1, "hyp"): 4912, (143, 1, "gli"): 4914, (143, 11, "hyp"): 1481, (143, 11, "gli"): 1481,
                     (143, 13, "hyp"): 1363, (143, 13, "gli"): 1363, (143, 143, "hyp"): 410, (143, 143, "gli"): 410}
    s, length, w = T.classes(143, 1, "hyp")
    assert s[0] == 3 and abs(length[0] - 2 * math.acosh(1.5)) < 1e-15 and length[-1] <= 17.0 and np.all(np.diff(length) > 0)
    assert abs(length[1] - 2 * math.acosh(2.0)) < 1e-15 and w[1] > 0             # trace 4: the systole, 2 arccosh 2
    assert T.classes(143, 1, "hyp", twist=True)[0].size == 4912
    with pytest.raises(ValueError):
        T.classes(13, 1)


def test_class_numbers_by_a_second_route():
    """Dirichlet's finite sum, which needs no PARI and no GRH, against closed forms and against the frozen weights."""
    assert abs(T.l1_dirichlet(5) - 2 * math.log((1 + math.sqrt(5)) / 2) / math.sqrt(5)) < 1e-13
    assert abs(T.l1_dirichlet(8) - math.log(1 + math.sqrt(2)) / math.sqrt(2)) < 1e-13
    assert abs(T.l1_dirichlet(12) - 2 * math.log(2 + math.sqrt(3)) / math.sqrt(12)) < 1e-13
    assert abs(T.l1_dirichlet(572) - 4 * math.log(12 + math.sqrt(143)) / math.sqrt(572)) < 1e-12       # the Pell unit 12 + sqrt 143
    assert T.core_discriminant(45) == (5, 3) and T.core_discriminant(12) == (12, 1) and T.core_discriminant(4 * 143) == (572, 1)
    assert list(T.kronecker_table(5)) == [0, 1, -1, -1, 1] and list(T.kronecker_table(8)) == [0, 1, 0, -1, 0, -1, 0, 1]
    assert [T.kronecker(5, p) for p in (2, 3, 5, 11)] == [-1, -1, 0, 1]
    worst, n = 0.0, 0
    for N, Q in ((1, 1), (11, 1), (11, 11), (143, 1), (143, 11), (143, 13), (143, 143)):
        for kind, sign in (("hyp", -1), ("gli", 1)):
            s, _, w = T.classes(N, Q, kind)
            for si, wi in list(zip(s, w))[:10]:
                v = 2 * math.sqrt(Q) * T.weight(Q * Q * int(si) ** 2 + sign * 4 * Q, N, Q)
                worst, n = max(worst, abs(v - wi) / max(1.0, abs(wi))), n + 1
    s, _, w = T.classes(143, 1, "hyp", twist=True)
    signed = 0
    for si, wi in list(zip(s, w))[:10]:
        v = 2 * T.weight(int(si) ** 2 - 4, 143, 1, twist_trace=int(si))
        worst, n = max(worst, abs(v - wi) / max(1.0, abs(wi))), n + 1
        signed += abs(v - 2 * T.weight(int(si) ** 2 - 4, 143, 1)) > 1e-9
    assert n == 150 and worst < 1e-9 and signed >= 1                # the twist changes some weights: not the same table twice
    lev = T.Level(143)
    assert lev.E == {11: 0.0, 13: 4.0, 143: 20.0} and {Q: T.fixed_points(143, Q) for Q in lev.Qs} == lev.E      # 20 = h(-572) + h(-143)
    assert T.class_number_negative(-143) == 10 and T.class_number_negative(-572) == 10
    assert all(abs(T.mirror_length(143, Q) - lev.L[Q]) < 1e-9 for Q in lev.Qs)
    assert abs(lev.L[143] - 4 * math.log(12 + math.sqrt(143))) < 1e-9          # the Fricke mirrors: twice the Pell geodesic's length
    assert lev.cusps == 4 and len(T.mirror_ends(143)) == 8 and abs(lev.sum_beta - 22.62396724327941) < 1e-12
    for bad in (lambda: T.l1_dirichlet(1), lambda: T.core_discriminant(-4), lambda: T.weight(9, 1, 1), lambda: T.kronecker_table(20)):
        with pytest.raises(ValueError):
            bad()


def test_integrals_by_second_routes():
    zs = np.array([1 + 0.3j, 0.5 + 2.7j, 1 + 9j, 0.5 + 0.01j, 1 + 0j])
    assert max(abs(complex(mp.digamma(complex(z))) - w) for z, w in zip(zs, T.digamma(zs))) < 1e-14
    tf = T.TestFunction(0.0, 0.5)
    I = tf.integrals()
    # Psi_1 from g instead of h:  -gamma g(0) + int_0^inf (g(0) - g(u)) / (e^u - 1) du
    g0 = float(tf.g(0.0))
    second = -T.EULER * g0 + float(mp.quad(lambda u: (g0 - float(tf.g(float(u)))) / mp.expm1(u), mp.linspace(mp.mpf("1e-30"), 60, 61)))
    assert abs(I["psi 1"] - second) < 1e-10
    assert abs(T.continuum_level1(tf, "series") - T.continuum_level1(tf, "quad")) < 1e-10        # the explicit formula against zeta'/zeta
    assert abs(tf.h_i_half - 2 * math.exp(0.25 / 0.5)) < 1e-14 and abs(float(tf.g(0.0)) - 2 * 0.5 / math.sqrt(2 * math.pi)) < 1e-15
    off = T.TestFunction(0.7, 0.45)
    num = float(mp.quad(lambda r: (mp.exp(-(r - 0.7) ** 2 / 0.405) + mp.exp(-(r + 0.7) ** 2 / 0.405)) * mp.cos(1.3 * r), mp.linspace(0, 8, 17))) / math.pi
    assert abs(num - float(off.g(1.3))) < 1e-12                        # g is the transform of h
    with pytest.raises(ValueError):
        T.continuum_level1(tf, "other")
    with pytest.raises(ValueError):
        T.TestFunction(0.0, 0.0)


def test_level_1_the_whole_formula_cancels():
    """PSL(2, Z) has no cusp form below R = 9.53, so seven pieces of size 0.06 to 3.7 must sum to zero."""
    for sigma, tol in ((0.45, 1e-9), (0.6, 1e-12)):
        p = T.level1_total(T.TestFunction(0.0, sigma))
        assert abs(p["total"]) < tol, (sigma, p["total"])
        assert min(abs(v) for k, v in p.items() if k != "total") > 0.05
    assert abs(T.level1_total(T.TestFunction(0.0, 0.6), cut=6.0)["total"]) > 1e-4       # without the long geodesics it does not
    with pytest.raises(ValueError):
        T.level1_total(T.TestFunction(0.0, 0.35))                     # needs geodesics longer than the frozen 17


def test_level_11_the_three_forms_of_the_record():
    lists = {(0, (1,)): [2.48359105955], (1, (1,)): [2.96820576382], (1, (-1,)): [2.03309099386], (0, (-1,)): []}
    r = T.check_lists(lists, T.TestFunction(0.0, 0.45), 11, 3.0)
    assert all(abs(v["residual"]) < 2e-9 for v in r.values())
    assert r[(1, (-1,))]["list"] > 7e-5 and r[(0, (1,))]["list"] > 4e-7           # the forms are seen: 1e-5 and 1e-3 of their size
    wide = T.check_lists(lists, T.TestFunction(0.0, 0.6), 11, 3.0)
    assert all(v["positive"] and 0 <= v["residual"] < 1e-5 and v["radius"] > 3.0 for v in wide.values())
    swapped = dict(lists)
    swapped[(1, (1,))], swapped[(1, (-1,))] = lists[(1, (-1,))], lists[(1, (1,))]   # the two odd forms in each other's sector
    bad = T.check_lists(swapped, T.TestFunction(0.0, 0.6), 11, 3.0)
    assert not bad[(1, (1,))]["positive"] and bad[(1, (-1,))]["residual"] > 5e-3
    with pytest.raises(ValueError):
        T.Level(13)


def test_level_143_every_sector_is_positive_and_complete_to_1p95():
    for sigma, floor in ((0.43, 1.95), (0.45, 1.95), (0.5, 1.93)):
        r = T.check_table(T.TestFunction(0.0, sigma))
        assert len(r) == 8 and all(v["positive"] for v in r.values())
        assert all(0 < v["residual"] < 2e-3 for v in r.values())       # the forms above R = 2, and nothing else
        assert min(v["radius"] for v in r.values()) > floor, (sigma, {k: round(v["radius"], 3) for k, v in r.items()})
    c = T.certificate()
    assert c["all_positive"] and c["no_exceptional_eigenvalue"] and c["lengths_to"] == 17.0
    assert c["untwisted_radius"] > 1.95 and c["twisted_radius"] > 1.6 and c["largest_residual"] < 0.03
    assert max(v["radius"] for v in c["sectors"].values()) > 2.15     # some sectors are certified past the list's top
    with pytest.raises(ValueError):
        T.check_table(T.TestFunction(0.0, 0.35))                      # the lab's width needs lengths to 23.5
    with pytest.raises(ValueError):
        T.completeness(1e-3, T.TestFunction(0.5, 0.5))


def test_controls_a_missing_or_a_spurious_form_is_seen():
    tf = T.TestFunction(0.0, 0.5)
    L = T.table_lists(False)
    key = (0, (1, 1))
    dropped = dict(L)
    dropped[key] = [R for R in L[key] if abs(R - 0.8992552254) > 1e-6]
    r = T.check_lists(dropped, tf)[key]
    assert len(dropped[key]) == len(L[key]) - 1 and r["residual"] > 0.39 and abs(r["radius"] - 0.8992552254) < 2e-3     # and it says where
    added = dict(L)
    added[(1, (1, 1))] = L[(1, (1, 1))] + [1.5]
    assert T.check_lists(added, tf)[(1, (1, 1))]["radius"] is None and not T.check_lists(added, tf)[(1, (1, 1))]["positive"]
    moved = dict(L)                                                    # the first form, put in the wrong parity
    moved[(1, (-1, -1))], moved[(0, (-1, -1))] = L[(1, (-1, -1))][1:], L[(0, (-1, -1))] + L[(1, (-1, -1))][:1]
    m = T.check_lists(moved, tf)
    assert m[(1, (-1, -1))]["residual"] > 1.0 and not m[(0, (-1, -1))]["positive"]
    short = T.check_table(tf, cut=8.0)                                 # and with the geodesics cut at length 8 the gate fails
    assert not all(v["positive"] and v["residual"] < 2e-3 for v in short.values())


def test_twisted_total():
    lev = T.Level(143)
    for sigma, floor in ((0.43, 1.6), (0.45, 1.57), (0.5, 1.48)):
        tw = T.check_twisted(T.TestFunction(0.0, sigma))
        assert tw["forms"] == 42 and tw["positive"] and 0 < tw["residual"] < 0.03 and tw["radius"] > floor
    parts = lev.twisted_total(T.TestFunction(0.0, 0.5))
    assert abs(parts["total"] - (parts["area"] + parts["geodesics"] + parts["cusps"])) < 1e-12 and parts["cusps"] < 0
    plain = float(np.sum(T.TestFunction(0.0, 0.5).g(T.classes(143, 1, "hyp")[1]) * T.classes(143, 1, "hyp")[2]))
    assert abs(parts["geodesics"] - plain) > 0.1                       # the twist's signs matter
    with pytest.raises(ValueError):
        T.Level(11).twisted_total(T.TestFunction(0.0, 0.5))
