"""TR-01: the two mirror classes (two routes), CM points on mirrors (two routes), the thirteen star-quotient curves, L-parity, the fold, the sealed I.4 test."""
import contextlib
import io
import json
from fractions import Fraction

import pytest

from mtft.matter import ledger as L
from mtft.matter import puzzle as Z
from mtft.matter import star_quotients as S
from mtft.matter import time_reversal as T
from mtft.matter.__main__ import main as matter_main

ONE = [57, 58, 65, 77, 91, 118, 123, 143, 155]
TWO = [82, 141, 142, 145]


def test_reflections_two_routes_and_conjugation():
    refl = T.reflections(12)
    assert len(refl) == 212 and all(T.is_reflection(m) for m in refl)
    assert all(T.mirror_class(m) == T.mirror_class(m, "farey") for m in refl)                        # content 2 <-> Farey neighbours
    assert all(T.content(T.fixed_form(m)) in (1, 2) and T.farey_det(m) in (1, 2) for m in refl)
    assert {T.mirror_class(m) for m in refl} == {T.THETA_0, T.THETA_PI}
    gs = (((1, 1), (0, 1)), ((0, -1), (1, 0)), ((1, 0), (0, -1)), ((2, 1), (1, 1)), ((3, -2), (-4, 3)))
    assert all(T.mirror_class(T.conjugate(m, g)) == T.mirror_class(m) for m in refl[:300] for g in gs)  # a conjugacy invariant
    assert all(T.is_reflection(T.conjugate(m, g)) for m in refl[:300] for g in gs)
    import cmath
    for th in (0.3, 1.2, 2.5):                                                                           # points of |tau + 1| = 1
        z = -1 + cmath.exp(1j * th)
        assert abs(T.reflect(T.STANDARD["|tau + 1| = 1 (the oloid move q -> -q/(1 + q))"], z) - z) < 1e-12
    assert abs(T.reflect(T.STANDARD["Re tau = -1/2 (the triangle mirror n -> -1 - n)"], -0.5 + 2j) - (-0.5 + 2j)) < 1e-12
    assert T.mirror_endpoints(((1, 0), (0, -1))) == (Fraction(0), None) and T.mirror_endpoints(((-1, 0), (1, 1))) == (Fraction(-2), Fraction(0))


def test_the_dictionary_and_the_standard_mirrors():
    classes = {k: T.mirror_class(m) for k, m in T.STANDARD.items()}
    assert classes["Re tau = 0 (theta = 0)"] == T.THETA_0
    assert all(v == T.THETA_PI for k, v in classes.items() if k != "Re tau = 0 (theta = 0)")       # the oloid and triangle mirrors: theta = pi
    assert T.physics_dictionary()["time_reversal"] == "tau -> -conj(tau)"


def test_cm_points_two_routes():
    for D in (-3, -4, -7, -8, -11, -12, -15, -20, -23, -24, -35, -40, -44, -52, -143, -572):
        forms = T.reduced_forms(D)
        assert len(forms) == S.class_number(D)                                                        # SQ-01's class numbers
        for f in forms:
            fixing = T.reflections_fixing(f)
            classes = {T.mirror_class(m) for m in fixing}
            route2 = None if not classes else ("both" if len(classes) == 2 else next(iter(classes)))
            assert T.cm_mirror_class(f) == route2
            assert all(abs(T.reflect(m, T.cm_point(f)) - T.cm_point(f)) < 1e-9 for m in fixing)
    assert T.cm_mirror_class((1, 0, 1)) == "both" and T.cm_mirror_class((1, 1, 1)) == T.THETA_PI
    red, g = T.reduce_form((11, 11, 6))
    assert red == (6, 1, 6) and abs(T.cm_point(red) - (g[0][0] * T.cm_point((11, 11, 6)) + g[0][1]) / (g[1][0] * T.cm_point((11, 11, 6)) + g[1][1])) < 1e-12
    import random
    random.seed(3)
    for _ in range(1500):                                                                                  # the reduction tracks the point and ends on the boundary cases
        A, B, C = random.randint(1, 40), random.randint(-80, 80), random.randint(1, 60)
        if B * B - 4 * A * C >= 0:
            continue
        (a, b, c), g = T.reduce_form((A, B, C))
        assert -a < b <= a <= c and (b >= 0 or (abs(b) != a and a != c)) and b * b - 4 * a * c == B * B - 4 * A * C
        z = T.cm_point((A, B, C))
        assert abs(T.cm_point((a, b, c)) - (g[0][0] * z + g[0][1]) / (g[1][0] * z + g[1][1])) < 1e-9
    assert T.reduce_form((7, -7, 7))[0] == (7, 7, 7) and T.reduce_form((2, -2, 7))[0] == (2, 2, 7)
    fp = T.fricke_fixed_point_classes()
    assert fp[-143]["real"] == [((1, 1, 36), T.THETA_PI), ((6, 1, 6), T.THETA_PI)]                 # W_143's real points: two of each class
    assert fp[-572]["real"] == [((1, 0, 143), T.THETA_0), ((11, 0, 13), T.THETA_0)]
    assert fp[-52]["real"] == [((1, 0, 13), T.THETA_0), ((2, 2, 7), T.THETA_PI)]                   # W_13's four points lie over these two classes
    with pytest.raises(ValueError):
        T.reduced_forms(5)
    with pytest.raises(ValueError):
        T.cm_mirror_class((1, 3, 1))


def test_curves():
    t = T.curve_table()
    assert T.discriminant((0, -1, 1, -1, -2)) == -1859 and T.discriminant((1, 0, 0, 4, 1)) == -65 ** 2 and T.discriminant((0, 0, 1, -1, 0)) == 37
    assert [N for N, r in t.items() if r["real_circles"] == 1] == ONE and [N for N, r in t.items() if r["real_circles"] == 2] == TWO
    assert all(r["lattice"] == (T.THETA_PI if N in ONE else T.THETA_0) for N, r in t.items())
    assert all(r["root_number"] == -1 and r["rank"] == 1 for r in t.values()) and sorted(t) == list(S.tori())
    assert all(list(r["class_moddeg"].values()).count(4) == 1 and r["class_moddeg"][r["curve"]] == 4 for r in t.values())
    assert all(d in (2, 4, 8) for r in t.values() for d in r["class_moddeg"].values())                  # no degree-1 curve: 4 = m^2 d forces m = 1
    assert t[65]["curve"] == "65a2" and t[65]["class_moddeg"] == {"65a1": 2, "65a2": 4}                # the quotient is not the optimal curve
    assert T.real_circles(T.COMPARISON_CURVES["11a1"]["a"]) == 1 and T.real_circles(T.COMPARISON_CURVES["37a1"]["a"]) == 2
    assert "odd" in T.l_parity(-1) and "even" in T.l_parity(1)
    with pytest.raises(ValueError):
        T.l_parity(0)
    assert T.circles_census() == {"one_circle": ONE, "two_circles": TWO, "all_odd": True, "moddeg_4_unique": True}


def test_fold_and_strong_cp():
    a, b = T.fold_preimages(2.0, 4.0)
    assert abs((a + b) - (-1)) < 1e-12 and abs(4 * a * (a + 1) - 2) < 1e-12 and T.fold_preimages(-5.0) == ()
    cp = T.strong_cp_test
    assert cp(("cm", (1, 0, 143)), det_phase=0.0)["verdict"] == "PASS"                                   # theta = 0 class, arg det 0
    assert cp(("cm", (1, 1, 36)), det_phase=T.PI)["verdict"] == "PASS"                                   # theta = pi class, arg det pi (CC-63)
    assert cp(("cm", (1, 0, 143)), det_phase=T.PI)["verdict"] == "FAIL: theta-bar = pi"
    assert cp(("cm", (1, 1, 36)), det_phase=0.0)["verdict"] == "FAIL: theta-bar = pi"
    assert cp(("cm", (1, 0, 1)), det_phase=0.0)["verdict"] == "PASS" and cp(("cm", (1, 0, 1)), det_phase=T.PI)["verdict"] == "PASS"   # tau = i: both
    assert cp(("cm", (1, 0, 143)))["verdict"] == "pending: class theta_0 needs arg det(M_u M_d) = 0"
    assert cp(("cm", (1, 1, 36)))["verdict"] == "pending: class theta_pi needs arg det(M_u M_d) = pi"
    assert cp(("cm", (1, 0, 143)), det_phase=0.3)["verdict"] == "FAIL: theta-bar = 0.3"
    assert cp(("cm", (2, 1, 18)), det_phase=0.0)["verdict"] == "FAIL: off every mirror (needs an axion)"
    assert cp(("cusp", None), det_phase=0.0)["verdict"].startswith("not a candidate")
    r = cp(("numeric", complex(-0.5, 1.0232745926964613)), tolerance=1e-9)
    assert r["mirror_class"] == T.THETA_PI and not r["exact"] and r["verdict"].startswith("pending")
    assert cp(("numeric", 2.0j), det_phase=0.0, tolerance=1e-9)["verdict"] == "pass to tolerance only"
    with pytest.raises(ValueError):
        cp(("numeric", 2.0j))
    with pytest.raises(ValueError):
        cp(("point", 1))
    for bad in (((1, 0), (0, 1)), ((2, 1), (1, -2)), ((1, 0, 0), (0, -1)), ((1.0, 0), (0, -1))):
        assert not T.is_reflection(bad)
    with pytest.raises(ValueError):
        T.fixed_form(((1, 0), (0, 1)))
    with pytest.raises(ValueError):
        T.conjugate(((1, 0), (0, -1)), ((2, 0), (0, 1)))
    with pytest.raises(ValueError):
        T.mirror_class(((1, 0), (0, -1)), route="guess")
    with pytest.raises(ValueError):
        T.reflect(((1, 0), (0, 1)), 1j)


def test_command_line():
    def run(*argv):
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            code = matter_main(["timerev", *argv])
        return code, buf.getvalue()
    code, out = run(); assert code == 0 and "65a2" in out and "theta_pi" in out and "Pre-registered test" in out
    code, out = run("--json"); d = json.loads(out); assert code == 0 and d["census"]["one_circle"] == ONE and d["sweep"]["routes_agree"] is True
    code, out = run("--candidate", "cm", "1,0,143", "--det-phase", "0"); assert code == 0 and "PASS" in out
    code, out = run("--candidate", "cm", "1,1,36"); assert code == 0 and "pending" in out
    code, out = run("--candidate", "cm", "1,1,36", "--det-phase", "0"); assert code == 0 and "theta-bar = pi" in out
    code, out = run("--candidate", "numeric", "0.5+1.2j"); assert code == 2


def test_the_record():
    E = {e["id"]: e for e in L.E}
    tr = [e for e in L.E if e["wave"] == "TR-01"]
    assert [e["id"] for e in tr] == [f"TR-{i}" for i in (1, 2, 3, 4, 5, 6, 8, 7)] and all(e["module"] == "mtft.matter.time_reversal" for e in tr)
    assert E["TR-6"]["status"] == "superseded" and "CC-63" in E["TR-6"]["note"] and E["TR-8"]["status"] == "active" and "PR-S2" in E["TR-8"]["statement"]
    a = Z.assignment()
    assert all(a[f"TR-{i}"] == "II" for i in range(1, 9))                                              # piece I keeps no ledger claims; I.4 cites TR-8
    subs = Z.subs()
    assert subs["I.4"][1].prereg == "T6" and {"TR-8", "CC-63"} <= set(subs["I.4"][1].claims) and subs["I.4"][1].status == "open"
    assert subs["II.10"][1].status == "answered" and "TR-3" in subs["II.10"][1].claims
    assert next(t for t in Z.PREREG if t["id"] == "T6")["outcome"] == "pending"


def test_pari_reproduces_the_frozen_curves():
    cypari2 = pytest.importorskip("cypari2")
    pari = cypari2.Pari(); pari.allocatemem(10 ** 9)
    for N in (65, 143, 145):
        c = T.STAR_QUOTIENT_CURVES[N]
        E = pari(f'ellinit("{c["curve"]}")')
        assert [int(E[i]) for i in range(5)] == list(c["a"]) and int(E[11]) == T.discriminant(c["a"])
        assert int(pari.ellmoddegree(E)) == 4 and int(pari.ellrootno(E)) == -1 and int(pari.elltors(E)[0]) == c["torsion"]
        lam = pari.lfunlambda(pari.lfuncreate(E), pari("1+I"))
        assert abs(float(pari.real(lam))) < 1e-12 and abs(float(pari.imag(lam))) > 0.1                   # odd under t -> -t
    E11 = pari('ellinit("11a1")'); lam = pari.lfunlambda(pari.lfuncreate(E11), pari("1+I"))
    assert abs(float(pari.imag(lam))) < 1e-12 and int(pari.ellrootno(E11)) == 1                          # even
