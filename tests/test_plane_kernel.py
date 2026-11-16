"""mtft.plane.kernel: one real 2 x 2 matrix read eight ways.  Expected numbers are the lab's pre-registered ones (record
MP06, 8 November 2026) or are derived in the test by a second route; nothing here is a value written from memory."""
import json
import math
import random
from fractions import Fraction as Fr
from pathlib import Path

import numpy as np
import pytest

from mtft import plane as PL
from mtft.plane import kernel as K

PELL = ((12, 1), (143, 12))
CAT = ((2, 1), (1, 1))
PHI = (1 + math.sqrt(5)) / 2
LAB_CORE = Path(__file__).resolve().parents[1] / "studies" / "plane_2026nov" / "lab_core_outputs.json"


def _flat(pairs):
    """A tuple of points as one flat list, for pytest.approx (which takes no nested sequences)."""
    return [v for pt in pairs for v in pt]


def _support(M, phi):
    """The support function of W(M) with no formula for the range: the top eigenvalue of the Hermitian part of e^{-i phi} M."""
    A = np.exp(-1j * phi) * np.array(M, dtype=complex)
    return float(np.linalg.eigvalsh((A + A.conj().T) / 2)[-1])


def test_the_pell_unit_read_eight_ways():
    r = PL.read(PELL, level=143)
    assert tuple(r["readings"]) == K.READINGS == ("motion", "form", "number", "flow", "field", "vacuum", "two", "range")
    assert all(r["readings"][k]["label"] == K.LABELS[k] for k in K.READINGS)
    inv = r["invariants"]
    assert (inv["tr"], inv["det"], inv["disc"], r["kind"], inv["unimodular"]) == (24, 1, 572, "hyp", True)
    m = r["readings"]["motion"]
    assert abs(m["rapidity"] - 3.176313) < 5e-7 and abs(m["length"] - 6.352626) < 5e-7 and abs(m["multiplier"] - (12 + math.sqrt(143))) < 1e-12
    assert m["psi"] == -21 and m["cycle_type"] == [1, 11, 13, 143] and m["lorenz_word"] == "x" + "y" * 22
    assert m["fixed"]["ends"] == pytest.approx((-1 / math.sqrt(143), 1 / math.sqrt(143)), abs=1e-15)      # the Fricke circle's ends
    R = r["readings"]["range"]                                                    # MP06: centre 12, semi-axes 72 and 71, w = 84
    assert (R["shape"], R["centre"], R["semi_real"], R["semi_imag"], R["minor_axis"], R["radius"]) == ("ellipse", 12, 72.0, 71.0, 142.0, 84.0)
    assert R["semi_real_sq"] - R["semi_imag_sq"] == 143 and abs(R["ratio"] - 1.7144) < 5e-5 and R["ratio"] <= 2
    assert _flat(R["foci"]) == pytest.approx([12 - math.sqrt(143), 0.0, 12 + math.sqrt(143), 0.0])
    f = r["readings"]["form"]
    assert f["form"] == (143, 0, -1) and f["disc"] == 572 and f["content"] == 1 and f["type"] == "indefinite" and f["momentum"]["kind"] == "spacelike"
    n = r["readings"]["number"]
    assert n["system"] == "split-complex" and n["kappa"] == 1 and n["norm"] == 1 and abs(n["angle"] - m["rapidity"]) < 1e-12
    assert abs(n["a"] ** 2 - n["b"] ** 2 - 1) < 1e-9                              # the norm a² − κb² is the determinant
    v = r["readings"]["vacuum"]
    assert v["symplectic"] and v["quanta"] == 5184 == 72 ** 2 and abs(math.sinh(v["distance"] / 2) ** 2 - 5184) < 1e-6
    assert r["readings"]["two"]["gap"] == 144.0 and r["readings"]["two"]["rotation"] == 71
    assert r["readings"]["field"]["bogoliubov"] == "a squeeze of a mode" and r["readings"]["flow"]["wins"].startswith("strain wins")
    j = K.jacobi_monodromy(PELL)                                                  # MP06: L = 6.352626, e^L = 573.998, 2 cosh L = tr M² = 574
    assert j["trace_M2"] == 574 and abs(j["two_cosh_length"] - 574) < 1e-9 and abs(j["exp_length"] - 573.998) < 1e-3 and abs(j["length"] - 6.352626) < 5e-7
    assert K.jacobi_monodromy(CAT)["trace_M2"] == 7 and K.jacobi_monodromy(((0, -1), (1, 0)))["defined"] is False


def test_numerical_ranges_of_the_lab_matrices():
    e = K.numerical_range(((0, 1), (0, 0)))                                       # the dual unit: Crouzeix's 2 is sharp
    assert (e["shape"], e["semi_real"], e["semi_imag"], e["radius"], e["norm"], e["ratio"]) == ("disk", 0.5, 0.5, 0.5, 1.0, 2.0)
    c = K.numerical_range(CAT)
    assert c["shape"] == "segment" and c["normal"] and _flat(c["foci"]) == pytest.approx([1 / PHI ** 2, 0.0, PHI ** 2, 0.0])
    assert c["radius"] == pytest.approx(PHI ** 2) and c["ratio"] == pytest.approx(1.0)
    for M in (((0, -1), (1, 0)), ((0, 1), (1, 0)), ((0, 1), (-1, 0)), ((1, 1), (1, 0))):      # S, G(1), G(−1), Fibonacci: normal
        assert K.numerical_range(M)["shape"] == "segment" and K.numerical_range(M)["ratio"] == pytest.approx(1.0)
    for kap in (-3.0, -1.0, -0.5, 0.0, 0.3, 1.0, 2.5):                            # G(κ) = [[0, 1], [κ, 0]]
        g = K.numerical_range(((0, 1), (kap, 0)))
        assert g["semi_major"] == pytest.approx((1 + abs(kap)) / 2) and g["semi_minor"] == pytest.approx(abs(1 - abs(kap)) / 2)
        assert g["norm"] == pytest.approx(max(1.0, abs(kap)))
        if abs(kap) <= 1:
            assert g["ratio"] == pytest.approx(2 / (1 + abs(kap)))
    assert K.numerical_range(((3, 0), (0, 3)))["shape"] == "point" and K.numerical_range(((3, 0), (0, 3)))["radius"] == 3.0
    rot = K.numerical_range(((1.2, -1.6), (1.6, 1.2)))                            # a scaled rotation, in floats: normal, a segment
    assert rot["shape"] == "segment" and _flat(rot["foci"]) == pytest.approx([1.2, -1.6, 1.2, 1.6]) and rot["radius"] == pytest.approx(2.0)
    assert len(K.range_boundary(PELL, 12)) == 13 and K.range_boundary(PELL, 12)[0] == pytest.approx((84.0, 0.0))


def test_the_range_is_strain_and_rotation_on_random_matrices():
    """Three routes: the closed forms here, the lab core's own computation, and numpy with no formula for the range."""
    rng = random.Random(8)
    worst = {"lab": 0.0, "support": 0.0, "radius": 0.0, "norm": 0.0}
    shapes = set()
    for trial in range(500):
        scale = 10 ** rng.uniform(-1, 2)
        M = tuple(tuple(rng.uniform(-1, 1) * scale for _ in range(2)) for _ in range(2))
        if trial % 5 == 0:                                                        # some with complex eigenvalues for certain
            M = ((M[0][0], -abs(M[0][1]) - scale), (abs(M[1][0]) + scale, M[1][1]))
        R, L, dec = K.numerical_range(M), K.lab_numerical_range(M), K.decompose(M)
        shapes.add((R["shape"], R["real_foci"]))
        assert R["shape"] == L["shape"] and R["real_foci"] == L["real_foci"]
        assert R["semi_real"] == pytest.approx(dec["strain"]) and R["semi_imag"] == pytest.approx(abs(dec["w"]))
        s = max(1.0, R["norm"])
        worst["lab"] = max(worst["lab"], abs(R["semi_major"] - L["semi_major"]) / s, abs(R["semi_minor"] - L["semi_minor"]) / s, abs(R["norm"] - L["norm"]) / s)
        worst["radius"] = max(worst["radius"], abs(R["radius"] - L["radius"]) / s)
        worst["norm"] = max(worst["norm"], abs(R["norm"] - float(np.linalg.norm(np.array(M), 2))) / s)
        alpha = float(dec["alpha"])
        for k in range(12):
            phi = 2 * math.pi * (k + 0.37) / 12
            mine = alpha * math.cos(phi) + math.sqrt(R["semi_real"] ** 2 * math.cos(phi) ** 2 + R["semi_imag"] ** 2 * math.sin(phi) ** 2)
            worst["support"] = max(worst["support"], abs(mine - _support(M, phi)) / s)
        grid = max(_support(M, 2 * math.pi * k / 720) for k in range(720))
        assert grid <= R["radius"] * (1 + 1e-12) and R["radius"] - grid < 2e-5 * s  # the radius is the top of the support function
        assert R["ratio"] <= 2 + 1e-12                                            # Crouzeix, for p(z) = z
    assert worst["lab"] < 1e-12 and worst["norm"] < 1e-12 and worst["support"] < 1e-9 and worst["radius"] < 1e-9, worst
    assert shapes == {("ellipse", True), ("ellipse", False)}                      # both orientations of the foci were met


def test_quanta_equal_strain_squared_for_determinant_one():
    rng = random.Random(3)
    for _ in range(200):
        L = ((rng.uniform(-1, 1), rng.uniform(-2, 2)), (rng.uniform(-2, 2), 0.0))
        L = ((L[0][0], L[0][1]), (L[1][0], -L[0][0]))                             # traceless: exp(tL) has determinant 1
        M = K.expm(L, rng.uniform(0.2, 2.0))
        assert abs(M[0][0] * M[1][1] - M[0][1] * M[1][0] - 1) < 1e-9
        r = PL.read(M)
        v, dec = r["readings"]["vacuum"], r["decomposition"]
        assert v["symplectic"] and v["quanta"] == pytest.approx(dec["strain_sq"], rel=1e-9, abs=1e-12)
        assert math.sinh(v["distance"] / 2) ** 2 == pytest.approx(dec["strain_sq"], rel=1e-7, abs=1e-9)
        z = complex(*v["image_of_i"])                                             # the image of i, by the Möbius map itself
        assert z == pytest.approx((M[0][0] * 1j + M[0][1]) / (M[1][0] * 1j + M[1][1]))
    exact = PL.read(((Fr(3, 2), Fr(1, 2)), (Fr(5, 2), Fr(3, 2))))                 # exact entries stay exact
    assert exact["invariants"]["det"] == 1 and exact["readings"]["vacuum"]["quanta"] == Fr(9, 4) == exact["decomposition"]["strain_sq"]


def test_motions_by_kind():
    S = PL.read(((0, -1), (1, 0)))
    assert S["kind"] == "ell" and S["readings"]["motion"]["order"] == 4 and S["readings"]["motion"]["fixed"]["tau"] == (0.0, 1.0)
    assert S["readings"]["motion"]["angle"] == pytest.approx(math.pi) and S["readings"]["number"]["system"] == "complex"
    assert S["readings"]["form"]["type"] == "definite" and S["readings"]["form"]["momentum"]["kind"] == "timelike"
    R3 = PL.read(((0, 1), (-1, -1)))
    assert R3["readings"]["motion"]["order"] == 3 and R3["readings"]["motion"]["fixed"]["tau"] == pytest.approx((-0.5, math.sqrt(3) / 2))
    assert PL.read(((1, -1), (1, 0)))["readings"]["motion"]["order"] == 6
    T = PL.read(((1, 1), (0, 1)), level=143)
    assert T["kind"] == "par" and T["readings"]["motion"]["fixed"]["cusp"] == math.inf and T["readings"]["motion"]["order"] == math.inf
    assert T["readings"]["motion"]["cycle_type"] == [1, 11, 13, 143] and T["readings"]["number"]["system"] == "dual"
    assert PL.read(((1, 0), (5, 1)))["readings"]["motion"]["fixed"]["cusp"] == 0
    one, minus = PL.read(((1, 0), (0, 1))), PL.read(((-1, 0), (0, -1)))
    assert one["invariants"]["scalar"] and one["readings"]["motion"]["order"] == 1 and minus["readings"]["motion"]["order"] == 2
    assert one["readings"]["number"]["system"] == "real" and one["readings"]["range"]["shape"] == "point"
    fib = PL.read(((1, 1), (1, 0)))                                               # determinant −1: not a motion of the plane
    assert fib["kind"] == "hyp" and fib["readings"]["motion"]["acts_on_the_plane"] is False and "rapidity" not in fib["readings"]["motion"]
    assert fib["readings"]["vacuum"]["symplectic"] is False and fib["readings"]["number"]["angle"] is None
    cat = PL.read(CAT, level=143)
    assert cat["readings"]["motion"]["cycle_type"] == [7] * 4 + [35] * 4 and cat["readings"]["motion"]["psi"] == 0
    assert cat["readings"]["motion"]["length"] == pytest.approx(2 * math.log(PHI ** 2))
    half = PL.read(((2.0, 0.0), (0.0, 0.5)))
    assert half["readings"]["motion"]["fixed"]["ends"] == (0.0, math.inf) and half["readings"]["motion"]["length"] == pytest.approx(math.log(4))
    with pytest.raises(ValueError):
        PL.read(((2.0, 0.0), (0.0, 0.5)), level=11)


def test_generators_and_their_motions():
    rng = random.Random(21)
    for _ in range(200):
        L = tuple(tuple(rng.uniform(-1.5, 1.5) for _ in range(2)) for _ in range(2))
        t = rng.uniform(-1.5, 1.5)
        E = np.array(K.expm(L, t))
        A = np.array(L) * t                                                       # a second route: the power series, by numpy
        series, term = np.eye(2), np.eye(2)
        for n in range(1, 40):
            term = term @ A / n
            series = series + term
        assert np.max(np.abs(E - series)) < 1e-10 * max(1.0, np.max(np.abs(series)))
        assert np.linalg.det(E) == pytest.approx(math.exp((L[0][0] + L[1][1]) * t))
    centre = PL.read(((0, 1), (-1, 0)), generator=True)
    assert centre["as"] == "generator" and centre["readings"]["flow"]["region"] == "centre" and centre["readings"]["motion"]["returns"]
    assert centre["readings"]["field"]["equation"] == "the wave equation" and centre["readings"]["vacuum"]["path"] == "the vacuum stands still"
    saddle = PL.read(((0, 1), (1, 0)), generator=True)
    assert saddle["readings"]["flow"]["region"] == "saddle" and saddle["readings"]["field"]["equation"] == "Laplace's equation"
    assert saddle["readings"]["field"]["time"] == "modular time, a thermometer" and saddle["readings"]["form"]["type"] == "indefinite"
    shear = PL.read(((0, 1), (0, 0)), generator=True)
    assert shear["readings"]["flow"]["region"] == "pure shear" and shear["readings"]["number"]["system"] == "dual"
    focus = PL.read(((-0.3, 1), (-1, -0.3)), generator=True)
    assert focus["readings"]["flow"]["region"] == "attracting focus" and "damped" in focus["readings"]["motion"]["words"]
    assert focus["readings"]["vacuum"]["on_the_plane"] is False and "only at α = 0" in focus["readings"]["form"]["preserves"]
    assert K.region(0.5, 0.1) == "repelling node" and K.region(-0.5, 0.25) == "det L = 0: one direction frozen" and K.region(0.2, 0.0) == "repelling star or improper node"
    assert K.kind(((0, 1), (-1, 0))) == "ell" and K.kind(((1, 1), (0, 1))) == "par" and K.kind(PELL) == "hyp"


def test_the_port_equals_the_lab_core():
    """The lab's own JavaScript core (version 8.3) was run on 52 matrices and its output frozen (studies/plane_2026nov); the package
    must return the same readings.  This is the port's own gate: no number in it was computed by the package."""
    lab = json.loads(LAB_CORE.read_text(encoding="utf-8"))
    assert len(lab["core_sha256"]) == 64 and len(lab["matrices"]) == len(lab["kernel"]) == 52
    worst = 0.0
    for M, J in zip(lab["matrices"], lab["kernel"], strict=True):
        Mt = tuple(tuple(row) for row in M)
        r, R, L = K.read(Mt), K.numerical_range(Mt), K.lab_numerical_range(Mt)
        dec, motion = r["decomposition"], r["readings"]["motion"]
        assert (R["shape"], R["real_foci"], r["kind"]) == (J["range"]["kind"], J["range"]["real"], J["kind"])
        assert K.region(float(dec["alpha"]), float(dec["beta"])) == J["region"]
        scale = max(1.0, J["range"]["norm2"])
        pairs = [(R["semi_major"], J["range"]["semiMajor"]), (R["semi_minor"], J["range"]["semiMinor"]), (R["radius"], J["range"]["w"]),
                 (R["norm"], J["range"]["norm2"]), (float(R["centre"]), J["range"]["centre"]), (L["radius"], J["range"]["w"]),
                 (L["semi_major"], J["range"]["semiMajor"]), (L["norm"], J["range"]["norm2"])]
        pairs += [(float(dec[k]), J["dec"][k]) for k in ("alpha", "s1", "s2", "w")]
        worst = max(worst, max(abs(a - b) for a, b in pairs) / scale, abs(float(dec["beta"]) - J["dec"]["beta"]) / scale ** 2)
        E, top = K.expm(Mt, lab["t"]), max(1.0, max(abs(v) for row in J["expm"] for v in row))
        worst = max(worst, max(abs(E[i][j] - J["expm"][i][j]) for i in range(2) for j in range(2)) / top)
        if J["length"] is not None and abs(J["det"] - 1) < 1e-9:
            assert motion["length"] == pytest.approx(J["length"], rel=1e-12)
        if J["psi"] is not None:
            assert motion["psi"] == round(J["psi"]) and abs(J["psi"] - round(J["psi"])) < 1e-9
        if J["order"] is not None:
            assert motion["order"] == (math.inf if J["order"] == "inf" else J["order"])
        if J["jac"] is not None:
            jac = K.jacobi_monodromy(Mt)
            if jac["defined"]:
                assert float(jac["trace_M2"]) == pytest.approx(J["jac"]["trM2"], rel=1e-12)
                assert jac["two_cosh_length"] == pytest.approx(J["jac"]["twoCoshL"], rel=1e-12)
    assert worst < 1e-12, worst
