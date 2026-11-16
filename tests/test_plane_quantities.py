"""mtft.plane.quantities: the lattice of physical quantities, its unit sublattices and the prospector."""
import itertools
import json
from pathlib import Path

import pytest

from mtft.plane import quantities as Q

LAB_CORE = Path(__file__).resolve().parents[1] / "studies" / "plane_2026nov" / "lab_core_outputs.json"


def _det(rows):
    """An integer determinant by cofactors (the second route to the Planck sublattice's index)."""
    if len(rows) == 1:
        return rows[0][0]
    return sum((-1) ** j * rows[0][j] * _det([r[:j] + r[j + 1:] for r in rows[1:]]) for j in range(len(rows)))


def test_the_unit_sublattices():
    lat = Q.lattice()
    assert lat["natural"] == {"smith": [1, 1, 2], "rank": 3, "quotient": "Z + Z/2"}
    assert lat["planck"] == {"smith": [1, 1, 2, 2], "rank": 4, "index": 4, "quotient": "Z/2 + Z/2"}
    assert abs(_det([list(Q.UNITS[k]) for k in ("c", "hbar", "eps0", "G")])) == 4              # the index, by a determinant
    assert all(Q.natural(Q.UNITS[k]) == (0, 0) for k in ("c", "hbar", "eps0")) and Q.natural(Q.UNITS["G"]) == (-2, 0)
    assert Q.natural(Q.MASS) == (1, 0) and Q.natural(Q.CHARGE) == (0, 1)                     # the two generators of the quotient
    box = list(itertools.product(range(-3, 4), repeat=4))
    assert {Q.planck_class(v) for v in box} == {(0, 0), (0, 1), (1, 0), (1, 1)}              # four classes in Planck units
    assert {Q.natural(v)[1] for v in box} == {0, 1} and len({Q.natural(v)[0] for v in box}) > 10


def test_the_kernel_of_natural_units_is_the_span_of_c_hbar_eps0():
    trivial, count = 0, 0
    for v in itertools.product(range(-3, 4), repeat=4):
        d = Q.planck_decompose(v)
        back = tuple(d["Delta"] * m + d["eps"] * q + d["n_c"] * c + d["n_hbar"] * h + d["n_eps0"] * e
                     for m, q, c, h, e in zip(Q.MASS, Q.CHARGE, Q.UNITS["c"], Q.UNITS["hbar"], Q.UNITS["eps0"]))
        assert back == v and (d["Delta"], d["eps"]) == Q.natural(v)
        count += 1
        if Q.natural(v) == (0, 0):                                                           # then v is a combination of the three units alone
            trivial += 1
            assert tuple(d["n_c"] * c + d["n_hbar"] * h + d["n_eps0"] * e for c, h, e in zip(Q.UNITS["c"], Q.UNITS["hbar"], Q.UNITS["eps0"])) == v
    assert count == 7 ** 4 and trivial > 50
    assert Q.planck_decompose((1, 2, -2, 0)) == {"Delta": 1, "eps": 0, "n_c": 2, "n_hbar": 0, "n_eps0": 0}      # E = m c²
    assert Q.planck_decompose(Q.UNITS["G"])["Delta"] == -2                                   # G is a mass to the power −2


def test_the_catalogue():
    assert (len(Q.CATALOGUE), len(Q.EXTENSION), len(Q.LAB)) == (33, 26, 13) and len(Q.everything()) == 72
    assert all(len(q) == 4 and len(q[1]) == 4 for q in Q.everything())
    names = {q[0]: q[1] for q in Q.everything()}
    assert names["energy · work · torque"] == (1, 2, -2, 0) and names["capacitance"] == (-1, -2, 4, 2) and names["speed of light c"] == Q.UNITS["c"]
    assert names["permittivity ε₀"] == Q.UNITS["eps0"] and names["gravitational constant"] == Q.UNITS["G"]
    # consistency of the catalogue with itself: E = F x, P = E/t, V = E/q, C = q/V, L = Φ/I, R = V/I, ε₀ μ₀ c² = 1
    add = lambda *vs: tuple(sum(c) for c in zip(*vs))
    neg = lambda v: tuple(-x for x in v)
    assert add(names["force"], names["position"]) == names["energy · work · torque"]
    assert add(names["energy · work · torque"], neg(names["time"])) == names["power"]
    assert add(names["energy · work · torque"], neg(names["charge"])) == names["potential φ"]
    assert add(names["charge"], neg(names["potential φ"])) == names["capacitance"]
    assert add(names["magnetic flux"], neg(names["current"])) == names["inductance"]
    assert add(names["potential φ"], neg(names["current"])) == names["resistance · impedance"]
    assert add(names["permittivity ε₀"], names["permeability μ₀"], Q.UNITS["c"], Q.UNITS["c"]) == (0, 0, 0, 0)
    assert Q.dim_str((1, 2, -1, 0)) == "M L² T⁻¹" and Q.dim_str((0, 0, 0, 0)) == "1" and Q.dim_str((-1, -3, 4, 2)) == "M⁻¹ L⁻³ T⁴ I²"


def test_the_prospector():
    visc = Q.site((1, -1, -1, 0))
    assert [q[0] for q in visc["here"]] == ["dynamic viscosity"] and visc["occupied"] and visc["Delta"] == 3 and visc["eps"] == 0
    assert "t · Π" in visc["products"] and visc["planck"] == {"Delta": 3, "eps": 0, "n_c": 3, "n_hbar": -2, "n_eps0": 0}
    row = Q.site((-1, -1, 4, 2))                                                             # the lab's preset "your row"
    assert row["dimension"] == "M⁻¹ L⁻¹ T⁴ I²" and not row["occupied"] and row["product_count"] >= 1
    assert "x · C" in row["products"]                                                        # a length times a capacitance: found, and it checks
    assert tuple(a + b for a, b in zip((0, 1, 0, 0), (-1, -2, 4, 2))) == (-1, -1, 4, 2)
    assert (row["Delta"], row["eps"]) == (-2, 0) == Q.natural(Q.UNITS["G"])                  # in natural units it sits in G's class
    assert all(Q.natural(q[1]) == (-2, 0) for q in row["mates"]) and any(q[3] == "G" for q in row["mates"])
    origin = Q.site((0, 0, 0, 0))
    assert [q[3] for q in origin["lab"]] == ["N, d, ψ"] and origin["planck"] == {"Delta": 0, "eps": 0, "n_c": 0, "n_hbar": 0, "n_eps0": 0}
    assert len(Q.names_at((0, 0, 1, 0))) == 4 and len(Q.names_at((0, 0, 1, 0), lab=False)) == 1
    with pytest.raises(ValueError):
        Q.site((1, 2, 3))


def test_a_change_of_units_moves_no_ratio():
    factors = {"T": 1e3, "L": 0.01, "M": 7.0}
    for v in ((1, 2, -2, 0), (0, 0, -1, 0), (-1, -3, 4, 2)):
        f = Q.rescale(v, factors)
        assert f == pytest.approx(7.0 ** -v[0] * 0.01 ** -v[1] * 1e3 ** -v[2])
        a, b = 3.2, 0.4                                                                      # two quantities at one site
        assert (a * f) / (b * f) == pytest.approx(a / b)
    assert Q.rescale((0, 0, 0, 0), factors) == 1.0 and Q.rescale((0, 0, -1, 0), {"T": 10.0}) == pytest.approx(10.0)


def test_the_prospector_equals_the_lab_core():
    """The lab's own JavaScript core was run on 49 sites and its output frozen (studies/plane_2026nov): what lives at each site, the
    pairs that build it, its class and its Planck coefficients must be the lab's, name for name."""
    lab = json.loads(LAB_CORE.read_text(encoding="utf-8"))
    assert lab["planck_index"] == Q.lattice()["planck"]["index"] == 4 and len(lab["sites"]) == len(lab["quantities"]) == 49
    occupied = 0
    for v, J in zip(lab["sites"], lab["quantities"]):                               # equal lengths, asserted above (no strict=: Python 3.9)
        s = Q.site(v)
        occupied += bool(J["here"])
        assert [q[0] for q in s["here"]] == J["here"] and s["products"] == J["products"] and len(s["mates"]) == J["mates"]
        assert (s["Delta"], s["eps"], s["dimension"]) == (J["Delta"], J["eps"], J["dim"])
        assert J["planck"]["ok"] and (s["planck"]["n_c"], s["planck"]["n_hbar"], s["planck"]["n_eps0"]) == (J["planck"]["nc"], J["planck"]["nh"], J["planck"]["ne"])
    assert 0 < occupied < 49                                                      # both occupied and empty sites were compared
