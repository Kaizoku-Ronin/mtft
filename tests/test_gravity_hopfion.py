"""mtft.gravity.hopfion: the squared hopfion's lines of force after t = 0 as an integrable system.

Expected numbers are the lab's recorded ones (MP22, 14 November 2026: studies/gravity_2026nov/lab/mp22) or are computed in the
test by another route, and the test says which.  The lab integrated its lines with scipy; here they are integrated by the
package's own routine, so the agreement is between two integrators as well."""
import json
import math
from pathlib import Path

import numpy as np
import pytest

from mtft.gravity import fields as FD
from mtft.gravity import hopfion as HP

MP22 = json.loads((Path(__file__).resolve().parents[1] / "studies" / "gravity_2026nov" / "lab" / "mp22" / "mp22_numbers.json").read_text(encoding="utf-8"))
POINT = np.array(HP.MP10_POINT)


def _sphere(n, seed):
    v = np.random.default_rng(seed).normal(size=(n, 4))
    v /= np.linalg.norm(v, axis=1)[:, None]
    return v[:, 0] + 1j * v[:, 1], v[:, 2] + 1j * v[:, 3]


def test_the_derivation_is_exact():
    checks = HP.exact_checks()
    assert len(checks) == 6 and all(checks.values()), checks


def test_the_model():
    for t in (0.1, 0.4, 0.8, 1.5):
        md = HP.Hopfion(t)
        assert md.m == pytest.approx(math.sqrt(1 + t * t)) and md.eps == pytest.approx(t / (2 * math.sqrt(1 + t * t))) and md.eps < 0.5
        assert md.sigma ** 2 == pytest.approx(1 - 1j * t) and md.focus == pytest.approx(2 * math.sqrt(md.m))
        assert md.leading_law() == pytest.approx(5 / 12 * md.eps ** 2)
        al, be = _sphere(20000, int(1000 * t))
        for kind in "BE":
            z, p = md.reduce(al, be, kind)
            assert np.max(np.abs(md.H(z, p) - 1)) < 1e-12                           # H is |α|² + |β|²
            a2, b2 = md.lift(z, p, kind)
            assert np.max(np.abs(a2 - al)) < 1e-12 and np.max(np.abs(b2 - be)) < 1e-12
            place = md.place(md.I2(z, p))                                            # MP22's P2: strictly between the two closed lines
            assert 0 < place.min() and place.max() < 1
            assert np.max(np.abs(md.place(md.I2(*md.reduce(*md.axial_line(np.linspace(0, 6, 50), kind), kind))))) < 1e-12
            assert md.place(md.I2(*md.reduce(*md.cross_point(kind), kind))) == pytest.approx(1.0, abs=1e-12)
        assert md.K_cross > md.K_axial == pytest.approx(1 / t ** 2) and md.torus(0) == md.K_axial and md.torus(1) == pytest.approx(md.K_cross)
    with pytest.raises(ValueError):
        HP.Hopfion(0.0)


def test_the_reduction_pointwise_in_space():
    """MP22's P0 with exact gradients: a step along the field moves (s, P) along the Hamiltonian flow (the lab, by differences: 2e-8)."""
    rng = np.random.default_rng(22)
    worst = 0.0
    for t in (0.1, 0.4, 0.8, 1.5):
        md = HP.Hopfion(t)
        X = rng.normal(0, 1.5, (60, 3))
        al, be = FD.labels(t, X)
        ga, gb = FD.grad_labels(t, X)
        for kind in "BE":
            V = FD.G(t, X)
            u = FD.unit(V.imag if kind == "B" else V.real)
            dal, dbe = np.sum(ga * u, -1), np.sum(gb * u, -1)
            s = np.sqrt(1 - 1j * t * (1 - al))
            P = np.conj(be) if kind == "B" else 1j * np.conj(be)
            ds, dP = (1j * t / (2 * s)) * dal, (np.conj(dbe) if kind == "B" else 1j * np.conj(dbe))
            gV = 4 * (s * s - md.sigma ** 2) * np.conj(s) / t ** 2
            v = np.stack([ds.real, ds.imag, dP.real, dP.imag], -1)
            w = np.stack([(2 * P).real, (2 * P).imag, (-gV).real, (-gV).imag], -1)
            across = np.linalg.norm(v - (np.sum(v * w, -1) / np.sum(w * w, -1))[:, None] * w, axis=1) / np.linalg.norm(v, axis=1)
            worst = max(worst, float(across.max()))
    assert worst < 1e-10, worst
    assert worst < MP22["P0"]["worst angle"]                                         # sharper than the lab's differences


def test_the_second_invariant_along_lines_of_force():
    """MP22's P1 and P6, in space on the field itself, by the package's own integrator."""
    assert np.max(np.abs(HP.line_field(0.4, "B")(POINT) - FD.line_field(0.4, "B", squared_field=True)(POINT))) < 1e-14   # two codes of one field
    for t, kind in ((0.4, "B"), (0.4, "E"), (0.8, "E")):
        run = HP.invariant_drift(t, POINT, 20, kind)
        assert run["drift"] < 1e-8 and run["H_off"] < 1e-12 and run["error"] < 1e-9, run
    run = HP.invariant_drift(0.8, POINT, 100, "B")                                   # P6: the line MP11 had lost
    rec = MP22["P6"]
    assert run["drift"] < 1e-8 and run["drift"] < rec["drift of I_2 / range"]
    assert run["place"] == pytest.approx(rec["place in the range"], abs=1e-10)
    assert 4.6 < run["farthest"] < 4.75 and rec["farthest from the origin"] == pytest.approx(4.697, abs=1e-3)   # MP22 forecast more than 10: a recorded miss
    md = HP.Hopfion(0.8)
    z, p = md.reduce(*FD.labels(0.8, POINT), "B")
    assert md.winding(K=float(-md.I2(z, p) / (4 * md.m))) == pytest.approx(rec["W"], abs=1e-10)


def test_the_table_of_the_winding():
    """MP22's P4 and its table of W − 1, to the plan's 1e-5 and to 1e-9 against the lab's number file."""
    table, lab = HP.table(), MP22["P4"]["table"]
    assert tuple(table) == HP.TABLE_TIMES
    for t, row in table.items():
        rec = lab[f"t={t:.1f}"]
        assert row["eps"] == pytest.approx(rec["eps"], abs=1e-12) and row["leading"] == pytest.approx(rec["leading law (5/12) eps^2"], abs=1e-12)
        assert row["axial"] == pytest.approx(rec["axial W - 1"], abs=1e-9)
        for fr in HP.TABLE_FRACTIONS:
            assert row[fr] == pytest.approx(rec[f"W - 1 at {fr:.3f} of the range"], abs=1e-9)
        assert row["leading"] < row[0.999] < row[0.5] < row["axial"]                 # the leading law is the floor, reached beside the cross line
    # the record's printed table at t = 0.4 and t = 3.0: axial line, middle torus, beside the cross line
    assert [round(table[0.4][k], 5) for k in ("axial", 0.5, 0.999)] == [0.01577, 0.01506, 0.01442]
    assert [round(table[3.0][k], 5) for k in ("axial", 0.5, 0.999)] == [0.29630, 0.13883, 0.09659]
    md = HP.Hopfion(0.4)
    assert md.winding() == pytest.approx(MP22["P4"]["axial W at t = 0.4"], abs=1e-12) and round(md.winding(), 6) == 1.015770
    small = HP.Hopfion(0.05)                                                         # the leading law at t = 0.05: within 2 % of 5/12 on five tori
    ratios = [(small.winding(fraction=f, dps=40) - 1) / small.eps ** 2 for f in (0.05, 0.25, 0.5, 0.75, 0.95)]
    assert ratios == pytest.approx(MP22["P4"]["(W - 1)/eps^2 at t = 0.05, at 0.05, 0.25, 0.5, 0.75, 0.95 of the range of K"], abs=1e-9)
    assert all(abs(r / (5 / 12) - 1) < 0.02 for r in ratios)
    # MP22's P7: at t = 0.4 the torus with W = 65/64 is a torus of closed lines, and so is the one with 66/65
    assert md.winding(K=MP22["P7"]["K with W = 65/64"]) == pytest.approx(65 / 64, abs=1e-12)
    assert md.winding(K=MP22["P7"]["K with W = 66/65"]) == pytest.approx(66 / 65, abs=1e-12)
    assert 1 / (table[0.4]["axial"]) == pytest.approx(63.4, abs=0.05) and 1 / table[0.4][0.999] == pytest.approx(69.4, abs=0.1)   # MP22: from 63.4 to 69.4 circuits


def test_the_winding_measured_along_an_orbit():
    """The second route: the two periods measured along an orbit of the plane system, against the two integrals."""
    md = HP.Hopfion(0.4)
    z, p = md.reduce(*FD.labels(0.4, POINT), "B")
    got, rec = md.measured_periods(complex(z), complex(p)), MP22["P5"]["b"]
    assert got["K"] == pytest.approx(rec["K"], abs=1e-10) and got["place"] == pytest.approx(rec["place in the range"], abs=1e-10)
    assert got["W"] == pytest.approx(rec["W"], abs=1e-10) and got["spread_mu"] < 1e-10 and got["spread_nu"] < 1e-10
    Tm, Tn, roots = md.periods(got["K"])
    assert got["T_mu"] == pytest.approx(float(Tm), abs=1e-10) and got["T_nu"] == pytest.approx(float(Tn), abs=1e-10)
    assert 0 < roots[0] < roots[1] < 1 < roots[2]
    assert round(1 / (got["W"] - 1), 1) == 63.9                                      # MP22: one turn in 63.9 circuits, not the older records' 29
    run = md.orbit(complex(z), complex(p), 3.0)                                      # H and I₂ are kept along the orbit
    zz, pp = run["y"][:, 0] + 1j * run["y"][:, 1], run["y"][:, 2] + 1j * run["y"][:, 3]
    assert np.ptp(md.H(zz, pp)) < 1e-12 and np.ptp(md.I2(zz, pp)) < 1e-12 * md.I2_range
    wide = HP.Hopfion(0.8)
    z, p = wide.reduce(*FD.labels(0.8, POINT), "B")
    assert wide.measured_periods(complex(z), complex(p))["W"] == pytest.approx(MP22["P6"]["W"], abs=1e-10)


def test_the_two_closed_lines_in_space():
    """MP22's P2: the axial great circle and the cross line are closed lines of force, linked once."""
    for t, kind, link in ((0.4, "B", -1), (0.4, "E", 1), (0.8, "B", -1)):
        md = HP.Hopfion(t)
        ph = np.linspace(0, 2 * math.pi, 400, endpoint=False) + 0.1
        al, be = md.axial_line(ph, kind)
        dal, dbe = md.axial_tangent(ph, kind)
        X = FD.chart(t, al, be)
        V = FD.G(t, X)
        u = FD.unit(V.imag if kind == "B" else V.real)
        ga, gb = FD.grad_labels(t, X)
        va, vb = np.sum(ga * u, -1), np.sum(gb * u, -1)                              # the field's velocity on the label sphere
        v = np.stack([va.real, va.imag, vb.real, vb.imag], -1)
        w = np.stack([dal.real, dal.imag, dbe.real, dbe.imag], -1)
        across = np.linalg.norm(v - (np.sum(v * w, -1) / np.sum(w * w, -1))[:, None] * w, axis=1) / np.linalg.norm(v, axis=1)
        assert across.max() < 1e-9                                                   # the circle is tangent to the field
        axial = HP.closed_line(t, X[3], kind)
        assert axial["closure"] < 1e-8 and axial["length"] == pytest.approx(2 * math.pi, abs=1e-9)   # a great circle of the label sphere
        cross = HP.closed_line(t, FD.chart(t, *md.cross_point(kind)), kind)
        assert cross["closure"] < 1e-8 and cross["length"] < 2 * math.pi
        zc, pc = md.reduce(*FD.labels(t, cross["points"]), kind)
        assert np.max(np.abs(md.elliptic(zc)[1] ** 2 - md.c0)) < 1e-10               # it stays on its hyperbola cos²ν = c₀
        assert np.max(np.abs(md.place(md.I2(zc, pc)) - 1)) < 1e-10
        circle = FD.chart(t, *md.axial_line(np.linspace(0, 2 * math.pi, 900, endpoint=False), kind))
        assert FD.gauss_linking(circle, cross["points"]) == pytest.approx(link, abs=0.02)   # the lab's signs: −1 magnetic, +1 electric
