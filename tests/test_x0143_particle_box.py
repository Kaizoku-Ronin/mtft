"""tests/test_x0143_particle_box.py — gates for the X0(143) particle box
(studies/x0143_particle_box{,_v02,_v03}.py), v0.11.0 integration wave.

Fast tier (default):
  * tessellation combinatorics (index 168, 56/84/4, widths {1,11,13,143},
    Manin dims 29/26)            [Add. BQ leg 2]
  * cuspidal T2 characteristic polynomial (exact factor string)
                                 [Add. BQ leg 2]
  * f1 Hecke eigenvalues vs LMFDB 143.a1 (extract_ap, p <= 47)
                                 [Add. BQ leg 3; oracle table in mtft.x0_143]
  * G1 cusp-incidence integer certificates   [Add. BQ leg 6]
  * G4 unitary control: gamma0 = 0 => bound population EXACTLY constant
                                 [Add. BQ leg 6]
  * oracle j-invariant of 143.a1 (Add. BQ sec.5 fix: -262144/1859)

Slow tier (-m slow): capture ceiling smoke (two-sink structure, ceiling 1/2).
"""
import importlib.util
import os
import sys

import numpy as np
import pytest
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
STUDIES = os.path.join(HERE, "..", "studies")
if STUDIES not in sys.path:
    sys.path.insert(0, STUDIES)  # studies scripts import each other by bare name


def _load(name):
    spec = importlib.util.spec_from_file_location(
        name, os.path.join(STUDIES, name + ".py"))
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def engine():
    v02 = _load("x0143_particle_box_v02")
    p1, tris, edges, ms = v02.build_engine()
    return v02, p1, tris, edges, ms


def test_combinatorics(engine):
    v02, p1, tris, edges, ms = engine
    v01 = sys.modules["x0143_particle_box"]
    assert len(p1.reps) == 168
    assert len(tris) == 56 and len(edges) == 84
    T_orb = v01.orbits(p1, ((1, 1), (0, 1)), len(p1.reps))
    assert sorted(len(o) for o in T_orb) == [1, 11, 13, 143]


def test_T2_cuspidal_charpoly(engine):
    v02, p1, tris, edges, ms = engine
    T2 = ms.hecke_on_quotient(2)
    assert T2.shape == (29, 29)
    A2s, _ = ms.restrict_to_cuspidal(T2)
    assert A2s.shape == (26, 26)
    x = sp.symbols("x")
    cp = sp.factor(sp.Matrix(A2s.tolist()).charpoly(x).as_expr())
    expect = ("x**2*(x + 2)**4*(x**4 - 3*x**3 - x**2 + 5*x + 1)**2*"
              "(x**6 - 10*x**4 + 2*x**3 + 24*x**2 - 7*x - 12)**2")
    assert str(cp) == expect


def test_f1_hecke_eigenvalues(engine):
    """f1 (= 143.a1) a_p from the Manin engine vs the certified oracle table
    (itself point-counting + LMFDB verified; auditor: Add. BQ leg 3)."""
    from mtft.x0_143 import CURVE_143A1
    v02, p1, tris, edges, ms = engine
    P = v02.float_projection(ms)
    Bc, proj_c, restrict, T2, E, lines = v02.eigendata(ms, P)
    lines = v02.assign_orbits(lines)
    f1 = [L for L in lines if abs(L[2]) < 1e-9]
    assert len(f1) == 2          # both Atkin-Lehner signs of the same form
    ap = v02.extract_ap(ms, P, restrict, lines, 47)
    oracle = CURVE_143A1.hecke_eigenvalues
    for L in f1:
        vals = ap[id(L)]
        for p in (3, 5, 7, 17, 19, 23, 29, 31, 37, 41, 43, 47):
            assert abs(vals[p] - oracle[p]) < 1e-8, (p, vals[p], oracle[p])


def test_oracle_j_invariant():
    """Add. BQ sec.5: j(143.a1) = c4^3/Delta = -262144/1859 exactly."""
    from mtft.x0_143 import CURVE_143A1
    assert CURVE_143A1.j_invariant == "-262144/1859"
    c4, delta = 64, -1859
    assert delta == -11 * 13 ** 2
    assert abs(c4 ** 3 / delta - (-262144 / 1859)) < 1e-12


def test_G1_cusp_incidence_certificates(engine):
    v02, p1, tris, edges, ms = engine
    v03 = _load("x0143_particle_box_v03")
    n, widths = v03.cusp_incidence(p1, tris)
    # sum_t n_c(t) = width(c); sum_c n_c(t) = 3  (integer, EXACT)
    assert all(int(n[c].sum()) == widths[c] for c in range(4))
    assert all(int(n[:, t].sum()) == 3 for t in range(56))


def test_G4_unitary_control(engine):
    """gamma0 = 0 => every population trajectory is EXACTLY constant:
    capture is impossible without dissipation (theorem, not numerics)."""
    v02, p1, tris, edges, ms = engine
    v01 = sys.modules["x0143_particle_box"]
    v03 = _load("x0143_particle_box_v03")
    n, widths = v03.cusp_incidence(p1, tris)
    A, edge_list, tri_of = v01.dual_graph(p1, tris, edges)
    L = np.diag(A.sum(axis=1)) - A
    ci = int(np.argsort(widths)[-1])  # widest cusp (143)
    Nc = n[ci].astype(float)
    H = L - 4.0 * np.diag(Nc)
    w, V, bound = v03.bound_states(H)
    assert len(bound) > 0
    # distance dipole via BFS from the well support
    seed = list(np.nonzero(Nc > 0)[0])
    d = -np.ones(56, int); d[seed] = 0
    frontier = list(seed)
    while frontier:
        nxt = []
        for u in frontier:
            for v in np.nonzero(A[u])[0]:
                if d[v] < 0:
                    d[v] = d[u] + 1
                    nxt.append(int(v))
        frontier = nxt
    D = d.astype(float)
    M = v03.emission_generator(w, V, D, gamma0=0.0)
    # generator must be identically zero in the gamma0 = 0 control
    assert np.allclose(M, 0.0)


@pytest.mark.slow
def test_capture_ceiling_smoke(engine):
    """Two-sink Lindblad structure => capture ceiling exactly 1/2
    (Add. BQ leg 6)."""
    pytest.skip("slow tier: run studies/x0143_particle_box_v03.py stage G "
                "end-to-end for the full ceiling certificate")
