"""
tests/test_hodge_tower.py
=========================
The Hodge endomorphism tower of the X0(143) particle box (v0.11.3).

Certifies, on the 22-dimensional NEW part of the cuspidal homology, the
three-story tower that places the generation decomposition inside Hodge
theory (session log 2026-08-05; Cert class):

    dim { X : [X, T2] = 0 }              = 44   (all Hecke-commuting maps)
    dim { X : [X, T2] = 0, [X, J] = 0 }  = 22   (the Hodge cut: one C
                                                 per embedding)
    dim_Q  Q[T2]                          = 11   (= 1 + 4 + 6 = the RM
                                                 fields K1 x K2 x K3
                                                 = End(J0(143)_new) x Q)

together with [T2, J] = 0 (Hecke operators preserve the Hodge
decomposition, as classes of algebraic correspondences must) and the
totally-real check on the Hecke fields (Albert type I: no Weil-class
pathology is possible for these isotypic factors).

Fast tier: the two J-free stories (44 and 11) plus totally-real roots.
Slow tier (-m slow): builds the complex structure J from the engine's own
period functionals (v03 cache or fresh extraction) and certifies the
middle story and the commutation residuals.
"""

import importlib.util
import os
import sys

import numpy as np
import pytest

_HERE = os.path.dirname(os.path.abspath(__file__))
_STUDIES = os.path.join(_HERE, "..", "studies")
if _STUDIES not in sys.path:
    sys.path.insert(0, _STUDIES)


def _load(name):
    spec = importlib.util.spec_from_file_location(
        name, os.path.join(_STUDIES, name + ".py"))
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def new_part():
    """T2 restricted to the 22-dim new part, plus the embedding data."""
    v02 = _load("x0143_particle_box_v02")
    p1, tris, edges, ms = v02.build_engine()
    P = v02.float_projection(ms)
    Bc, proj_c, restrict, T2, E, lines = v02.eigendata(ms, P)
    lines = v02.assign_orbits(lines)
    cols = [np.asarray(L[3], float) for L in lines]      # 22 eigenvectors
    Bn = np.array(cols).T                                # 26 x 22
    assert np.linalg.matrix_rank(Bn, tol=1e-8) == 22
    T2n = np.linalg.pinv(Bn) @ T2 @ Bn
    return v02, ms, P, restrict, lines, Bn, T2n


def _nullity(M):
    s = np.linalg.svd(M, compute_uv=False)
    return int(np.sum(s < 1e-7 * s.max()))


def test_hecke_commutant_dim_44(new_part):
    """First story: sum of squared multiplicities = 4 * (1+4+6) = 44."""
    *_, T2n = new_part
    I22 = np.eye(22)
    K_T = np.kron(I22, T2n) - np.kron(T2n.T, I22)
    assert _nullity(K_T) == 44


def test_rational_hecke_algebra_dim_11(new_part):
    """Third story: dim_Q Q[T2] = 11 = K1 + K2 + K3 = End(J0_new) x Q."""
    *_, T2n = new_part
    pows, Mk = [np.eye(22).ravel()], np.eye(22)
    for _ in range(21):
        Mk = Mk @ T2n
        pows.append(Mk.ravel())
    assert np.linalg.matrix_rank(np.array(pows), tol=1e-7) == 11


def test_hecke_fields_totally_real():
    """Albert type I: every root of q2 and q3 is real, so the isotypic
    abelian varieties carry real multiplication and no Weil classes."""
    r2 = np.roots([1, -3, -1, 5, 1])
    r3 = np.roots([1, 0, -10, 2, 24, -7, -12])
    assert np.all(np.abs(r2.imag) < 1e-9)
    assert np.all(np.abs(r3.imag) < 1e-9)


@pytest.mark.slow
def test_hodge_cut_dim_22_and_commutation(new_part):
    """Middle story: adding the complex structure J (built from the
    engine's own period functionals) cuts the commutant to exactly 22,
    with J^2 = -I and [T2, J] = 0 to machine precision."""
    from collections import defaultdict
    v02, ms, P, restrict, lines, Bn, T2n = new_part
    v03 = _load("x0143_particle_box_v03")
    _, _, _, _, _, lines3, Pis, meta = v03.get_period_data()
    emb = defaultdict(list)
    for L, Pi in zip(lines3, Pis):
        emb[(L[1], round(L[2], 8))].append((L, np.asarray(Pi)))
    assert len(emb) == 11
    cols, J_blocks = [], []
    for (orbit, lam), pair in sorted(emb.items(), key=lambda kv: kv[0][1]):
        (L1, Pi1), (L2, _) = pair
        u1 = np.asarray(L1[3], float)
        u2 = np.asarray(L2[3], float)
        z1, z2 = complex(Pi1 @ u1), complex(Pi1 @ u2)
        M = np.array([[z1.real, z2.real], [z1.imag, z2.imag]])
        a1 = np.linalg.solve(M, [(1j * z1).real, (1j * z1).imag])
        a2 = np.linalg.solve(M, [(1j * z2).real, (1j * z2).imag])
        cols += [u1, u2]
        J_blocks.append(np.array([[a1[0], a2[0]], [a1[1], a2[1]]]))
    Bn2 = np.array(cols).T
    T2n2 = np.linalg.pinv(Bn2) @ (Bn @ T2n @ np.linalg.pinv(Bn)) @ Bn2 \
        if False else np.linalg.pinv(Bn2) @ (restrict(
            v02.hecke_float(ms, P, 2))) @ Bn2
    Jn = np.zeros((22, 22))
    for k, Jb in enumerate(J_blocks):
        Jn[2 * k:2 * k + 2, 2 * k:2 * k + 2] = Jb
    assert np.linalg.norm(Jn @ Jn + np.eye(22)) < 1e-10
    assert (np.linalg.norm(T2n2 @ Jn - Jn @ T2n2)
            / np.linalg.norm(T2n2)) < 1e-10
    I22 = np.eye(22)
    K_T = np.kron(I22, T2n2) - np.kron(T2n2.T, I22)
    K_J = np.kron(I22, Jn) - np.kron(Jn.T, I22)
    assert _nullity(np.vstack([K_T, K_J])) == 22
