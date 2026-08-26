"""Core period and Hodge geometry for X_0(143).

The frozen v6 artifact supplies the 13x26 cuspidal period matrix, the
intersection form, and the exact symplectic change of basis.  The shipped
``Omega_symplectic_13x26`` field is retained for provenance only: its source
JSON has a known scientific-notation serialization defect.  This module
*always* reconstructs

    Omega_sym = Omega_cusp @ S

at call time before forming tau = A^{-1}B.

Epistemic classes
-----------------
* exact integer forms/change-of-basis: EXACT
* period/tau identities: CERTIFIED numerical replay of the v6 artifact
* Hodge metric/complex structure: exact formulas evaluated from certified tau
"""
from __future__ import annotations

from functools import lru_cache
import json
from pathlib import Path
import re

import mpmath as mp

GENUS = 13
LEVEL = 143
_DATA = Path(__file__).resolve().parent / "_data"


def data_path(name: str) -> Path:
    p = _DATA / name
    if not p.exists():
        raise FileNotFoundError(f"period data file missing: {name}")
    return p


@lru_cache(maxsize=1)
def period_record() -> dict:
    """Raw frozen v6 period record (deep-copy if mutating)."""
    return json.loads(data_path("X0_143_period_data_v6.json").read_text())


def _real_token(x) -> mp.mpf:
    """Parse a frozen real token, including the legacy bare ``E-57`` zero.

    The old GP->JSON exporter occasionally tokenized ``0.E-57 7.9E-59`` as
    ``['0.E-57', '7.9']`` plus a dropped exponent in the *derived*
    Omega_symplectic field.  Bare exponent tokens therefore mean numerical
    zero.  We never use that defective derived field for reconstruction.
    """
    if isinstance(x, (int, float)):
        return mp.mpf(x)
    s = str(x).strip()
    if re.fullmatch(r"[+-]?E[+-]?\d+", s, re.I):
        return mp.mpf("0")
    if re.fullmatch(r"[+-]?0(?:\.0*)?E[+-]?\d+", s, re.I):
        return mp.mpf("0")
    return mp.mpf(s)


def _complex_pair(z) -> mp.mpc:
    return mp.mpc(_real_token(z[0]), _real_token(z[1]))


def _complex_matrix(rows) -> mp.matrix:
    return mp.matrix([[_complex_pair(z) for z in row] for row in rows])


def _integer_matrix(rows) -> mp.matrix:
    return mp.matrix([[int(x) for x in row] for row in rows])


def intersection_inverse() -> tuple[tuple[int, ...], ...]:
    """Q = E^{-1} in the v6 cuspidal K basis (26x26), EXACT."""
    return tuple(tuple(int(x) for x in r)
                 for r in period_record()["Q_intersection_inverse"])


def intersection_form() -> tuple[tuple[int, ...], ...]:
    """E in the v6 cuspidal K basis (26x26), EXACT, unimodular."""
    return tuple(tuple(int(x) for x in r)
                 for r in period_record()["E_intersection"])


def symplectic_change() -> tuple[tuple[int, ...], ...]:
    """S with S^T E S = J, EXACT unimodular."""
    return tuple(tuple(int(x) for x in r)
                 for r in period_record()["S_symplectic"])


def symplectic_form(g: int = GENUS) -> mp.matrix:
    """J = [[0,I],[-I,0]] in the package's pinned symplectic orientation."""
    J = mp.zeros(2*g)
    for i in range(g):
        J[i, g+i] = 1
        J[g+i, i] = -1
    return J


def omega_cusp(dps: int = 60) -> mp.matrix:
    """13x26 periods in the v6 cuspidal K basis."""
    with mp.workdps(dps):
        return +_complex_matrix(period_record()["Omega_cusp_13x26"])


def legacy_omega_symplectic(dps: int = 60) -> mp.matrix:
    """Forensic only: the defective derived JSON field.  Do not use in math."""
    with mp.workdps(dps):
        return +_complex_matrix(period_record()["Omega_symplectic_13x26"])


def omega_symplectic(dps: int = 60) -> mp.matrix:
    """13x26 symplectic period matrix, reconstructed as Omega_cusp @ S."""
    with mp.workdps(dps):
        return +(omega_cusp(dps) * _integer_matrix(symplectic_change()))


def frozen_riemann_matrix(dps: int = 60) -> mp.matrix:
    """The v6 frozen tau record, for E2/replay comparisons only."""
    with mp.workdps(dps):
        return +_complex_matrix(period_record()["tau_13x13"])


def riemann_matrix(dps: int = 60) -> mp.matrix:
    """tau=A^{-1}B in H_13, recomputed from Omega_cusp and exact S."""
    with mp.workdps(dps):
        O = omega_symplectic(dps)
        A, B = O[:, :GENUS], O[:, GENUS:]
        return +(A ** -1 * B)


def normalized_periods(dps: int = 60) -> mp.matrix:
    """[I | tau], periods of the alpha-normalized holomorphic basis."""
    with mp.workdps(dps):
        out = mp.zeros(GENUS, 2*GENUS)
        t = riemann_matrix(dps)
        for i in range(GENUS):
            out[i, i] = 1
            for j in range(GENUS):
                out[i, GENUS+j] = t[i, j]
        return +out


def hodge_complex_structure(dps: int = 60) -> mp.matrix:
    r"""Complex structure on real homology in the symplectic cycle frame.

    It is pinned by

        [I tau] J_H = i [I tau],   J_H^2 = -I.

    For tau=X+iY,

      J_H = [[-X Y^-1, -Y-X Y^-1 X],
             [ Y^-1,     Y^-1 X    ]].
    """
    with mp.workdps(dps):
        t = riemann_matrix(dps)
        X = mp.matrix([[mp.re(t[i, j]) for j in range(GENUS)]
                       for i in range(GENUS)])
        Y = mp.matrix([[mp.im(t[i, j]) for j in range(GENUS)]
                       for i in range(GENUS)])
        Yi = Y ** -1
        a, b, c, d = -X * Yi, -Y - X * Yi * X, Yi, Yi * X
        out = mp.zeros(2*GENUS)
        for i in range(GENUS):
            for j in range(GENUS):
                out[i, j] = a[i, j]
                out[i, GENUS+j] = b[i, j]
                out[GENUS+i, j] = c[i, j]
                out[GENUS+i, GENUS+j] = d[i, j]
        return +out


def hodge_metric(dps: int = 60) -> mp.matrix:
    r"""Positive symplectic metric G=E J_H in the symplectic frame.

    In block form it is

      [[Y^-1, Y^-1 X],
       [X Y^-1, Y + X Y^-1 X]],

    so for gamma=(n,m), gamma^T G gamma =
    (n+Xm)^T Y^-1(n+Xm)+m^T Y m.
    """
    with mp.workdps(dps):
        return +(symplectic_form() * hodge_complex_structure(dps))


def charge_energy(n, m, dps: int = 60) -> mp.mpf:
    """Dimensionless positive quadratic energy in the symplectic frame.

    This is a mathematical Hodge/polarization invariant for a chosen integral
    charge vector.  Interpreting it as a physical mass/energy requires an
    external model and scale.
    """
    if len(n) != GENUS or len(m) != GENUS:
        raise ValueError("n and m must each have length 13")
    with mp.workdps(dps):
        v = mp.matrix([*[mp.mpf(x) for x in n], *[mp.mpf(x) for x in m]])
        G = hodge_metric(dps)
        return +(v.T * G * v)[0]


__all__ = [
    "GENUS", "LEVEL", "data_path", "period_record", "intersection_inverse",
    "intersection_form", "symplectic_change", "symplectic_form", "omega_cusp",
    "omega_symplectic", "frozen_riemann_matrix", "riemann_matrix",
    "normalized_periods", "hodge_complex_structure", "hodge_metric",
    "charge_energy",
]
