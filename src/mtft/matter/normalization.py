"""Conditional U(1) dimensional reduction, not a selected MTFT vacuum.

The zero-mode photon is assumed to exist and commute with the internal flux.
Its coupling is a parent-action INPUT. Gram matrices and charges are supplied;
this module does not identify an electron or generate its mass.
"""
from __future__ import annotations

import math
import numbers
import numpy as np


def _positive(value, name):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, numbers.Real):
        raise ValueError(f"{name} must be a finite positive real number")
    value = float(value)
    if not math.isfinite(value) or value <= 0:
        raise ValueError(f"{name} must be a finite positive real number")
    return value


def canonical_u1_vertex(gram, charge, gauge_coupling):
    """Normalize a same-charge mode block and its constant-photon vertex.

    With psi_raw=C psi_phys, C† K C=I and C†(g q K)C=g q I.
    Returns numerical diagnostics; matrix identities are conditional algebra.
    This is a vector-like charged block, not a chiral SM anomaly calculation.
    """
    coupling = _positive(gauge_coupling, "gauge_coupling")
    if isinstance(charge, (bool, np.bool_)) or not isinstance(charge, numbers.Real):
        raise ValueError("charge must be finite and real")
    charge = float(charge)
    if not math.isfinite(charge):
        raise ValueError("charge must be finite and real")
    k = np.asarray(gram, dtype=complex)
    if k.ndim != 2 or k.shape[0] != k.shape[1] or k.shape[0] == 0:
        raise ValueError("gram must be a nonempty square matrix")
    if not np.isfinite(k).all():
        raise ValueError("gram must be finite")
    scale = float(np.linalg.norm(k, ord=2))
    if scale == 0 or np.linalg.norm(k-k.conj().T, ord=2) > 1e-12*scale:
        raise ValueError("gram must be Hermitian")
    k = (k+k.conj().T)/2
    try:
        lower = np.linalg.cholesky(k)
    except np.linalg.LinAlgError as exc:
        raise ValueError("gram must be positive definite") from exc
    # K=L L†, hence C=(L†)^-1, not L^-1.
    c = np.linalg.solve(lower.conj().T, np.eye(len(k), dtype=complex))
    kinetic = c.conj().T @ k @ c
    vertex = c.conj().T @ (charge*coupling*k) @ c
    return {
        "basis_change": c, "kinetic": kinetic, "vertex": vertex,
        "kinetic_residual": float(np.linalg.norm(kinetic-np.eye(len(k)), ord=2)),
        "vertex_residual": float(np.linalg.norm(vertex-charge*coupling*np.eye(len(k)), ord=2)),
        "condition_number": float(np.linalg.cond(k)),
        "status": "DIAGNOSTIC evaluation of an EXACT conditional identity",
        "photon_zero_mode": "ASSUMED constant, unbroken, massless U(1) generator",
        "gauge_coupling": "INPUT", "electron_identification": "OPEN",
    }


def four_dimensional_coupling(parent_coupling, internal_volume):
    """g4=gD/sqrt(Vol) for the unwarped constant U(1) mode.

    Caller must use compatible natural units: Vol has dimension mass^-d and
    gD dimension mass^-d/2. The output is dimensionless. Both values are INPUT.
    """
    return _positive(parent_coupling, "parent_coupling") / math.sqrt(
        _positive(internal_volume, "internal_volume")
    )


def electroweak_matching(g2, gy):
    """Conditional Q=T3+Y matching, with no kinetic mixing or GS photon mass.

    Higgs vacuum and neutral gauge-boson kinetic action are assumed standard.
    U(1)_Y is not by itself U(1)_em. Couplings are supplied, not predicted.
    """
    g2, gy = _positive(g2, "g2"), _positive(gy, "gy")
    norm = math.hypot(g2, gy)
    e = g2 * (gy/norm)
    return {"e": e, "sin2_theta_w": (gy/norm)**2,
            "alpha": e*e/(4*math.pi), "generator": "Q = T3 + Y",
            "status": "EXACT conditional tree-level matching; couplings INPUT"}


def breathing_gauge_kinetics(chi_over_mp, internal_dimension=2):
    """Gauge kinetic prefactor for uniform Einstein-frame breathing.

    ds²_(4+d)=exp(-d*sigma)ds²_4E+exp(2*sigma)ds²_int,
    chi=Mp*sqrt(d(d+2)/2)*sigma. Z=exp(d*sigma) multiplies -F²/4.
    Thus Z=exp(sqrt(2d/(d+2))*chi/Mp); g_eff/g0=Z^-1/2 for
    a CONSTANT background chi. For varying chi retain Z F² in the action;
    a field-dependent canonical rescaling introduces derivative terms.
    """
    if isinstance(internal_dimension, bool) or not isinstance(internal_dimension, numbers.Integral) or internal_dimension < 1:
        raise ValueError("internal_dimension must be a positive integer")
    if isinstance(chi_over_mp, (bool, np.bool_)) or not isinstance(chi_over_mp, numbers.Real):
        raise ValueError("chi_over_mp must be finite and real")
    x = float(chi_over_mp)
    if not math.isfinite(x):
        raise ValueError("chi_over_mp must be finite and real")
    a = math.sqrt(2*internal_dimension/(internal_dimension+2))
    try:
        z, relative_g = math.exp(a*x), math.exp(-a*x/2)
    except OverflowError as exc:
        raise ValueError("breathing prefactor overflows floating point") from exc
    if z == 0 or relative_g == 0 or not math.isfinite(z*relative_g):
        raise ValueError("breathing prefactor outside numerical domain")
    return {"Z": z, "g_eff_over_g0": relative_g,
            "mp_dZ_dchi": a*z, "canonical_exponent": a,
            "internal_dimension": int(internal_dimension),
            "status": "EXACT conditional reduction; numerical evaluation DIAGNOSTIC",
            "implemented_evolution": False}


def electromagnetic_invariants(electric_field, magnetic_field):
    """Flat +--- invariants in natural units for supplied classical E, B.

    F_mu_nu F^mu_nu=2(B²-E²); canonical Maxwell energy=(E²+B²)/2
    and flux=E cross B use Z=1. For an unrescaled field in -Z F²/4,
    multiply this energy and flux by Z; this function does not apply Z.
    A null plane wave has zero F² yet nonzero stress-energy: the direct
    chi F² coupling does not replace gravitational coupling to T_mu_nu.
    """
    e, b = np.asarray(electric_field, dtype=float), np.asarray(magnetic_field, dtype=float)
    if e.shape != (3,) or b.shape != (3,) or not np.isfinite(e).all() or not np.isfinite(b).all():
        raise ValueError("electric and magnetic fields must be finite real 3-vectors")
    return {"F_squared": float(2*(b@b-e@e)),
            "energy_density": float((e@e+b@b)/2),
            "E_dot_B": float(e@b), "poynting": np.cross(e,b).tolist(),
            "status": "classical diagnostic in declared units; canonical Z=1"}
