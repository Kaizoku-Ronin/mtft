"""Matter-03 field-coordinate normalization with supplied physical inputs.

This module uses the lower-index superpotential convention
W = kappa * sum(Y[a,b,c] * left[a] * right[b] * higgs[c]).
The raw fields equal C times the canonical fields, with C^dagger G C = I.
All three holomorphic tensor slots therefore transform with C, without
complex conjugation. This convention is distinct from the historical
section-basis/target-index API in ``surface.hym.normalise_yukawa``.

The metrics, complete tensor, selected light spaces, Higgs vacuum direction,
vev, and overall parent coefficient are INPUT. The functions perform finite
numerical linear algebra; they do not derive a vacuum or predict MTFT masses.
"""
from __future__ import annotations

import math
import numbers

import numpy as np

_HERMITIAN_RTOL = 1e-12
_UNIT_ATOL = 1e-10


def _complex_array(value, name):
    """Reject nonnumerical and Boolean entries before NumPy coerces them."""
    try:
        entries = np.asarray(value, dtype=object)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be a finite numerical array") from exc
    if any(isinstance(x, (bool, np.bool_)) or not isinstance(x, numbers.Complex)
           for x in entries.flat):
        raise ValueError(f"{name} must contain numerical entries, excluding booleans")
    try:
        array = np.asarray(value, dtype=complex)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(f"{name} must be a finite numerical array") from exc
    if not np.isfinite(array).all():
        raise ValueError(f"{name} must be finite")
    return array


def _kappa(value):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, numbers.Complex):
        raise ValueError("kappa must be an explicit finite numerical scalar")
    try:
        result = complex(value)
    except (ValueError, OverflowError) as exc:
        raise ValueError("kappa must be finite") from exc
    if not math.isfinite(result.real) or not math.isfinite(result.imag):
        raise ValueError("kappa must be finite")
    return result


def _vev(value):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, numbers.Real):
        raise ValueError("vev must be an explicit finite nonnegative real scalar")
    try:
        result = float(value)
    except (ValueError, OverflowError) as exc:
        raise ValueError("vev must be finite and nonnegative") from exc
    if not math.isfinite(result) or result < 0:
        raise ValueError("vev must be finite and nonnegative")
    return result


def _finite_result(result, name):
    if not np.isfinite(result).all():
        raise ValueError(f"{name} exceeds the finite numerical domain")
    return result


def canonical_frame(gram):
    """Return C with C^dagger G C = I for a supplied positive Hermitian G.

    Hermiticity is checked to relative entry tolerance 1e-12 before harmless
    roundoff symmetrization. No null directions are dropped or regularized.
    The result is a numerical diagnostic, not an error-certified mode metric.
    """
    g = _complex_array(gram, "gram")
    if g.ndim != 2 or g.shape[0] != g.shape[1] or g.shape[0] == 0:
        raise ValueError("gram must be a nonempty square matrix")
    with np.errstate(over="ignore", invalid="ignore"):
        scale = float(np.max(np.abs(g)))
    if not math.isfinite(scale) or scale <= 0:
        raise ValueError("gram must be positive definite in the finite numerical domain")
    scaled = g / scale
    if float(np.max(np.abs(scaled-scaled.conj().T))) > _HERMITIAN_RTOL:
        raise ValueError("gram must be Hermitian")
    # Scale before symmetrization to avoid overflow on finite large inputs.
    scaled = (scaled + scaled.conj().T) / 2
    try:
        lower = np.linalg.cholesky(scaled)
        frame = np.linalg.solve(lower.conj().T, np.eye(len(g), dtype=complex))
    except np.linalg.LinAlgError as exc:
        raise ValueError("gram must be positive definite") from exc
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        frame = frame / math.sqrt(scale)
    return _finite_result(frame, "canonical frame")


def canonical_yukawa(tensor, gram_left, gram_right, gram_higgs, *, kappa):
    """Return the lower-index Yukawa tensor in canonical field coordinates.

    The supplied tensor shape is (n_left, n_right, n_higgs). This returns
    kappa * sum(C_left[a,i] C_right[b,j] C_higgs[c,k] Y[a,b,c]).
    ``kappa`` is required and may be any finite complex scalar, including
    zero. It contains the overall parent superpotential normalization.
    Rank preservation statements require nonzero kappa; mass magnitudes
    remain conditional on these INPUT metrics and parent coefficient.
    """
    y = _complex_array(tensor, "tensor")
    if y.ndim != 3 or any(size == 0 for size in y.shape):
        raise ValueError("tensor must have three nonempty field indices")
    coefficient = _kappa(kappa)
    frames = tuple(canonical_frame(g) for g in (gram_left, gram_right, gram_higgs))
    if tuple(len(c) for c in frames) != y.shape:
        raise ValueError("Gram dimensions must match the three tensor field indices")
    with np.errstate(over="ignore", invalid="ignore"):
        normalized = coefficient * np.einsum("ai,bj,ck,abc->ijk", *frames, y, optimize=True)
    return _finite_result(normalized, "canonical Yukawa tensor")


def canonical_mass_matrix(tensor, gram_left, gram_right, gram_higgs,
                          higgs_direction, *, vev, kappa):
    """Return M=(vev/sqrt(2)) sum_k Y_canonical[:,:,k] n[k].

    ``higgs_direction`` n is required to be a unit vector in the canonical
    Higgs frame: n^dagger n=1 to absolute tolerance 1e-10. It is never
    silently normalized. The neutral expectation value is (vev/sqrt(2))*n.
    ``vev`` is explicit, finite and nonnegative in the caller's declared
    physical units. With multiple Higgs sectors it is this sector's scale,
    not automatically the total electroweak scale. Both vev and kappa are
    required INPUT; no Standard Model value or observed mass is assumed.
    """
    scale = _vev(vev)
    normalized = canonical_yukawa(tensor, gram_left, gram_right, gram_higgs, kappa=kappa)
    direction = _complex_array(higgs_direction, "higgs_direction")
    if direction.ndim != 1 or direction.shape != (normalized.shape[2],):
        raise ValueError("higgs_direction must match the canonical Higgs field index")
    with np.errstate(over="ignore", invalid="ignore"):
        norm_squared = float(np.vdot(direction, direction).real)
    if not math.isfinite(norm_squared) or abs(norm_squared-1) > _UNIT_ATOL:
        raise ValueError("higgs_direction must be unit normalized in the canonical frame")
    with np.errstate(over="ignore", invalid="ignore"):
        matrix = (scale/math.sqrt(2))*np.einsum("ijk,k->ij", normalized, direction)
    return _finite_result(matrix, "canonical mass matrix")


def mass_singular_values(mass_matrix):
    """Return descending tree-level Dirac mass singular values of canonical M.

    The supplied finite matrix may be rectangular. There are min(m,n)
    singular values; this does not identify chiral unpaired modes or choose
    a light-family sector. No metric, Higgs scale, or unit convention is
    inferred from an unnormalized matrix.
    """
    matrix = _complex_array(mass_matrix, "mass_matrix")
    if matrix.ndim != 2 or any(size == 0 for size in matrix.shape):
        raise ValueError("mass_matrix must have two nonempty field indices")
    try:
        singular = np.linalg.svd(matrix, compute_uv=False)
    except np.linalg.LinAlgError as exc:
        raise ValueError("mass singular-value decomposition did not converge") from exc
    return _finite_result(singular, "mass singular values")
