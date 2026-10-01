"""Matter-01: a massive tree-level QED reference process.

This module implements e- e+ -> gamma gamma in natural GeV units, with metric
diag(+1,-1,-1,-1), spinor normalization ubar u=2m, vbar v=-2m, and
e**2=4*pi*alpha.  It averages the two initial spins independently (factor 1/4)
and sums the two transverse polarizations of each final photon.  Both electron
exchange diagrams are included.  An overall phase of the amplitude is omitted.

The electron mass and alpha below are frozen *reference inputs*, not outputs of
MTFT, and not a declaration that these are the latest measured constants.
The module is ordinary QED.  No geometric mode, overlap tensor, quantum gravity,
running coupling, radiative correction, or real-time many-body evolution is
silently substituted into this calculation.

The angular cross sections describe EVENTS over the full photon sphere.  They
include the identical-photon factor 1/2.  An inclusive photon-number histogram
has twice that normalization.  The mass is retained at all angles, including
the endpoints; the massless approximation must not be integrated to them.

Internally extended precision is used where NumPy/platform support it.  Results
are floating point diagnostics, not interval-certified bounds.
"""

from dataclasses import dataclass
from numbers import Real
from typing import Dict, Iterable, Optional, Tuple

import numpy as np


ELECTRON_MASS_GEV = 0.00051099895069
ALPHA_REFERENCE = 1.0 / 137.035999177
# Rounded table conversion: 1 GeV^-2 = 0.3893793721 mb = 389379.3721 nb.
GEV2_TO_NB = 389379.3721

REFERENCE_INPUTS = {
    "electron_mass_gev": ELECTRON_MASS_GEV,
    "alpha": ALPHA_REFERENCE,
    "alpha_inverse": 137.035999177,
    "gev_inverse_squared_to_nb": GEV2_TO_NB,
    "status": "INPUT: frozen reference constants; not MTFT predictions",
    "constant_provenance": "PDG 2026 physical-constants review",
    "source_url": "https://pdg.lbl.gov/2026/reviews/rpp2026-rev-phys-constants.pdf",
}

_R = np.longdouble
_C = np.clongdouble
_PI = np.arccos(_R(-1))
_I2 = np.eye(2, dtype=_C)
_Z2 = np.zeros((2, 2), dtype=_C)
_PAULI = (
    np.array([[0, 1], [1, 0]], dtype=_C),
    np.array([[0, -1j], [1j, 0]], dtype=_C),
    np.array([[1, 0], [0, -1]], dtype=_C),
)
_GAMMA = (np.block([[_I2, _Z2], [_Z2, -_I2]]),) + tuple(
    np.block([[_Z2, sigma], [-sigma, _Z2]]) for sigma in _PAULI
)
_I4 = np.eye(4, dtype=_C)


def _finite_real(value: float, name: str) -> np.longdouble:
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Real):
        raise ValueError("{} must be a finite real number".format(name))
    result = _R(value)
    if not np.isfinite(result):
        raise ValueError("{} must be a finite real number".format(name))
    return result


def _positive(value: float, name: str) -> np.longdouble:
    result = _finite_real(value, name)
    if result <= 0:
        raise ValueError("{} must be finite and positive".format(name))
    return result


def _four_vector(value: Iterable[float], name: str = "four-vector") -> np.ndarray:
    uncast = np.asarray(value)
    if uncast.shape != (4,) or uncast.dtype.kind not in "iuf":
        raise ValueError("{} must contain four finite real components".format(name))
    result = np.asarray(uncast, dtype=_R)
    if not np.all(np.isfinite(result)):
        raise ValueError("{} must contain four finite real components".format(name))
    return result


def _dot(a: np.ndarray, b: np.ndarray) -> np.longdouble:
    return a[0] * b[0] - np.dot(a[1:], b[1:])


def minkowski_dot(a: Iterable[float], b: Iterable[float]) -> float:
    """Return a.b for the (+---) metric; order is (energy, px, py, pz)."""
    return float(_dot(_four_vector(a), _four_vector(b)))


def gamma_matrices() -> Tuple[np.ndarray, ...]:
    """Return copies of the four Dirac-representation gamma matrices."""
    return tuple(g.copy() for g in _GAMMA)


def slash(momentum: Iterable[float]) -> np.ndarray:
    """Return gamma^mu p_mu = gamma^0 E - gamma^i p_i."""
    p = _four_vector(momentum)
    return _GAMMA[0] * p[0] - sum(_GAMMA[i] * p[i] for i in range(1, 4))


def _validate_mass_shell(p: np.ndarray, mass: np.longdouble, name: str) -> None:
    if p[0] <= 0:
        raise ValueError("{} energy must be positive".format(name))
    scale = max(p[0] ** 2, mass ** 2, np.finfo(_R).tiny)
    if abs(_dot(p, p) - mass ** 2) > _R(2e-11) * scale:
        raise ValueError("{} is not on its stated mass shell".format(name))


def _spin_label(spin: int) -> int:
    if (isinstance(spin, (bool, np.bool_)) or
            not isinstance(spin, (int, np.integer)) or spin not in (0, 1)):
        raise ValueError("spin must be the integer 0 or 1")
    return int(spin)


def _spinor(p: np.ndarray, mass: np.longdouble, spin: int, antiparticle: bool):
    xi = np.zeros(2, dtype=_C)
    xi[spin] = 1
    lower = sum(_PAULI[i] * p[i + 1] for i in range(3)) @ xi / (p[0] + mass)
    return np.sqrt(p[0] + mass) * (
        np.concatenate((lower, xi)) if antiparticle else np.concatenate((xi, lower))
    )


def u_spinor(momentum: Iterable[float], spin: int, mass_gev=ELECTRON_MASS_GEV):
    """Positive-energy electron spinor; ubar u=2m, u^dagger u=2E."""
    p, m = _four_vector(momentum), _positive(mass_gev, "mass_gev")
    _validate_mass_shell(p, m, "electron")
    return _spinor(p, m, _spin_label(spin), False)


def v_spinor(momentum: Iterable[float], spin: int, mass_gev=ELECTRON_MASS_GEV):
    """Positive-energy positron spinor; vbar v=-2m, v^dagger v=2E."""
    p, m = _four_vector(momentum), _positive(mass_gev, "mass_gev")
    _validate_mass_shell(p, m, "positron")
    return _spinor(p, m, _spin_label(spin), True)


def dirac_adjoint(spinor: Iterable[complex]) -> np.ndarray:
    """Return the row spinor psi^dagger gamma^0."""
    psi = np.asarray(spinor, dtype=_C)
    if psi.shape != (4,) or not np.all(np.isfinite(psi)):
        raise ValueError("spinor must have four finite complex components")
    return psi.conjugate() @ _GAMMA[0]


def boost(vector: Iterable[float], velocity: Iterable[float]) -> np.ndarray:
    """Active Lorentz boost by a three-velocity, in units c=1.

    E' = gamma (E + velocity.p).  Apply this same transformation to all external
    momenta and polarization four-vectors.  Recomputing the unpolarized spin sum
    in the new frame checks covariance without a spin-basis convention problem.
    """
    p = _four_vector(vector)
    uncast = np.asarray(velocity)
    if uncast.shape != (3,) or uncast.dtype.kind not in "iuf":
        raise ValueError("velocity must contain three finite real components")
    b = np.asarray(uncast, dtype=_R)
    if not np.all(np.isfinite(b)):
        raise ValueError("velocity must contain three finite real components")
    b2 = np.dot(b, b)
    if b2 >= 1:
        raise ValueError("boost speed must be strictly smaller than 1")
    if b2 == 0:
        return p.copy()
    gamma = 1 / np.sqrt(1 - b2)
    bp = np.dot(b, p[1:])
    result = np.empty(4, dtype=_R)
    result[0] = gamma * (p[0] + bp)
    result[1:] = p[1:] + (((gamma - 1) * bp / b2) + gamma * p[0]) * b
    return result


@dataclass(frozen=True)
class Kinematics:
    """External momenta (GeV); p_minus,p_plus incoming; k1,k2 outgoing."""

    p_minus: np.ndarray
    p_plus: np.ndarray
    k1: np.ndarray
    k2: np.ndarray
    mass_gev: float = ELECTRON_MASS_GEV

    def __post_init__(self):
        mass = _positive(self.mass_gev, "mass_gev")
        for name in ("p_minus", "p_plus", "k1", "k2"):
            value = _four_vector(getattr(self, name), name).copy()
            value.flags.writeable = False
            object.__setattr__(self, name, value)
        object.__setattr__(self, "mass_gev", float(mass))
        validate_kinematics(self)

    @property
    def s_gev2(self) -> float:
        p = self.p_minus + self.p_plus
        return float(_dot(p, p))

    def boosted(self, velocity: Iterable[float]):
        """Return a new configuration with all four momenta boosted together."""
        return Kinematics(
            *(boost(getattr(self, name), velocity) for name in
              ("p_minus", "p_plus", "k1", "k2")), mass_gev=self.mass_gev
        )

    def as_dict(self) -> Dict[str, object]:
        return {name: [float(x) for x in getattr(self, name)] for name in
                ("p_minus", "p_plus", "k1", "k2")}


def validate_kinematics(kinematics: Kinematics) -> None:
    """Reject nonfinite, off-shell, nonconserving, or threshold configurations."""
    kin = kinematics
    m = _positive(kin.mass_gev, "mass_gev")
    for name in ("p_minus", "p_plus"):
        _validate_mass_shell(_four_vector(getattr(kin, name), name), m, name)
    for name in ("k1", "k2"):
        _validate_mass_shell(_four_vector(getattr(kin, name), name), _R(0), name)
    incoming, outgoing = kin.p_minus + kin.p_plus, kin.k1 + kin.k2
    if np.max(np.abs(incoming - outgoing)) > _R(2e-11) * max(incoming[0], outgoing[0]):
        raise ValueError("external momenta do not conserve four-momentum")
    s = _dot(incoming, incoming)
    if s <= 4 * m * m:
        raise ValueError("incoming invariant energy must be strictly above 2m")
    # If the mass is below the working precision in E^2-p^2, a putatively
    # 'massive' forward calculation is not trustworthy.  Reject, never replace
    # it silently by a massless amplitude.
    if m * m < 64 * np.finfo(_R).eps * max(kin.p_minus[0] ** 2, kin.p_plus[0] ** 2):
        raise ValueError("massive kinematics exceed this platform's floating-point range")


def centre_of_mass(sqrt_s_gev: float, cos_theta: float, phi: float = 0,
                   mass_gev: float = ELECTRON_MASS_GEV) -> Kinematics:
    """Construct CM momenta; electron points +z and cos_theta describes k1."""
    energy = _positive(sqrt_s_gev, "sqrt_s_gev") / 2
    m = _positive(mass_gev, "mass_gev")
    c, angle = _finite_real(cos_theta, "cos_theta"), _finite_real(phi, "phi")
    if abs(c) > 1:
        raise ValueError("cos_theta must be finite and in [-1,1]")
    if energy <= m:
        raise ValueError("sqrt_s_gev must be strictly greater than 2*mass_gev")
    p = np.sqrt((energy - m) * (energy + m))
    sin_theta = np.sqrt((1 - c) * (1 + c))
    direction = np.array([sin_theta * np.cos(angle), sin_theta * np.sin(angle), c], dtype=_R)
    return Kinematics(
        np.array([energy, 0, 0, p], dtype=_R),
        np.array([energy, 0, 0, -p], dtype=_R),
        np.concatenate((np.array([energy], dtype=_R), energy * direction)),
        np.concatenate((np.array([energy], dtype=_R), -energy * direction)),
        mass_gev=float(m),
    )


def transverse_polarizations(photon_momentum: Iterable[float]) -> Tuple[np.ndarray, np.ndarray]:
    """Two real orthonormal transverse polarizations, epsilon^0=0 in this frame.

    Basis selection uses the least parallel Cartesian axis; no polar endpoint
    singularity is introduced.  Only their summed probabilities are basis
    independent; a particular polarized amplitude depends on that basis.
    """
    k = _four_vector(photon_momentum)
    _validate_mass_shell(k, _R(0), "photon")
    direction = k[1:] / np.sqrt(np.dot(k[1:], k[1:]))
    pivot = np.zeros(3, dtype=_R)
    pivot[int(np.argmin(np.abs(direction)))] = 1
    e1 = np.cross(pivot, direction)
    e1 /= np.sqrt(np.dot(e1, e1))
    e2 = np.cross(direction, e1)
    return (np.concatenate((np.zeros(1, dtype=_R), e1)),
            np.concatenate((np.zeros(1, dtype=_R), e2)))


def _amplitude_operator(kin: Kinematics, e1: np.ndarray, e2: np.ndarray,
                        alpha: np.longdouble):
    m = _R(kin.mass_gev)
    # On-shell identity avoids squaring and subtracting nearly equal p and m.
    t_den = -2 * _dot(kin.p_minus, kin.k1)
    u_den = -2 * _dot(kin.p_minus, kin.k2)
    if t_den >= 0 or u_den >= 0:
        raise ValueError("invalid electron propagator denominators")
    e_squared = 4 * _PI * alpha
    t = e_squared * (slash(e2) @ (slash(kin.p_minus - kin.k1) + m * _I4)
                     @ slash(e1)) / t_den
    u = e_squared * (slash(e1) @ (slash(kin.p_minus - kin.k2) + m * _I4)
                     @ slash(e2)) / u_den
    return t, u


def diagram_amplitudes(p_minus: Iterable[float], p_plus: Iterable[float],
                       k1: Iterable[float], k2: Iterable[float],
                       epsilon1: Iterable[float], epsilon2: Iterable[float],
                       spin_minus: int, spin_plus: int,
                       mass_gev: float = ELECTRON_MASS_GEV,
                       alpha: float = ALPHA_REFERENCE) -> Dict[str, complex]:
    """Return named complex t, u, and total amplitudes for one spin/polarization.

    The t diagram emits photon k1 first along the incoming electron line; the u
    diagram emits k2 first.  Polarizations may be replaced by momenta for Ward
    checks, or shifted by c*k.  They are deliberately not required to be unit
    transverse vectors here.  Real linear-polarization four-vectors are used.
    """
    kin = Kinematics(p_minus, p_plus, k1, k2, mass_gev=mass_gev)
    a = _positive(alpha, "alpha")
    t_op, u_op = _amplitude_operator(kin, _four_vector(epsilon1), _four_vector(epsilon2), a)
    u = _spinor(kin.p_minus, _R(mass_gev), _spin_label(spin_minus), False)
    vbar = dirac_adjoint(_spinor(kin.p_plus, _R(mass_gev), _spin_label(spin_plus), True))
    t, crossed = vbar @ t_op @ u, vbar @ u_op @ u
    return {"t": complex(t), "u": complex(crossed), "total": complex(t + crossed)}


def unpolarized_amplitude_squared(kinematics: Kinematics,
                                  alpha: float = ALPHA_REFERENCE,
                                  polarizations1: Optional[Tuple[np.ndarray, np.ndarray]] = None,
                                  polarizations2: Optional[Tuple[np.ndarray, np.ndarray]] = None) -> float:
    """Spin-averaged, physical-polarization-summed |M|^2 (dimensionless).

    Custom bases must contain two real four-vectors each.  Boost a CM transverse
    basis together with the momenta to test covariance; adding longitudinal
    gauge shifts is allowed.  The function does not silently renormalize them.
    """
    validate_kinematics(kinematics)
    a = _positive(alpha, "alpha")
    kin, m = kinematics, _R(kinematics.mass_gev)
    eps1 = transverse_polarizations(kin.k1) if polarizations1 is None else polarizations1
    eps2 = transverse_polarizations(kin.k2) if polarizations2 is None else polarizations2
    if len(eps1) != 2 or len(eps2) != 2:
        raise ValueError("each physical polarization basis must have two vectors")
    for photon, basis in ((kin.k1, eps1), (kin.k2, eps2)):
        basis = tuple(_four_vector(e) for e in basis)
        for i, e in enumerate(basis):
            scale = max(np.max(np.abs(e)) * np.max(np.abs(photon)), _R(1e-30))
            if abs(_dot(e, photon)) > _R(2e-10) * scale:
                raise ValueError("physical polarization must be transverse to its photon")
            for j, other in enumerate(basis):
                expected = _R(-1 if i == j else 0)
                scale = max(np.max(np.abs(e)) * np.max(np.abs(other)), _R(1))
                if abs(_dot(e, other) - expected) > _R(2e-10) * scale:
                    raise ValueError("physical polarizations must be orthonormal in (+---)")
    electrons = tuple(_spinor(kin.p_minus, m, spin, False) for spin in (0, 1))
    positrons = tuple(dirac_adjoint(_spinor(kin.p_plus, m, spin, True)) for spin in (0, 1))
    total = _R(0)
    for e1 in eps1:
        for e2 in eps2:
            t, crossed = _amplitude_operator(kin, _four_vector(e1), _four_vector(e2), a)
            op = t + crossed
            for u in electrons:
                for vbar in positrons:
                    amp = vbar @ op @ u
                    total += amp.real ** 2 + amp.imag ** 2
    return float(total / 4)


def differential_cross_section(sqrt_s_gev: float, cos_theta: float, phi: float = 0,
                               mass_gev: float = ELECTRON_MASS_GEV,
                               alpha: float = ALPHA_REFERENCE) -> Dict[str, float]:
    """CM event cross sections including the full-sphere photon factor 1/2.

    d_sigma/d_Omega=|M|^2/(128*pi^2*s*beta), and
    d_sigma/d_cos(theta)=2*pi*d_sigma/d_Omega.  Individual labelled-photon
    differential counts are twice these event cross sections.
    """
    kin = centre_of_mass(sqrt_s_gev, cos_theta, phi=phi, mass_gev=mass_gev)
    alpha = _positive(alpha, "alpha")
    energy = _R(sqrt_s_gev)
    ratio = 2 * _R(mass_gev) / energy
    beta = np.sqrt((1 - ratio) * (1 + ratio))
    msq = unpolarized_amplitude_squared(kin, alpha=alpha)
    per_solid_angle = _R(msq) / (128 * _PI ** 2 * energy ** 2 * beta)
    per_costheta = 2 * _PI * per_solid_angle
    if not all(np.isfinite(float(value)) for value in
               (msq, per_solid_angle, per_costheta, per_costheta * GEV2_TO_NB)):
        raise ValueError("cross section exceeds the output floating-point range")
    return {
        "sqrt_s_gev": float(energy), "mass_gev": float(mass_gev),
        "alpha": float(alpha), "cos_theta": float(cos_theta), "phi": float(phi),
        "beta": float(beta), "msq": msq,
        "event_dsigma_domega_gev2": float(per_solid_angle),
        "event_dsigma_dcostheta_gev2": float(per_costheta),
        "event_dsigma_dcostheta_nb": float(per_costheta * GEV2_TO_NB),
    }


def integrate_cross_section(sqrt_s_gev: float, mass_gev: float = ELECTRON_MASS_GEV,
                            alpha: float = ALPHA_REFERENCE, order: int = 96,
                            cos_theta_max: float = 1.0) -> Dict[str, object]:
    """Integrate event d_sigma/d_cos(theta) over [-c_max,+c_max].

    Gauss-Legendre integration uses y=atanh(beta*cos(theta)), resolving the
    narrow massive forward/backward peaks without discarding their endpoint
    neighborhoods.  Convergence still needs comparison with a higher order or
    an independent oracle.  At threshold the cross section diverges as 1/beta;
    beta>0 is required and the physically finite sigma*v limit is separate.
    """
    energy = _positive(sqrt_s_gev, "sqrt_s_gev")
    mass = _positive(mass_gev, "mass_gev")
    _positive(alpha, "alpha")
    if isinstance(order, (bool, np.bool_)) or not isinstance(order, (int, np.integer)) or order < 8:
        raise ValueError("quadrature order must be an integer >=8")
    cut = _finite_real(cos_theta_max, "cos_theta_max")
    if cut <= 0 or cut > 1:
        raise ValueError("cos_theta_max must be in (0,1]")
    # Checks threshold and platform precision before any quadrature work.
    centre_of_mass(energy, 0, mass_gev=mass)
    ratio = 2 * mass / energy
    beta = np.sqrt((1 - ratio) * (1 + ratio))
    ymax = np.arctanh(beta * cut)
    nodes, weights = np.polynomial.legendre.leggauss(int(order))
    summands = []
    for node, weight in zip(nodes, weights):
        y = _R(node) * ymax
        c = np.tanh(y) / beta
        jacobian = ymax / (beta * np.cosh(y) ** 2)
        point = differential_cross_section(energy, c, mass_gev=mass, alpha=alpha)
        summands.append(_R(weight) * jacobian * _R(point["event_dsigma_dcostheta_gev2"]))
    sigma = sum(summands, _R(0))
    return {"sqrt_s_gev": float(energy), "mass_gev": float(mass), "alpha": float(alpha),
            "beta": float(beta), "cos_theta_max": float(cut), "quadrature_order": int(order),
            "event_sigma_gev2": float(sigma), "event_sigma_nb": float(sigma * GEV2_TO_NB),
            "sigma_times_2beta_gev2": float(2 * beta * sigma),
            "integration_variable": "y=atanh(beta*cos(theta))",
            "convention": "events; identical-photon 1/2; full azimuth"}


def amplitude_sample(kinematics: Kinematics, alpha: float = ALPHA_REFERENCE) -> Dict[str, object]:
    """Serialize all 16 physical spin/polarization amplitudes, retaining diagrams."""
    records = []
    eps1 = transverse_polarizations(kinematics.k1)
    eps2 = transverse_polarizations(kinematics.k2)
    for s1 in (0, 1):
        for s2 in (0, 1):
            for pol1 in (0, 1):
                for pol2 in (0, 1):
                    amps = diagram_amplitudes(
                        kinematics.p_minus, kinematics.p_plus, kinematics.k1, kinematics.k2,
                        eps1[pol1], eps2[pol2], s1, s2,
                        mass_gev=kinematics.mass_gev, alpha=alpha)
                    records.append({"spin_minus": s1, "spin_plus": s2,
                                    "polarization1": pol1, "polarization2": pol2,
                                    **{name: {"real": float(value.real), "imag": float(value.imag)}
                                       for name, value in amps.items()}})
    return {"polarizations1": [[float(x) for x in e] for e in eps1],
            "polarizations2": [[float(x) for x in e] for e in eps2], "amplitudes": records}
