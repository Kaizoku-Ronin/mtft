"""Matter-01: explicitly labelled reference quantum-field calculations.

The QED kernel is a benchmark for future MTFT reductions.  It is not a claim
that the Standard Model spectrum or physical parameters have been derived.
"""

from .qed import (
    ALPHA_REFERENCE, ELECTRON_MASS_GEV, GEV2_TO_NB, REFERENCE_INPUTS,
    Kinematics, amplitude_sample, boost, centre_of_mass, diagram_amplitudes,
    differential_cross_section, dirac_adjoint, gamma_matrices,
    integrate_cross_section, minkowski_dot, slash, transverse_polarizations,
    u_spinor, unpolarized_amplitude_squared, v_spinor, validate_kinematics,
)

__all__ = [
    "ALPHA_REFERENCE", "ELECTRON_MASS_GEV", "GEV2_TO_NB", "REFERENCE_INPUTS",
    "Kinematics", "amplitude_sample", "boost", "centre_of_mass", "diagram_amplitudes",
    "differential_cross_section", "dirac_adjoint", "gamma_matrices",
    "integrate_cross_section", "minkowski_dot", "slash", "transverse_polarizations",
    "u_spinor", "unpolarized_amplitude_squared", "v_spinor", "validate_kinematics",
]
