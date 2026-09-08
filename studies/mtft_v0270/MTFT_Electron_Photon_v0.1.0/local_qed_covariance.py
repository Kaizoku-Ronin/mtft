"""Finite Euclidean Wilson–Dirac covariance on the MTFT internal space.

Flat 4D Euclidean spacetime, Dirac spin, local U(1), and a Wilson regulator
are supplied. The 13-dimensional internal matrices come from the stored
Hecke experiment; their trial masses and charges are not SM assignments.
Run after charge_study.py to create local_qed_results.json.
"""
import hashlib
import numpy as np
from mtft import GAUGE
from common import ROOT, dump, gamma_matrices, herm, comm


SHAPE = (3, 3, 3, 3)
SPACING = 1.0
WILSON_R = 1.0
SEEDS = range(2026090716, 2026090720)
TOLERANCE = 1e-10


def norm(value):
    """Euclidean norm over every site, spin, and flavour entry."""
    return float(np.linalg.norm(value.ravel()))


def flavour_action(matrix, field):
    return np.einsum('...ij,...sj->...si', matrix, field, optimize=True)


def spin_action(matrix, field):
    return np.einsum('st,...tj->...sj', matrix, field, optimize=True)


def adjoint(matrix):
    return np.swapaxes(matrix.conj(), -2, -1)


def unitary_phases(Q, phases):
    """Apply the spectral theorem to exp(i phases(x) Q), without projection."""
    eigenvalues, vectors = np.linalg.eigh(Q)
    factors = np.exp(1j * phases[..., None] * eigenvalues)
    return np.einsum('ij,...j,kj->...ik', vectors, factors, vectors.conj(), optimize=True)


def wilson_dirac(field, links, mass, gamma, spacing=SPACING, r=WILSON_R):
    """Matrix-free D; mass may be constant or a site-dependent spurion.

    D = M + sum_mu [gamma_mu (U psi_+ - U_-^dagger psi_-)/(2a)
                     + r (2 psi - U psi_+ - U_-^dagger psi_-)/(2a)].
    """
    result = flavour_action(mass, field)
    for mu in range(4):
        link = links[..., mu, :, :]
        forward = flavour_action(link, np.roll(field, -1, axis=mu))
        backward_link = adjoint(np.roll(link, 1, axis=mu))
        backward = flavour_action(backward_link, np.roll(field, 1, axis=mu))
        result = result + spin_action(gamma[mu], forward-backward)/(2*spacing)
        result = result + r*(2*field-forward-backward)/(2*spacing)
    return result


def transform_links(links, omega):
    transformed = np.empty_like(links)
    for mu in range(4):
        transformed[..., mu, :, :] = (
            omega @ links[..., mu, :, :] @ adjoint(np.roll(omega, -1, axis=mu))
        )
    return transformed


def covariance_record(field, links, mass, gamma, omega):
    before = flavour_action(omega, wilson_dirac(field, links, mass, gamma))
    transformed_links = transform_links(links, omega)
    transformed_field = flavour_action(omega, field)
    after = wilson_dirac(transformed_field, transformed_links, mass, gamma)
    defect = after-before
    expected = flavour_action(mass @ omega - omega @ mass, field)
    scale = max(norm(before), norm(after), np.finfo(float).tiny)
    spurion = omega @ mass @ adjoint(omega)
    repaired = wilson_dirac(transformed_field, transformed_links, spurion, gamma)
    repaired_scale = max(norm(before), norm(repaired), np.finfo(float).tiny)
    expected_norm = norm(expected)
    return {
        'normalization': scale,
        'transformed_D_field_norm': norm(after),
        'omega_D_field_norm': norm(before),
        'absolute_covariance_defect': norm(defect),
        'relative_covariance_defect': norm(defect)/scale,
        'absolute_predicted_mass_defect': expected_norm,
        'relative_predicted_mass_defect': expected_norm/scale,
        'relative_defect_identity_error': norm(defect-expected)/scale,
        'defect_relative_identity_error': (
            norm(defect-expected)/expected_norm if expected_norm > 1e-8 else None
        ),
        'spurion_covariance_normalization': repaired_scale,
        'relative_spurion_covariance_defect': norm(repaired-before)/repaired_scale,
        'transformed_link_unitarity_max': float(np.max(np.linalg.norm(
            transformed_links @ adjoint(transformed_links)-np.eye(13), axis=(-2,-1)
        ))),
    }


def main():
    frozen_path = ROOT/'charge_operators.npz'
    with np.load(frozen_path) as operators:
        Ps = operators['projectors'].copy()
        mass = herm(operators['M0'] + 0.2*operators['balanced_V'])
    gamma = gamma_matrices()
    gamma[1:] *= 1j
    gamma5 = gamma[0] @ gamma[1] @ gamma[2] @ gamma[3]
    clifford = max(norm(gamma[i]@gamma[j]+gamma[j]@gamma[i]
                       -2*(i == j)*np.eye(4)) for i in range(4) for j in range(4))
    structural = {
        'Euclidean_Clifford_absolute_residual': clifford,
        'gamma_Hermiticity_absolute_residual': max(norm(g-adjoint(g)) for g in gamma),
        'gamma5_Hermiticity_absolute_residual': norm(gamma5-adjoint(gamma5)),
        'gamma5_squared_identity_absolute_residual': norm(gamma5@gamma5-np.eye(4)),
        'mass_Hermiticity_absolute_residual': norm(mass-adjoint(mass)),
    }
    coupling = float(np.sqrt(4*np.pi*GAUGE.alpha))
    rows = []
    assignments = [
        ('common_charge', [-1,-1,-1,-1], -np.eye(13)),
        ('equal_mixed_charges', [0,-1,-1,2], herm(-Ps[1]-Ps[2]+2*Ps[3])),
        ('different_mixed_charges', [0,-1,0,2], herm(-Ps[1]+2*Ps[3])),
    ]
    for seed in SEEDS:
        rng = np.random.default_rng(seed)
        potential = rng.normal(size=SHAPE+(4,))
        chi = rng.normal(size=SHAPE)
        field = rng.normal(size=SHAPE+(4,13))+1j*rng.normal(size=SHAPE+(4,13))
        probe = rng.normal(size=SHAPE+(4,13))+1j*rng.normal(size=SHAPE+(4,13))
        for label, charges, Q in assignments:
            links = unitary_phases(Q, coupling*SPACING*potential)
            omega = unitary_phases(Q, chi)
            record = covariance_record(field, links, mass, gamma, omega)
            D_field = wilson_dirac(field, links, mass, gamma)
            Ddag_probe = spin_action(gamma5, wilson_dirac(
                spin_action(gamma5, probe), links, mass, gamma))
            lhs = np.vdot(probe, D_field)
            rhs = np.vdot(Ddag_probe, field)
            adjoint_scale = max(norm(probe)*norm(D_field),
                                norm(Ddag_probe)*norm(field), np.finfo(float).tiny)
            record.update({
                'seed': seed, 'case': label, 'sector_charges': charges,
                'mass_charge_commutator_absolute_norm': norm(comm(mass,Q)),
                'link_unitarity_max': float(np.max(np.linalg.norm(
                    links @ adjoint(links)-np.eye(13), axis=(-2,-1)))),
                'omega_unitarity_max': float(np.max(np.linalg.norm(
                    omega @ adjoint(omega)-np.eye(13), axis=(-2,-1)))),
                'gamma5_adjoint_relative_residual': float(abs(lhs-rhs)/adjoint_scale),
                'gamma5_adjoint_normalization': adjoint_scale,
            })
            rows.append(record)

    commuting = [r for r in rows if r['case'] != 'different_mixed_charges']
    forbidden = [r for r in rows if r['case'] == 'different_mixed_charges']
    summary = {
        'commuting_covariance_max': max(r['relative_covariance_defect'] for r in commuting),
        'different_charge_covariance_min': min(r['relative_covariance_defect'] for r in forbidden),
        'different_charge_covariance_max': max(r['relative_covariance_defect'] for r in forbidden),
        'defect_identity_error_max': max(r['relative_defect_identity_error'] for r in rows),
        'spurion_covariance_max': max(r['relative_spurion_covariance_defect'] for r in rows),
        'gamma5_adjoint_error_max': max(r['gamma5_adjoint_relative_residual'] for r in rows),
        'link_unitarity_max': max(r['link_unitarity_max'] for r in rows),
        'transformed_link_unitarity_max': max(r['transformed_link_unitarity_max'] for r in rows),
        'omega_unitarity_max': max(r['omega_unitarity_max'] for r in rows),
    }
    checks = {
        'structural_checks': max(structural.values()) < TOLERANCE,
        'commuting_cases_covariant': summary['commuting_covariance_max'] < TOLERANCE,
        'distinct_charges_fail_covariance': summary['different_charge_covariance_min'] > 1e-4,
        'mass_defect_identity': summary['defect_identity_error_max'] < TOLERANCE,
        'spurion_restores_covariance': summary['spurion_covariance_max'] < TOLERANCE,
        'gamma5_adjoint_check': summary['gamma5_adjoint_error_max'] < TOLERANCE,
        'unitarity': max(summary[key] for key in [
            'link_unitarity_max','transformed_link_unitarity_max','omega_unitarity_max'
        ]) < TOLERANCE,
    }
    assert all(checks.values()), checks
    result = {
        'status': 'DIAGNOSTIC finite matrix-free gauge-covariance experiment; all gates passed',
        'configuration': {
            'lattice_shape': list(SHAPE), 'boundary_conditions': 'periodic in all four directions',
            'lattice_spacing': SPACING, 'Wilson_r': WILSON_R, 'spin_components': 4,
            'internal_complex_dimensions': 13, 'total_complex_field_entries': int(np.prod(SHAPE)*4*13),
            'seeds': list(SEEDS), 'configuration_count': len(rows),
            'mass_operator': 'M0 + 0.2 balanced_V from charge_operators.npz',
            'mass_units': 'dimensionless trial values, not observed particle masses',
            'gauge_coupling': coupling, 'alpha_inv_MTFT_three_term': float(GAUGE.alpha_inv),
            'inputs_distribution': 'A_mu, chi, and each real/imaginary spinor entry: independent N(0,1)',
            'numerical_tolerance': TOLERANCE,
        },
        'normalizations': {
            'covariance': '||D_Uprime(Omega psi)-Omega D_U psi||_2 / max(||D_Uprime(Omega psi)||_2,||Omega D_U psi||_2)',
            'defect_identity': '||covariance_defect-[M,Omega]psi||_2 / same covariance denominator',
            'spurion_covariance': 'Same ratio with Mprime(x)=Omega(x) M Omega(x)^dagger in D_Uprime',
            'gamma5_adjoint': '|<phi,Dpsi>-<gamma5 D gamma5 phi,psi>| / max(||phi||_2 ||Dpsi||_2,||gamma5 D gamma5 phi||_2 ||psi||_2)',
            'unitarity': 'Maximum per-site/link Frobenius norm ||U U^dagger-I||_F',
        },
        'input_sha256': hashlib.sha256(frozen_path.read_bytes()).hexdigest(),
        'structural_residuals': structural, 'summary': summary, 'checks': checks, 'runs': rows,
        'interpretation': [
            'The supplied local U(1) lattice action is compatible with the Hecke mixing exactly when the fixed mass commutes with charge.',
            'Different old and q4 charges give a nonzero local mass covariance defect on these backgrounds.',
            'Spurion covariance only illustrates the necessary transformation law of a charged mass source; no dynamical scalar or Higgs mechanism is derived.',
            'This is a finite Euclidean regulator test. No continuum spectrum, Lorentz reconstruction, electron identity, charge quantization, or Standard Model embedding is established.',
        ],
    }
    dump('local_qed_results.json', result)
    print('Finite local gauge covariance:', len(rows), 'configurations; all gates passed.')
    for key, value in summary.items():
        print(key + ':', value)


if __name__ == '__main__':
    main()
