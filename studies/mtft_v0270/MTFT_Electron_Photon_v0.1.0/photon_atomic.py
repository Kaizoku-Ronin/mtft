"""Conditional Maxwell and Coulomb controls; no spacetime or spin derivation.

Run from any directory. Outputs photon_atomic_results.json beside this script.
Spacetime dimension, Maxwell dynamics, and the Schroedinger Coulomb Hamiltonian
are supplied physical assumptions. The outputs here use dimensionless units;
the common electromagnetic coupling is propagated separately by the main study.
"""
import numpy as np
from scipy.linalg import eigh_tridiagonal

from common import dump


def curl_symbol(d):
    """Forward-difference curl: C(d) A = d cross A."""
    dx, dy, dz = d
    return np.array([[0, -dz, dy], [dz, 0, -dx], [-dy, dx, 0]], complex)


def maxwell_controls():
    """Analyze selected Fourier symbols of a periodic spatial discretization.

    Time remains continuous, c is set to 1, and box length is 2*pi. This is not
    an N**3 real-space simulation. Translation invariance reduces each sampled
    Fourier momentum to a 3 by 3 matrix exactly at the discrete-operator level.
    """
    records, ratios = [], []
    for mode in [(1, 0, 0), (1, 1, 0), (1, 2, 3)]:
        k = np.asarray(mode, float)
        continuum_omega = np.linalg.norm(k)
        errors = []
        for N in [16, 32, 64, 128]:
            a = 2 * np.pi / N
            d = np.expm1(1j * k * a) / a
            C = curl_symbol(d)
            K = C.conj().T @ C
            d2 = float(np.vdot(d, d).real)
            eig, vectors = np.linalg.eigh(K)
            expected = np.array([0, d2, d2])
            numerical_omega = float(np.sqrt(np.mean(eig[1:])))
            speed_ratio = numerical_omega / continuum_omega
            error = 1 - speed_ratio
            errors.append(float(error))
            records.append({
                'mode': list(mode), 'points_per_axis': N, 'spacing': a,
                'omega_continuum_c1': float(continuum_omega),
                'omega_discrete': numerical_omega,
                'phase_speed_over_assumed_c': speed_ratio,
                'relative_dispersion_error': float(error),
                'eigenvalues_omega_squared': eig.tolist(),
                'expected_eigenvalues': expected.tolist(),
                'transverse_mode_count': int(np.count_nonzero(eig > 1e-10*d2)),
                'gauge_nullity': int(np.count_nonzero(np.abs(eig) < 1e-10*d2)),
                'spectral_identity_relative_residual': float(np.max(np.abs(eig-expected))/d2),
                'matrix_identity_relative_residual': float(np.linalg.norm(K-(d2*np.eye(3)-np.outer(d,d.conj())))/d2),
                'gauge_kernel_relative_residual': float(np.linalg.norm(K@d)/(d2*np.linalg.norm(d))),
                'transversality_relative_residual': float(np.linalg.norm(d.conj()@vectors[:,1:])/np.linalg.norm(d)),
                'eigenpair_relative_residual': float(np.linalg.norm(K@vectors-vectors*eig)/d2),
            })
        ratios.append({'mode': list(mode), 'coarse_to_fine_error_ratios':
                       [errors[i]/errors[i+1] for i in range(len(errors)-1)]})
    zero = curl_symbol(np.zeros(3))
    all_ratios = [v for record in ratios for v in record['coarse_to_fine_error_ratios']]
    worst = max(r[field] for r in records for field in [
        'spectral_identity_relative_residual', 'matrix_identity_relative_residual',
        'gauge_kernel_relative_residual', 'transversality_relative_residual',
        'eigenpair_relative_residual'])
    gates = {
        'two_transverse_one_gauge_at_each_sampled_nonzero_momentum': all(
            r['transverse_mode_count'] == 2 and r['gauge_nullity'] == 1 for r in records),
        'relative_algebra_residual_below_1e_minus_10': worst < 1e-10,
        'dispersion_refinement_ratios_between_3_and_5': all(3 < q < 5 for q in all_ratios),
        'zero_momentum_curl_exactly_zero_in_float_arithmetic': bool(np.all(zero == 0)),
    }
    return {
        'classification': 'DIAGNOSTIC: imported Maxwell theory and spatial geometry',
        'units': 'c=1; periodic spatial box length 2*pi; continuous time',
        'operator': 'd_i=(exp(i*k_i*a)-1)/a; K=C(d)^dagger C(d)=|d|^2 I-d d^dagger',
        'mode_selection': 'Three fixed physical Fourier momenta; not a scan over every grid momentum.',
        'zero_momentum': {
            'K': (zero.conj().T@zero).real.tolist(), 'kernel_dimension': 3,
            'interpretation': 'Separate homogeneous zero-frequency sector. The nonzero-momentum two-polarization count is not applied to k=0.'},
        'records': records, 'refinement': ratios,
        'worst_relative_algebra_residual': worst, 'gates': gates,
        'all_gates_passed': all(gates.values()),
        'limitations': [
            'The assumed continuum c=1 is recovered under refinement, not predicted by MTFT.',
            'No time discretization, matter coupling, photon quantization, or Lorentz restoration theorem is tested.',
            'Two transverse modes follow from the supplied Maxwell operator in three spatial dimensions.',
        ],
    }


def radial_solve(ell, spacing, radius):
    intervals = int(round(radius / spacing))
    assert abs(intervals*spacing-radius) < 1e-12
    r = spacing*np.arange(1, intervals)
    potential = ell*(ell+1)/(2*r*r)-1/r
    diagonal = 1/spacing**2 + potential
    off_diagonal = np.full(intervals-2, -0.5/spacing**2)
    # For fixed ell the first radial state has principal n=ell+1.
    count = 3-ell
    energies, vectors = eigh_tridiagonal(
        diagonal, off_diagonal, select='i', select_range=(0,count-1),
        eigvals_only=False, check_finite=True, lapack_driver='stebz', tol=0)
    Hv = diagonal[:,None]*vectors
    Hv[:-1] += off_diagonal[:,None]*vectors[1:]
    Hv[1:] += off_diagonal[:,None]*vectors[:-1]
    residuals = np.linalg.norm(Hv-vectors*energies, axis=0)
    records = []
    for idx, energy in enumerate(energies):
        n = ell+1+idx
        exact = -1/(2*n*n)
        records.append({
            'n': n, 'ell': ell, 'spacing_bohr': spacing, 'box_radius_bohr': radius,
            'interior_grid_points': intervals-1,
            'energy_hartree': float(energy), 'continuum_energy_hartree': exact,
            'signed_continuum_error_hartree': float(energy-exact),
            'absolute_continuum_error_hartree': float(abs(energy-exact)),
            'relative_continuum_error': float(abs((energy-exact)/exact)),
            'eigenpair_l2_residual_hartree': float(residuals[idx]),
            'eigenpair_residual_relative_to_abs_energy': float(residuals[idx]/abs(energy)),
        })
    return records


def coulomb_controls():
    spacings = [.1, .05, .025, .0125]
    records, refinements, boundary = [], [], []
    for ell in [0,1,2]:
        local = []
        for h in spacings:
            local.extend(radial_solve(ell,h,100.0))
        records.extend(local)
        for n in range(ell+1,4):
            states = [r for r in local if r['n'] == n]
            errors = [r['absolute_continuum_error_hartree'] for r in states]
            refinements.append({'n': n, 'ell': ell,
                'coarse_to_fine_error_ratios': [errors[i]/errors[i+1] for i in range(3)]})
        large_box = radial_solve(ell,spacings[-1],150.0)
        for large in large_box:
            small = next(r for r in local if r['n'] == large['n'] and r['spacing_bohr'] == spacings[-1])
            delta = large['energy_hartree']-small['energy_hartree']
            boundary.append({
                'n': large['n'], 'ell': ell, 'spacing_bohr': spacings[-1],
                'energy_R100_hartree': small['energy_hartree'],
                'energy_R150_hartree': large['energy_hartree'],
                'signed_boundary_shift_hartree': delta,
                'absolute_boundary_shift_hartree': abs(delta),
                'relative_boundary_shift': abs(delta/large['continuum_energy_hartree']),
                'R150_eigenpair_l2_residual_hartree': large['eigenpair_l2_residual_hartree'],
            })
    finest = [r for r in records if r['spacing_bohr'] == spacings[-1]]
    ground = next(r for r in finest if r['n'] == 1 and r['ell'] == 0)
    ratios = [q for row in refinements for q in row['coarse_to_fine_error_ratios']]
    # Independent closed form for the infinite radial grid's ell=0 ground state.
    # Stable rationalization avoids cancellation in (1-sqrt(1+h*h))/(h*h).
    closed = []
    for row in records:
        if row['n'] == 1 and row['ell'] == 0:
            h = row['spacing_bohr']
            expected = -1/(1+np.sqrt(1+h*h))
            closed.append({'spacing_bohr': h, 'infinite_grid_ground_energy_hartree': float(expected),
                'finite_box_solver_minus_infinite_grid_hartree': float(row['energy_hartree']-expected)})
    gates = {
        'finest_ground_relative_continuum_error_below_1e_minus_4': ground['relative_continuum_error'] < 1e-4,
        'all_refinement_ratios_between_3_and_5': all(3 < q < 5 for q in ratios),
        # These are numerical consistency diagnostics, not interval certification.
        'max_relative_boundary_shift_below_1e_minus_8': max(r['relative_boundary_shift'] for r in boundary) < 1e-8,
        'max_absolute_eigenpair_residual_below_1e_minus_8_hartree': max(
            [r['eigenpair_l2_residual_hartree'] for r in records] +
            [r['R150_eigenpair_l2_residual_hartree'] for r in boundary]) < 1e-8,
        'ground_infinite_grid_control_agrees_below_1e_minus_8_hartree': max(
            abs(r['finite_box_solver_minus_infinite_grid_hartree']) for r in closed) < 1e-8,
    }
    return {
        'classification': 'DIAGNOSTIC: imported nonrelativistic Coulomb Hamiltonian',
        'units': 'Bohr radius a0=1/(m_e*alpha), Hartree E_h=m_e*alpha^2, hbar=c=1',
        'hamiltonian': '-(1/2)u_second_derivative+[ell*(ell+1)/(2*r^2)-1/r]*u',
        'boundary_conditions': 'u(0)=u(R)=0; positive interior grid; point nucleus with infinite mass',
        'solver': 'scipy.linalg.eigh_tridiagonal; selected bottom 3-ell eigenpairs; LAPACK stebz, tol=0',
        'normalization': 'Eigenvectors have unit discrete Euclidean norm; residuals refer to this same norm.',
        'records': records, 'finest_grid': finest, 'refinement': refinements,
        'boundary_check': boundary, 'independent_ground_control': closed,
        'independent_ground_derivation': (
            'On the infinite positive radial grid, u_j=j*exp(-j*asinh(h)) solves the ell=0 '
            'difference equation: sinh(asinh(h))/h^2=1/h cancels the 1/j Coulomb term. '
            'Its eigenvalue is (1-sqrt(1+h^2))/h^2=-1/(1+sqrt(1+h^2)). '
            'This verifies the discretized ground state independently of the eigensolver; '
            'the continuum value remains -1/2.'),
        'gates': gates, 'all_gates_passed': all(gates.values()),
        'limitations': [
            'Finite-difference convergence is distinct from a small matrix eigensolver residual.',
            'No reduced-mass correction, relativity, radiative shift, nuclear structure, or hyperfine effects.',
            'The Coulomb law and electron mass are supplied; this is not an independent MTFT spectrum prediction.',
            'The boundary comparison is numerical evidence, not a rigorous tail bound.',
        ],
    }


def main():
    result = {'maxwell': maxwell_controls(), 'coulomb': coulomb_controls()}
    result['all_gates_passed'] = all(result[k]['all_gates_passed'] for k in ['maxwell','coulomb'])
    dump('photon_atomic_results.json',result)
    print('All photon/Coulomb gates passed:',result['all_gates_passed'])
    print('Worst Maxwell algebra relative residual:',result['maxwell']['worst_relative_algebra_residual'])
    print('Finest Coulomb states:')
    for row in result['coulomb']['finest_grid']:
        print('n, ell, energy, relative continuum error:',row['n'],row['ell'],row['energy_hartree'],row['relative_continuum_error'])
    assert result['all_gates_passed']


if __name__ == '__main__':
    main()
