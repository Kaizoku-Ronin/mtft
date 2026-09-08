"""Conditional QED coupling propagation and a Pauli-operator freedom test.

The loop expression is an imported standard QED parameter integral, not a
new amplitude derivation. The added Pauli coefficient is defined at low
energy so it shifts F2(0). No matching from MTFT to this coefficient exists
in this experiment. Corrections of order alpha*kappa and higher are omitted.
"""
import math
import numpy as np
import sympy as sp
from scipy.integrate import quad
from mtft.constants import GAUGE, LEPTONS
from common import (dump, gamma_matrices, ALPHA_INV_REFERENCE,
                    ALPHA_INV_REFERENCE_UNCERTAINTY, ELECTRON_REST_ENERGY_EV,
                    REFERENCE_URL)


def exact_pauli_transversality():
    gs = [sp.Matrix(g).applyfunc(sp.nsimplify) for g in gamma_matrices()]
    q = sp.symbols('q0:4')  # Covariant momentum components.
    out = sp.zeros(4)
    for mu in range(4):
        for nu in range(4):
            sigma = sp.I*(gs[mu]*gs[nu]-gs[nu]*gs[mu])/2
            out += q[mu]*q[nu]*sigma
    out = out.applyfunc(sp.expand)
    assert out == sp.zeros(4)
    return {'status':'EXACT identity in the supplied Dirac representation',
            'contraction':'q_mu sigma^(mu nu) q_nu = 0 for arbitrary q',
            'zero_matrix':[[str(v) for v in row] for row in out.tolist()]}


def main():
    inputs = {'NIST_CODATA_2022':ALPHA_INV_REFERENCE,
              'MTFT_three_term':GAUGE.alpha_inv,
              'MTFT_four_term':GAUGE.alpha_inv_4term,
              'MTFT_Monster':GAUGE.alpha_inv_monster}
    alpha_ref = 1/ALPHA_INV_REFERENCE
    reference_energy = -.5*ELECTRON_REST_ENERGY_EV*alpha_ref**2
    rows = []
    # Zero photon-regulator-mass limit; the electron remains massive:
    # alpha/(2pi) int_0^1 2*x*(1-x)^2 / (1-x)^2 dx.
    # quad does not evaluate at the endpoints; its error is an estimate.
    factor, estimated_error = quad(lambda x: 2*x*(1-x)**2/(1-x)**2,
                                   0.,1.,epsabs=1e-13,epsrel=1e-13)
    assert abs(factor-1) < 1e-13
    for name,inverse in inputs.items():
        alpha = 1/inverse
        E1 = -.5*ELECTRON_REST_ENERGY_EV*alpha**2
        anomaly = alpha*factor/(2*math.pi)
        rows.append({'coupling':name,'alpha_inv':inverse,'alpha':alpha,
                     'e_natural_units':math.sqrt(4*math.pi*alpha),
                     'inverse_alpha_offset_in_reference_sigma':
                         (inverse-ALPHA_INV_REFERENCE)/ALPHA_INV_REFERENCE_UNCERTAINTY,
                     'alpha_relative_offset':alpha/alpha_ref-1,
                     'coulomb_E1_infinite_nuclear_mass_eV':E1,
                     'coulomb_E1_relative_shift':E1/reference_energy-1,
                     'coulomb_E1_over_electron_rest_energy':-.5*alpha**2,
                     'one_loop_anomaly_only':anomaly,
                     'g_magnitude_one_loop_only':2*(1+anomaly)})
    alpha = GAUGE.alpha
    a_qed = alpha*factor/(2*math.pi)
    magnetic = []
    for kappa in [-.01,0.,.01]:
        magnetic.append({'kappa_P':kappa,'F1_at_zero':1.,'F2_at_zero':a_qed+kappa,
                         'g_magnitude_to_retained_order':2*(1+a_qed+kappa),
                         'spin_splitting_over_Bohr_magneton_times_B':2*(1+a_qed+kappa),
                         'mass_anchor_eV':ELECTRON_REST_ENERGY_EV,
                         'electron_charge_in_units_of_e':-1,
                         'leading_coulomb_E1_eV':-.5*ELECTRON_REST_ENERGY_EV*alpha**2})
    gamma = gamma_matrices()
    eta = [1.,-1.,-1.,-1.]
    clifford = max(np.linalg.norm(gamma[mu]@gamma[nu]+gamma[nu]@gamma[mu]
                    -(2*eta[mu]*np.eye(4) if mu==nu else np.zeros((4,4))))
                   for mu in range(4) for nu in range(4))
    residual = 0.
    rng = np.random.default_rng(2026090730)
    for _ in range(32):
        qcov = rng.normal(size=4)*.2
        contraction = np.zeros((4,4),complex)
        for mu in range(4):
            for nu in range(4):
                sigma = .5j*(gamma[mu]@gamma[nu]-gamma[nu]@gamma[mu])
                contraction += qcov[mu]*(.01j/2)*sigma*qcov[nu]  # m=1 units
        residual = max(residual,float(np.linalg.norm(contraction)))
    assert clifford == 0
    assert residual < 1e-14
    catalog_mass = LEPTONS.e*1e9
    result = {'status':'CONDITIONAL QED baseline plus a deliberately free low-energy Pauli coefficient',
              'reference_url':REFERENCE_URL,
              'electron_rest_energy_eV':ELECTRON_REST_ENERGY_EV,
              'reference_inverse_alpha_uncertainty':ALPHA_INV_REFERENCE_UNCERTAINTY,
              'sigma_scope':'Reference offsets only; no MTFT theory-error distribution or joint likelihood supplied.',
              'alpha_variants':rows,
              'loop_parameter_integral':factor,'quadrature_estimated_error':estimated_error,
              'pauli_cases':magnetic,'exact_pauli_identity':exact_pauli_transversality(),
              'gamma_Clifford_residual':float(clifford),'pauli_transversality_residual':residual,
              'pauli_random_momenta':32,
              'catalog_mass_resolution':{'mtft_catalog_electron_rest_energy_eV':catalog_mass,
                                         'reference_electron_rest_energy_eV':ELECTRON_REST_ENERGY_EV,
                                         'relative_difference':catalog_mass/ELECTRON_REST_ENERGY_EV-1},
              'excluded':'Higher QED loops, recoil, relativistic/nuclear effects, Pauli-loop terms, and UV matching.',
              'independence':'Coulomb and magnetic numbers obtained by reusing alpha are correlated consequences, not independent alpha predictions.'}
    dump('magnetic_coupling_results.json',result)
    print('Pauli vertex is exactly transverse; numerical residual',residual)
    print('MTFT default: alpha^-1 =',GAUGE.alpha_inv,'E1 =',rows[1]['coulomb_E1_infinite_nuclear_mass_eV'],'eV')
    print('Same mass and charge; |g| values:',[row['g_magnitude_to_retained_order'] for row in magnetic])


if __name__ == '__main__':
    main()
