"""Read-only MTFT 0.26.0 compatibility audit; records findings without editing MTFT.

Run with MTFT 0.26.0 importable. This is a dependency/consistency audit, not
a statistical fit to current particle data or a proof of Standard Model matching.
"""
import hashlib
import inspect
import json
from pathlib import Path
import platform
from unittest.mock import patch
import numpy as np
import mtft
from mtft import constants, hosotani, lattice, particles, koide, falsify
from mtft.surface import frozen

ROOT = Path(__file__).resolve().parent

def snapshot(h):
    return dict(h.gauge_masses(), lambda_from_mass=h.higgs_self_coupling(),
                fermion_mass_for_kappa_1=h.fermion_mass(1.0))

def main():
    assert mtft.__version__ == '0.26.0', mtft.__version__
    base = snapshot(hosotani.HosotaniMTFT())
    twice_local_vev = snapshot(hosotani.HosotaniMTFT(v_ew=492.44))
    with patch.object(constants.SM, 'v_ew', 492.44), patch.object(constants.HIGGS, 'V_EW', 492.44):
        twice_all_vevs = snapshot(hosotani.HosotaniMTFT(v_ew=492.44))
    assert snapshot(hosotani.HosotaniMTFT()) == base  # temporary overrides restored

    # A reference-data sentinel: mass calculation should not consult constants.PDG.
    class RefuseReference:
        def __getattr__(self, name):
            raise AssertionError('Mass calculation attempted a PDG lookup: '+name)
    with patch.object(constants, 'PDG', RefuseReference()):
        without_pdg = snapshot(hosotani.HosotaniMTFT())
    assert without_pdg == base

    catalog = particles.StandardModel()
    electron = catalog.by_name('Electron')
    theta, radius = base['theta_0'], base['R_tau_GeV_inv']
    kappa = electron.compute_kappa(theta, radius)
    reconstructed = kappa*abs(np.sin(theta))/radius
    assert abs(reconstructed-electron.mass_GeV) < 1e-15
    su3 = {}
    for name, fn in [('random_su_n', lattice.random_su_n), ('su_n_near_identity', lattice.su_n_near_identity)]:
        records = []
        for seed in range(8):
            U = fn(3, rng=np.random.default_rng(seed))
            records.append({'seed':seed,'determinant_error':float(abs(np.linalg.det(U)-1)),
                            'unitarity_error':float(np.linalg.norm(U.conj().T@U-np.eye(3)))})
        su3[name] = {'samples':records,
                    'determinant_one_failures_at_1e_10':sum(r['determinant_error']>1e-10 for r in records),
                    'max_determinant_error':max(r['determinant_error'] for r in records),
                    'max_unitarity_error':max(r['unitarity_error'] for r in records)}

    source_root = Path(inspect.getfile(mtft)).parent
    names = ['constants.py','hosotani.py','particles.py','koide.py','lattice.py','quantum.py','falsify.py',
             'surface/gauge.py','surface/bimodule.py','surface/ising.py','periods/physics.py']
    hashes = {name:hashlib.sha256((source_root/name).read_bytes()).hexdigest() for name in names}
    data = frozen.x0143()
    ratios = {'mass_routine_mW_over_mZ':base['m_W']/base['m_Z'],
              'separate_GAUGE_W_Z_ratio':constants.GAUGE.W_Z_ratio,
              'relative_disagreement':constants.GAUGE.W_Z_ratio/(base['m_W']/base['m_Z'])-1,
              'public_HIGGS_lambda':constants.HIGGS.lambda_quartic,
              'mass_based_lambda':base['lambda_from_mass'],
              'lambda_factor':constants.HIGGS.lambda_quartic/base['lambda_from_mass']}
    result = {
        'date':'2026-09-07','mtft_version':mtft.__version__,'python':platform.python_version(),'numpy':np.__version__,
        'scope':'Read-only source audit and bounded numerical checks, not full-suite validation or an experimental exclusion.',
        'default_mass_pipeline':base,
        'stored_package_references_not_a_current_PDG_fit':{'mW':constants.PDG.m_W,'mZ':constants.PDG.m_Z,'mH':constants.PDG.m_H},
        'mass_output_independent_of_PDG_object':True,
        'VEV_dependency_check':{'default':base,'only_instance_v_ew_doubled':twice_local_vev,
                                'instance_and_two_global_vev_anchors_doubled':twice_all_vevs,
                                'interpretation':'The public instance VEV is not a single global scale control: mW is read from constants.SM and mH is calibrated to constants.HIGGS.'},
        'internal_relations':ratios,
        'electron_catalog':{'stored_mass_GeV':electron.mass_GeV,'stored_kappa':electron.kappa,
                            'representation_string':electron.representation,
                            'inferred_kappa':float(kappa),'mass_reconstructed_from_inferred_kappa_GeV':float(reconstructed),
                            'interpretation':'Inverting an observed mass to obtain kappa and then reconstructing the mass is an identity, not an independent prediction.'},
        'koide_tau':{'result':koide.predict_tau_mass(),
                    'electron_muon_defaults_MeV':list(inspect.signature(koide.predict_tau_mass).parameters[p].default for p in ['m_e_MeV','m_mu_MeV']),
                    'interpretation':'Conditional on two measured masses, Koide=2/3, and the larger-root choice.'},
        'su3_group_membership':su3,
        'frozen_surface_gates':{k:bool(v) for k,v in data['gates'].items()},
        'source_sha256':hashes,
        'findings':[
            {'id':'SM-AUDIT-01','status':'reproduced implementation inconsistency',
             'finding':'HIGGS.lambda_quartic is 4 times HosotaniMTFT.higgs_self_coupling under the latter explicit lambda=mH^2/(2v^2) convention. The prediction test checks the lower formula directly and misses the public property.',
             'class':'EXACT factor-of-four formula comparison, numerical reproduction'},
            {'id':'SM-AUDIT-02','status':'reproduced implementation defect',
             'finding':'random_su_n normalizes one row by det(Q)^(1/N), leaving determinant det(Q)^(1-1/N). All eight N=3 samples fail determinant one while remaining unitary.',
             'class':'EXACT determinant identity, finite numerical reproduction'},
            {'id':'SM-AUDIT-03','status':'requires a common physical prescription',
             'finding':'The two MW/MZ formulas differ. A running-angle versus pole-mass interpretation would require an explicit scale, scheme, and radiative matching calculation; none is supplied by these evaluated routines.',
             'class':'DIAGNOSTIC consistency issue, not a statistical exclusion'},
            {'id':'SM-AUDIT-04','status':'dependency distinction established',
             'finding':'HosotaniMTFT reuses arithmetic mass formulas and calibrates potential coefficients. The separate particle catalog supplies observed masses, spin/representation labels, and kappas. These layers must not be counted as independent derivations.',
             'class':'source inspection plus dependency checks'}
        ]}
    (ROOT/'audit_results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    (ROOT/'package_honest_report.txt').write_text(falsify.honest_report()+'\n')
    print(json.dumps({'default_mass_pipeline':base,'internal_relations':ratios,
                      'random_SU3_failures':su3['random_su_n']['determinant_one_failures_at_1e_10'],
                      'near_identity_SU3_failures':su3['su_n_near_identity']['determinant_one_failures_at_1e_10'],
                      'frozen_gates_passed':sum(bool(v) for v in data['gates'].values())},indent=2))

if __name__ == '__main__':
    main()
