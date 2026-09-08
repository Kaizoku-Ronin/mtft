"""Read-only v0.26.2 checks: conventions, Mellin normalization, and row counts."""
from pathlib import Path
import hashlib
import json
import math
import mpmath as mp
import numpy as np
import mtft
from mtft.constants import HIGGS
from mtft.hosotani import HosotaniMTFT
from mtft.lattice import random_su_n, random_su_n_defective_v0261
from mtft.falsify import prediction_table, honest_report

ROOT = Path(__file__).resolve().parent


def main():
    assert mtft.__version__ == '0.26.2'
    conventions = {'original_lambda': HIGGS.lambda_quartic,
                   'canonical_lambda': HIGGS.lambda_quartic_pdg,
                   'hosotani_lambda': HosotaniMTFT().higgs_self_coupling(),
                   'original_mass_relation_relative_error':
                       abs(HIGGS.m_H**2-HIGGS.lambda_quartic*HIGGS.V_EW**2/2)/HIGGS.m_H**2,
                   'canonical_mass_relation_relative_error':
                       abs(HIGGS.m_H**2-2*HIGGS.lambda_quartic_pdg*HIGGS.V_EW**2)/HIGGS.m_H**2,
                   'interpretation':'Two conventions for the same potential, not two Higgs predictions.'}
    assert conventions['canonical_lambda'] == conventions['original_lambda']/4
    assert abs(conventions['canonical_lambda']-conventions['hosotani_lambda']) < 1e-14
    group = []
    for N in [2,3,5]:
        for seed in range(8):
            U = random_su_n(N,np.random.default_rng(seed))
            old = random_su_n_defective_v0261(N,np.random.default_rng(seed))
            group.append({'N':N,'seed':seed,'det_error':float(abs(np.linalg.det(U)-1)),
                          'unitarity_error':float(np.linalg.norm(U.conj().T@U-np.eye(N))),
                          'retired_det_error':float(abs(np.linalg.det(old)-1))})
    assert max(r['det_error'] for r in group) < 1e-13
    assert max(r['unitarity_error'] for r in group) < 1e-13

    # Direct finite-weight integral checks. The infinite identity follows by
    # absolute convergence and exchanging sum/integral for Re(s)>1.
    # This finite quadrature is a diagnostic, not a proof of analytic continuation.
    mellin = []
    with mp.workdps(50):
        n_max = 24
        weights = [mp.mpf(0)]*(n_max+1)
        for n in range(2,n_max+1):
            weights[n] = mp.fsum(mp.log(d)/d for d in range(2,n+1) if n%d==0)
        def ZL(y):
            return mp.fsum(weights[n]*mp.exp(-2*mp.pi*y*n) for n in range(2,n_max+1))
        for s in [mp.mpf(2),mp.mpf(3),mp.mpc('2.5','0.4')]:
            direct = mp.fsum(weights[n]*mp.power(n,-s) for n in range(2,n_max+1))
            integral = mp.quad(lambda y: mp.power(y,s-1)*ZL(y),[0,.01,.1,1,mp.inf])
            recovered = (2*mp.pi)**s/mp.gamma(s)*integral
            error = abs((recovered-direct)/direct)
            assert error < mp.mpf('1e-40')
            mellin.append({'s':str(s),'finite_n_max':n_max,'direct_Dirichlet_sum':str(direct),
                           'Mellin_recovery':str(recovered),'relative_error':str(error),
                           'quadrature_integral':str(integral)})
        # A common scalar rescaling cannot make e^-a*n and n^-beta the same
        # weight kernel at n=2,3,4: the log's second finite difference is nonzero.
        log_second_difference = 2*mp.log(3)-mp.log(2)-mp.log(4)
        kernel_obstruction = {'expression':'2 log 3 - log 2 - log 4 = log(9/8)',
                              'value':str(log_second_difference),
                              'scope':'Rules out a pointwise common exponential kernel at beta>0, not every conceivable analytic-continuation construction.'}

    predictions = prediction_table()
    groups = {}
    for p in predictions:
        key = p.group or 'row_'+str(p.number)
        groups.setdefault(key,[]).append(p.number)
    non_definition = {k:v for k,v in groups.items() if v != [18]}
    # Row17 is the squared/scaled version of the same mass constraint as row8.
    by_number = {p.number:p for p in predictions}
    from mtft.constants import PDG
    p8,p17 = by_number[8],by_number[17]
    dependence = {'predicted_lambda_from_row8':p8.predicted**2/(2*PDG.v_ew**2),
                  'row17_predicted':p17.predicted,
                  'observed_lambda_from_row8':p8.observed**2/(2*PDG.v_ew**2),
                  'row17_observed':p17.observed}
    assert abs(dependence['predicted_lambda_from_row8']-p17.predicted) < 1e-14
    assert abs(dependence['observed_lambda_from_row8']-p17.observed) < 1e-14
    accounting = {'raw_rows':len(predictions),'metadata_groups':groups,
                  'metadata_group_count':len(groups),'after_removing_definition_row18':len(non_definition),
                  'after_also_merging_row17_into_higgs_family':len(non_definition)-1,
                  'row8_row17_algebraic_dependency':dependence,
                  'additional_accounting_notes':[
                      'Rows14 and15 compare arithmetic to rational/numerical mathematical targets, not independent measured particle observables.',
                      'Row13 predicts tau conditionally on measured electron/muon input defaults; row19 uses PDG.alpha_inv as an input.',
                      'Shared measured masses, scheme choices and tolerance choices create further dependencies.',
                      'None of the counts is a statistically established number of independent tests. No joint null probability calculated.'
                  ]}

    trace3 = [d for d in range(143) if (d*d-3*d+1)%143==0]
    trace4 = [d for d in range(143) if (d*d-4*d+1)%143==0]
    d = trace4[0]; a=4-d; b=(a*d-1)//143
    assert not trace3 and a*d-143*b==1
    systole = {'trace3_residues':trace3,'trace4_residues':trace4,
               'trace4_matrix':[[a,b],[143,d]],'length':2*math.acosh(2),
               'metric':'Curvature -1 upper-half-plane quotient Gamma0(143)\\H, with cusps.',
               'quotient_area': '56*pi',
               'compact_uniformizing_area': '48*pi',
               'normalization_scope':'A positive scalar ruler changes scales; it cannot change inertia/signature or remove non-scalar anticommutator defects.'}
    source_root = Path(mtft.__file__).resolve().parent
    source_hashes = {name:hashlib.sha256((source_root/name).read_bytes()).hexdigest()
                     for name in ['constants.py','lattice.py','arithmetic_wick.py','falsify.py']}
    result = {'version':mtft.__version__,'source_hashes':source_hashes,'Higgs_conventions':conventions,
              'SU_draws':group,'mellin_finite_checks':mellin,'mellin_working_dps':50,
              'kernel_obstruction':kernel_obstruction,'prediction_accounting':accounting,'systole':systole}
    (ROOT/'release_mellin_accounting_results.json').write_text(json.dumps(result,indent=2)+'\n')
    (ROOT/'v0262_honest_report.txt').write_text(honest_report()+'\n')
    print('Version',mtft.__version__,'convention and SU checks passed.')
    print('Mellin max relative error:',max(float(r['relative_error']) for r in mellin))
    print('Accounting: raw',len(predictions),'metadata',len(groups),'without definition',len(non_definition),
          'with additional exact Higgs dependency',len(non_definition)-1)


if __name__ == '__main__':
    main()
