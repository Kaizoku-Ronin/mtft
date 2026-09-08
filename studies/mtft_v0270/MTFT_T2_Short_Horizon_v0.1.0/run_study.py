#!/usr/bin/env python3
"""Run exact certificate, sampled diagnostics, and independent expm controls."""
import hashlib
import json
from pathlib import Path
from fractions import Fraction as F
import numpy as np
from scipy.linalg import expm
from dyadic_bounds import encode, certify, rational

ROOT = Path(__file__).resolve().parent


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    frozen = np.load(ROOT/'inputs/frozen_models.npz')
    h, hr = frozen['H'], frozen['Hr']
    parent = np.load(ROOT/'inputs/T2_blocks.npz')
    q = np.load(ROOT/'inputs/selected_models.npz')['Q_svd_r4']
    bq = parent['B']@q
    dr = q.conj().T@parent['D']@q; dr = (dr+dr.conj().T)/2
    expected_h = np.block([[parent['A'],parent['B']],[parent['B'].conj().T,parent['D']]])
    expected_hr = np.block([[parent['A'],bq],[bq.conj().T,dr]])
    # Different BLAS builds may differ in final rounding: the stored matrices
    # remain the certificate target; log the numerical reconstruction error.
    reconstruction_error = max(float(np.linalg.norm(h-expected_h)),float(np.linalg.norm(hr-expected_hr)))
    if reconstruction_error > 1e-12:
        raise RuntimeError('Parent-model reconstruction mismatch')
    encoded = {'schema': 1, 'H': encode(h), 'Hr': encode(hr),
               'interpretation': 'exact dyadic block-coordinate finite matrices',
               'npz_sha256': sha(ROOT/'inputs/frozen_models.npz')}
    (ROOT/'FROZEN_MODELS.json').write_text(json.dumps(encoded, indent=2)+'\n')
    horizons = [F(1,8), F(1,4), F(1,2), F(3,4), F(1)]
    orders = [4,6,8,10,12]
    cert = certify(encoded['H'], encoded['Hr'], 8, horizons, orders, F(1,100))
    cert['provenance'] = {str(p.relative_to(ROOT)): sha(p) for p in
                          [ROOT/'PROTOCOL.md', ROOT/'FROZEN_MODELS.json',
                           ROOT/'dyadic_bounds.py', Path(__file__)]}
    (ROOT/'BOUND_CERTIFICATE.json').write_text(json.dumps(cert, indent=2)+'\n')
    times = np.unique(np.r_[np.linspace(0,1,1001), [float(t) for t in horizons]])
    def prop(mat):
        w, v = np.linalg.eigh(mat); v = v[:8]
        return np.einsum('ai,ti,bi->tab',v,np.exp(-1j*times[:,None]*w),v.conj())
    up, ur = prop(h), prop(hr)
    err = np.linalg.norm(up-ur,axis=(1,2))/np.linalg.norm(up,axis=(1,2))
    controls = []
    for t in [0,.125,.25,.5,.75,1]:
        k = int(np.argmin(abs(times-t)))
        direct = expm(-1j*t*h)[:8,:8]; direct_r = expm(-1j*t*hr)[:8,:8]
        controls.append({'t':t, 'full_expm_error':float(np.linalg.norm(direct-up[k])),
                         'reduced_expm_error':float(np.linalg.norm(direct_r-ur[k]))})
    table = []
    for t in horizons:
        entries = [r for r in cert['rows'] if rational(r['horizon']) == t]
        best = min(entries,key=lambda r:rational(r['relative_bound']))
        table.append({'horizon':float(t),'best_order':best['order'],
                      'relative_uniform_bound':float(rational(best['relative_bound'])),
                      'certifies_one_percent':best['certifies_budget'],
                      'sampled_maximum_relative_error':float(max(err[times<=float(t)]))})
    summary = {'status':'EXACT_RATIONAL_BOUND_ON_FROZEN_FINITE_MODEL',
               'physical_time_identified':False, 'exact_period_geometry_certified':False,
               'sample_count':len(times),'horizons':table,'expm_controls':controls,
               'parent_model_reconstruction_error':reconstruction_error,
               'minimum_sampled_full_norm':float(min(np.linalg.norm(up,axis=(1,2)))),
               'checks':{'exact_hermitian_full':bool(np.array_equal(h,h.conj().T)),
                         'exact_hermitian_reduced':bool(np.array_equal(hr,hr.conj().T)),
                         'moments_zero_and_one_exact':all(rational(x)==0 for x in cert['difference_norm_uppers'][:2]),
                         'expm_agreement':max(max(c['full_expm_error'],c['reduced_expm_error']) for c in controls)<1e-12,
                         'all_sampled_errors_below_corresponding_bounds':all(x['sampled_maximum_relative_error']<=x['relative_uniform_bound'] for x in table)}}
    np.savez_compressed(ROOT/'short_horizon_pointwise.npz',times=times,relative_error=err,
                        full_norm=np.linalg.norm(up,axis=(1,2)))
    (ROOT/'short_horizon_results.json').write_text(json.dumps(summary,indent=2)+'\n')
    if not all(summary['checks'].values()):
        raise RuntimeError('Study validation failed')
    print(json.dumps(summary,indent=2))


if __name__ == '__main__':
    main()
