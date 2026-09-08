"""Numerical group-membership and local gauge-action checks of the patch."""
import numpy as np
from mtft.constants import HIGGS
from mtft.hosotani import HosotaniMTFT
from mtft.lattice import (LatticeConfig, MTFTAction, random_su_n,
                         random_su_n_RETIRED_SM_AUDIT_02)
from common import dump


def main():
    rows = []
    for N in [2,3,5]:
        for seed in range(16):
            U = random_su_n(N,np.random.default_rng(seed))
            old = random_su_n_RETIRED_SM_AUDIT_02(N,np.random.default_rng(seed))
            rows.append({'N':N,'seed':seed,'unitarity_error':float(np.linalg.norm(U.conj().T@U-np.eye(N))),
                         'determinant_error':float(abs(np.linalg.det(U)-1)),
                         'retired_determinant_error':float(abs(np.linalg.det(old)-1))})
    rng = np.random.default_rng(20260907)
    cfg = LatticeConfig(N=3,L=2); cfg.hot_start(rng)
    transformed = LatticeConfig(N=3,L=2)
    G = np.empty((2,2,2,2,3,3),complex)
    for site in np.ndindex(2,2,2,2):
        G[site] = random_su_n(3,rng)
    for site in np.ndindex(2,2,2,2):
        for mu in range(4):
            neighbor = cfg._shift(site,mu)
            transformed.set_link(site,mu,G[site]@cfg.get_link(site,mu)@G[neighbor].conj().T)
    action = MTFTAction(beta=6.,kappa=1.,y=.18174,n_max=12)
    action_results = {}
    for name in ['wilson_action','arithmetic_action','total_action']:
        f = getattr(action,name); a,b = f(cfg),f(transformed)
        action_results[name] = {'before':float(a),'after':float(b),'absolute_difference':float(abs(a-b)),
                                'relative_difference':float(abs(a-b)/max(1.,abs(a)))}
    data = {'quartic_current':HIGGS.lambda_quartic,
            'quartic_retired':HIGGS.lambda_quartic_RETIRED_SM_AUDIT_01,
            'quartic_hosotani':HosotaniMTFT().higgs_self_coupling(),
            'mass_relation_residual':float(abs(2*HIGGS.lambda_quartic*HIGGS.V_EW**2-HIGGS.m_H**2)),
            'group_samples':rows,'local_gauge_action':action_results,
            'scope':'One 2^4 SU(3) hot lattice; finite n_max=12 action invariance, not continuum physics.'}
    assert max(row['determinant_error'] for row in rows) < 4e-14
    assert max(row['unitarity_error'] for row in rows) < 4e-14
    assert max(row['relative_difference'] for row in action_results.values()) < 1e-12
    dump('correction_results.json',data)
    print('Corrections: 48 group draws and local SU(3) action invariance passed.')


if __name__ == '__main__':
    main()
