"""Check CC-21/CC-22 without changing either public quartic convention.

This is a finite numerical release-contract check, not a physics validation.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    source = args.source.resolve()
    if not (source / 'src/mtft/__init__.py').is_file():
        raise SystemExit('Expected an unpacked MTFT source tree with src/mtft/__init__.py')
    sys.path.insert(0, str(source / 'src'))
    import numpy as np
    import mtft
    from mtft.constants import HIGGS
    from mtft.hosotani import HosotaniMTFT
    from mtft.lattice import LatticeConfig, MTFTAction, random_su_n, random_su_n_defective_v0261

    imported = Path(mtft.__file__).resolve()
    assert imported.is_relative_to(source / 'src')
    quartic = float(HIGGS.lambda_quartic)
    pdg = float(HIGGS.lambda_quartic_pdg)
    hosotani = float(HosotaniMTFT().higgs_self_coupling())
    scale = float(HIGGS.m_H ** 2)
    conventions = {
        'lambda_quartic_original_convention': quartic,
        'lambda_quartic_pdg': pdg,
        'lambda_hosotani_at_default_vev': hosotani,
        'conversion_absolute_residual': abs(pdg - quartic / 4),
        'original_mass_relation_relative_residual': abs(quartic * HIGGS.V_EW**2 / 2 - scale) / scale,
        'pdg_mass_relation_relative_residual': abs(2 * pdg * HIGGS.V_EW**2 - scale) / scale,
        'hosotani_absolute_residual': abs(pdg - hosotani),
    }
    draws = []
    for n in (2, 3, 5):
        for seed in range(16):
            u = random_su_n(n, np.random.default_rng(seed))
            draws.append({'N': n, 'seed': seed,
                          'unitarity_frobenius_residual': float(np.linalg.norm(u.conj().T @ u - np.eye(n))),
                          'determinant_absolute_residual': float(abs(np.linalg.det(u) - 1))})
    retired = random_su_n_defective_v0261(3, np.random.default_rng(0))
    retired_defect = float(abs(np.linalg.det(retired) - 1))
    rng = np.random.default_rng(20260907)
    original = LatticeConfig(N=3, L=2)
    original.hot_start(rng)
    transformed = LatticeConfig(N=3, L=2)
    gauges = np.empty((2, 2, 2, 2, 3, 3), complex)
    for site in np.ndindex(2, 2, 2, 2):
        gauges[site] = random_su_n(3, rng)
    for site in np.ndindex(2, 2, 2, 2):
        for mu in range(4):
            neighbor = original._shift(site, mu)
            transformed.set_link(site, mu, gauges[site] @ original.get_link(site, mu) @ gauges[neighbor].conj().T)
    action = MTFTAction(beta=6., kappa=1., y=.18174, n_max=12)
    action_rows = {}
    for name in ('wilson_action', 'arithmetic_action', 'total_action'):
        function = getattr(action, name)
        before, after = float(function(original)), float(function(transformed))
        action_rows[name] = {'before': before, 'after': after,
                             'relative_residual': abs(before - after) / max(1., abs(before))}
    gates = {
        'quartic_conversion': conventions['conversion_absolute_residual'] < 1e-14,
        'original_mass_convention': conventions['original_mass_relation_relative_residual'] < 1e-12,
        'pdg_mass_convention': conventions['pdg_mass_relation_relative_residual'] < 1e-12,
        'hosotani_default_convention': conventions['hosotani_absolute_residual'] < 1e-12,
        'all_48_SU_N_draws': max(row['determinant_absolute_residual'] for row in draws) < 4e-14
                           and max(row['unitarity_frobenius_residual'] for row in draws) < 4e-14,
        'historical_defect_reproduced': retired_defect > 1e-3,
        'finite_local_action_covariance': max(row['relative_residual'] for row in action_rows.values()) < 1e-12,
    }
    result = {'schema_version': 1, 'version': mtft.__version__, 'imported_from': str(imported),
              'scope': 'Finite CC-21/CC-22 API checks plus one periodic 2^4 SU(3) action covariance check; no continuum or Standard Model claim.',
              'conventions': conventions, 'group_draws': draws, 'historical_determinant_defect': retired_defect,
              'local_gauge_action': action_rows, 'gates': gates, 'all_pass': all(gates.values())}
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    print(json.dumps({'version': mtft.__version__, 'all_pass': result['all_pass'], 'gates': len(gates)}))
    if not result['all_pass']:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
