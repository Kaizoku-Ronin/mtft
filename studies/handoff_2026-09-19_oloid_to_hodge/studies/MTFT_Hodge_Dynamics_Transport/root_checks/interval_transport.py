#!/usr/bin/env python3
"""Chosen control: a chaotic base with a flat, bounded coboundary cocycle.

This does not select an MTFT parameter law. It tests the exact telescoping
obstruction on the compact nonsingular interval c in [-3,-2]. Float results
are diagnostics; the report supplies the proof for arbitrary base itineraries.
"""
from pathlib import Path
import csv
import json
import math
import numpy as np
from scipy.integrate import solve_ivp

OUT = Path(__file__).resolve().parent

def connection(c):
    return np.array([[-1/(4*c), 1/(4*c*(c+1))], [-.25, 1/(4*c)]])

def fundamental(rtol):
    def rhs(c, flattened):
        return (connection(c) @ flattened.reshape(2,2)).ravel()
    options = dict(method='DOP853', rtol=rtol, atol=rtol/100,
                   dense_output=True, max_step=.025)
    lower = solve_ivp(rhs, (-2.5,-3), np.eye(2).ravel(), **options)
    upper = solve_ivp(rhs, (-2.5,-2), np.eye(2).ravel(), **options)
    assert lower.success and upper.success
    def y(c):
        return (lower if c < -2.5 else upper).sol(c).reshape(2,2)
    return y

def experiment(rtol):
    y = fundamental(rtol)
    rows = []
    seeds = [.123456789, math.sqrt(2)-1, math.pi/10]
    for seed in seeds:
        x = seed
        c0 = x-3
        y0inv = np.linalg.inv(y(c0))
        product = np.eye(2)
        base_log = 0.
        for n in range(1,8193):
            derivative = abs(4*(1-2*x))
            assert derivative > 0
            base_log += math.log(derivative)
            xn = 4*x*(1-x)
            c, cn = x-3, xn-3
            transfer = y(cn) @ np.linalg.inv(y(c))
            product = transfer @ product
            x = xn
            if n in [256,1024,4096,8192]:
                exact_form = y(cn) @ y0inv
                smax = np.linalg.svd(exact_form, compute_uv=False)[0]
                rows.append(dict(seed=seed, steps=n, solver_rtol=rtol,
                    base_finite_exponent=base_log/n,
                    transport_finite_exponent=math.log(smax)/n,
                    endpoint_transport_norm=smax,
                    product_telescope_defect=float(np.linalg.norm(product-exact_form)),
                    determinant_defect=abs(float(np.linalg.det(exact_form))-1)))
    grid = np.linspace(-3,-2,1001)
    ymax = max(np.linalg.norm(y(c),ord=2) for c in grid)
    invmax = max(np.linalg.norm(np.linalg.inv(y(c)),ord=2) for c in grid)
    return rows, dict(sampled_Y_norm_max=ymax, sampled_inverse_norm_max=invmax,
                      sampled_uniform_bound=ymax*invmax,
                      qualification='sampled diagnostic, not certified supremum')

if __name__ == '__main__':
    coarse, _ = experiment(1e-9)
    fine, bounds = experiment(2e-13)
    differences = {key:max(abs(a[key]-b[key]) for a,b in zip(coarse,fine))
        for key in ['transport_finite_exponent','endpoint_transport_norm']}
    with (OUT/'interval_results.csv').open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(fine[0]))
        writer.writeheader(); writer.writerows(fine)
    result={'status':'DIAGNOSTIC', 'base':'c_next=-3+4(c+3)(-c-2)',
        'domain':[-3,-2], 'transport':'Y(c_next) inverse(Y(c))',
        'exact_theorem':'Products telescope and all cocycle exponents vanish on a compact nonsingular interval.',
        'base_reference_exponent':math.log(2), 'tolerance_comparison':differences,
        'sampled_bounds':bounds, 'rows':fine}
    (OUT/'interval_results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'final_rows':[r for r in fine if r['steps']==8192],
        'tolerance_comparison':differences,'sampled_bounds':bounds},indent=2))
