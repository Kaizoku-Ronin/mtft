"""Explore the precomputed MTFT Ising companion data, or draw exact samples.

Examples:
  python playground.py field --beta .64 --eta .03
  python playground.py twist --basis 25 --beta .64
  python playground.py defect --edge 4
  python playground.py bilayer --N 35 --beta .64 --kappa .1
  python playground.py sample --a 4 --d 1 --count 5 --seed 143
  python playground.py zeros
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
from ising_extensions import IntegerSampler, joint_observables
from ising_exact import thermodynamics

ROOT=Path(__file__).resolve().parent


def read(name):
    return json.loads((ROOT/name).read_text())


def bilayer_observables(row, beta, kappa):
    n,E=row['spins_per_layer'],row['E_per_layer']
    D=np.asarray(row['counts'],dtype=float)
    k,j=np.nonzero(D)
    S,Q=2*E-2*k,n-2*j
    logw=np.log(D[k,j])+beta*S+kappa*Q
    shift=float(logw.max());w=np.exp(logw-shift);p=w/w.sum()
    meanQ=float(p@Q)
    return {'N':row['N'],'spins_per_layer':n,'beta':float(beta),'kappa':float(kappa),
            'log_Z':float(shift+np.log(w.sum())),
            'overlap_per_site':meanQ/n,
            'overlap_response_to_kappa_per_site':float(p@((Q-meanQ)**2)/n),
            'intralayer_interaction_energy_per_spin':float(-(p@S)/(2*n)),
            'evaluation':'Numerical finite sums from exact integer counts.'}


def main():
    parser=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    sub=parser.add_subparsers(dest='command',required=True)
    field=sub.add_parser('field',help='Field response and Fisher geometry at N=143')
    field.add_argument('--beta',type=float,default=.64)
    field.add_argument('--eta',type=float,default=0.)
    twist=sub.add_parser('twist',help='A recorded flat sign twist on N=143')
    twist.add_argument('--basis',type=int,choices=range(26),default=0)
    twist.add_argument('--beta',type=float,default=.64)
    defect=sub.add_parser('defect',help='Exact ground-state effect of reversing one AF bond')
    defect.add_argument('--edge',type=int,choices=range(84),default=4)
    layer=sub.add_parser('bilayer',help='Two coupled spin fields on a small graph')
    layer.add_argument('--N',type=int,choices=(6,11,15,35),default=35)
    layer.add_argument('--beta',type=float,default=.64)
    layer.add_argument('--kappa',type=float,default=.1)
    sample=sub.add_parser('sample',help='Independent configurations using exact integer conditional weights')
    for key,default in [('a',4),('d',1),('w0',1),('w1',1),('count',5),('seed',143)]:
        sample.add_argument('--'+key,type=int,default=default)
    sub.add_parser('zeros',help='Summary of the finite-polynomial zero calculations')
    args=parser.parse_args()
    for name in ('beta','eta','kappa'):
        if hasattr(args,name) and not np.isfinite(getattr(args,name)):
            parser.error(name+' must be finite')
    if args.command=='field':
        r=read('joint_143.json')
        out=joint_observables(r['counts'],r['spins'],r['E'],args.beta,args.eta,higher=True)
        out['convention']='Natural coordinates (beta, eta); full Fisher metric; Gaussian curvature K, scalar R=2K.'
        out['evaluation']='Numerical finite sums from exact counts. Very large parameters may underflow; null curvature means a determinant guard was reached.'
    elif args.command=='twist':
        if args.beta<0:parser.error('twist free-energy display expects beta >= 0')
        r=next(r for r in read('surface_twists.json')['levels'] if r['N']==143)
        t=r['basis_twists'][args.basis]
        z0=thermodynamics(r['baseline_counts'],r['spins'],r['E'],args.beta)['log_Z']
        z=thermodynamics(t['counts'],r['spins'],r['E'],args.beta)['log_Z']
        out={k:t[k] for k in ('basis_index','negative_edge_indices','min_frustrated','ground_degeneracy')}
        out.update(beta=args.beta,dimensionless_twist_cost=z0-z,
                   free_energy_cost_J1=(z0-z)/args.beta if args.beta else 0.,
                   convention=r['basis_warning'])
    elif args.command=='defect':
        r=read('defect_census.json')
        out={k:v for k,v in r['rows'][args.edge].items() if k!='counts'}
        out['baseline_spin_bits']=r['baseline_ground']['spin_bits']
        out['distance_definition']=r['distance_definition']
    elif args.command=='bilayer':
        r=next(r for r in read('bilayer.json')['levels'] if r['N']==args.N)
        out=bilayer_observables(r,args.beta,args.kappa)
    elif args.command=='sample':
        if min(args.a,args.d,args.w0,args.w1)<=0 or not 1<=args.count<=100000:
            parser.error('weights must be positive; count must be between 1 and 100000')
        r=read('joint_143.json')
        sampler=IntegerSampler(r['spins'],r['edges'],args.a,args.d,(args.w0,args.w1))
        out={'beta':sampler.beta,'eta':sampler.eta,'seed':args.seed,
             'spin_convention':'0=-1, 1=+1. Independent equilibrium draws, not a time trajectory.',
             'exact_integer_partition':str(sampler.Z_integer),
             'spin_bits':sampler.sample(args.count,args.seed).tolist()}
    else:
        f=read('fisher_zeros.json')
        out={'Lee_Yang':[{'a':r['satisfied_weight'],'beta':r['beta'],'closest_angle':r['closest_angle'],
                          'gates':r['gates']} for r in read('lee_yang.json')['rows']],
             'Fisher':{k:f[k] for k in ('class','factorization','max_normalized_residual','max_40_vs_65_digit_root_shift')}}
    print(json.dumps(out,indent=2,allow_nan=False))


if __name__=='__main__':main()
