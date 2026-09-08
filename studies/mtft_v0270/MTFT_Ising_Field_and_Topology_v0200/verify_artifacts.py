"""Verify stored exact identities, rational root certificates, numerical gates, and CLI examples.

This verifier reads the independent C++ enumeration. Re-run that program to
regenerate its input-bound output; README.md gives the complete commands.
"""
from __future__ import annotations
from collections import Counter
import hashlib
import itertools
import json
import math
from pathlib import Path
import subprocess
import sys
import numpy as np
import sympy as sp
from ising_exact import thermodynamics
from run_extensions import signed_even_transform

ROOT=Path(__file__).resolve().parent


def read(name):return json.loads((ROOT/name).read_text())


def main():
    gates=[]
    def check(name,condition,kind='EXACT'):
        if not condition:raise AssertionError(name)
        gates.append({'name':name,'class':kind,'passed':True})
    baseline={r['N']:r for r in read('baseline_results.json')['levels']}
    joints={N:read(f'joint_{N}.json') for N in (6,11,15,35,143)}
    for N,r in joints.items():
        n,E,D=r['spins'],r['E'],r['counts']
        check(f'joint {N}: integer domain and state total',len(D)==E+1 and all(len(row)==n+1 and all(type(c) is int and c>=0 for c in row) for row in D) and sum(map(sum,D))==1<<n)
        check(f'joint {N}: energy marginal',list(map(sum,D))==baseline[N]['density_of_states'])
        check(f'joint {N}: binomial field marginal',[sum(row[m] for row in D) for m in range(n+1)]==[math.comb(n,m) for m in range(n+1)])
        check(f'joint {N}: reversal and cubic parity',all(row==row[::-1] and all(c==0 or (k-m)%2==0 for m,c in enumerate(row)) for k,row in enumerate(D)))
        check(f'joint {N}: recorded construction gates',all(r['gates'].values()),'Recorded EXACT construction checks')
    twists={r['N']:r for r in read('surface_twists.json')['levels']}
    for N,r in twists.items():
        check(f'twists {N}: basis size',len(r['basis_twists'])==2*r['genus'])
        check(f'twists {N}: each census and minimum',all(sum(t['counts'])==1<<r['spins'] and t['min_frustrated']==next(k for k,d in enumerate(t['counts']) if d) and t['min_frustrated']>0 for t in r['basis_twists']))
        check(f'twists {N}: recorded flatness, quotient, gauge checks',all(r['gates'].values()),'Recorded EXACT construction checks')
    exact=twists[143]['two_cycle_projection']['class_even_counts']
    independent=read('twisted_cycle_counts143.json')
    check('two homology coordinates: all 4 x 85 C++ coefficients',independent==exact)
    check('two homology coordinates: all 2^29 cycles accounted for',sum(map(sum,independent))==1<<29 and all(sum(row)==1<<27 for row in independent))
    transformed=signed_even_transform(baseline[143]['density_of_states'],56,84)
    check('independent cycle marginal vs exact signed transform',[sum(row[k] for row in independent) for k in range(85)]==transformed)
    tokens=list(map(int,(ROOT/'twist_cycle_input143.txt').read_text().split()))
    expected=[56,84]+list(itertools.chain.from_iterable(joints[143]['edges']))+list(itertools.chain.from_iterable(twists[143]['basis_representatives'][:2]))
    check('C++ input matches graph and two cocycles',tokens==expected)
    check('N6/N11 untwisted graph control',Counter(tuple(sorted(e)) for e in joints[6]['edges'])==Counter(tuple(sorted(e)) for e in joints[11]['edges']) and joints[6]['counts']==joints[11]['counts'])
    for r in read('bilayer.json')['levels']:
        N,n,E,B=r['N'],r['spins_per_layer'],r['E_per_layer'],r['counts'];single=baseline[N]['density_of_states']
        convolution=[sum(single[i]*single[k-i] for i in range(E+1) if 0<=k-i<=E) for k in range(2*E+1)]
        check(f'bilayer {N}: independent layers',list(map(sum,B))==convolution)
        check(f'bilayer {N}: overlap binomial',[sum(row[j] for row in B) for j in range(n+1)]==[(1<<n)*math.comb(n,j) for j in range(n+1)])
        check(f'bilayer {N}: locked layers and reversal',all(row==row[::-1] and row[0]==(single[k//2] if k%2==0 else 0) for k,row in enumerate(B)))
    defects=read('defect_census.json');ref=defects['baseline_ground']['spin_bits'];edges=joints[143]['edges']
    for r in defects['rows']:
        s=r['witness_bits'];bits=[1]*84;bits[r['edge']]=0
        cost=sum((s[u]^s[v])^b for (u,v),b in zip(edges,bits))
        check(f'defect {r["edge"]}: witness and exact ground count',cost==r['min_frustrated']==next(k for k,c in enumerate(r['counts']) if c) and sum(a!=b for a,b in zip(s,ref))==r['minimum_rearrangement'] and r['counts'][cost]==r['ground_degeneracy'])
    dc=defects['largest_rearrangement_counting_certificate']
    check('edge 4: independent ground-distance count',[ (m,c) for m,c in enumerate(dc['ground_state_counts_by_hamming_distance']) if c]==[(9,1),(47,1)] and dc['edge']==4)
    y,u,x=sp.symbols('y u x')
    for r in read('lee_yang.json')['rows']:
        a=r['satisfied_weight'];D=joints[143]['counts']
        field=sp.Poly.from_list(list(reversed([sum(D[k][m]*a**(84-k) for k in range(85)) for m in range(57)])),y)
        reduced=sp.Poly.from_list([int(c) for c in r['reduced_u_coefficients_descending']],u)
        substitution=sp.Poly(sp.expand(y**28*reduced.as_expr().subs(u,y+1/y)),y)
        check(f'Lee-Yang a={a}: exact reciprocal reduction',field.primitive()[1]==substitution.primitive()[1] and reduced.degree()==28)
        intervals=[tuple(sp.Rational(t) for t in z['u_interval']) for z in r['positive_imaginary_roots']]
        check(f'Lee-Yang a={a}: 28 disjoint sign-changing rational intervals',len(intervals)==28 and all(-2<lo<hi<2 and reduced.eval(lo)*reduced.eval(hi)<0 for lo,hi in intervals) and all(intervals[i][1]<intervals[i+1][0] for i in range(27)))
    fisher=read('fisher_zeros.json');Q=sp.Poly.from_list([int(c) for c in fisher['Q45_coefficients_descending']],x)
    D=sp.Poly.from_list(list(reversed(baseline[143]['density_of_states'])),x)
    check('Fisher exact factorization and multiplicities',D==sp.Poly(2*(1+x)**28*(1+x*x),x)*Q and Q.degree()==45 and Q.eval(-1)!=0 and Q.rem(sp.Poly(x*x+1,x))!=0)
    check('Fisher numerical residual and precision agreement',float(fisher['max_normalized_residual'])<1e-55 and float(fisher['max_40_vs_65_digit_root_shift'])<1e-30,'DIAGNOSTIC')
    geometry=read('geometry_checks.json')
    check('Fisher geometry: sign convention and three-route checks',geometry['categorical_gate']['K_determinant']=='0.25' and all(float(c['max_absolute_disagreement'])<1e-45 and c['double_absolute_error']<1e-9 for c in geometry['cross_checks']),'DIAGNOSTIC')
    grid=np.load(ROOT/'field_geometry_grid.npz')
    check('Fisher geometry: positive metric determinants and spin reversal',np.all(grid['fisher_determinant']>0) and np.max(abs(grid['magnetization_per_spin']+grid['magnetization_per_spin'][:,::-1]))<1e-10 and np.max(abs(grid['gaussian_curvature']-grid['gaussian_curvature'][:,::-1]))<1e-9,'DIAGNOSTIC')
    sampling=read('sampling.json');samples=np.load(ROOT/'equilibrium_samples.npz')
    check('Sampler: recorded all-configuration probability gate',sampling['small_graph_all_configuration_probabilities_exact'],'Recorded EXACT construction check')
    for r in sampling['rows']:
        a=r['a'];count=r['sample_count']
        check(f'Sampler a={a}: stored histograms and bit samples',sum(samples[f'q{a}_cut_histogram'])==count and sum(samples[f'q{a}_magnetization_histogram'])==count and samples[f'q{a}_spin_bits'].shape==(256,56) and np.all(samples[f'q{a}_spin_bits']<=1))
        check(f'Sampler a={a}: finite-sample mean diagnostic',max(r['S_standardized_error'],r['M_standardized_error'])<6,'DIAGNOSTIC')
    def cli(*args):
        return json.loads(subprocess.check_output([sys.executable,str(ROOT/'playground.py'),*args],cwd=ROOT,text=True))
    check('CLI field example',abs(cli('field','--beta','.64','--eta','.03')['gaussian_curvature']+.055307893621000476)<1e-10,'Integration')
    check('CLI twist example',cli('twist','--basis','25')['min_frustrated']==3,'Integration')
    check('CLI defect example',cli('defect','--edge','4')['minimum_rearrangement']==9,'Integration')
    bil=cli('bilayer','--N','35','--beta','.64','--kappa','0')
    single=thermodynamics(baseline[35]['density_of_states'],16,24,.64)['log_Z']
    check('CLI bilayer factorization',abs(bil['log_Z']-2*single)<1e-12 and abs(bil['overlap_per_site'])<1e-12,'Integration')
    check('CLI sampling example',len(cli('sample','--count','3')['spin_bits'])==3,'Integration')
    check('CLI zeros example',len(cli('zeros')['Lee_Yang'])==3,'Integration')
    source=ROOT/'mtft-0.26.0.tar.gz'
    if source.exists():
        check('Source archive SHA256',hashlib.sha256(source.read_bytes()).hexdigest()=='758f29a7a0c2ca08896a15964b3a607c7143ac1fd335bcc8f4d84a4880e53cdf')
    result={'all_passed':True,'number_of_gates':len(gates),'class_totals':dict(Counter(g['class'] for g in gates)),
            'independent_cycle_enumeration':{'cycles':1<<29,'coefficient_checks':4*85,
                'log':(ROOT/'twisted_cycle_enumeration.log').read_text().strip()},'gates':gates}
    (ROOT/'verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='gates'},indent=2))


if __name__=='__main__':main()
