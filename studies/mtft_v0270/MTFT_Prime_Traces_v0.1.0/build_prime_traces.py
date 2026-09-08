"""Exact finite prime traces, grouped by mathematical infinitude status."""
import argparse
import json
from math import isqrt
from pathlib import Path
import platform
import time
import numpy as np
from prime_families import RESIDUES,MARKS,LABELS,catalogue

ROOT=Path(__file__).resolve().parent


def sieve(limit):
    out=np.ones(limit+1,dtype=bool);out[:2]=False
    for d in range(2,isqrt(limit)+1):
        if out[d]:out[d*d::d]=False
    return out


def trial_prime(n):
    if n<2:return False
    return all(n%d for d in range(2,isqrt(n)+1))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--n-max',type=int,default=1000)
    args=parser.parse_args()
    if not 60<=args.n_max<=3000:parser.error('n-max must be between 60 and 3000 for this in-memory experiment')
    started=time.perf_counter();N=args.n_max;limit=(N+1)**2
    flags=sieve(2*limit+2)
    primes=np.flatnonzero(flags[:limit+1]).astype(np.int64)
    bins=np.array([isqrt(int(p)) for p in primes],dtype=np.int64)
    assert bins.min()==1 and bins.max()<=N
    properties=np.zeros((len(primes),4),dtype=bool)
    properties[:,0]=flags[primes-2]|flags[primes+2]
    properties[:,1]=flags[2*primes+1]
    properties[:,2]=(primes>2)&flags[(primes-1)//2]
    quotients=np.full(len(primes),-1,dtype=np.int64)
    for i,p0 in enumerate(primes):
        p=int(p0)
        if p==2:continue
        remainder=pow(2,p-1,p*p)
        assert (remainder-1)%p==0
        quotients[i]=(remainder-1)//p
    properties[:,3]=(primes>2)&(quotients==0)
    for j,threshold,allowed in ((0,7,[1,11,13,17,19,29]),(1,5,[11,23,29]),(2,11,[17,23,29])):
        assert np.all(np.isin(primes[(primes>threshold)&properties[:,j]]%30,allowed))
    n=np.arange(1,N+1,dtype=np.int64)
    cuts=(n+1)**2
    interval_total=np.bincount(bins,minlength=N+1)[1:]
    interval_residue=np.column_stack([np.bincount(bins[primes%30==r],minlength=N+1)[1:] for r in RESIDUES])
    interval_types=np.column_stack([np.bincount(bins[properties[:,j]],minlength=N+1)[1:] for j in range(4)])
    exceptions=np.bincount(bins[primes<=5],minlength=N+1)[1:]
    assert np.array_equal(interval_residue.sum(axis=1)+exceptions,interval_total)
    assert np.all(interval_types<=interval_total[:,None])
    assert interval_total.min()>0
    first_indices=np.searchsorted(primes,n*n,side='right')
    first_prime=primes[first_indices]
    assert np.all(first_prime<(n+1)**2)
    previous_prime=np.r_[-1,primes[:-1]]
    same_family_gap=np.full((len(primes),4),-1,dtype=np.int64)
    for j in range(4):
        ids=np.flatnonzero(properties[:,j]);same_family_gap[ids[1:],j]=np.diff(primes[ids])
    data={'p':primes,'square_n':bins,'position_numerator':primes-bins*bins,
          'position_denominator':2*bins+1,'residue_mod30':primes%30,
          'properties':properties,'fermat_quotient_base2':quotients,
          'previous_prime':previous_prime,'previous_same_family_gap':same_family_gap,
          'n':n,'upper_square':cuts,'interval_total':interval_total,
          'interval_residue':interval_residue,'interval_types':interval_types,
          'interval_exceptions':exceptions,'first_prime':first_prime,
          'cumulative_total':np.cumsum(interval_total),
          'cumulative_residue':np.cumsum(interval_residue,axis=0),
          'cumulative_types':np.cumsum(interval_types,axis=0)}
    np.savez_compressed(ROOT/'prime_traces.npz',**data)
    # Independent trial division checks the sieve and family tests, including
    # partner primality beyond the displayed square interval.
    assert all(bool(flags[k])==trial_prime(k) for k in range(5001))
    ids=np.unique(np.r_[np.arange(min(50,len(primes))),np.linspace(0,len(primes)-1,160,dtype=int)])
    for i in ids:
        p=int(primes[i]);assert trial_prime(p)
        expected=[trial_prime(p-2) or trial_prime(p+2),trial_prime(2*p+1),p>2 and trial_prime((p-1)//2)]
        assert properties[i,:3].tolist()==expected
    known_pi={10:4,100:25,1000:168,10000:1229,100000:9592,1000000:78498}
    anchors={x:k for x,k in known_pi.items() if x<=limit}
    for x,k in anchors.items():assert int(np.count_nonzero(flags[:x+1]))==k
    wieferich=primes[properties[:,3]].tolist();assert wieferich==[1093,3511]
    for p in wieferich:
        # An exact unreduced power gives a second route for the two positives.
        assert ((2**(p-1)-1)//p)%p==0
    # Exact transport: safe primes up to x correspond to Sophie Germain primes
    # up to (x-1)/2, not to SG primes with the same cutoff x.
    for x in cuts:
        safe=int(np.count_nonzero(properties[:,2]&(primes<=x)))
        sg=int(np.count_nonzero(properties[:,1]&(primes<=(int(x)-1)//2)))
        assert safe==sg
    intersections=np.zeros((8,4),dtype=np.int64)
    for i,r in enumerate(RESIDUES):
        for j in range(4):intersections[i,j]=np.count_nonzero((primes%30==r)&properties[:,j])
    examples=[]
    for p in (2,3,5,7,11,23,1091,1093,3511):
        i=int(np.searchsorted(primes,p))
        examples.append({'p':p,'square_n':int(bins[i]),'residue_mod30':int(primes[i]%30),
                         'position':[int(data['position_numerator'][i]),int(data['position_denominator'][i])],
                         'properties':dict(zip(MARKS,properties[i].tolist())),
                         'fermat_quotient_base2':None if p==2 else int(quotients[i])})
    summary={'scope':{'n_min':1,'n_max':N,'prime_upper_bound':limit,'partner_sieve_bound':2*limit+2},
             'prime_count':len(primes),'residue_order':list(RESIDUES),'property_order':list(MARKS),
             'property_labels':list(LABELS),
             'total_primes_by_residue':data['cumulative_residue'][-1].tolist(),
             'property_counts':dict(zip(MARKS,data['cumulative_types'][-1].tolist())),
             'residue_property_intersections':intersections.tolist(),'wieferich_base2':wieferich,
             'legendre_empty_intervals':(np.flatnonzero(interval_total==0)+1).tolist(),
             'examples':examples,'gates':{'all_5001_small_integers_checked_by_trial_division':True,
             'sampled_family_rows_checked_independently':len(ids),'prime_count_anchors':anchors,
             'residue_partition_conserves_all_counts':True,'Wieferich_positives_two_exact_routes':True,
             'necessary_congruence_restrictions_checked':True,
             'SG_safe_cutoff_transport_all_square_endpoints':True,'all_square_intervals_nonempty':bool(np.all(interval_total>0))},
             'arithmetic_class':'EXACT finite integer sieve, modular arithmetic, counts, and positions; no infinitude inferred from data.',
             'runtime_seconds':time.perf_counter()-started,'python':platform.python_version(),'numpy':np.__version__}
    (ROOT/'trace_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    (ROOT/'family_catalogue.json').write_text(json.dumps(catalogue(),indent=2)+'\n')
    print(json.dumps({k:summary[k] for k in ('scope','prime_count','property_counts','total_primes_by_residue','wieferich_base2','gates','runtime_seconds')},indent=2))


if __name__=='__main__':main()
