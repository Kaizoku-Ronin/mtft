"""Inspect a prime or a square interval with family infinitude labels."""
import argparse
import json
from pathlib import Path
import numpy as np
from prime_families import RESIDUES,MARKS,catalogue

ROOT=Path(__file__).resolve().parent


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    pick=parser.add_mutually_exclusive_group(required=True)
    pick.add_argument('--prime',type=int)
    pick.add_argument('--interval',type=int,help='n for (n^2,(n+1)^2)')
    pick.add_argument('--catalogue',action='store_true')
    args=parser.parse_args()
    if args.catalogue:print(json.dumps(catalogue(),indent=2));return
    data=np.load(ROOT/'prime_traces.npz');p=data['p']
    if args.interval is not None:
        n=args.interval
        if not 1<=n<=len(data['n']):parser.error('interval outside computed range')
        ids=np.flatnonzero(data['square_n']==n)
        out={'n':n,'open_interval':[n*n,(n+1)**2],'prime_count':len(ids),
             'primes':p[ids].tolist(),'proved_infinite_residue_family_counts':dict(zip(map(str,RESIDUES),data['interval_residue'][n-1].tolist())),
             'open_infinitude_family_counts':dict(zip(MARKS,data['interval_types'][n-1].tolist())),
             'membership_arithmetic':'EXACT finite calculation'}
    else:
        value=args.prime;i=int(np.searchsorted(p,value))
        if i==len(p) or p[i]!=value:parser.error('not a prime in the computed range')
        r=int(data['residue_mod30'][i]);q=int(data['fermat_quotient_base2'][i])
        out={'p':value,'square_n':int(data['square_n'][i]),
             'position_fraction':[int(data['position_numerator'][i]),int(data['position_denominator'][i])],
             'proved_infinite_families':['all_primes']+([f'residue30_{r}'] if value>5 else []),
             'finite_mod30_exception':value<=5,
             'open_infinitude_properties':dict(zip(MARKS,data['properties'][i].tolist())),
             'fermat_quotient_base2':None if q<0 else q,
             'previous_same_family_gap':{name:(None if int(g)<0 else int(g)) for name,g in zip(MARKS,data['previous_same_family_gap'][i])},
             'membership_arithmetic':'EXACT; OPEN describes infinitude, not uncertainty of these membership tests'}
    print(json.dumps(out,indent=2))


if __name__=='__main__':main()
