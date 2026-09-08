"""Inspect a prime's splitting/regularity trace or the symmetry census."""
import argparse
import json
from common import read

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--prime',type=int)
    group.add_argument('--interval',type=int)
    group.add_argument('--sectors',action='store_true')
    args=parser.parse_args()
    if args.sectors:
        c=read('sector_certificate.json')
        result={k:c[k] for k in ['block_order','real_dimensions','AL_sign_order','AL_intersection_real_dimensions','AL_allowed_direct_coupling']}
    else:
        rows=read('prime_atlas.json')
        if args.prime is not None:
            matches=[r for r in rows if r['p']==args.prime]
            if not matches:
                parser.error('The requested integer is not a prime in the census 2 <= p < 10201.')
            result=matches[0]
        else:
            if not 1<=args.interval<=100:
                parser.error('The interval index must be between 1 and 100.')
            result=[r for r in rows if r['square_n']==args.interval]
    print(json.dumps(result,indent=2))

if __name__=='__main__':
    main()
