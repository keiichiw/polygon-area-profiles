#!/usr/bin/env python3
"""Exact proportionality tests for region areas in Q(zeta_(2n))."""
import argparse
import json
import sys
from cyclotomic import Field, profile


def run(n_max):
    count, hits = 0, []
    pairs = ((0,1),(0,2),(1,2))
    for n in range(5,n_max+1):
        f = Field(2*n)
        for m in range(2,(n+1)//2):
            for k in range(1,m):
                count += 1
                areas = profile(n,k,m,f)
                for i,j in pairs:
                    q = areas[i].rational_ratio(areas[j])
                    if q is not None:
                        hits.append({'n':n,'k':k,'m':m,'ratio':f'a{i}/a{j}','value':str(q)})
        print(f'Exact ratios: n={n}, configurations={count}',file=sys.stderr,flush=True)
    return {'method':'exact coefficient-vector proportionality in Q(zeta_(2n))',
            'n_min':5,'n_max':n_max,'configurations':count,'rational_quotients':hits}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--n-max',type=int,default=100)
    args = parser.parse_args()
    if args.n_max < 5:
        parser.error('--n-max must be at least 5')
    print(json.dumps(run(args.n_max),indent=2))
