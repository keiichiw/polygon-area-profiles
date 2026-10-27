#!/usr/bin/env python3
"""Compare rotational quotient with independent full expansion, including areas."""
import argparse
from collections import Counter
import json
from pathlib import Path
import sys

from explore_cells import build_arrangement, enumerate_faces, digest
from explore_cell_orbits import calculate


def compare(n):
    orbit,_ = calculate(n)
    f,points,chords = build_arrangement(n)
    histogram = Counter()
    reference = enumerate_faces(f,points,chords,60,histogram)
    for key in ('vertices','edges','cells','rational_cells','rational_area_counts',
                'exact_area_sum','euler_characteristic'):
        if orbit[key]!=reference[key]:
            raise ArithmeticError(f'n={n}, mismatch: {key}')
    checksum = digest(sorted(histogram.items()))
    if orbit['exact_area_multiset_sha256']!=checksum:
        raise ArithmeticError(f'n={n}, exact area multiset differs')
    return {'n':n,'exact_area_multiset_sha256':checksum,'reference_report':reference}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--n',type=int,nargs='+',default=list(range(3,31))+[36,40])
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    if any(n<3 for n in args.n) or len(set(args.n))!=len(args.n):
        parser.error('--n requires distinct integers at least 3')
    results=[]
    for n in args.n:
        results.append(compare(n))
        print(f'PASS: full expansion and rotational quotient agree, n={n}',file=sys.stderr,flush=True)
    text=json.dumps(results,indent=2)+'\n'
    if args.output:
        args.output.write_text(text)
    else:
        print(text,end='')


if __name__=='__main__':
    main()
