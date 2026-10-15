#!/usr/bin/env python3
"""Exact PR Table 3 matching and all 65 Appendix d_j values (stdlib only)."""
import csv
import json
from fractions import Fraction as F
from itertools import permutations
from math import lcm
from pathlib import Path

HERE = Path(__file__).resolve().parent


def target(n, b, c):
    return (F(b - 2, 2*n), F(b, 2*n), F(n - c - 2, 2*n),
            F(n - b + 2, 2*n), F(1, n), F(c - b, 2*n))


def admissible(n, b, c):
    return (n >= 5 and 2 < b < c < n - 2 and c - b >= 2
            and b % 2 == c % 2 == n % 2)


def match_row(left, right):
    """Check all 72 permitted orderings, reconstructing integer n,B,C."""
    matches = set()
    for a, b in ((left, right), (right, left)):
        for u in permutations(a):
            for v in permutations(b):
                if u[0] + v[0] != F(1, 2) or u[1] - u[0] != v[1]:
                    continue
                if v[1] <= 0:
                    continue
                n = 1 / v[1]
                big_b = 2*n*u[1]
                big_c = big_b + 2*n*v[2]
                if any(x.denominator != 1 for x in (n, big_b, big_c)):
                    continue
                integers = tuple(map(int, (n, big_b, big_c)))
                if admissible(*integers) and target(*integers) == u + v:
                    matches.add(integers)
    return sorted(matches)


def affine(a=0, b=0):
    """a*s+b, represented by a pair of exact rational numbers."""
    return F(a), F(b)


def families():
    a = affine
    return [
        ('I', (a(0,'1/6'), a(1), a(-2,'1/3')),
         (a(1,'1/3'), a(1), a(-1,'1/6')), F(1,6)),
        ('II', (a(0,'1/6'), a(-3,'1/2'), a(1)),
         (a(-1,'1/6'), a(2), a(1,'1/6')), F(1,6)),
        ('III', (a(0,'1/6'), a(-2,'1/6'), a(2)),
         (a(-2,'1/6'), a(1), a(1,'1/2')), F(1,12)),
        ('IV', (a(-4,'1/3'), a(1), a(1,'1/3')),
         (a(-2,'1/6'), a(3), a(1,'1/6')), F(1,12)),
    ]


def check_families():
    expected = {'I':['1/12'], 'II':['1/24','1/12'],
                'III':[], 'IV':['1/24']}
    report, all_matches = {}, []
    for name, left, right, upper in families():
        candidates = set()
        for a, b in left:
            for c, d in right:
                slope, constant = a+c, b+d
                if slope == 0:
                    assert constant != F(1,2), 'An entire family needs checking'
                    continue
                s = (F(1,2) - constant) / slope
                if 0 < s < upper:
                    candidates.add(s)
        strings = [str(s) for s in sorted(candidates)]
        assert strings == expected[name]
        matches = []
        for s in sorted(candidates):
            lhs = tuple(a*s+b for a,b in left)
            rhs = tuple(a*s+b for a,b in right)
            assert all(x > 0 for x in lhs+rhs)
            assert sum(lhs+rhs) == 1
            for n,b,c in match_row(lhs,rhs):
                matches.append([str(s),n,b,c])
                all_matches.append([name,str(s),n,b,c])
        report[name] = {'cross_pair_parameters':strings, 'matches':matches}
    assert all_matches == [['I','1/12',12,4,6]]
    return report


def check_sporadic():
    records = list(csv.DictReader((HERE/'pr-table4.csv').open(newline='')))
    assert len(records) == 65
    expected = [1,1,1,1,6,1,1,1,7,1,1,1,2,1,1,1,1,2,6,
                3,1,1,3,2,2,2,2,2,3,2,2,9,5,1,1,1,1,5,1,1,3,
                7,2,2,11,1,1,1,5,1,3,3,7,6,4,4,19,1,1,1,7,3,2,7,27]
    distances = []
    for j,row in enumerate(records,1):
        assert int(row['row']) == j
        values = tuple(F(row[key]) for key in ('U','V','W','X','Y','Z'))
        assert all(x > 0 for x in values) and sum(values) == 1
        denominator = int(row['denominator'])
        assert denominator == lcm(*(x.denominator for x in values))
        assert denominator % 2 == 0
        distance = min(abs(denominator*(a+b-F(1,2)))
                       for a in values[:3] for b in values[3:])
        assert distance.denominator == 1
        assert distance == int(row['d_j']) == expected[j-1]
        assert distance > 0
        assert match_row(values[:3],values[3:]) == []
        distances.append(int(distance))
    return {'rows':len(records), 'd_j':distances, 'matches':[]}


def main():
    print(json.dumps({'arithmetic':'fractions.Fraction (exact)',
                      'table3':check_families(),
                      'table4':check_sporadic(), 'status':'PASS'}, indent=2))


if __name__ == '__main__':
    main()
