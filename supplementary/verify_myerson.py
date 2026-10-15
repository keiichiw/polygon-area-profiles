#!/usr/bin/env python3
"""Exact verification of the Myerson matching used in the manuscript.

Python 3.9+; standard library only. Run from the repository root:
    python3 supplementary/verify_myerson.py
Optionally use --manuscript /path/to/main.tex to also check a local manuscript table.

The 15 rows are transcribed from Luk--Cheung (2024), Table 5, p. 503:
https://doi.org/10.1007/s00026-023-00685-9
Myerson's original classification is Theorem 4 of:
https://doi.org/10.1007/BF01844426

This optionally checks a local manuscript table against these rows, checks all eight
equation-preserving orderings, solves each affine
parameter equation exactly, and verifies every tabulated trigonometric
identity by polynomial reduction modulo a cyclotomic polynomial. It does
not re-prove Myerson's classification, the manuscript's norm restrictions,
or its Poonen--Rubinstein argument. No bound on n and no floating-point
rationality test are used.
"""
from __future__ import annotations

import argparse
import re
from fractions import Fraction as F
from functools import cache
from math import lcm
from pathlib import Path

Quad = tuple[F, F, F, F]
Affine = tuple[F, F]  # (coefficient of sigma, constant)

RAW_ROWS = (
    ('1/21','8/21','1/14','3/14'),
    ('1/14','5/14','2/21','5/21'),
    ('4/21','10/21','3/14','5/14'),
    ('1/20','9/20','1/15','4/15'),
    ('2/15','7/15','3/20','7/20'),
    ('1/30','3/10','1/15','2/15'),
    ('1/15','7/15','1/10','7/30'),
    ('1/10','13/30','2/15','4/15'),
    ('4/15','7/15','3/10','11/30'),
    ('1/30','11/30','1/10','1/10'),
    ('7/30','13/30','3/10','3/10'),
    ('1/15','4/15','1/10','1/6'),
    ('2/15','7/15','1/6','3/10'),
    ('1/12','5/12','1/10','3/10'),
    ('1/10','3/10','1/6','1/6'),
)
ROWS = tuple(tuple(F(x) for x in row) for row in RAW_ROWS)
ORDERS = (
    (0,1,2,3), (0,1,3,2), (1,0,2,3), (1,0,3,2),
    (2,3,0,1), (2,3,1,0), (3,2,0,1), (3,2,1,0),
)
FAMILY: tuple[Affine, ...] = (
    (F(0), F(1,6)), (F(1), F(0)),
    (F(1,2), F(0)), (F(-1,2), F(1,2)),
)
HALF = F(1,2)


def check(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def check_manuscript(path: Path) -> None:
    """Compare the numbered, four-column LaTeX table with all 15 rows."""
    source = path.read_text(encoding='utf-8')
    label = r'\label{tab:myerson}'
    check(source.count(label) == 1, 'Expected one manuscript Myerson table')
    table = source.split(label, 1)[1].split(r'\end{tabular}', 1)[0]
    row_pattern = (r'^\s*(\d+)\s*' +
                   r'&\s*\$(\d+/\d+)\$\s*' * 4 + r'\\\\\s*$')
    records = re.findall(row_pattern, table, flags=re.MULTILINE)
    check(len(records) == len(ROWS), 'Wrong manuscript Myerson table size')
    for j, record in enumerate(records, 1):
        check(int(record[0]) == j, 'Wrong manuscript Myerson row number')
        check(tuple(F(x) for x in record[1:]) == ROWS[j-1],
              f'Manuscript Myerson row {j} differs from verification data')


def check_orders() -> None:
    """Generate the group independently from pair swaps and pair interchange."""
    expected = set()
    for left, right in (((0,1), (2,3)), ((2,3), (0,1))):
        for first in (left, left[::-1]):
            for second in (right, right[::-1]):
                expected.add(first + second)
    check(len(ORDERS) == 8 and set(ORDERS) == expected,
          'Expected all eight equation-preserving orderings')


def trim(p: list[int]) -> list[int]:
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return p


def div_monic(p: list[int], q: list[int]) -> tuple[list[int], list[int]]:
    """Polynomial division in Z[x]; q must be monic. Coefficients ascend."""
    check(bool(q) and q[-1] == 1, 'Divisor must be monic')
    p = trim(p.copy())
    quotient = [0] * max(1, len(p)-len(q)+1)
    while len(p) >= len(q) and p != [0]:
        degree = len(p)-len(q)
        coeff = p[-1]
        quotient[degree] += coeff
        for j, value in enumerate(q):
            p[degree+j] -= coeff*value
        trim(p)
    return trim(quotient), trim(p)


@cache
def cyclotomic(n: int) -> tuple[int, ...]:
    check(n >= 1, 'Cyclotomic order must be positive')
    p = [-1] + [0]*(n-1) + [1]
    for d in range(1,n):
        if n % d == 0:
            p, remainder = div_monic(p, list(cyclotomic(d)))
            check(remainder == [0], 'Cyclotomic construction not exact')
    return tuple(p)


def sine_identity(row: Quad) -> bool:
    """Test sin(pi*a)sin(pi*b)=sin(pi*c)sin(pi*d) exactly."""
    common = lcm(*(x.denominator for x in row))
    order = 2*common
    exponents = [int(common*x) for x in row]
    polynomial = [0]*order
    # sin(pi*x)=(z^(common*x)-z^(-common*x))/(2*i), z primitive order.
    for sign, (u,v) in ((1, exponents[:2]), (-1, exponents[2:])):
        for exponent, coeff in ((u+v,1), (u-v,-1), (-u+v,-1), (-u-v,1)):
            polynomial[exponent % order] += sign*coeff
    _, remainder = div_monic(polynomial, list(cyclotomic(order)))
    return remainder == [0]


def lincomb(terms: list[tuple[int, Affine]], constant: int = 0) -> Affine:
    return (sum((k*x[0] for k,x in terms), F(0)),
            sum((k*x[1] for k,x in terms), F(constant)))


def roots(equations: tuple[Affine, ...]) -> list[F]:
    result: set[F] = set()
    for coefficient, constant in equations:
        if coefficient == 0:
            check(constant != 0, 'Unexpected free parameter: handle explicitly')
            continue
        root = -constant/coefficient
        if 0 < root <= HALF:
            result.add(root)
    return sorted(result)


def evaluate(assignment: tuple[Affine, ...], value: F) -> Quad:
    return tuple(a*value+b for a,b in assignment)


def s_parameters(q: Quad) -> tuple[F,F,F]:
    a,_,c,d = q
    n = 2/c
    return n, 2*n*a, n*(HALF-d)


def admissible(n: F, l: F, t: F) -> bool:
    return (all(x.denominator == 1 for x in (n,l,t))
            and n >= 5 and l >= 1 and t >= 1
            and l+2*t <= n-2 and (n-l) % 2 == 0)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manuscript', type=Path,
                        help='Optional local manuscript for table comparison')
    args = parser.parse_args()
    if args.manuscript is not None:
        check_manuscript(args.manuscript)
        print('PASS: all 15 rows agree with the manuscript Myerson table.')
    check_orders()
    for j,row in enumerate(ROWS,1):
        check(all(0 < x <= HALF for x in row), f'Invalid angle in row {j}')
        check(sine_identity(row), f'Incorrect sine identity in row {j}')
    print('PASS: all 15 tabulated identities, by exact cyclotomic reduction.')

    print('\nS=1: parameter-family roots in 0 < sigma <= 1/2:')
    expected_roots = (
        [F(4,9)], [F(1,3)], [], [F(1,9)],
        [F(1,4), F(1,2)], [F(1,3)], [F(1,2)], [],
    )
    s_geometric: set[tuple[F,F,F]] = set()
    s_admissible: set[tuple[F,F,F]] = set()
    for j, perm in enumerate(ORDERS,1):
        assignment = tuple(FAMILY[i] for i in perm)
        a,b,c,d = assignment
        equations = (
            lincomb([(1,b),(-1,a),(-1,d)]),
            lincomb([(1,a),(1,b),(1,d)], -1),
        )
        values = roots(equations)
        check(values == expected_roots[j-1], f'S family row {j} mismatch')
        print(f'  ordering {j}: {", ".join(map(str, values)) or "none"}')
        for value in values:
            q = evaluate(assignment,value)
            aa,bb,cc,dd = q
            check(sine_identity(q), 'Parameter-family identity failed')
            if 0 < aa < dd < HALF:
                params = s_parameters(q)
                s_geometric.add(params)
                if admissible(*params):
                    s_admissible.add(params)
    check(s_geometric == {(F(9),F(3),F(2)), (F(12),F(3),F(3))},
          'Wrong S-family candidates before parity')
    check(s_admissible == {(F(9),F(3),F(2))}, 'Wrong admissible S family')
    print('PASS: candidates (9,3,2), (12,3,3); parity leaves (9,3,2).')

    s_sporadic = set()
    r_sporadic = set()
    for j,row in enumerate(ROWS,1):
        for perm in ORDERS:
            q = tuple(row[i] for i in perm)
            a,b,c,d = q
            if b == a+d or a+b+d == 1:
                s_sporadic.add((j,q))
            if d == a+b and d < HALF:
                r_sporadic.add((j,q))
    expected_s = {
        (12,(F(1,10),F(1,6),F(4,15),F(1,15))),
        (13,(F(1,6),F(3,10),F(7,15),F(2,15))),
    }
    check(s_sporadic == expected_s, 'Wrong S sporadic matching')
    print('\nS=1: sporadic matching (all 15 x 8 orderings):')
    for j,q in sorted(s_sporadic):
        a,_,_,d = q
        check(a > d, 'Unexpected admissible angle ordering')
        print(f'  row {j}: ({", ".join(map(str,q))}), a={a} > d={d}')
    print('PASS: both candidates contradict a < d.')

    print('\nR=1: parameter-family roots in 0 < sigma <= 1/2:')
    expected_r_roots = (
        [F(2,9)], [], [F(2,9)], [], [HALF], [], [HALF], [],
    )
    r_family = set()
    for j, perm in enumerate(ORDERS,1):
        assignment = tuple(FAMILY[i] for i in perm)
        a,b,c,d = assignment
        values = roots((lincomb([(1,d),(-1,a),(-1,b)]),))
        check(values == expected_r_roots[j-1], f'R family row {j} mismatch')
        print(f'  ordering {j}: {", ".join(map(str, values)) or "none"}')
        for value in values:
            q = evaluate(assignment,value)
            check(sine_identity(q), 'Parameter-family identity failed')
            aa,bb,cc,dd = q
            if all(0 < x < HALF for x in q) and dd == aa+bb:
                n = 2/cc
                r_family.add((n, 2*n*aa, 2*n*bb))
    check(r_family == {(F(18),F(6),F(8)), (F(18),F(8),F(6))},
          'Wrong R-family result')
    expected_r = {
        (12,(F(1,10),F(1,6),F(1,15),F(4,15))),
        (12,(F(1,6),F(1,10),F(1,15),F(4,15))),
        (13,(F(1,6),F(3,10),F(2,15),F(7,15))),
        (13,(F(3,10),F(1,6),F(2,15),F(7,15))),
    }
    check(r_sporadic == expected_r, 'Wrong R sporadic matching')
    print('\nR=1: family gives n=18, {A,B}={6,8}.')
    r_all = set(r_family)
    for j,q in sorted(r_sporadic):
        a,b,c,d = q
        n = 2/c
        r_all.add((n,2*n*a,2*n*b))
    good_r = set()
    for n,A,B in r_all:
        k = A/2
        l = B
        t = (n-A-B)/2
        if k.denominator == 1 and k >= 1 and admissible(n,l,t):
            good_r.add((int(n),int(k),int(l)))
    check(good_r == {(18,3,8),(18,4,6),(30,3,10),(30,5,6)},
          'Wrong final R classification')
    print('R=1: sporadic row 12 gives n=30, {A,B}={6,10}.')
    print('R=1: sporadic row 13 gives n=15, {A,B}={5,9}; A is not even.')
    print('PASS: final (n,k,l) =', sorted(good_r))
    print('\nALL CHECKS PASSED. External classification theorems are assumed.')


if __name__ == '__main__':
    main()
