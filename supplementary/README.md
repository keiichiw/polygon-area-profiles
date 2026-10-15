# Verification material

Article: *Rational Area Profiles of Regular Polygons Cut by Two Congruent Diagonals*

Journal: *Journal of Geometry*

Paper authors: Shoichi Inoue and Keiichi Watanabe

Software author: Keiichi Watanabe

Author e-mail addresses: shoichi31415@gmail.com; keiichi@keiichiw.dev

Caption: Exact verification of the Myerson and Poonen--Rubinstein finite
matching and numerical checks of the polygon area formulas and classified
values. The repository contains the verification programs, input data, instructions, and recorded outputs.

This directory accompanies the manuscript on rational areas obtained by cutting
a regular polygon with two congruent diagonals. It contains exact checks of
the finite Myerson and Poonen--Rubinstein matching and separate numerical
consistency checks.
The proofs remain in the manuscript.

## Run

Run these commands from the repository root. No manuscript or TeX installation
is required.

Python 3.9 or later is sufficient; the three main scripts use only the standard
library. Run from the extracted archive root or repository root (or invoke each
script by its full path):

```sh
python3 supplementary/verify_exact.py
python3 supplementary/verify_myerson.py
python3 supplementary/verify_geometry.py
```

All three commands fail with an assertion error on a discrepancy. The Myerson
check prints a text report; the other two print JSON reports. Run without
Python's `-O` option, which disables assertions in the latter two scripts.
Recorded outputs are `myerson-results.txt`, `exact-results.json`, and
`numerical-results.json`.

## Myerson matching

`verify_myerson.py` uses 15 rows transcribed from Luk--Cheung, *Rational Angles and Tilings of the Sphere by Congruent
Quadrilaterals*, **Annals of Combinatorics 28** (2024), 485--527,
[Theorem 2 and Table 5, p. 503](https://doi.org/10.1007/s00026-023-00685-9).
This restates Myerson's Theorem 4 and corrects the second entry of row 13 to
`7/15`. The row numbering and order of entries follow that published table.

The script verifies all 15 sine-product identities by integer polynomial
reduction modulo cyclotomic polynomials. Using `fractions.Fraction`, it then
checks the eight permutations preserving the two pairs, solves the affine
parameter equations for both S and R, and compares every sporadic assignment
with the linear angle conditions in the proofs. It checks the angle ranges,
integrality, parity, and admissibility of the surviving candidates. For S,
the two sporadic assignments violate `a < d`; for R, row 13 is excluded by
the parity of A. No finite bound on n or numerical rationality test is used.

The default run checks the mathematical identities and matching without a manuscript.
Optionally pass `--manuscript /path/to/main.tex` to additionally compare all 15 rows
with a local manuscript table. The manuscript is not distributed here; this
optional run adds one table-comparison PASS line to the recorded default output.
This verifies the manuscript's finite deductions, not the completeness of
Myerson's external classification or the separate norm arguments.

## Poonen--Rubinstein matching

`verify_exact.py` uses `fractions.Fraction` throughout. For each of the four
Table 3 families, it solves all nine affine cross-pair equations
`u1 + v1 = 1/2`, retains the permitted parameter range, and checks all
permutations within the triples and both orders of the triples. It reconstructs
integer `n,B,C`, checks admissibility, and compares all six coordinates to the
target in the manuscript. The only match is family I at `s=1/12`, with
`(n,B,C)=(12,4,6)`.

`pr-table4.csv` records all 65 Table 4 rows, their least common denominators, and
the manuscript's `d_j`. The script recomputes all 585 cross-sums exactly,
checks all 65 values against the expected sequence embedded in the script (transcribed from the manuscript, not parsed from it at runtime), verifies that
each is positive, and independently runs the full ordered-tuple matching.
No sporadic row matches.

The source is B. Poonen and M. Rubinstein, *The number of intersection points
made by the diagonals of a regular polygon*, **SIAM Journal on Discrete
Mathematics 11** (1998), no. 1, 135--156,
[DOI 10.1137/S0895480195281246](https://doi.org/10.1137/S0895480195281246).
The table transcription was checked against the
[author-hosted version dated November 18, 1997](https://math.mit.edu/~poonen/papers/ngon.pdf):
Table 3 on printed page 13 and Table 4 on printed page 14.
Rows are numbered 1--65 from top to bottom, preserving the printed `U,V,W;X,Y,Z`
column order. Denominator groups are 30, 42, 60, 84, 90, 120, 210.
This explicit version fixes the row convention without assuming identical
pagination in the journal version.

The checked PDF's SHA-256 is
`1e70b304b7441b91b8c37dd2de731009b15f344bdaeeb1b118967426f2b98936`.
The PDF itself is not redistributed. An optional source-transcription check
requires `pypdf` and a locally downloaded copy:

```sh
python3 supplementary/verify_source.py /path/to/ngon.pdf
```

It compares all 65 six-tuples to PDF page 14, prints the source hash, and performs
no network access. `verify_exact.py` does not need the PDF or `pypdf`.

## Numerical checks and their limits

`verify_geometry.py` uses 60 significant decimal digits and an absolute
tolerance of `1e-50`. Pi is computed with Gauss--Legendre iteration; sine and
cosine are evaluated by convergent Taylor series using `decimal.Decimal`.
No binary floating-point arithmetic or external numerical library is used.

For every admissible `(n,k,m)` with `5 <= n <= 80` and `1 <= k < m < n/2`,
the script constructs Cartesian vertices on the unit circle, intersects the
two straight lines, and obtains the four region areas by the shoelace formula.
Those measurements are compared to all the manuscript's area formulas,
the equality of the two expressions for `T(n,l,t)`, the symmetry `a1=a3`,
and total area 1. It checks all stated rational values in this range and also checks the reflection description of the `E2`
family. The latter verifies that `X = V_m + V_(n-m+k)` in Cartesian coordinates,
the reflection of the origin in their chord.

For configurations outside each classified set, it additionally checks that
the quantity recovered from measured areas is separated from that theorem's
finite set of possible rational values: `R` in `{1/3,1}`, `S` in `{2/3,1}`,
`T` in `{1/3,1}`, and `U` in `{1,2}`. This is a consistency check of the stated
classification; it is not an independent test for arbitrary rational numbers.

The recorded default run covered **19,760 configurations**. The largest
area-formula residual was `1.8e-59`; the largest residual for a stated rational
area value was `8e-60`. The counts were 5 configurations in `E0`, 39 in `E2`,
2 in `E1`, and 743 in `E_adj` within the tested range. Complete error and
separation statistics are in `numerical-results.json`.

Numerical checks do not prove rationality or irrationality and do not replace
the exact arguments. They neither classify cases with `n>80` nor establish the
completeness of the external Poonen--Rubinstein theorem.
