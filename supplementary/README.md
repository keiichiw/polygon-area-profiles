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


## Relative-area and cell experiments (v1.1.0)

These independently implemented additions support the manuscript's concluding
remarks. They are not included in archived v1.0.0. They use only the Python 3.9+
standard library and run offline:

```sh
python3 supplementary/verify_arithmetic.py
python3 supplementary/verify_relative_areas.py
python3 supplementary/explore_cells.py
python3 scripts/check_reproduction.py --extended
```

`cyclotomic.py` represents Q(zeta_N) in the power basis modulo Phi_N using integer
coefficient vectors with a common rational denominator. Polynomial Euclidean
inversion, equality, intersections, shoelace areas and rationality tests are exact.
It reuses the cyclotomic-polynomial construction in `verify_myerson.py`.

`verify_relative_areas.py` constructs four areas directly from boundary vertices
and the diagonal intersection in Q(zeta_(2n)). A quotient is rational exactly
when the two reduced coefficient vectors are rationally proportional. The saved
`relative-area-results.json` covers 39,200 configurations, 5 <= n <= 100.
Only (n,k,m)=(6,1,2) has rational quotients a0/a1=1/2, a0/a2=1/13 and a1/a2=2/13.
An optional `--n-max N` restricts or extends this finite computation.

`explore_cells.py` builds the arrangement of every chord, including polygon
sides and all diagonals (diameters included). Every alternating quadruple of
vertices supplies an exact crossing; exact field equality merges concurrency.
At n=3 the polygon is one cell. Distinct intersections on a chord and cyclic
orders of incident edges are determined numerically using Decimal embeddings.
Half-edge traversal then enumerates the faces. This numerical ordering is not
interval-certified, so the computation remains exploratory even though each
reported area and its rationality test are exact conditional on the arrangement.

Both 60- and 120-digit orderings must yield identical face cycles and exact area
vectors. Gaps at or below 10^(-(precision-20)) cause an error, except that
components used to select an angular half-plane are treated as zero within
that threshold. Connectivity, simple face boundaries, Euler characteristic 2,
one exterior face of normalized signed area -1, positive bounded areas and
an exact bounded-area sum of 1 are checked. These checks do not certify the
numerical ordering or prove the conjecture for untested polygons.

The saved `cell-results.json` records all n=3,...,26 and 30. Rational cells occur
only for n=3,4,6, where every cell is rational (1, 4 and 24 cells respectively).
For n=30 there are 21,480 cells, none with rational normalized area.
`--n 5 12 30` selects particular polygons. Reports include exact rational-area
histograms, graph counts, and hashes of canonical face cycles and exact area
vectors. Progress goes to stderr, JSON reports to stdout.

`verify_arithmetic.py` checks field inverses, conjugation, nonrational controls,
zero-division failures, the two unequal-diagonal dodecagon areas, exhaustive
reflection representatives in that dodecagon, agreement with independently
evaluated trigonometric formulas, and small arrangements with concurrency.
The existing n<=80 numerical profile check is separate from these experiments.


## Extension to n <= 100 by rotational symmetry

`explore_cell_orbits.py` uses the same exact cyclotomic arithmetic but stores a
rotational quotient instead of expanding all intersections. Its saved report is
`cell-orbit-results.json`; the earlier full-expansion `cell-results.json` remains
as an independent reference. Run from an extracted Online Resource 1 root:

```sh
python3 supplementary/explore_cell_orbits.py --output cells.json --benchmark-output timings.json
python3 supplementary/explore_cell_orbits.py --n 40 50 60 100 --checkpoint-dir checkpoints --output selected.json
python3 supplementary/verify_cell_orbits.py
```

The default range is every integer 3 <= n <= 100, including cases already proved
in Appendix B. `--checkpoint-dir` writes one completed n atomically; rerunning
the same command reuses these records. A fingerprint of the algorithm, exact
arithmetic, reference geometry and public helper sources rejects stale checkpoints.
The output and benchmark files are written only after the selected run completes;
completed checkpoints survive an interrupted run. Benchmark times and process
peak RSS are separate from the deterministic mathematical reports. RSS uses bytes
and is a process-lifetime high-water mark, so use a fresh process for each n when
comparing individual memory requirements. Checkpoints and timings are optional.

For four vertices, the four positive cyclic gaps sum to n. Keeping only the
lexicographically smallest cyclic gap tuple selects one representative of each
rotation orbit of crossing pairs, including configurations with a stabilizer.
A nonzero intersection is assigned to the half-open angular sector [0,2*pi/n),
and its canonical coordinate is reduced exactly in Q(zeta_n). The sector is
checked independently at both precisions; near-boundary collinearity is tested
by an exact wedge product. The fixed origin is handled separately. Original
vertices form one orbit. No numerical closeness test merges exact points.

Every chord is a rotation of (0,s), 1 <= s <= floor(n/2). Crossings are accumulated
on these representative chords, with both endpoint choices for a diameter.
The two precision runs order their cuts and incident edges independently.
Directed edge orbits and a rotation offset determine a face walk. At the origin,
the cyclic order includes every actual incident ray. A face's multiplicity is
n divided by its stabilizer, computed from the repeated directed edge orbits in
its closed boundary. This covers the fixed exterior face, the odd-n central
cell, and the even-n fixed central vertex. Edge orbits must have size n.

The polygon boundary cycle and connectivity of the quotient imply connectivity
of the lifted graph. Weighted vertex, edge and face counts satisfy Euler's
formula. Quotient face-cycle signatures, including rotation offsets, must agree at 60 and 120 digits. Areas are
then computed exactly once per representative cell, since rotations preserve
area. Every representative has positive real area, the exterior has signed area
-1, and multiplicity-weighted bounded areas sum exactly to 1. Rationality and
rational-area histograms are exact. As with full expansion, angular ordering is
not interval-certified; this remains an exploratory computation.

`verify_cell_orbits.py` compares against independent full expansion for every
3 <= n <= 30 and also n=36,40. It checks vertex/edge/cell counts, rational-area
histograms, and a hash of the complete exact area multiset, with orbit
multiplicities. The recorded comparison is `cell-orbit-reference-results.json`.
The original full-expansion reports are unchanged. Its face walk now visits
directed edges once rather than repeatedly scanning the entire remaining set.
Small-range comparisons and rotational-sector tests also run in the normal
regression suite. `../scripts/check_reproduction.py --extended` additionally reproduces the
full n <= 100 report and the full reference comparison.

`cell-orbit-benchmarks.json` records indicative local measurements for n=40,50,60,100
in separate processes, including the Python version and architecture. These
measurements are not correctness tests or portable performance guarantees.

The exact arithmetic also skips zero terms in conjugation, uses cyclic coefficient
shifts for rotations, and folds products modulo x^n-1 before reduction modulo
Phi_n. Dense products use exact integer Kronecker substitution with a checked
signed-digit carry bound; sparse products retain direct convolution. Regression
tests compare these operations with naive convolution and long division,
including negative coefficients, 30-digit coefficients, and unreduced polynomials
longer than 3n, through n=100. These changes preserve the field and its exact
rationality test and require no additional libraries.

For a resumable complete verification, use
`python3 scripts/check_reproduction.py --extended --cell-checkpoint-dir checkpoints`.
Without that optional flag, the checker recomputes every report from scratch.
Fingerprint-matched checkpoints only replace the n <= 100 orbit calculation;
all report comparisons, helper hashes, other experiments and regression tests
still run. A fresh full n <= 100 computation is what produced the saved report.

The completed report covers all 98 integers 3 <= n <= 100. Rational cells occur
only at n=3,4,6. For n=100 it records 3,731,301 vertices, 7,561,800 edges and
3,830,500 bounded cells, represented by 38,305 bounded-cell rotation orbits.
The separately measured n=100 run took approximately 30 seconds and 279 MB peak
RSS on the recorded machine; these are illustrative measurements, not bounds.
