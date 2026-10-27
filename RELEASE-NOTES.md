# v1.1.0 — Relative-area and all-diagonal experiments

- Exact relative-area quotients for all 39,200 admissible configurations through n=100.
- Exploratory all-diagonal cell calculations for every 3 <= n <= 100 using rotational quotient arrangements, with exact intersections and weighted areas and numerical ordering checked at 60 and 120 digits.
- Independent full-expansion comparisons through n=30 and at n=36,40, including exact area multisets; the earlier full-expansion report is retained.
- Arithmetic and geometric regression tests, recorded results, and reproduction instructions.
- Standard-library arithmetic optimizations and resumable checkpoints with source-fingerprint invalidation.

All calculation material is distributed together. The default `python3 scripts/check_reproduction.py` checks every payload hash, runs arithmetic regression tests, and reproduces the original three reports. `--extended` also recomputes the additional full-range experiments. Full-range results were validated on 2026-10-23; the default suite and current manuscript table comparison were rerun for this release. Numerical planar ordering remains exploratory, not interval-certified.

The public archive excludes manuscript sources and figure-generation code. Software authorship and license scope are unchanged.

# v1.0.0 — Initial release

Verification code and reproducibility data accompanying *Rational Area Profiles of Regular Polygons Cut by Two Congruent Diagonals*.

- Exact Myerson sine-product identities and finite matching.
- Exact Poonen–Rubinstein matching for four families and 65 sporadic rows.
- Numerical area checks for 19,760 configurations with 5 <= n <= 80.
- Input data, recorded results, source provenance and an offline reproduction command.

Run `python3 scripts/check_reproduction.py` with Python 3.9 or later. The three main checks use only the standard library. Manuscript comparison is optional; manuscript sources and figures are not included.

Software author: Keiichi Watanabe. Paper authors: Shoichi Inoue and Keiichi Watanabe.
Code and software documentation are MIT-licensed; see README for scope and data provenance.
