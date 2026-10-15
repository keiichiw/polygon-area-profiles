# Polygon area profiles
[![DOI](https://zenodo.org/badge/1420863632.svg)](https://doi.org/10.5281/zenodo.23399092)

Verification code and reproducibility data for *Rational Area Profiles of Regular
Polygons Cut by Two Congruent Diagonals*, by Shoichi Inoue and Keiichi Watanabe.

## Reproduce the calculations

Python 3.9 or later is required. The three main checks use only the standard
library; no package installation is needed. Run from this directory:

```sh
python3 scripts/check_reproduction.py
```

This runs the checks and compares their output with the recorded results. It
also checks the SHA-256 hashes in `manifest.json`. Individual commands are:

```sh
python3 supplementary/verify_exact.py
python3 supplementary/verify_myerson.py
python3 supplementary/verify_geometry.py
```

Do not enable Python optimization (`-O` or `PYTHONOPTIMIZE`): some checks rely
on assertions. See [verification details](supplementary/README.md) for inputs,
source citations, optional PDF transcription checking, and the scope of each check.

## Relation to the paper

| Program | Manuscript labels | Scope |
| --- | --- | --- |
| `verify_myerson.py` | `tab:myerson`, `lem:S-one`, `lem:R-one` | Exact identities and finite matching for S=1 and R=1 |
| `verify_exact.py` | `lem:T-one`, `app:finite-matching` | Four PR families and 65 sporadic rows |
| `verify_geometry.py` | `prop:area-reduction`, `thm:profile`, `rem:Q2-geometry` | Area formulas, classified values, symmetry and reflection |
| `verify_source.py` | PR Table 4 used in the appendix | Optional transcription comparison with a local source PDF |

The exact checks assume the external classification theorems. They do not
prove those theorems or the manuscript's norm arguments. Numerical checks
cover 19,760 configurations with 5 <= n <= 80 at 60-digit precision and absolute
tolerance 1e-50; they do not prove irrationality or classify larger polygons.
The PR expected sequence is embedded in code; the Myerson table is compared with a local manuscript only
when the optional `--manuscript` argument is provided.

## Provenance and citation

The manuscript, bibliography, paper PDFs, figure artwork and figure-generation
scripts remain in the private manuscript repository. This repository contains
only verification code, required inputs, recorded results and instructions.
`manifest.json` records the source manuscript commit and file hashes for provenance.

The main verification commands work without the manuscript. Authors can also
check a local manuscript's Myerson table explicitly:

```sh
python3 supplementary/verify_myerson.py --manuscript /path/to/private/main.tex
```

`CITATION.cff` identifies **Keiichi Watanabe** as the sole software author.
The accompanying paper is authored by **Shoichi Inoue and Keiichi Watanabe**.
Please cite the software when using the code and also cite the paper for its
mathematical results. A release DOI will be added after a Zenodo deposit.
The repository is hosted at <https://github.com/keiichiw/polygon-area-profiles>.

## License and contributions

The Python source code and software documentation are licensed under the
[MIT License](LICENSE), copyright 2026 Keiichi Watanabe. Keiichi Watanabe is
the sole contributor to the code; paper authorship is separate.

The MIT grant covers the `.py` files in `scripts/` and `supplementary/`, and
the software documentation (`README.md`, `AGENTS.md` and
`supplementary/README.md`). It does not apply to transcribed source data
or recorded calculation outputs. No additional license for those
materials is granted here; retain their provenance and source citations.

To update the payload, intentionally update the input data and expected
outputs, rerun verification, and regenerate the manifest hashes after reviewing
the changes. Do not silently update a published release.
