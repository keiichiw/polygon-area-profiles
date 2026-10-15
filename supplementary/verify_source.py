#!/usr/bin/env python3
"""Optional local-PDF provenance check. Requires pypdf; performs no network I/O."""
import argparse
import csv
import hashlib
import re
from pathlib import Path
from pypdf import PdfReader

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('pdf', type=Path, help='author version from math.mit.edu/~poonen/papers/ngon.pdf')
args = parser.parse_args()
reader = PdfReader(args.pdf)
text = reader.pages[13].extract_text()
assert 'Table 4.' in text and '65 sporadic' in text
rows = []
for line in text.splitlines():
    fractions = re.findall(r'\b\d+/\d+\b',line)
    if len(fractions) == 6:
        rows.append(fractions)
with (Path(__file__).resolve().parent/'pr-table4.csv').open(newline='') as stream:
    csv_rows = [[row[key] for key in ('U','V','W','X','Y','Z')]
                for row in csv.DictReader(stream)]
assert len(rows) == 65 and rows == csv_rows
print('PASS: all 65 Table 4 rows match the local author PDF, in order.')
print('Source SHA-256:',hashlib.sha256(args.pdf.read_bytes()).hexdigest())
