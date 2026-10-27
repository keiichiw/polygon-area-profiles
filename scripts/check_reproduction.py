#!/usr/bin/env python3
"""Verify payload integrity and reproduce recorded reports offline (new full-range checks: --extended)."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parent.parent

def main():
    if sys.flags.optimize or os.environ.get("PYTHONOPTIMIZE"):
        raise SystemExit("Run without -O or PYTHONOPTIMIZE; verification uses assertions.")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--extended', action='store_true',
                        help='also rerun all exact ratios and exploratory cell arrangements')
    parser.add_argument('--cell-checkpoint-dir',type=Path,
                        help='reuse matching cell checkpoints during --extended')
    args = parser.parse_args()
    if args.cell_checkpoint_dir and not args.extended:
        parser.error('--cell-checkpoint-dir requires --extended')
    manifest = json.loads((ROOT / "manifest.json").read_text())
    for name, expected in manifest["payload_sha256"].items():
        actual = hashlib.sha256((ROOT / name).read_bytes()).hexdigest()
        if actual != expected:
            raise SystemExit("Payload hash mismatch: " + name)
    print("PASS: payload SHA-256 hashes", flush=True)
    checks = [("verify_exact.py", "exact-results.json"),
                           ("verify_myerson.py", "myerson-results.txt"),
                           ("verify_geometry.py", "numerical-results.json")]
    if args.extended:
        checks += [("verify_relative_areas.py", "relative-area-results.json"),
                   ("explore_cells.py", "cell-results.json"),
                   ("explore_cell_orbits.py", "cell-orbit-results.json"),
                   ("verify_cell_orbits.py", "cell-orbit-reference-results.json")]
    subprocess.run([sys.executable, str(ROOT / "supplementary/verify_arithmetic.py")],
                   cwd=ROOT, check=True)
    for script, report in checks:
        command = [sys.executable, str(ROOT / "supplementary" / script)]
        if script=="explore_cell_orbits.py" and args.cell_checkpoint_dir:
            command += ["--checkpoint-dir",str(args.cell_checkpoint_dir.resolve())]
        result = subprocess.run(command,cwd=ROOT,text=True,capture_output=True,check=True)
        expected = (ROOT / "supplementary" / report).read_text()
        matches = (json.loads(result.stdout) == json.loads(expected)
                   if report.endswith(".json") else result.stdout == expected)
        if not matches:
            raise SystemExit("Recorded output differs: " + script)
        print("PASS: " + script + " matches " + report, flush=True)

if __name__ == "__main__":
    main()
