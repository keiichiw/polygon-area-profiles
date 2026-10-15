#!/usr/bin/env python3
"""Verify payload integrity and reproduce the three recorded reports offline."""
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
    manifest = json.loads((ROOT / "manifest.json").read_text())
    for name, expected in manifest["payload_sha256"].items():
        actual = hashlib.sha256((ROOT / name).read_bytes()).hexdigest()
        if actual != expected:
            raise SystemExit("Payload hash mismatch: " + name)
    print("PASS: payload SHA-256 hashes", flush=True)
    for script, report in [("verify_exact.py", "exact-results.json"),
                           ("verify_myerson.py", "myerson-results.txt"),
                           ("verify_geometry.py", "numerical-results.json")]:
        result = subprocess.run([sys.executable, str(ROOT / "supplementary" / script)],
                                cwd=ROOT, text=True, capture_output=True, check=True)
        expected = (ROOT / "supplementary" / report).read_text()
        matches = (json.loads(result.stdout) == json.loads(expected)
                   if report.endswith(".json") else result.stdout == expected)
        if not matches:
            raise SystemExit("Recorded output differs: " + script)
        print("PASS: " + script + " matches " + report, flush=True)

if __name__ == "__main__":
    main()
