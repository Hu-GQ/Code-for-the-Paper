"""Verify the deposit against DEPOSIT.sha256.

Run from the paper-2 directory:  python python/verify_deposit.py
Exits non-zero if any listed file is missing or has changed.
"""
import hashlib
import sys
from pathlib import Path

root = Path(__file__).resolve().parents[1]
bad = 0
n = 0
for line in (root / "DEPOSIT.sha256").read_text(encoding="utf-8").splitlines():
    if not line.strip():
        continue
    want, rel = line.split("  ", 1)
    p = root / rel
    if not p.exists():
        print("MISSING", rel)
        bad += 1
        continue
    got = hashlib.sha256(p.read_bytes()).hexdigest()
    if got != want:
        print("CHANGED", rel)
        bad += 1
    n += 1
print(f"checked {n} files, {bad} problem(s)")
sys.exit(1 if bad else 0)
