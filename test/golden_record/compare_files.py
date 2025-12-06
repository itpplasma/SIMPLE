#!/usr/bin/env python3
"""Compare two numerical files for bit-identical results.

With deterministic floating-point builds (-ffp-contract=off, no -ffast-math),
golden record tests should produce bit-identical results. This script requires
exact equality, not approximate tolerance.
"""
import sys
import numpy as np

OLD_FILE = sys.argv[1]
NEW_FILE = sys.argv[2]

old_data = np.loadtxt(OLD_FILE)
new_data = np.loadtxt(NEW_FILE)

# Handle 1D arrays
if old_data.ndim == 1:
    old_data = old_data.reshape(-1, 1)
if new_data.ndim == 1:
    new_data = new_data.reshape(-1, 1)

# Check shapes match
if old_data.shape != new_data.shape:
    print(f"Shape mismatch: {OLD_FILE} has {old_data.shape}, {NEW_FILE} has {new_data.shape}")
    sys.exit(1)

# Require bit-identical results (deterministic FP builds)
mismatches = []
for i, (old_row, new_row) in enumerate(zip(old_data, new_data)):
    for j, (old_d, new_d) in enumerate(zip(old_row, new_row)):
        if old_d != new_d:  # Exact comparison, not np.isclose
            rel_diff = abs(old_d - new_d) / max(abs(old_d), 1e-300)
            mismatches.append((i, j, old_d, new_d, rel_diff))

if len(mismatches) == 0:
    print(f"Files {OLD_FILE} and {NEW_FILE} match (bit-identical).")
    sys.exit(0)

print(f"Found {len(mismatches)} non-matching entries comparing {OLD_FILE} and {NEW_FILE}.")
print("First 5 mismatches (row, col, ref_val, cur_val, rel_diff):")
for row, col, ref, cur, rel in mismatches[:5]:
    print(f"  [{row},{col}]: ref={ref:.16e}, cur={cur:.16e}, rel_diff={rel:.2e}")
if len(mismatches) > 5:
    print(f"  ... and {len(mismatches) - 5} more")
sys.exit(1)
