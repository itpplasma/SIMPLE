#!/usr/bin/env python3
"""Acceptance gate on retained orbit-study measurements; failure is deliberate."""
import json
from pathlib import Path

root=Path(__file__).resolve().parents[1]
summary=json.loads((root/'results/orbit_summary.json').read_text())
failed=[]
for row in summary['fits']:
    if row['status']=='FAILED':
        failed.append(f"{row['field']}: selector {row['n']} phase-space exponent "
                      f"{row['phase_space_rms_exponent']:.3f}, requested {row['expected']}")
if summary['n2']['status']!='available':
    failed.append('n=2 transform/RHS and its convergence measurements are unavailable')
if failed:
    print('FAIL: requested GC orbit-order acceptance')
    for item in failed:
        print('  '+item)
    raise SystemExit(1)
print('PASS: requested GC orbit-order acceptance')
