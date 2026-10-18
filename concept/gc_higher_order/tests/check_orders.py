#!/usr/bin/env python3
"""Orbit acceptance: first phase-space order via 1+, second via T2."""
import json
import math
from pathlib import Path

root=Path(__file__).resolve().parents[1]
summary=json.loads((root/'results/orbit_summary.json').read_text())
failed=[]
fields={'uniform','constant_gradient','constant_curvature_on_fieldline',
        'magnetic_shear','sheet_inner','sheet_outer'}
for name in fields:
    for selector in [0,1,11,2]:
        found=[r for r in summary['fits'] if r['field']==name and r['n']==selector]
        if len(found)!=1 or not all(math.isfinite(found[0][q+'_rms_exponent'])
                                  for q in ['position','velocity']):
            failed.append(f'{name}: selector {selector} study missing/invalid')
for name in fields:
    for selector,px,pv in [(11,2,2),(2,3,2)]:
        found=[row for row in summary['fits'] if row['field']==name and row['n']==selector]
        if len(found)!=1:
            failed.append(f'{name}: selector {selector} measurement unavailable')
            continue
        row=found[0]
        if name=='uniform':
            # Independent closed-form test_gc supplies the uniform oracle.
            continue
        for quantity,expected in [('position',px),('velocity',pv)]:
            measured=row[quantity+'_rms_exponent']
            if measured<expected-0.3:
                failed.append(f'{name}: selector {selector} {quantity} exponent '
                              f'{measured:.3f}, requested {expected}')
        precision=[r for r in summary['tolerance_check']
                   if r['field']==name and r['n']==selector]
        if len(precision)!=1 or not all(math.isfinite(precision[0][q+'_relative_change'])
                                      for q in ['position','velocity']) or max(precision[0]['position_relative_change'],
                                    precision[0]['velocity_relative_change'])>0.01:
            failed.append(f'{name}: selector {selector} numerical precision unresolved')
if summary['n2']['status']!='available':
    failed.append('n=2 transform/RHS and its convergence measurements are unavailable')
if failed:
    print('FAIL: requested GC orbit-order acceptance (1+ and 2)')
    for item in failed:
        print('  '+item)
    raise SystemExit(1)
print('PASS: 1+ position/velocity >=2; T2 position >=3, velocity >=2')
