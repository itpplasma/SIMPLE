#!/usr/bin/env python3
"""Acceptance gates for retained measurements; never substitutes for a study."""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
def read(name):
    return json.loads((ROOT/'results'/name).read_text())

symbolic = read('canonical_symbolic.json')
assert symbolic['status'] == 'PASS'
checks = read('canonical_checks.json')
assert checks['status'] == 'PASS'
fields = {'uniform', 'constant_gradient', 'constant_curvature_on_fieldline',
          'magnetic_shear', 'sheet_inner', 'sheet_outer'}
methods = {'midpoint', 'gauss2'}
expected = {(f, n, m) for f in fields for n in (1, 2) for m in methods}
assert {(z['field'], z['order'], z['method']) for z in checks['rows']} == expected
assert max(z['symplectic_defect'] for z in checks['rows']) < 2e-12
assert max(z['tangent_fd_error'] for z in checks['rows']) < 3e-7

orders = read('canonical_orders.json')
rk = {(z['field'], z['n']): z for z in read('orbit_summary.json')['fits']}
expected_orders = {(f, n, m) for f in fields for n in (1, 11, 2) for m in methods}
assert {(z['field'], z['order'], z['method']) for z in orders['fits']} == expected_orders
for z in orders['fits']:
    if z['field'] == 'uniform':
        continue
    for key in ('position_rms_exponent', 'velocity_rms_exponent'):
        assert abs(z[key]-rk[z['field'], z['order']][key]) < .04, z
    assert z['position_rms_exponent'] >= (3 if z['order'] == 2 else 2)-.3
    assert z['velocity_rms_exponent'] >= (2 if z['order'] in (11, 2) else 1)-.3
assert max(z['relative_change'] for z in orders['step_halving']) < .01
assert {(z['field'], z['order'], z['method']) for z in orders['step_halving']} == {
    (f, n, m) for f, n, m in expected_orders if f.startswith('sheet')}
observed = {(z['field'], z['order'], z['method'], z['epsilon']) for z in orders['rows']}
assert observed == {(f, n, m, e) for f, n, m in expected_orders for e in orders['epsilon']}
assert max(z['rk_relative_error_change'] for z in orders['rows'] if z['field'] != 'uniform') < .01

energy = read('canonical_energy.json')
assert energy['status'] == 'PASS'
assert {(z['field'], z['order'], z['method'], z['dt']) for z in energy['rows']} == {
    (f, n, m, h) for f, n, m in expected for h in (1., 2.)}
for z in energy['rows']:
    assert z['time'] == 4000
    assert z['J_error'] == 0
    assert z['p_phi_error'] in (0, None)
    assert z['energy_error_blocks'][-1] < 1.25*max(z['energy_error_blocks'][0], 2e-12)
for f in ('sheet_inner', 'sheet_outer'):
    for n in (1, 2):
        for m, power in [('midpoint', 2), ('gauss2', 4)]:
            pair = sorted([z for z in energy['rows'] if (z['field'], z['order'], z['method']) == (f, n, m)], key=lambda z: z['dt'])
            ratio = pair[1]['energy_error_max']/pair[0]['energy_error_max']
            assert 2**(power-.5) < ratio < 2**(power+.5)
cost = read('canonical_cost.json')
assert {(z['field'], z['order'], z['method']) for z in cost['rows']} == {
    (f, n, m) for f in ('magnetic_shear', 'sheet_inner', 'sheet_outer') for n in (1, 2)
    for m in ('midpoint', 'gauss2', 'RK45_canonical', 'RK45_D2')}
print('PASS: canonical derivatives/maps, bounded energy, cyclic momenta, D2 orbit exponents and cost matrix')
