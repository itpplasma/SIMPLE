#!/usr/bin/env python3
"""Canonical midpoint/GL2 long-time and full-orbit studies, without alignment."""
import argparse
import csv
import compact_json
import sys
import time
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'gen'))
from gc_canonical import CanonicalGC, step
from orbit_study import Field, fields


def initial(field, eps, order):
    seed = np.r_[field.initial, .35, .25, .7]
    particle = field.evaluate(0, 0, eps, seed)
    gc = field.inverse(order, eps, particle)
    model = CanonicalGC(field, eps, order)
    return model, model.to_canonical(gc), particle


def integrate(model, y, times, method, maxstep):
    trajectory = [y.copy()]
    count = 0
    start = time.perf_counter()
    for t0, t1 in zip(times[:-1], times[1:]):
        n = int(np.ceil((t1-t0)/maxstep-1e-10))
        dt = (t1-t0)/n
        for _ in range(n):
            y = step(model, y, dt, method)
        count += n
        trajectory.append(y.copy())
    return np.array(trajectory), count, time.perf_counter()-start


def rms(v):
    return float(np.sqrt(np.mean(np.sum(v*v, axis=1))))


def order_study(order_select=(1, 11, 2), output=None):
    reference = {(row['field'], int(row['n']), float(row['epsilon'])): row
                 for row in csv.DictReader((ROOT/'results/orbit_errors.csv').open())}
    epsilons = [.08, .04, .02, .01, .005]
    times = np.linspace(0, 1, 401)
    rows, fits, precision = [], [], []
    for name, entry in fields.items():
        field = Field(name, entry)
        full = {}
        for eps in epsilons:
            _, _, particle = initial(field, eps, 2)
            fo = solve_ivp(lambda t, z: np.r_[z[3:], np.cross(z[3:], field.b(z[:3]))/eps],
                           (0, 1), particle, t_eval=times, rtol=2e-12, atol=2e-14, max_step=eps/3)
            assert fo.success
            full[eps] = fo.y.T
        for order in order_select:
            for method in ('midpoint', 'gauss2'):
                sample = []
                for eps in epsilons:
                    model, y, _ = initial(field, eps, order)
                    dt = (.25*eps**1.5 if method == 'midpoint' and name.startswith('sheet') else .0025)
                    for refinement in range(5):
                        canonical, steps, seconds = integrate(model, y, times, method, dt)
                        reconstructed = np.array([field.evaluate(0, order, eps, model.to_gc(z)) for z in canonical])
                        error = reconstructed-full[eps]
                        px, pv = rms(error[:, :3]), rms(error[:, 3:])
                        rk = reference[name, order, eps]
                        change = max(abs(px/float(rk['position_rms'])-1), abs(pv/float(rk['velocity_rms'])-1))
                        if name == 'uniform' or change < .01:
                            break
                        dt /= 2
                    assert name == 'uniform' or change < .01, (name, order, method, eps, change)
                    if name.startswith('sheet'):
                        rr = np.hypot(np.hypot(full[eps][:, 0], full[eps][:, 1])-10, full[eps][:, 2])
                        rg = np.array([model.inverse(z)[0][0] for z in canonical])
                        assert (max(np.max(rr), np.max(rg)) < .5 if name == 'sheet_inner'
                                else min(np.min(rr), np.min(rg)) > .5)
                    row = dict(field=name, order=order, method=method, epsilon=eps,
                        position_rms=px, velocity_rms=pv, steps=steps, seconds=seconds,
                        maxstep=dt, rk_relative_error_change=float(change), refinements=refinement)
                    rows.append(row)
                    sample.append(row)
                    print(name, order, method, eps, 'steps', steps, 'RK difference', round(change, 5), flush=True)
                fit = dict(field=name, order=order, method=method)
                for key in ('position_rms', 'velocity_rms'):
                    fit[key+'_exponent'] = float(np.polyfit(np.log(epsilons[1:]), np.log([z[key] for z in sample[1:]]), 1)[0])
                fit['status'] = 'solver_floor' if name == 'uniform' else 'PASS'
                if name != 'uniform':
                    assert fit['position_rms_exponent'] >= (3 if order == 2 else 2)-.3
                    assert fit['velocity_rms_exponent'] >= (2 if order in (11, 2) else 1)-.3
                fits.append(fit)
                # Independent step-halving check at the finest scale.
                if name.startswith('sheet'):
                    model, y, _ = initial(field, epsilons[-1], order)
                    fine, _, _ = integrate(model, y, times, method, sample[-1]['maxstep']/2)
                    reconstructed = np.array([field.evaluate(0, order, epsilons[-1], model.to_gc(z)) for z in fine])
                    error = reconstructed-full[epsilons[-1]]
                    change = max(abs(rms(error[:, :3])/sample[-1]['position_rms']-1),
                                 abs(rms(error[:, 3:])/sample[-1]['velocity_rms']-1))
                    assert change < .01, (name, order, method, 'step-halving', change)
                    precision.append(dict(field=name, order=order, method=method, relative_change=float(change)))
    result = dict(status='PASS', times=[0, 1], samples=401, epsilon=epsilons,
        fit='last four scales; raw RMS position/velocity; matched physical initial state',
        temporal_gate='1% agreement with D2 RK45 errors and finest sheet step halving',
        rows=rows, fits=fits, step_halving=precision)
    target = ROOT/'results/canonical_orders.json' if output is None else output
    target.write_text(compact_json.dumps(result))


def energy_study():
    rows = []
    for name, entry in fields.items():
        field = Field(name, entry)
        for order in (1, 2):
            for method in ('midpoint', 'gauss2'):
                for dt in (2., 1.):
                    model, y, _ = initial(field, .04, order)
                    times = np.arange(0, 4000+dt/2, dt)
                    trajectory, steps, seconds = integrate(model, y, times, method, dt)
                    H = np.array([model.evaluate(z, hessian=False)[0] for z in trajectory])
                    error = H-H[0]
                    blocks = [float(np.max(np.abs(z))) for z in np.array_split(error, 8)]
                    pp_index = 4 if name.startswith('sheet') else 3
                    pp = float(np.max(np.abs(trajectory[:, pp_index]-y[pp_index])))
                    action = float(np.max(np.abs(trajectory[:, 5]-y[5])))
                    assert action == 0
                    if name.startswith('sheet') or name == 'constant_curvature_on_fieldline':
                        assert pp == 0
                    assert blocks[-1] < 1.25*max(blocks[0], 2e-12), (name, method, blocks)
                    assert np.max(np.abs(error)) < 1e-3
                    radii = np.array([model.inverse(z)[0][0] for z in trajectory])
                    if name.startswith('sheet'):
                        assert (np.max(radii) < .5 if name == 'sheet_inner' else np.min(radii) > .5)
                    rows.append(dict(field=name, order=order, method=method, dt=dt,
                        time=4000, steps=steps, seconds=seconds,
                        energy_error_max=float(np.max(np.abs(error))),
                        energy_error_blocks=blocks, energy_endpoint=float(error[-1]),
                        p_phi_error=pp if name.startswith('sheet') or name == 'constant_curvature_on_fieldline' else None,
                        J_error=action, radial_range=[float(np.min(radii)), float(np.max(radii))]))
                    print('long', name, order, method, dt, 'energy', rows[-1]['energy_error_max'], flush=True)
    for name in ('sheet_inner', 'sheet_outer'):
        for order in (1, 2):
            for method, power in [('midpoint', 2), ('gauss2', 4)]:
                pair = [z for z in rows if z['field'] == name and z['order'] == order and z['method'] == method]
                ratio = pair[0]['energy_error_max']/pair[1]['energy_error_max']
                assert 2**(power-.5) < ratio < 2**(power+.5), (name, order, method, ratio)
    (ROOT/'results/canonical_energy.json').write_text(compact_json.dumps(dict(status='PASS',
        criterion='eight equal time blocks: final envelope <=1.25 initial envelope plus roundoff floor; step order on sheets',
        rows=rows)))


def cost_study():
    rows = []
    for name in ('magnetic_shear', 'sheet_inner', 'sheet_outer'):
        field = Field(name, fields[name])
        for order in (1, 2):
            model, y, _ = initial(field, .04, order)
            times = np.linspace(0, 20, 201)
            t0 = time.perf_counter()
            reference = solve_ivp(lambda t, z: model.evaluate(z, hessian=False)[1], (0, 20), y,
                                  t_eval=times, rtol=2e-13, atol=2e-15)
            assert reference.success
            reference_seconds = time.perf_counter()-t0
            for method in ('midpoint', 'gauss2', 'RK45_canonical', 'RK45_D2'):
                model, y, _ = initial(field, .04, order)
                gc_initial = model.to_gc(y) if method == 'RK45_D2' else None
                calls0 = model.calls
                if method in ('midpoint', 'gauss2'):
                    out, count, seconds = integrate(model, y, times, method, .1)
                else:
                    t0 = time.perf_counter()
                    if method == 'RK45_canonical':
                        sol = solve_ivp(lambda t, z: model.evaluate(z, hessian=False)[1], (0, 20), y,
                                        t_eval=times, rtol=1e-9, atol=1e-11)
                        out = sol.y.T
                        seconds = time.perf_counter()-t0
                    else:
                        sol = solve_ivp(lambda t, z: field.evaluate(2, order, .04, z), (0, 20), gc_initial,
                                        t_eval=times, rtol=1e-9, atol=1e-11)
                        seconds = time.perf_counter()-t0
                        calls = model.calls-calls0
                        out = np.array([model.to_canonical(z) for z in sol.y.T])
                    assert sol.success
                    count = sol.nfev
                if method != 'RK45_D2':
                    calls = model.calls-calls0
                # Phase is compared without wrapping or aligning; a static
                # chart angle differs by multiples of 2pi only beyond this run.
                gc_out = np.array([model.to_gc(z) for z in out])
                gc_ref = np.array([model.to_gc(z) for z in reference.y.T])
                rows.append(dict(field=name, order=order, method=method, seconds=seconds,
                    steps_or_nfev=count, analytic_kernel_calls=calls,
                    slow_position_rms=rms(gc_out[:, :3]-gc_ref[:, :3]),
                    carried_phase_rms=rms((gc_out[:, 5]-gc_ref[:, 5])[:, None]),
                    reference_seconds=reference_seconds))
                print('cost', name, order, method, seconds, flush=True)
    (ROOT/'results/canonical_cost.json').write_text(compact_json.dumps(dict(rows=rows,
        scope='single thread Python canonical kernels versus same-kernel RK45 and Fortran D2 RHS; excludes initial inversion and particle reconstruction',
        interval=[0, 20], dt=.1, rk_tolerance=dict(rtol=1e-9, atol=1e-11))))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('study', choices=['orders', 'firstplus', 'energy', 'cost'])
    args = parser.parse_args()
    {'orders': order_study, 'firstplus': lambda: order_study((11,), ROOT/'build/canonical_firstplus.json'),
     'energy': energy_study, 'cost': cost_study}[args.study]()
