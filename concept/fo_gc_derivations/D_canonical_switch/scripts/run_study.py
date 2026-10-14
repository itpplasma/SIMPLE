#!/usr/bin/env python3
"""Run the item-D switch study on the SPECTRE fixtures and fit eps exponents.

Usage: python3 scripts/run_study.py   (after scripts/build.sh)
Writes results/switch_study.jsonl (raw) and results/switch_study.json (fits).
"""
import json
import math
import subprocess
import sys
from pathlib import Path

D = Path(__file__).resolve().parents[1]
DATA = D.parents[2] / 'test' / 'tests' / 'data'
EXE = D / 'build' / 'switch_study'
NPER, SPP, SPER = 5, 1000, 40
CASES = {
    'tok2vol_sheet': (DATA / 'tok2vol_sheet.h5', 0.02,
                      [(1, (0.5, 0.3, 0.2), 0.3, 0.7), (1, (0.6, 2.5, 1.0), -0.6, 2.0),
                       (2, (1.5, 1.2, 0.5), 0.3, 0.7), (2, (1.4, 4.0, 2.0), -0.5, 4.0)]),
    'ncsx_freeboundary': (DATA / 'ncsx_freeboundary.h5', 0.01,
                          [(1, (0.5, 0.4, 0.3), 0.3, 0.7), (2, (1.5, 2.0, 0.1), -0.4, 2.0),
                           (3, (2.5, 4.0, 0.7), 0.2, 5.0)]),
}
NEPS = 7


def lines():
    out = []
    for label, (path, eps0, seeds) in CASES.items():
        for owner, u, lam, phi in seeds:
            for k in range(NEPS):
                eps = eps0 / 2**k
                out.append((label, path, owner, u, lam, phi, eps, SPP))
                if k >= NEPS - 2:  # time-resolution check at the two smallest eps
                    out.append((label, path, owner, u, lam, phi, eps, 2 * SPP))
    return out


def fit(eps, val):
    xs = [math.log(e) for e in eps]
    ys = [math.log(v) for v in val]
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    return sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sum((x - mx)**2 for x in xs)


def main():
    raw = D / 'results' / 'switch_study.jsonl'
    if '--fit-only' not in sys.argv:
        text = '\n'.join(
            f"{lab} '{p}' {o} {u[0]} {u[1]} {u[2]} {lam} {phi} {eps!r} {NPER} {spp} {SPER}"
            for lab, p, o, u, lam, phi, eps, spp in lines()) + '\n'
        res = subprocess.run([str(EXE)], input=text, capture_output=True, text=True, check=True)
        raw.write_text(res.stdout)
    rows = [json.loads(l) for l in raw.read_text().splitlines() if l.strip()]
    base = [r for r in rows if r['steps_per_period'] == SPP]
    fine = {(r['label'], r['owner'], tuple(r['seed']), r['eps']): r
            for r in rows if r['steps_per_period'] == 2 * SPP}
    summary = {'protocol': dict(nperiods=NPER, steps_per_period=SPP, samples_per_period=SPER,
                                n_eps=NEPS, eps_ratio=2, fit='least squares of log defect vs '
                                'log eps; all points and the smallest four'),
               'cases': [], 'max_roundtrip': {}}
    worst = dict(z=0.0, w=0.0, simple=0.0, chart=0.0, chart_rel=0.0)
    time_check = 0.0
    for r in base:
        worst['z'] = max(worst['z'], *r['roundtrip_z'])
        worst['w'] = max(worst['w'], *r['roundtrip_w'])
        worst['simple'] = max(worst['simple'], *r['roundtrip_simple'])
        worst['chart'] = max(worst['chart'], r['roundtrip_fo_chart'])
        worst['chart_rel'] = max(worst['chart_rel'], r['roundtrip_fo_chart_rel'])
        f = fine.get((r['label'], r['owner'], tuple(r['seed']), r['eps']))
        if f:
            for key in ('order0', 'order1', 'energy_matching'):
                a, b = r['mu_defect'][key], f['mu_defect'][key]
                time_check = max(time_check, abs(a - b) / b)
    summary['max_roundtrip'] = worst
    summary['max_relative_change_mu_defect_on_step_halving'] = time_check
    groups = {}
    for r in base:
        groups.setdefault((r['label'], r['owner'], tuple(r['seed']), r['lambda']), []).append(r)
    for (lab, owner, seed, lam), rs in groups.items():
        rs.sort(key=lambda r: -r['eps'])
        eps = [r['eps'] for r in rs]
        case = dict(field=lab, owner=owner, seed=list(seed), pitch=lam, eps=eps,
                    rho_larmor=[r['rho_larmor'] for r in rs], fits={})
        for key in ('order0', 'order1', 'energy_matching'):
            v = [r['mu_defect'][key] for r in rs]
            case[key] = v
            case['fits'][key] = dict(all=fit(eps, v), last4=fit(eps[-4:], v[-4:]))
        summary['cases'].append(case)
    (D / 'results' / 'switch_study.json').write_text(json.dumps(summary, indent=1) + '\n')
    print(json.dumps(worst), 'time-check', time_check)
    for c in summary['cases']:
        print(c['field'], c['owner'], c['pitch'],
              ' '.join(f"{k}:{c['fits'][k]['all']:.3f}/{c['fits'][k]['last4']:.3f}"
                       for k in c['fits']))


if __name__ == '__main__':
    main()
