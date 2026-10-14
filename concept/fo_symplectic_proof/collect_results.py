"""Collect bounded numerical evidence from verbose registered make-test logs.

Usage: collect_results.py OUT.json LOG... (the output is not committed).
"""
import hashlib
import json
import re
import sys
from pathlib import Path

root = Path(__file__).resolve().parent
result = {
    'kernel_proof_commit': 'a27d1ec50d32e853b3ea0806e4662ed8d1be5b0c',
    'source_base_commit': '45c0598',
    'reference_commit': '89babb30c0c42504cccd52ec8cfd7ab2d78d38e0',
    'source_sha256': hashlib.sha256((root / 'spectre_fo_gauss2.f90').read_bytes()).hexdigest(),
    'reference_sha256': hashlib.sha256((root / 'reference_gauss2.f90').read_bytes()).hexdigest(),
    'norm': 'Frobenius in normalized canonical q,p coordinates',
    'parameters_frozen': ['E0', 'width', 'ro0', 'length', 'gmin', 'h', 'm_ref'],
    'slope_units': 'relative H error per fictitious step',
    'slope_standard_error': 'OLS on all per-step samples; correlated samples, descriptive only',
    'tokamak': '1000-step sanity smoke only; excluded from stellarator acceptance evidence',
    'symbolic': json.loads((root / 'symbolic.json').read_text()),
    'fixtures': json.loads((root / 'fixtures.json').read_text()),
    'rho_monitor': json.loads((root / 'monitor.json').read_text()),
    'tests': {},
}
float_pattern = r'[-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?'

def pairs(line):
    for label in ('tangent defect', 'FD defect', 'FD disagreement'):
        line = line.replace(label, label.replace(' ', '_'))
    return {key: float(value) for key, value in re.findall(
        r'([^\s=]+)=\s*(' + float_pattern + r')', line)}

for name in sys.argv[2:]:
    state = {}
    for line in Path(name).read_text().splitlines():
        match = re.search(r'Start\s+(\d+):\s+((?:test_fo_symplectic|test_fo_gauss2_(?:resolution_|original_oracle))\S*)', line)
        if match:
            tid, test = match.groups()
            if test == 'test_fo_symplectic_stellarator_fixture':
                continue
            state[tid] = {'test': test, 'fixture': None, 'velocity': None,
                          'rho': None, 'gmin': None}
            result['tests'][test] = {'maps': [], 'energies': [], 'diagnostics': [],
                                     'monitor_derivatives': []}
        match = re.match(r'(\d+): (.*)', line)
        if not match or match[1] not in state:
            continue
        context = state[match[1]]
        line = match[2]
        target = result['tests'][context['test']]
        if line.startswith('velocity multiplier='):
            context['velocity'] = pairs(line)['multiplier']
        elif line.startswith('CASE '):
            context['fixture'] = line.split(maxsplit=2)[2]
        elif line.startswith('MAP configuration'):
            context.update(pairs(line))
        elif line.startswith('MAP reference='):
            record = pairs(line)
            record.update(fixture=context['fixture'], velocity=context['velocity'],
                          rho=context['rho'], gmin=context['gmin'],
                          reference='reference=T' in line)
            target['maps'].append(record)
        elif line.startswith('START marker='):
            context['marker'] = int(pairs(line)['marker'])
            context['start'] = [float(x) for x in re.findall(float_pattern,line.split(' z=')[1])]
        elif line.startswith('PARAM '):
            context['parameters'] = dict(zip(['width','ro0','h','length','gmin','m_ref'],
                map(float,re.findall(float_pattern,line.split('=')[1]))))
        elif line.startswith('ENVELOPE '):
            record = pairs(line)
            record.update(fixture=context['fixture'],velocity=context['velocity'],
                          start=context['start'],parameters=context['parameters'])
            target['energies'].append(record)
        elif line.startswith('MONITOR '):
            record = pairs(line)
            record.update(fixture=context['fixture'],velocity=context['velocity'])
            target['monitor_derivatives'].append(record)
        elif line.startswith(('NEWTON_', 'ENERGY_FAILURE', 'FAILURE_KERNEL',
                              'LOOSE_TOLERANCE_MAP', 'ERROR STOP')):
            target['diagnostics'].append({'fixture':context['fixture'],
                                          'velocity':context['velocity'],'text':line})
    for line in Path(name).read_text().splitlines():
        match = re.search(r'Test\s+#\d+:\s+((?:test_fo_symplectic|test_fo_gauss2_(?:resolution_|original_oracle))\S*)\s+(.+?)(\d+\.\d+) sec',line)
        if match and match[1] in result['tests']:
            result['tests'][match[1]]['outcome'] = match[2].lstrip('.').strip()
            result['tests'][match[1]]['wall_seconds'] = float(match[3])

Path(sys.argv[1]).write_text(json.dumps(result,indent=2) + '\n')
print('Collected',len(result['tests']),'registered tests,',
      sum(len(t['maps']) for t in result['tests'].values()),'maps,',
      sum(len(t['energies']) for t in result['tests'].values()),'energy runs')
