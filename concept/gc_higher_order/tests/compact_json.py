"""Write measurement JSON with one record per line to keep diffs reviewable."""
import json


def dumps(obj):
    """Top-level keys on separate lines; each dict in a list of dicts on one line."""
    lines = ['{']
    items = list(obj.items())
    for i, (key, value) in enumerate(items):
        sep = ',' if i < len(items)-1 else ''
        if isinstance(value, list) and value and all(isinstance(x, dict) for x in value):
            lines.append(f'  {json.dumps(key)}: [')
            lines += ['    '+json.dumps(x)+(',' if j < len(value)-1 else '')
                      for j, x in enumerate(value)]
            lines.append('  ]'+sep)
        else:
            lines.append(f'  {json.dumps(key)}: {json.dumps(value)}{sep}')
    lines.append('}')
    return '\n'.join(lines)+'\n'
