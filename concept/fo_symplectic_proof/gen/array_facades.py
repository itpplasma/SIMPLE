"""Bind the checked FortSym scalar DAGs directly to the consumer's array ABI.

Only symbol indexing, declarations and continuation layout change. Arithmetic
and statement order remain the emitted FortSym order. Normal builds need no CAS.
"""
from pathlib import Path
import re
import sys
import textwrap

root = Path(sys.argv[1])
source = (root / 'fo_kernel.f90').read_text()
indices = {f'p{i}': f'p({i})' for i in range(1, 4)}
for prefix, count, target in [('a', 3, 'aj'), ('g', 6, 'gj')]:
    for i in range(1, count + 1):
        for j in range(1, 11):
            indices[f'{prefix}{i}_{j}'] = f'{target}({i},{j})'
for j in range(1, 11):
    indices[f't1_{j}'] = f'tj({j})'
for prefix, count in [('gh', 6), ('f', 6), ('inverse', 6), ('velocity', 3)]:
    for i in range(1, count + 1):
        indices[f'{prefix}{i}'] = f'{prefix}({i})'
for prefix in ['hh', 'df']:
    for i in range(1, 7):
        for j in range(1, 7):
            indices[f'{prefix}{i}_{j}'] = f'{prefix}({i},{j})'


def body(leaf):
    scalar = (root / f'{leaf}.inc').read_text()
    result = scalar[scalar.index('    real(dp) :: work'):scalar.rindex('end subroutine')]
    return re.sub(r'\b\w+\b', lambda m: indices.get(m[0], m[0]), result)


def interface(name):
    return f'''subroutine {name}(p,aj,gj,mass,charge,e0,tj,energy,kenergy,gh,hh,f,df)
real(dp),intent(in)::p(3),aj(3,10),gj(6,10),mass,charge,e0,tj(10)
real(dp),intent(out)::energy,kenergy,gh(6),hh(6,6),f(6),df(6,6)
'''


for name, leaf in [('fo_eval_first', 'fo_scalar_first')]:
    replacement = body(leaf)
    if name == 'fo_eval_first':
        index = replacement.index('\n    work1 =')
        replacement = replacement[:index] + '\n    hh=0; df=0\n' + replacement[index:]
    replacement = interface(name) + replacement + f'end subroutine {name}'
    pattern = rf'^\s*subroutine {name}\(.*?^\s*end subroutine {name}'
    source, count = re.subn(pattern, lambda _: replacement, source, flags=re.M | re.S)
    if count != 1:
        raise ValueError(f'expected one {name} facade, found {count}')

# The monitor is differentiated after H: Hess(T(H-E0)) follows the product rule.
# The velocity/Hessian/flow leaves and their readback are proved by FortSym.
factored = interface('fo_eval_factored') + '''real(dp)::inverse(6),velocity(3)
call fo_eval_first(p,aj,gj,mass,charge,e0,tj,energy,kenergy,gh,hh,f,df)
'''
for leaf in ['fo_factor_velocity', 'fo_factor_hessian', 'fo_factor_flow']:
    # Flow's scalar t1..t10 are the monitor inputs, not CSE temporaries.
    if leaf == 'fo_factor_flow':
        indices.update({f't{i}': f'tj({i})' for i in range(1, 11)})
    factored += 'block\n' + body(leaf) + 'end block\n'
factored += 'end subroutine fo_eval_factored\n'
source = source.replace('end module fo_kernel', factored + 'end module fo_kernel')

# Array indexing lengthens tokens. Reflow complete logical lines before formatting.
logical = []
pending = ''
for line in source.splitlines():
    text = line.strip()
    pending += (' ' if pending else '') + text.removesuffix('&').lstrip('&').strip()
    if not text.endswith('&'):
        logical.append(pending)
        pending = ''
if pending:
    raise ValueError('unterminated generated continuation')
lines = []
for line in logical:
    if line.startswith('!'):
        lines.append(line)
        continue
    wrapped = textwrap.wrap(line, width=68, break_long_words=False, break_on_hyphens=False)
    wrapped = wrapped or ['']
    lines.extend(text + (' &' if i < len(wrapped) - 1 else '')
                 for i, text in enumerate(wrapped))
(root / 'fo_kernel.f90').write_text('\n'.join(lines) + '\n')
