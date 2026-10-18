#!/usr/bin/env python3
"""Exact local one-form and inverse-map identities, independent of time stepping."""
import json
from pathlib import Path
import sympy as s
from derive_canonical import physics, r, u, mu, th, ph, eps

checks = []
def prove(label, expression):
    assert s.simplify(expression) == 0, label
    checks.append(label)

# Arbitrary smooth local chart, no assumed axisymmetry or magnetic geometry.
zeta, J = s.symbols('zeta J', real=True)
q = [r, th, ph]
A1, A2, b1, b2, Rr, R1, R2 = [s.Function(n)(*q)
    for n in ['A1', 'A2', 'b1', 'b2', 'Rr', 'R1', 'R2']]
K = s.Function('K')(*q)
phi = -zeta-K
P = [A1/eps+u*b1-J*(R1-s.diff(K, th)),
     A2/eps+u*b2-J*(R2-s.diff(K, ph))]
# Pull back alpha_gc to (r,th,ph,u,J,zeta), retaining variable action.
original = s.Matrix([-J*Rr, A1/eps+u*b1-J*R1,
                      A2/eps+u*b2-J*R2, 0, 0, 0])
variables = [r, th, ph, u, J, zeta]
original -= J*s.Matrix([s.diff(phi, v) for v in variables])
canonical = s.Matrix([0, P[0], P[1], 0, 0, J])
for i, v in enumerate(variables):
    prove('one-form coefficient '+str(v), (original-canonical)[i].subs(s.diff(K, r), Rr))
for i in range(6):
    for j in range(i):
        difference = s.diff((original-canonical)[i], variables[j])-s.diff((original-canonical)[j], variables[i])
        prove(f'extended two-form {i},{j}', difference.subs(Rr, s.diff(K, r)))

# Leading implicit near-identity change at fixed (q,P,J), not fixed u.
a1, a2, h1, h2, c1, c2 = s.symbols('a1 a2 h1 h2 c1 c2')
D = a1*h2-a2*h1
delta = s.Matrix([eps**2*mu*(c1*h2-c2*h1)/D,
                   eps*mu*(a1*c2-a2*c1)/D])
M = s.Matrix([[a1/eps, h1], [a2/eps, h2]])
for i, value in enumerate(M*delta-eps*mu*s.Matrix([c1, c2])):
    prove('fixed canonical momentum displacement '+str(i), value)

# A concrete curvature obstruction: constant magnetic shear, smooth on R^3.
k = s.symbols('k', real=True, nonzero=True)
B = s.sqrt(1+k*k*r*r)
rs_y, rs_z = k*k*r/(2*B**3), k/(2*B**3)
curl_rs_z = s.diff(rs_y, r)
prove('connection curvature at r=0 is nonzero', curl_rs_z.subs(r, 0)-k*k/2)
assert curl_rs_z.subs(r, 0) != 0
checks.append('scalar-gauge removal of full connection is obstructed')

# Exact sheet primitive: substitute v=tanh(asinh(r/c)/2) to rationalize.
c, v, t = s.symbols('c v t', positive=True)
ct = s.symbols('ct', real=True)
d = s.sqrt(100+c*c*ct*ct)
primitive = (s.atanh((10*v-c*ct)/d)-s.atanh(-c*ct/d))/d
prove('rational radial primitive', s.diff(primitive, v)-1/(10+2*c*ct*v-10*v*v))
rr = 2*c*v/(1-v*v)
prove('sheet radial phase primitive', c*s.diff(primitive, v)/s.diff(rr, v) -
      c/(2*(10+rr*ct)*(c*(1+v*v)/(1-v*v))))

# Analytic A_theta uses 20 terms of a geometric series in r cos(theta)/10.
# For r<=1, every omitted derivative through order three has a bound <1e-15.
remainder = sum((n+2)**3*s.Rational(1, 10)**n for n in range(20, 100))
assert remainder < s.Rational(1, 10**15)
checks.append('sheet vector potential first three derivative remainder <1e-15')

# Slab independently known from the exact action/energy derivation of D2.
_, _, _, _, slab_h2 = physics('constant_gradient')
prove('canonical Hamiltonian retains exact slab H2', slab_h2+
      3*s.Rational(2, 5)**2*mu**2/(4*(1+s.Rational(2, 5)*r)**2))

result = dict(status='PASS', exact_identities=checks,
    scope='local chart; nonzero magnetic Jacobian; mu>0; sheet 0<r<=1',
    no_claim='global Darboux chart, arbitrary-order GC or uniform orbit bounds')
target = Path(__file__).resolve().parents[1]/'results/canonical_symbolic.json'
target.write_text(json.dumps(result, indent=2)+'\n')
print('PASS:', len(checks), 'canonical identities/obstructions')
