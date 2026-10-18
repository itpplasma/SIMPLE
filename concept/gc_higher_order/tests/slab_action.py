#!/usr/bin/env python3
"""Independent second-order slab magnetic action, derived in prove_identities.py."""
import compact_json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp

ROOT=Path(__file__).resolve().parents[1]
g=0.4
scales=[0.08,0.04,0.02,0.01,0.005]
rows=[]
for ep in scales:
    initial=np.array([0.2,0.1,0.1,0.5,0.4,0.35])
    def rhs(t,z):
        B=1+g*z[0]
        return [*z[3:],z[4]*B/ep,-z[3]*B/ep,0]
    sol=solve_ivp(rhs,(0,1),initial,t_eval=np.linspace(0,1,801),
                  rtol=2e-12,atol=2e-14,max_step=ep/4)
    assert sol.success
    B=1+g*sol.y[0]
    E=(sol.y[3]**2+sol.y[4]**2)/2
    mu0=E/B
    mu1=-g*E*sol.y[4]/B**3
    mu2=g*g/B**5*(1.5*E*sol.y[4]**2+0.75*E**2)
    vals=[mu0,mu0+ep*mu1,mu0+ep*mu1+ep**2*mu2]
    for n,value in enumerate(vals):
        rows.append(dict(epsilon=ep,n=n,rms=float(np.sqrt(np.mean((value-value[0])**2)))))
fits=[]
for n in range(3):
    selected=[r for r in rows if r['n']==n][1:]
    exponent=float(np.polyfit(np.log([r['epsilon'] for r in selected]),
                             np.log([r['rms'] for r in selected]),1)[0])
    assert exponent>n+0.8,(n,exponent)
    fits.append(dict(n=n,exponent=exponent,expected=n+1))
    print('slab magnetic action',n,'exponent',round(exponent,4))
(ROOT/'results/slab_action.json').write_text(compact_json.dumps(dict(
    scope='constant-direction linear-magnitude field only; magnetic action, not orbit tracking',
    formula='J=E/B - epsilon*g*E*vy/B**3 + epsilon**2*g**2*(3*E*vy**2/2+3*E**2/4)/B**5',
    rows=rows,fits=fits)))
