#!/usr/bin/env python3
"""Independent exact slab checks of BT 2016 fourth-order phase closure."""
import json
from pathlib import Path
import sympy as s

X, Y, Z, u, mu, phi, ep = s.symbols('X Y Z u mu phi epsilon', real=True)
B0, g = s.symbols('B0 g', real=True)
B = B0 + g*X
w = s.sqrt(2*mu*B)
z = [X,Y,Z,u,mu,phi]
b = s.Matrix([0,0,1])
vt = s.Matrix([w*s.cos(phi),w*s.sin(phi),0])
vel = u*b+vt
rho = b.cross(vt)/B
# TB (41)-(45), specialized before differentiating; frame connection=0.
G1 = s.Matrix([*(-rho),0,mu*g*rho[0]/B,-g*vt[0]/B**2])
G2x = g*vt[0]*vt/(2*B**3)
H2 = -3*g**2*mu**2/(4*B**2)
checks = []
def simp(v):
    return s.trigsimp(s.expand_trig(s.simplify(v)))
def prove(name, v):
    entries = list(v) if isinstance(v,s.MatrixBase) else [v]
    residuals = [simp(a) for a in entries]
    assert all(a == 0 for a in residuals),(name,residuals)
    print('proved:',name,flush=True)
    checks.append(dict(identity=name,status='proved',method='SymPy exact'))
def L(v):
    if isinstance(v,s.MatrixBase):
        return v.applyfunc(L)
    return simp(sum(G1[j]*s.diff(v,z[j]) for j in range(6)))
def curl(v):
    return s.Matrix([s.diff(v[2],Y)-s.diff(v[1],Z),
                     s.diff(v[0],Z)-s.diff(v[2],X),
                     s.diff(v[1],X)-s.diff(v[0],Y)])
def D1(v):
    return (sum((v.diff(z[j])*G1[j] for j in range(3,6)),s.zeros(3,1))
            + rho.cross(curl(v))).applyfunc(simp)
def da(v,a):
    out = s.diff(v,z[a+3])
    return -out if a == 2 else out

d1p2 = D1(u*b+vt/2)
d1p4 = D1(u*b/3+vt/4)
d11p3 = D1(D1(u*b/2+vt/3))
d11p4 = D1(d1p4)
G2u = simp(b.dot(d11p3))
G2mu = simp(-mu*g*G2x[0]/B-u*G2u/B-H2/B)
F = [simp(da(d1p4,a).dot(rho)) for a in range(3)]
Gpaper = [G1[3],G1[4],-G1[5]]
C = [simp(L(F[a])-sum(Gpaper[m]*da(F[m],a) for m in range(3)))
     for a in range(3)]
gt = -G1[5]
Q = simp(G2mu+(G2x.dot(da(d1p2,2))+rho.dot(da(d11p4,2))+C[2])/2
         -(G1[4]*da(gt,2)-gt*da(G1[4],2)+L(G1[4]))/4)
prove('slab fourth-order sigma4 solvability',s.integrate(Q,(phi,0,2*s.pi)))
sigma = simp(-s.integrate(s.expand_trig(Q),phi))
sigma -= s.integrate(sigma,(phi,0,2*s.pi))/(2*s.pi)
prove('slab periodic homological equation partial_theta sigma4=Q',
      -s.diff(sigma,phi)-Q)
G2theta = simp(-s.diff(sigma,mu)
    +(G2x.dot(da(d1p2,1))+rho.dot(da(d11p4,1))+C[1])/2
    -(G1[4]*da(gt,1)-gt*da(G1[4],1)-L(gt))/4)
G3par = simp(-s.diff(sigma,u)
    +(G2x.dot(da(d1p2,0))+rho.dot(da(d11p4,0))+C[0])/2
    -(G1[4]*da(gt,0)-gt*da(G1[4],0))/4)
G3x = G3par*b-b.cross(d11p3)/B
G2 = s.Matrix([*G2x,G2u,G2mu,-G2theta])
L2 = lambda v: sum((v.diff(z[j])*G2[j] for j in range(6)),s.zeros(3,1))
rho1 = (-G2x-L(rho)/2).applyfunc(simp)
rho2 = (-G3x-L2(rho)+L(L(rho))/6).applyfunc(simp)
v1 = -L(vel)
v2 = (-L2(vel)+L(L(vel))/2).applyfunc(simp)
xp = s.Matrix([X,Y,Z])+ep*rho+ep**2*rho1+ep**3*rho2
vp = vel+ep*v1+ep**2*v2
# Oracles derive from the original particle equations, not generator text.
energy = s.expand(vp.dot(vp)/2-u**2/2-mu*B-ep**2*H2)
for k in range(3):
    prove(f'slab particle energy pullback coefficient {k}',energy.coeff(ep,k))
Ay = lambda q: B0*q+g*q**2/2
momentum = s.expand(Ay(xp[0])+ep*vp[1]-Ay(X))
for k in range(4):
    prove(f'slab canonical momentum pullback coefficient {k}',momentum.coeff(ep,k))
# Independent particle action, derived from energy and canonical momentum.
E0 = simp((vel[0]**2+vel[1]**2)/2)
E1 = simp(vel[0]*v1[0]+vel[1]*v1[1])
E2 = simp(vel[0]*v2[0]+vel[1]*v2[1]+(v1[0]**2+v1[1]**2)/2)
dB1, dB2 = g*rho[0], g*rho1[0]
Jcoef = E0/B + ep*(E1/B-E0*dB1/B**2-g*E0*vel[1]/B**3) + ep**2*(
    E2/B-E1*dB1/B**2+E0*(dB1**2/B**3-dB2/B**2)
    -g*(E1*vel[1]+E0*v1[1])/B**3+3*g*E0*vel[1]*dB1/B**4
    +g**2*(3*E0*vel[1]**2/2+3*E0**2/4)/B**5)
Jcoef = s.expand(Jcoef)
for k in range(3):
    prove(f'independent slab magnetic action pullback coefficient {k}',
          Jcoef.coeff(ep,k)-(mu if k == 0 else 0))
# Symbolic proof of the finite Fourier primitive used by gc_second_order.
# Orthogonality of sampled modes follows from the finite geometric sum:
# sum_{j=0}^{15} r^(m*j)=0 for 0<|m|<=12, r=exp(2*pi*i/16).
r = s.symbols('r')
cyclotomic = s.cyclotomic_poly(16,r)
for m in range(1,13):
    prove(f'16-point Fourier orthogonality frequency {m}',
          s.rem(sum(r**((m*j)%16) for j in range(16)),cyclotomic,r))
# It then recovers every product of modes <=6 (sum/difference <=12).
for k in range(1,7):
    aa,bb=s.symbols('a b')
    primitive=(bb*s.cos(k*phi)-aa*s.sin(k*phi))/k
    prove(f'negative Fourier primitive mode {k}',
          s.diff(primitive,phi)+aa*s.cos(k*phi)+bb*s.sin(k*phi))
root = Path(__file__).resolve().parents[1]
(root/'results/second_order_symbolic.json').write_text(json.dumps(dict(
    sympy=s.__version__,checks=checks,slab=dict(
        sigma4=str(simp(sigma)),G2_theta=str(G2theta),G3_parallel=str(G3par),
        rho2=[str(v) for v in rho2],v2=[str(v) for v in v2])),indent=2)+'\n')
print(f'{len(checks)} exact second-order slab/Fourier identities passed')
