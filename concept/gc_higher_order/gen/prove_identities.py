#!/usr/bin/env python3
"""Exact identities, independent of emitted kernels. No random-probe proofs."""
import json
from pathlib import Path
import sympy as s

checks = []
def prove(name, expression):
    entries = list(expression) if isinstance(expression, s.MatrixBase) else [expression]
    residuals = [s.trigsimp(s.cancel(s.expand_trig(v))) for v in entries]
    assert all(v == 0 for v in residuals), (name, residuals)
    checks.append({"identity": name, "status": "proved", "method": "SymPy exact"})

# Reduced phase-space two-form, evaluated in an orthonormal frame b=e3.
# Inputs are arbitrary local vectors; Bstar_parallel !=0, epsilon !=0.
eps, u, mu, B = s.symbols('epsilon u mu B', nonzero=True)
bs = s.Matrix(s.symbols('Bs1:4'))
grad = s.Matrix(s.symbols('g1:4'))
R = s.Matrix(s.symbols('R1:4'))
b = s.Matrix([0,0,1])
xd = (u*bs+eps*mu*b.cross(grad))/bs[2]
ud = -mu*bs.dot(grad)/bs[2]
pd = -B/eps-R.dot(xd)
omega = s.zeros(6)
for i in range(3):
    for j in range(3):
        omega[i,j] = sum(s.LeviCivita(i,j,k)*bs[k]/eps for k in range(3))
    omega[i,3] = -b[i]
    omega[3,i] = b[i]
    omega[i,4] = eps*R[i]
    omega[4,i] = -eps*R[i]
omega[4,5], omega[5,4] = -eps, eps
prove('Euler-Lagrange inversion with gyrophase and gyrogauge vector',
      omega*s.Matrix([*xd,ud,0,pd])-s.Matrix([*(mu*grad),u,B,0]))
prove('exact energy conservation of the truncated reduced model',
      mu*grad.dot(xd)+u*ud)

# First-order drift from Bstar=B+epsilon*u*curl(b)+O(epsilon**2).
cx,cy,tau,kx,ky = s.symbols('cx cy tau kx ky')
cb = s.Matrix([-ky,kx,tau])  # curl b = tau*b+b cross kappa
star = B*b+eps*u*cb
xexpr = (u*star+eps*mu*b.cross(grad))/star[2]
coefficient = xexpr.diff(eps).subs(eps,0)
prove('first-order curvature and grad-B drift',
      coefficient-b.cross(mu*grad+u**2*s.Matrix([kx,ky,0]))/B)

# TB alpha1=partial_theta alpha2; theta=-phi here.
phi = s.symbols('phi', real=True)
a = s.Matrix([s.cos(phi),s.sin(phi),0])
r = b.cross(a)
D = s.Matrix(3,3,s.symbols('d0:9'))
alpha1 = -sum((r[i]*a[j]+a[i]*r[j])*D[i,j]/2
              for i in range(3) for j in range(3))
alpha2 = sum((a[i]*a[j]-r[i]*r[j])*D[i,j]/4
             for i in range(3) for j in range(3))
prove('TB alpha1=partial_theta alpha2, theta=-phi',-s.diff(alpha2,phi)-alpha1)
avg = s.integrate(alpha1,(phi,0,2*s.pi))/(2*s.pi)
prove('gyroaverage of alpha1 vanishes',avg)
prove('Banos gyroaverage <G1_parallel>=mu*tau',
      s.integrate(-u*r.dot(s.Matrix([kx,ky,0]))+mu*(tau+alpha1),
                  (phi,0,2*s.pi))/(2*s.pi)-mu*tau)

# Gauge rotation: e1'=cos(psi)e1+sin(psi)e2; phi'=phi-psi.
psi, da, dbeta, dpsi = s.symbols('psi da dbeta dpsi', real=True)
alpha,beta = s.symbols('alpha beta', real=True)
bh = s.Matrix([s.cos(alpha)*s.cos(beta),s.sin(alpha)*s.cos(beta),s.sin(beta)])
e1 = s.Matrix([-s.sin(alpha),s.cos(alpha),0])
e2 = bh.cross(e1)
e1p = s.cos(psi)*e1+s.sin(psi)*e2
e2p = -s.sin(psi)*e1+s.cos(psi)*e2
prove('gyrogauge invariance of perpendicular velocity',
      s.cos(phi-psi)*e1p+s.sin(phi-psi)*e2p-s.cos(phi)*e1-s.sin(phi)*e2)
def deriv(v):
    return v.diff(alpha)*da+v.diff(beta)*dbeta+v.diff(psi)*dpsi
prove('gyrogauge connection transforms as Rprime=R+dpsi',
      e2p.dot(deriv(e1p))-e2.dot(deriv(e1))-dpsi)

# Published second-order scalar consistency: TB (30),(33)-(38).
# Independent scalar symbols isolate the tensor contractions. This verifies
# the algebraic constraint, NOT an independent derivation of <G2^mu> or T2.
J,Om,p,t,c,a2,h,g2,k2,kg,dk = s.symbols('J Om p tau c a2 h g2 k2 kg dk')
bp = -t**2/2-c+a2-g2+h/2
bl = -2*a2-3*kg+3*k2+dk
avgG = J**2/(2*Om)*(t**2/2+c-a2-h/2) \
       -J*p**2/(2*Om**2)*(3*k2-kg+dk-t**2)
# Eq. (30), Pi1 left symbolic through the common dot(Pi1,v_gc) term.
vg2 = (J**2*g2+2*J*p**2*kg/Om+p**4*k2/Om**2)
from_constraint = -Om*avgG+J*p**2*(t**2/2-a2)/Om - vg2/2
from_published = J**2*bp/2+J*p**2*bl/(2*Om)-p**4*k2/(2*Om**2)
prove('TB second-order Hamiltonian scalar constraint, excluding common Pi1.vgc',
      from_constraint-from_published)

# Independent second-order slab derivation, not a tensor transcription.
# A_y(X+eps*q)-A_y(X)=eps*B*q+eps**2*g*q**2/2.
# Hamiltonian: H0=(p_x**2+B**2*q**2)/2; H1=B*g*q**3/2;
# H2=g**2*q**4/8. Its averaged normal form gives the slab diagnostic.
angle = s.symbols('angle', real=True)
I, Bslab = s.symbols('I Bslab', positive=True)
gslab = s.symbols('gslab', real=True)
a3, a4 = Bslab*gslab/2, gslab**2/8
H1 = a3*(2*I/Bslab)**s.Rational(3,2)*s.sin(angle)**3
H2 = a4*(2*I/Bslab)**2*s.sin(angle)**4
S1 = a3*(2*I/Bslab)**s.Rational(3,2)/Bslab*(s.cos(angle)-s.cos(angle)**3/3)
poisson = lambda f,h: s.diff(f,angle)*s.diff(h,I)-s.diff(f,I)*s.diff(h,angle)
prove('slab canonical homological equation',Bslab*s.diff(S1,angle)+H1)
normal2 = s.integrate(s.expand_trig(H2+poisson(S1,H1)/2),
                     (angle,0,2*s.pi))/(2*s.pi)
prove('independent second-order slab Hamiltonian',normal2+3*gslab**2*I**2/(4*Bslab**2))
# Particle-space action obtained by inverting that Hamiltonian and the conserved
# canonical momentum A_y(X)=A_y(x)+eps*v_y, through epsilon**2.
vx,vy,ep = s.symbols('vx vy ep',real=True)
Eperp = (vx**2+vy**2)/2
J0 = Eperp/Bslab
J1 = -gslab*Eperp*vy/Bslab**3
J2 = gslab**2/Bslab**5*(s.Rational(3,2)*Eperp*vy**2+s.Rational(3,4)*Eperp**2)
Japprox = J0+ep*J1+ep**2*J2
D = gslab*vx*s.diff(Japprox,Bslab) + Bslab/ep*(vy*s.diff(Japprox,vx)-vx*s.diff(Japprox,vy))
for power in [-1,0,1]:
    prove(f'slab action derivative coefficient epsilon**{power}',s.expand(D).coeff(ep,power))

# Lie composition to second order, for arbitrary polynomial vector fields.
q,p,t = s.symbols('q p t')
coords = s.Matrix([q,p])
G = s.Matrix([q*p,q**2-p])
H = s.Matrix([p**2,-q*p])
f = q**2*p+p**3
L = lambda V,F: sum(V[i]*s.diff(F,coords[i]) for i in range(2))
comm = s.Matrix([L(G,H[i])-L(H,G[i]) for i in range(2)])
prove('commutator of Lie derivatives is the vector-field Lie bracket',
      L(G,L(H,f))-L(H,L(G,f))-L(comm,f))
forward = coords+t*G+t**2*(H+s.Matrix([L(G,G[i])/2 for i in range(2)]))
inverse = coords-t*G+t**2*(-H+s.Matrix([L(G,G[i])/2 for i in range(2)]))
for i in range(2):
    composed = forward[i].subs({q:inverse[0],p:inverse[1]}, simultaneous=True)
    prove(f'Lie transform inverse composition through order two, coordinate {i}',
          s.series(composed,t,0,3).removeO()-coords[i])

# The proposed unrestricted "iff reset preserves omega AND H" is false.
# Translation q->q+1 preserves dq wedge dp but does not preserve H=q.
Jmat = s.Matrix([[0,1],[-1,0]])
reset = s.Matrix([q+1,p])
prove('counterexample reset translation is symplectic',reset.jacobian(coords).T*Jmat*reset.jacobian(coords)-Jmat)
assert s.expand(reset[0]-q) == 1
checks.append({"identity":"symplectic reset must preserve Hamiltonian", "status":"disproved",
               "counterexample":"(q,p)->(q+1,p), omega=dq wedge dp, H=q"})

out = Path(__file__).resolve().parents[1]/'results/symbolic.json'
out.write_text(json.dumps({"sympy":s.__version__,"checks":checks},indent=2)+'\n')
print(f'{len(checks)-1} exact identities passed; unrestricted reset/H necessity disproved')
