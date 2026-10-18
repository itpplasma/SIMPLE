"""Exact symbolic checks of the note's local formulas; no repository edits."""
import sympy as s,json
from pathlib import Path
out={}
def zero(name,expr):
 if isinstance(expr,s.MatrixBase):r=expr.applyfunc(s.simplify);ok=all(a==0 for a in r)
 else:r=s.simplify(expr);ok=r==0
 out[name]={'passed':bool(ok),'residual':str(r)}
 if not ok:raise AssertionError((name,r))
m,B,q=s.symbols('m B q',nonzero=True,real=True)
x,y,z,vx,vy,vz=s.symbols('x y z vx vy vz',real=True)
Z=s.Matrix([x,y,z,vx,vy,vz]);I=s.eye(3);O=s.zeros(3)
Q=s.Matrix([[0,q*B/m**2,0],[-q*B/m**2,0,0],[0,0,0]])
Pi=O.row_join(I/m).col_join((-I/m).row_join(Q))
C=s.Matrix([x+m*vy/(q*B),-q*B*(y-m*vx/(q*B)),z,m*vz,s.atan2(-vy,vx),m**2*(vx**2+vy**2)/(2*q*B)])
J=s.diag(*[s.Matrix([[0,1],[-1,0]])]*3)
D=C.jacobian(Z)
zero('uniform_field_canonical_poisson',D*Pi*D.T-J)
u,w,c,a,r,Bm,Bp=s.symbols('u w c a r Bm Bp',real=True,nonzero=True)
ct=s.symbols('ct',real=True)
up=u*c+w*a*ct;wp2=(w*c*ct-u*a)**2+w**2*(1-ct**2)
zero('direction_reset_speed',s.expand(up**2+wp2-u**2-w**2).subs(c**2,1-a**2))
dmu=m*(u**2+w**2-up**2)/(2*Bp)-m*w**2/(2*Bm)
def avg(poly,flux=False):
 moments={0:1,1:0,2:s.Rational(1,3) if flux else s.Rational(1,2),3:0,4:s.Rational(1,5) if flux else s.Rational(3,8)}
 return s.expand(sum(coeff*moments[int(power[0])] for power,coeff in s.Poly(s.expand(poly),ct).terms()))
mean=m*(u**2*a**2+w**2*(1-a**2/2))/(2*Bp)-m*w**2/(2*Bm)
var=m**2*u**2*w**2*c**2*a**2/(2*Bp**2)+m**2*w**4*a**4/(32*Bp**2)
zero('uniform_phase_mean',s.expand(avg(dmu)-mean).subs(c**2,1-a**2))
zero('uniform_phase_variance',avg((dmu-avg(dmu))**2)-var)
fmean=m*(u**2*a**2+w**2*(1-a**2/s.Integer(3)))/(2*Bp)-m*w**2/(2*Bm)
fvar=m**2*u**2*w**2*c**2*a**2/(3*Bp**2)+m**2*w**4*a**4/(45*Bp**2)
zero('flux_phase_mean',s.expand(avg(dmu,True)-fmean).subs(c**2,1-a**2))
zero('flux_phase_variance',avg((dmu-avg(dmu,True))**2,True)-fvar)
# Saltation in physical coordinates with one-sided magnetic symplectic forms.
bx0,by0,bx1,by1=s.symbols('bx0 by0 bx1 by1',real=True)
v=s.Matrix([vx,vy,vz]);n=s.Matrix([0,0,1]);bm=s.Matrix([bx0,by0,0]);bp=s.Matrix([bx1,by1,0])
def poiss(b):
 bbx,bby,bbz=b
 V=q/m**2*s.Matrix([[0,bbz,-bby],[-bbz,0,bbx],[bby,-bbx,0]])
 return O.row_join(I/m).col_join((-I/m).row_join(V))
omegam=-poiss(bm).inv();omegap=-poiss(bp).inv()
S=s.eye(6);S[3:6,0:3]=q/m*v.cross(bp-bm)*n.T/vz
zero('sheet_saltation_two_form',S.T*omegap*S-omegam)
# Helix and oscillator equations, signed frequency.
t,Om=s.symbols('t Om',real=True,nonzero=True);bb=s.Matrix([0,0,1]);vp=s.Matrix([vx,vy,0]);vv=vz*bb+s.cos(Om*t)*vp+s.sin(Om*t)*vp.cross(bb)
xx=s.Matrix([x,y,z])+vz*t*bb+s.sin(Om*t)/Om*vp+(1-s.cos(Om*t))/Om*vp.cross(bb)
zero('helix_position_derivative',xx.diff(t)-vv)
zero('helix_lorentz_equation',vv.diff(t)-Om*vv.cross(bb))
zero('helix_constant_speed',vv.dot(vv)-v.dot(v))
aeq,vn,omega=s.symbols('aeq vn omega',real=True,nonzero=True)
ss=aeq*(1-s.cos(omega*t))+vn/omega*s.sin(omega*t)
integ=aeq*(t-s.sin(omega*t)/omega)+vn/omega**2*(1-s.cos(omega*t))
zero('normal_oscillator',ss.diff(t,2)+omega**2*(ss-aeq))
zero('normal_integral',integ.diff(t)-ss)
# Ring fraction derivative and FLR small-argument expansions.
d,rho=s.symbols('d rho',real=True,positive=True)
zero('ring_fraction_derivative',s.diff(s.Rational(1,2)+s.asin(d/rho)/s.pi,d)-1/(s.pi*s.sqrt(rho**2-d**2)))
Path(__file__).with_name('algebra_results.json').write_text(json.dumps(out,indent=2))
print('\n'.join(f'{k}: PASS' for k in out));print('Exact checks:',len(out))
