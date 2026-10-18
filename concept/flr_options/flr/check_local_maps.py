"""Independent smooth-side DOP853 and finite-difference checks of exact sheet maps."""
import numpy as np,json
from scipy.integrate import solve_ivp
from pathlib import Path
bm=np.array([1.,0.,0.]);bp=1.3*np.array([.8,.6,0.]);n=np.array([0.,0.,1.]);m=q=1.
def helix(z,B,t):
 x,v=z[:3],z[3:];b=B/np.linalg.norm(B);om=q*np.linalg.norm(B)/m
 u=np.dot(v,b);w=v-u*b;v1=u*b+np.cos(om*t)*w+np.sin(om*t)*np.cross(w,b)
 x1=x+u*b*t+np.sin(om*t)*w/om+(1-np.cos(om*t))*np.cross(w,b)/om
 return np.r_[x1,v1]
def next_time(z,B):
 b=B/np.linalg.norm(B);om=q*np.linalg.norm(B)/m;w=abs(om)
 a=z[5]/w;k=np.dot(n,np.cross(z[3:],b))/om;c=z[2]+k;r=np.hypot(a,k)
 if r==0 or abs(c)>r*(1+1e-14):return np.inf
 delta=np.arctan2(a,-k);alpha=np.arccos(np.clip(-c/r,-1,1));times=[]
 for th in [delta+alpha,delta-alpha]:
  th=th%(2*np.pi)
  if th<1e-10:th=2*np.pi
  times.append(th/w)
 return min(times)
def flow(z,T):
 z=z.copy();remain=T;count=0
 side=1 if z[2]>1e-12 else -1 if z[2]<-1e-12 else np.sign(z[5])
 while remain>1e-13:
  B=bp if side>0 else bm;t=next_time(z,B)
  if t>=remain-1e-12:return helix(z,B,remain),count
  z=helix(z,B,t);z[2]=0.;remain-=t;side=-side;count+=1
  if count>10000:raise RuntimeError('event count')
 return z,count
def form(B):
 bx,by,bz=B;C=q*np.array([[0,bz,-by],[-bz,0,bx],[by,-bx,0.]])
 return np.block([[-C,m*np.eye(3)],[-m*np.eye(3),np.zeros((3,3))]])
z0=np.array([.1,-.2,0.,.4,.25,.6]);z=z0.copy();errs=[];dur=[]
for B in [bp,bm]:
 t=next_time(z,B);ze=helix(z,B,t)
 def rhs(t,z):return np.r_[z[3:],q/m*np.cross(z[3:],B)]
 sol=solve_ivp(rhs,[0,t],z,method='DOP853',rtol=2e-13,atol=2e-14,max_step=.035)
 errs.append(float(np.max(np.abs(sol.y[:,-1]-ze))));dur.append(float(t));z=ze;z[2]=0.
cycle_disp=z[:2]-z0[:2];Tcycle=sum(dur)
zseq=z0.copy()
for i in range(200):zseq,_=flow(zseq,Tcycle)
zskip=z0.copy();zskip[:2]+=200*cycle_disp
# Each cycle returns the same velocity, so the direct composition is an exact oracle.
cycle_err=float(np.max(np.abs(zseq-zskip)))
za=helix(z0,bp,.17);T=3.3;zf,count=flow(za,T);h=2e-6
D=np.column_stack([(flow(za+h*np.eye(6)[j],T)[0]-flow(za-h*np.eye(6)[j],T)[0])/(2*h) for j in range(6)])
Bfinal=bp if zf[2]>0 else bm
symp=float(np.max(np.abs(D.T@form(Bfinal)@D-form(bp))))
mu0=(np.dot(z0[3:],z0[3:])-np.dot(z0[3:],bm)**2)/(2*np.linalg.norm(bm))
v=z0[3:];b=bp/np.linalg.norm(bp);mup=(np.dot(v,v)-np.dot(v,b)**2)/(2*np.linalg.norm(bp))
out={'smooth_side_DOP853_max_errors':errs,'cycle_times':dur,'cycle_tangent_displacement':cycle_disp.tolist(),'cycle_velocity_return_error':float(np.max(np.abs(z[3:]-z0[3:]))),'200_cycles_direct_vs_compressed_max_error':cycle_err,'fixed_time_map_event_count':count,'finite_difference_symplectic_defect':symp,'energy_relative_error':float(abs(np.dot(zf[3:],zf[3:])-np.dot(za[3:],za[3:]))/np.dot(za[3:],za[3:])),'example_mu_minus':float(mu0),'example_mu_plus':float(mup)}
Path(__file__).with_name('local_map_results.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
assert max(errs)<1e-10 and cycle_err<1e-9 and symp<2e-7
