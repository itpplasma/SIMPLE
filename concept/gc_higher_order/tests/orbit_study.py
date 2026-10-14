#!/usr/bin/env python3
"""RK45 full-orbit comparison with phase carried, without phase alignment."""
import ctypes as ct
import csv
import compact_json
import os
import platform
from pathlib import Path
import numpy as np
import scipy
from scipy.integrate import solve_ivp
import sympy as s

ROOT = Path(__file__).resolve().parents[1]
lib = ct.CDLL(str(os.environ.get('GC_LIB',ROOT/'build/libgc.so')))
ptr = np.ctypeslib.ndpointer(dtype=np.float64, flags='F_CONTIGUOUS')
lib.gc_eval_c.argtypes = [ct.c_int,ct.c_int,ct.c_double,ct.c_int,
                         ptr,ptr,ptr,ptr,ptr,ptr,ptr,ct.POINTER(ct.c_int)]
x,y,z = s.symbols('x y z',real=True)
coords = [x,y,z]
R = s.sqrt(x*x+y*y)
fields = {
    'uniform': (s.Matrix([0,0,1]),[0.2,0.1,0.1],1),
    'constant_gradient': (s.Matrix([0,0,1+s.Rational(2,5)*x]),[0.2,0.1,0.1],1),
    'constant_curvature_on_fieldline': (s.Matrix([-y/R,x/R,0]),[1.4,0.0,0.1],3),
    'magnetic_shear': (s.Matrix([0,s.Rational(2,5)*x,1]),[0.2,0.1,0.1],1),
}
# Original sheet_field.f90: B_R=-2*a*Z/R, B_phi=F/R, B_Z=2*a*(R-R0)/R.
# Only smooth branch continuations, never a Taylor series across the sheet.
for label,aval,rminor in [('sheet_inner',s.Rational(1),0.22),
                         ('sheet_outer',s.Rational(13,10),0.78)]:
    br = -2*aval*z/R
    bt = 10/R
    bz = 2*aval*(R-10)/R
    fields[label] = (s.Matrix([br*x/R-bt*y/R,br*y/R+bt*x/R,bz]),
                     [10+rminor,0.0,0.04],3)

class Field:
    def __init__(self,name,entry):
        self.name = name
        expr,self.initial,self.axis = entry
        assert s.simplify(sum(s.diff(expr[i],coords[i]) for i in range(3))) == 0
        self.bfun = s.lambdify(coords,expr,'numpy',cse=True)
        self.dfun = s.lambdify(coords,expr.jacobian(coords),'numpy',cse=True)
        h = [s.diff(expr[i],coords[j],coords[k])
             for i in range(3) for j in range(3) for k in range(3)]
        self.hfun = s.lambdify(coords,h,'numpy',cse=True)
        third = [s.diff(expr[i],coords[j],coords[k],coords[l])
                 for i in range(3) for j in range(3) for k in range(3) for l in range(3)]
        self.tfun = s.lambdify(coords,third,'numpy',cse=True)
    def b(self,q):
        return np.asarray(self.bfun(*q),dtype=float).reshape(3)
    def jet(self,q):
        return (np.asfortranarray(self.b(q)),
                np.asfortranarray(self.dfun(*q),dtype=float),
                np.asfortranarray(np.asarray(self.hfun(*q),dtype=float).reshape(3,3,3)),
                np.asfortranarray(np.asarray(self.tfun(*q),dtype=float).reshape(3,3,3,3)))
    def evaluate(self,mode,n,eps,zz,anchor=None):
        zz = np.asfortranarray(zz,dtype=float)
        q = zz[:3] if anchor is None else np.asarray(anchor,dtype=float)
        b,d,h,third = self.jet(q)
        origin = np.asfortranarray(q)
        out = np.zeros(6,order='F')
        status = ct.c_int()
        lib.gc_eval_c(mode,n,eps,self.axis,origin,b,d,h,third,zz,out,ct.byref(status))
        if status.value:
            raise RuntimeError((self.name,mode,n,status.value,zz))
        return out
    def inverse(self,n,eps,particle):
        anchor = particle[:3].copy()
        for _ in range(6):
            result = self.evaluate(1,n,eps,particle,anchor)
            if np.linalg.norm(result[:3]-anchor) < 1e-12:
                break
            anchor = result[:3].copy()
        residual = np.linalg.norm(self.evaluate(0,n,eps,result)-particle)
        assert residual < 2e-10,(self.name,'refreshed inverse residual',residual)
        return result

# The identical physical reference orbit is reused across selectors. Its
# tolerance and initial state are unchanged; reconstruction remains separate.
full_orbits = {}

def run(field,eps,n,rtol=2e-12,atol=2e-14):
    # Identical physical initial state for every selector, from the n=0 chart.
    seed = np.array([*field.initial,0.35,0.25,0.7])
    particle = field.evaluate(0,0,eps,seed)
    initial = field.inverse(n,eps,particle)
    times = np.linspace(0,1.0,401)
    def full_rhs(t,zz):
        return np.r_[zz[3:],np.cross(zz[3:],field.b(zz[:3]))/eps]
    key = (field,eps,rtol,atol)
    if key not in full_orbits:
        full_orbits[key] = solve_ivp(full_rhs,(0,1),particle,method='RK45',t_eval=times,
                                    rtol=rtol,atol=atol,max_step=eps/3)
    fo = full_orbits[key]
    assert np.array_equal(fo.y[:,0],particle)
    gc = solve_ivp(lambda t,zz:field.evaluate(2,n,eps,zz),(0,1),initial,
                   method='RK45',t_eval=times,rtol=rtol,atol=atol,max_step=0.04)
    assert fo.success and gc.success
    reconstructed = np.array([field.evaluate(0,n,eps,zz) for zz in gc.y.T])
    error = reconstructed-fo.y.T
    def rms(v):
        return float(np.sqrt(np.mean(np.sum(v*v,axis=1))))
    position = rms(error[:,0:3])
    velocity = rms(error[:,3:6])
    energy = np.sum(fo.y[3:]**2,axis=0)/2
    hgc = np.array([field.evaluate(3,n,eps,zz)[0] for zz in gc.y.T])
    if field.name.startswith('sheet'):
        radius = np.hypot(np.hypot(fo.y[0],fo.y[1])-10,fo.y[2])
        if field.name == 'sheet_inner':
            assert np.max(radius) < .5,'inner test crossed sheet'
        else:
            assert np.min(radius) > .5,'outer test crossed sheet'
    return dict(field=field.name,n=n,epsilon=eps,position_rms=position,
                velocity_rms=velocity,phase_space_rms=rms(error),
                fo_energy_range=float(np.ptp(energy)),gc_energy_range=float(np.ptp(hgc)),
                fo_nfev=fo.nfev,gc_nfev=gc.nfev)

def main():
    epsilons = [0.08,0.04,0.02,0.01,0.005]
    rows, fits, convergence = [], [], []
    for name,entry in fields.items():
        field = Field(name,entry)
        for n in [0,1,11,2]:
            sample = [run(field,eps,n) for eps in epsilons]
            rows += sample
            fit = dict(field=name,n=n,expected_position=(3 if n==2 else 2 if n in (1,11) else 1),
                       expected_velocity=(2 if n in (2,11) else n+1))
            for key in ['position_rms','velocity_rms','phase_space_rms']:
                # The last four scales; uniform fields measure solver error.
                fit[key+'_exponent'] = float(np.polyfit(np.log(epsilons[1:]),
                    np.log([row[key] for row in sample[1:]]),1)[0])
            fit['status'] = ('solver_floor' if name=='uniform' else
                ('observed' if fit['position_rms_exponent'] > fit['expected_position']-0.3
                 and fit['velocity_rms_exponent'] > fit['expected_velocity']-0.3 else 'FAILED'))
            fit['selector'] = '1+' if n == 11 else str(n)
            if n in (0,1) and fit['status']=='FAILED':
                fit['status']='baseline_phase_limit'
            fits.append(fit)
            # Independent precision check on the finest scale: RK tolerance / 10.
            tighter = run(field,epsilons[-1],n,rtol=2e-13,atol=2e-15)
            convergence.append(dict(field=name,n=n,
                position_relative_change=abs(tighter['position_rms']/sample[-1]['position_rms']-1),
                velocity_relative_change=abs(tighter['velocity_rms']/sample[-1]['velocity_rms']-1)))
            print(name,fit['selector'],'p(x)=',round(fit['position_rms_exponent'],3),
                  'p(v)=',round(fit['velocity_rms_exponent'],3),fit['status'],flush=True)
    with (ROOT/'results/orbit_errors.csv').open('w') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]),lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    (ROOT/'results/orbit_summary.json').write_text(compact_json.dumps(dict(
        python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,
        epsilon_convention='dv/dt=v cross B/epsilon; normalized speed and gradient scales O(1)',
        time_interval=[0,1],samples=401,rk_tolerance=dict(rtol=2e-12,atol=2e-14),
        finest_tolerance_factor=0.1,fit='log-log last four epsilon values',
        norms='time-sampled RMS Euclidean; no gyroaveraging or phase alignment',
        fits=fits,tolerance_check=convergence,
        n2=dict(status='available',representation='Pi1=-mu curl(b)/2, Pi2=0; zero-mean sigma4'))))

if __name__ == '__main__':
    main()
