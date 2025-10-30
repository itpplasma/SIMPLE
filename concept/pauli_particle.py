# %%
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import root
from scipy.interpolate import lagrange
import common
from common import (f,r0,th0,ph0,pph0,timesteps,get_val,get_der,mu,qe,m,c,)
from plotting import plot_orbit, plot_cost_function

dt, nt = timesteps(steps_per_bounce=8, nbounce=100)

z = np.zeros([3, nt + 1])
z[:, 0] = [r0, th0, ph0]
[H, pth, vpar] = get_val(np.array([r0, th0, ph0, pph0]))

'''
def newton(z, vpar, v):
    r = z[0]
    th = z[1]
    e = 1

    while e >= 1e-11:
        f.evaluate(r,th,z[2])
        g = m*vpar*f.hth + qe/c*f.Ath
        dg = m*vpar*f.dhth[0] + qe/c*f.dAth[0]
        rnew = r - g/dg
        r = rnew
        e = g
        return rnew
'''

calc_fields = lambda vpar: {
    "A2": m * vpar * f.hth + qe * f.Ath / c,
    "dA2": m * vpar * f.dhth + qe * f.dAth / c,
    "d2A2": m * vpar * f.d2hth + qe * f.d2Ath / c,
    "A3": m * vpar * f.hph + qe * f.Aph / c,
    "dA3": m * vpar * f.dhph + qe * f.dAph / c,
    "d2A3": m * vpar * f.d2hph + qe * f.d2Aph / c,
    "dB": mu * f.dB,
    "d2B": mu * f.d2B,
    "dPhie": f.dPhie,
    "d2Phie": f.d2Phie,
}

def F(z, vpar, v):
    ret = np.zeros(2)

    f.evaluate(z[0],z[1],0)

    ret[0] = m*vpar*f.hth + qe/c *f.Ath - dt*v[0]
    ret[1] = m*vpar*f.hph + qe/c *f.Aph - dt*v[1]
    return ret



from time import time
tic = time()
for kt in range(nt):
    vparold = vpar
    f.evaluate(z[0, kt], z[1, kt], z[2, kt])
    cov = calc_fields(vparold)
    v1 = (qe*cov['dPhie'][2] - mu*cov['dB'][2] + vparold*(m * f.hph)) / (cov['dA2'][2] - cov['dA3'][1])
    v2 = (qe*cov['dPhie'][1] - mu*cov['dB'][1] + vparold*(m * f.hth)) / (cov['dA3'][1] - cov['dA2'][2])
   # v2 = (qe*f.dPhie[0] - f.dAph[0]*vparold - mu*f.dB[0]) / (f.dAth[0])
   # v1 = ((f.dAph[1] - f.dAth[2])*vparold - qe*f.dPhie[1] + mu*f.dB[1]) / (f.dAth[0])
    
    #vpar = vparold + dt*(-f.dAph[0]*v1 + (f.dAth[2] - f.dAph[1])*v2)
    vpar = vparold + dt*(m*f.hth*v1 + m*f.hph*v2)
    v = np.array([v1,v2,vpar])
    
    x0 = z[:2,kt]
    sol = root(F, x0, method='hybr',tol=1e-12,args=(vpar, v))
    z[:2,kt+1] = sol.x
    #z[0,kt+1] = newton(z[:,kt], vpar, v)

plot_orbit(z)
plt.show()




#implicit





