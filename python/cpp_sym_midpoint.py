# %%
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import root
from scipy.interpolate import lagrange
import common
from common import (r0,th0,ph0,pph0,timesteps,get_val,get_der,mu,qe,m,c,)
from field_correct_test import field
from plotting import plot_orbit, plot_cost_function

f = field()

dt = 1400
nt = 1000



metric = lambda z: {
    "_11": 1,
    "_22": z[0]**2,
    "_33": (1 + z[0]*np.cos(z[1]))**2,
    "^11": 1,
    "^22": 1/z[0]**2,
    "^33": 1/(1 + z[0]*np.cos(z[1]))**2,
    "d_11": [0,0,0],
    "d_22": [2*z[0], 0,0],
    "d_33": [2*(1+z[0]*np.cos(z[1]))*np.cos(z[1]), -2*(1+z[0]*np.cos(z[1]))*np.sin(z[1]), 0],
}

def ctrv(x, xold):
    v = np.zeros(3)
    v[0] = (x[0] - xold[0])/dt
    v[1] = (x[1] - xold[1])/dt
    v[2] = (x[2] - xold[2])/dt
    return v

def pdot(v, dAth, dAph, dB, g):
    pdot = np.zeros(3)
    pdot[0] = (m/2*(g['d_11'][0]*v[0]**2 + g['d_22'][0]*v[1]**2 + g['d_33'][0]*v[2]**2) + qe/c*(dAth[0]*v[1] + dAph[0]*v[2]) - mu*dB[0])
    pdot[1] = (m/2*(g['d_11'][1]*v[0]**2 + g['d_22'][1]*v[1]**2 + g['d_33'][1]*v[2]**2) + qe/c*(dAth[1]*v[1] + dAph[1]*v[2]) - mu*dB[1])
    pdot[2] = (m/2*(g['d_11'][2]*v[0]**2 + g['d_22'][2]*v[1]**2 + g['d_33'][2]*v[2]**2) + qe/c*(dAth[2]*v[1] + dAph[2]*v[2]) - mu*dB[2])
    return pdot



def F(x, xold, pold):
    global p
    ret = np.zeros(3)

    xmid = np.zeros(3)
    xmid = (x + xold)/2

    vmid = ctrv(x, xold)

    f.evaluate(xmid[0], xmid[1], xmid[2])
    g = metric(xmid)
    pmid = pold + dt/2 * pdot(vmid, f.co_dAth, f.co_dAph, f.dB, g)
    p = pold + dt * pdot(vmid, f.co_dAth, f.co_dAph, f.dB, g)

    ret[0] = x[0] - xold[0] - dt/m*g['^11']*(pmid[0])
    ret[1] = x[1] - xold[1] - dt/m*g['^22']*(pmid[1] - qe/c*f.co_Ath)
    ret[2] = x[2] - xold[2] - dt/m*g['^33']*(pmid[2] - qe/c*f.co_Aph)
    return ret


# Initial Conditions
z = np.zeros([3, nt + 1])
z[:, 0] = [r0, th0, ph0]

f.evaluate(r0, th0, ph0)
g = metric(z[:,0])
v = np.zeros(3)
p = np.zeros(3)
p[0] = v[0] 
p[1] = v[1] + qe/c*f.co_Ath
p[2] = v[2] + qe/c*f.co_Aph



from time import time
tic = time()
for kt in range(nt):
    pold = p

    sol = root(F, z[:,kt], method='hybr',tol=1e-12,args=(z[:,kt], pold))
    z[:,kt+1] = sol.x





plot_orbit(z)
plt.plot(z[0,:5]*np.cos(z[1,:5]), z[0,:5]*np.sin(z[1,:5]), "o")
plt.show()






# %%
