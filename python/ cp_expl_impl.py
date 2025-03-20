# %%
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import root
from scipy.interpolate import lagrange
import common
from common import (r0,th0,ph0,pph0,timesteps,get_val,get_der,mu,qe,m,c,)
from field_correct_test import field
from plotting import plot_orbit, plot_cost_function

dt, nt = timesteps(steps_per_bounce=8, nbounce=100)
print(dt)
print(nt)

f = field()
dt = 1
nt = 1000

z = np.zeros([3, nt + 1])
z[:, 0] = [r0, th0, ph0]

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

def implicit_p(p, pold, Ath, Aph, dAth, dAph, g):
    ret = np.zeros(3)

    ctrv = np.zeros(3)
    ctrv[0] = 1/m *g['^11'] *(p[0])
    ctrv[1] = 1/m *g['^22'] *(p[1] - qe/c*Ath)
    ctrv[2] = 1/m *g['^33'] *(p[2] - qe/c*Aph)

    ret[0] = p[0] - pold[0] - dt*(qe/c*(ctrv[1]*dAth[0] + ctrv[2]*dAph[0]) + m/2*(g['d_11'][0]*ctrv[0]**2 + g['d_22'][0]*ctrv[1]**2 + g['d_33'][0]*ctrv[2]**2))
    ret[1] = p[1] - pold[1] - dt*(qe/c*(ctrv[1]*dAth[1] + ctrv[2]*dAph[1]) + m/2*(g['d_11'][1]*ctrv[0]**2 + g['d_22'][1]*ctrv[1]**2 + g['d_33'][1]*ctrv[2]**2))
    ret[2] = p[2] - pold[2] - dt*(qe/c*(ctrv[1]*dAth[2] + ctrv[2]*dAph[2]) + m/2*(g['d_11'][2]*ctrv[0]**2 + g['d_22'][2]*ctrv[1]**2 + g['d_33'][2]*ctrv[2]**2))

    return ret

#Initial Conditions 
f.evaluate(r0, th0, ph0)
g = metric(z[:,0])
ctrv = np.zeros(3)
ctrv[0] = np.sqrt(g['^11']*mu*2*f.B)
ctrv[1] = 0 #-np.sqrt(g['^22']*mu*2*f.B/3)
ctrv[2] = 0 #-np.sqrt(g['^33']*mu*2*f.B/3)
p = np.zeros(3)
p[0] = g['_11']*ctrv[0]
p[1] = g['_22']*ctrv[1] + qe/c * f.co_Ath
p[2] = g['_33']*ctrv[2] + qe/c * f.co_Aph
pold = p

from time import time
tic = time()
for kt in range(nt):

    sol = root(implicit_p, p, method='hybr',tol=1e-12,args=(pold, f.co_Ath, f.co_Aph, f.co_dAth, f.co_dAph, g))
    p = sol.x
    pold = p 

    z[0,kt+1] = z[0,kt] + dt/m * g['^11']*(p[0])
    z[1,kt+1] = z[1,kt] + dt/m * g['^22']*(p[1] - qe/c*f.co_Ath)
    z[2,kt+1] = z[2,kt] + dt/m * g['^33']*(p[2] - qe/c*f.co_Aph)

    f.evaluate(z[0, kt+1], z[1, kt+1], z[2, kt+1])
    g = metric(z[:,kt+1])

    p = np.zeros(3)
    p[0] = g['_11']*m*(z[0,kt+1]-z[0,kt])/dt
    p[1] = g['_22']*m*(z[1,kt+1]-z[1,kt])/dt + qe/c*f.co_Ath
    p[2] = g['_33']*m*(z[2,kt+1]-z[2,kt])/dt + qe/c*f.co_Aph


plot_orbit(z)
plt.plot(z[0,:5]*np.cos(z[1,:5]), z[0,:5]*np.sin(z[1,:5]), "o")
plt.show()




# %%
