# %%
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import root
from scipy.interpolate import lagrange
import common
from common import (f,r0,th0,ph0,pph0,timesteps,get_val,get_der,mu,qe,m,c,)
from plotting import plot_orbit, plot_cost_function

dt, nt = timesteps(steps_per_bounce=8, nbounce=100)

print(dt)
print(nt)
dt = 8
nt = 5000

z = np.zeros([3, nt + 1])
z[:, 0] = [r0, th0, ph0]

f.evaluate(r0, th0, ph0)
v = np.zeros(3)
vold = v

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
'''
metric = lambda z: {
    "_11": 1,
    "_22": z[0]**2,
    "_33": z[0]**2 * np.sin(z[1])**2,
    "^11": 1,
    "^22": 1/z[0]**2,
    "^33": 1/(z[0]**2 * np.sin(z[1])**2),
    "d_11": [0,0,0],
    "d_22": [2*z[0], 0,0],
    "d_33": [2*z[0]*np.sin(z[1])**2, 2*z[0]**2 * np.sin(z[1])*np.cos(z[1]), 0],
}
'''


def implicit_v(v, vold, dAth, dAph, dB, g):
    ret = np.zeros(3)


    ret[0] = g['_11']*(v[0] - vold[0]) - dt*m*0.5*(v[0]**2 *g['d_11'][0] + v[1]**2 *g['d_22'][0] + v[2]**2 *g['d_33'][0]) + dt*m*v[0]*(g['d_11'][0]*v[0] + g['d_11'][1]*v[1] + g['d_11'][2]*v[2]) - dt*qe/c *(v[1]*(dAth[0]) +  v[2]*(dAph[0])) + dt*mu*dB[0]
    
    ret[1] = g['_22']*(v[1] - vold[1]) - dt*m*0.5*(v[0]**2 *g['d_11'][1] + v[1]**2 *g['d_22'][1] + v[2]**2 *g['d_33'][1]) + dt*m*v[1]*(g['d_22'][0]*v[0] + g['d_22'][1]*v[1] + g['d_22'][2]*v[2]) - dt*qe/c *(v[0]*(-dAth[0]) + v[2]*(dAph[1] - dAth[2])) + dt*mu*dB[1]
    
    ret[2] = g['_33']*(v[2] - vold[2]) - dt*m*0.5*(v[0]**2 *g['d_11'][2] + v[1]**2 *g['d_22'][2] + v[2]**2 *g['d_33'][2]) + dt*m*v[2]*(g['d_33'][0]*v[0] + g['d_33'][1]*v[1] + g['d_33'][2]*v[2]) - dt*qe/c *(v[0]*(-dAph[0]) + v[1]*(dAth[2] - dAph[1])) + dt*mu*dB[2]

    return ret



from time import time
tic = time()
for kt in range(nt):

    f.evaluate(z[0, kt], z[1, kt], z[2, kt])
    g = metric(z[:,kt])



    sol = root(implicit_v, v, method='hybr',tol=1e-12,args=(vold, f.dAth, f.dAph, f.dB, g))
    vold = v
    v = sol.x

    z[0,kt+1] = z[0,kt] + dt*v[0]
    z[1,kt+1] = z[1,kt] + dt*v[1]
    z[2,kt+1] = z[2,kt] + dt*v[2]



plot_orbit(z)
plt.plot(z[0,:5]*np.cos(z[1,:5]), z[0,:5]*np.sin(z[1,:5]), "o")
plt.show()




#implicit






# %%
