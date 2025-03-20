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

dt = 1050
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

def pdot(v, dAth, dAph, dB, g):
    pdot = np.zeros(3)
    pdot[0] = (m/2*(g['d_11'][0]*v[0]**2 + g['d_22'][0]*v[1]**2 + g['d_33'][0]*v[2]**2) + qe/c*(dAth[0]*v[1] + dAph[0]*v[2]) - mu*dB[0])
    pdot[1] = (m/2*(g['d_11'][1]*v[0]**2 + g['d_22'][1]*v[1]**2 + g['d_33'][1]*v[2]**2) + qe/c*(dAth[1]*v[1] + dAph[1]*v[2]) - mu*dB[1])
    pdot[2] = (m/2*(g['d_11'][2]*v[0]**2 + g['d_22'][2]*v[1]**2 + g['d_33'][2]*v[2]**2) + qe/c*(dAth[2]*v[1] + dAph[2]*v[2]) - mu*dB[2])
    return pdot

def p(v, Ath, Aph, g):
    p = np.zeros(3)
    p[0] = m*g['_11']*v[0] + 0
    p[1] = m*g['_22']*v[1] + qe/c * Ath
    p[2] = m*g['_33']*v[2] + qe/c * Aph
    return p


def F(x, xold, dLdxold, dLdxdotold):
    ret = np.zeros(3)

    xmid = (x + xold)/2
    vmid = (x - xold)/dt

    f.evaluate(xmid[0], xmid[1], xmid[2])
    g = metric(xmid)

    dLdx = pdot(vmid, f.co_dAth, f.co_dAph, f.dB, g)
    dLdxdot = p(vmid, f.co_Ath, f.co_Aph, g)

    ret = (dLdx + dLdxold)*dt/2 - (dLdxdot - dLdxdotold)
    return ret


# Initial Conditions
z = np.zeros([3, nt + 1])
z[:, 0] = [r0, th0, ph0]

f.evaluate(r0, th0, ph0)
g = metric(z[:,0])
vmid = np.zeros(3)


from time import time
tic = time()
for kt in range(nt):


    dLdx = pdot(vmid, f.co_dAth, f.co_dAph, f.dB, g)
    dLdxdot = p(vmid, f.co_Ath, f.co_Aph, g)

    sol = root(F, z[:,kt], method='hybr',tol=1e-12,args=(z[:,kt], dLdx, dLdxdot))
    z[:,kt+1] = sol.x

    x = z[:,kt+1]
    xold = z[:,kt]
    xmid = (x + xold)/2
    vmid = (x - xold)/dt
    f.evaluate(xmid[0], xmid[1], xmid[2])
    g = metric(xmid)



plot_orbit(z)
plt.plot(z[0,:5]*np.cos(z[1,:5]), z[0,:5]*np.sin(z[1,:5]), "o")
plt.show()






# %%
