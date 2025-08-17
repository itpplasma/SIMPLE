# %%
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import root
from scipy.interpolate import lagrange
import common
from common import (r0,th0,ph0,pph0,timesteps,get_val,get_der,mu,qe,m,c,)
from field_correct_test import field
from plotale import plot_orbit, plot_mani, plot_cost_function

f = field()

dt = 1.9

nt = 4000

metric = lambda z: {
    "_ii":  np.array([1, z[0]**2, (1 + z[0]*np.cos(z[1]))**2]),
    "^ii":  np.array([1, 1/z[0]**2, 1/(1 + z[0]*np.cos(z[1]))**2]),
    "d_11": np.array([0,0,0]),
    "d_22": np.array([2*z[0], 0,0]),
    "d_33": np.array([2*(1+z[0]*np.cos(z[1]))*np.cos(z[1]), -2*(1+z[0]*np.cos(z[1]))*np.sin(z[1]), 0]),
}

def dLdx(v, dAth, dAph, dB, g):
    ret = np.zeros(3)
    ret = (m/2*(g['d_11']*v[0]**2 + g['d_22']*v[1]**2 + g['d_33']*v[2]**2) + qe/c*(dAth*v[1] + dAph*v[2]) - mu*dB)
    return ret


def F(p, xold, pold):
    global A

    ret = np.zeros(3)
    f.evaluate(xold[0], xold[1], xold[2])
    g = metric(xold)
    A = np.array([0, f.co_Ath, f.co_Aph])
    v = 1/m * g['^ii'] *(p - qe/c * A)
    
    ret = p - pold - dt * dLdx(v, f.co_dAth, f.co_dAph, f.dB, g)
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

    sol = root(F, p, method='hybr',tol=1e-12,args=(z[:,kt], pold))
    p = sol.x

    z[:,kt+1] = z[:,kt] + dt/m * g['^ii']*(p - qe/c * A)




plot_orbit(z)
plt.plot(z[0,:5]*np.cos(z[1,:5]), z[0,:5]*np.sin(z[1,:5]), "o")
plt.show()

plot_mani(z)
plt.show()


# %%