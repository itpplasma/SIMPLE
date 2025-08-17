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

qe = 1.0; m = 1.0; c = 1.0 # particle charge, mass and speed of light
mu = 1e-5 # magnetic moment

# Initial conditions
r0    = 0.1
th0   = 1.5
ph0   = 1.0
vpar0 = 0.0

dt = 800
nt = 2000

metric = lambda z: {
    "_ii":  np.array([1, z[0]**2, (1 + z[0]*np.cos(z[1]))**2]),
    "^ii":  np.array([1, 1/z[0]**2, 1/(1 + z[0]*np.cos(z[1]))**2]),
    "d_11": np.array([0,0,0]),
    "d_22": np.array([2*z[0], 0,0]),
    "d_33": np.array([2*(1+z[0]*np.cos(z[1]))*np.cos(z[1]), -2*(1+z[0]*np.cos(z[1]))*np.sin(z[1]), 0]),
}



def F(x, xold, dLdxold, dLdxdotold):
    global p, dpdt, vmid
    ret = np.zeros(3)

    xmid = (x + xold)/2
    vmid = (x - xold)/dt

    f.evaluate(xmid[0], xmid[1], xmid[2])
    Amid = np.array([0, f.co_Ath, f.co_Aph])
    g = metric(xmid)

    dpdt =  m/2*(g['d_11']*vmid[0]**2 + g['d_22']*vmid[1]**2 + g['d_33']*vmid[2]**2) + qe/c*(f.co_dAth*vmid[1] + f.co_dAph*vmid[2]) - mu*f.dB
    p = m*g['_ii']*vmid + qe/c*Amid

    ret = (dpdt + dLdxold)*dt/2 - (p - dLdxdotold)
    return ret

def action(x, xold):
    X = np.array([x, xold])

    for i in 2:
        print(i)
        f.evaluate(X[1], X[2], X[3])
        g = metric(X)
        



    xmid = (x + xold)/2
    vmid = (x - xold)/dt
    f.evaluate(xmid[0], xmid[1], xmid[2])
    g = metric(xmid)

    Amid = np.array([0, f.co_Ath, f.co_Aph])
    hmid = np.array([f.co_hr, f.co_hth, f.co_hph])
    pmid = m*g['_ii']*vmid + qe/c*Amid
    vpar = m* np.sqrt(np.sum(pmid**2))

    S =  np.sum((qe/c * Amid + m*vpar*hmid)*vmid) - (m/2*vpar**2 - mu*f.B)* dt
    h = 6.6260755*1e-34

    return S


# Initial Conditions
z = np.zeros([3, nt + 1])
z[:, 0] = [r0, th0, ph0]
L = np.zeros(nt)

f.evaluate(r0, th0, ph0)
g = metric(z[:,0])
Amid = np.array([0, f.co_Ath, f.co_Aph])
vmid = np.zeros(3)
p = m*g['_ii']*vmid + qe/c*Amid
dpdt = m/2*(g['d_11']*vmid[0]**2 + g['d_22']*vmid[1]**2 + g['d_33']*vmid[2]**2) + qe/c*(f.co_dAth*vmid[1] + f.co_dAph*vmid[2]) - mu*f.dB


from time import time
tic = time()
for kt in range(nt):

    sol = root(F, z[:,kt], method='hybr',tol=1e-12,args=(z[:,kt], dpdt, p))
    z[:,kt+1] = sol.x




plot_orbit(z)
plt.plot(z[0,:5]*np.cos(z[1,:5]), z[0,:5]*np.sin(z[1,:5]), "o")
plt.show()






# %%
