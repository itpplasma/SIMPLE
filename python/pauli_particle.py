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
dt = 0.01
nt = 5000

z = np.zeros([3, nt + 1])
z[:, 0] = [r0, th0, ph0]

f.evaluate(r0, th0, ph0)
p = np.zeros(3)
p[0] = 0
p[1] = qe/(m*c) * f.Ath
p[2] = qe/(m*c) * f.Aph

pold = p

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




def implicit_p(p, pold, Ath, Aph, dB, dAth, dAph, g):
    ret = np.zeros(3)

    ctrv = np.zeros(3)
    ctrv[0] = 1/m *g['^11'] *(p[0])
    ctrv[1] = 1/m *g['^22'] *(p[1] - qe/c*Ath)
    ctrv[2] = 1/m *g['^33'] *(p[2] - qe/c*Aph)

    ret[0] = p[0] - pold[0] - dt*(qe/c*(ctrv[1]*dAth[0] + ctrv[2]*dAph[0]) - mu*dB[0] + m/2*(g['d_11'][0]*ctrv[0]**2 + g['d_22'][0]*ctrv[1]**2 + g['d_33'][0]*ctrv[2]**2))
    ret[1] = p[1] - pold[1] - dt*(qe/c*(ctrv[1]*dAth[1] + ctrv[2]*dAph[1]) - mu*dB[1] + m/2*(g['d_11'][1]*ctrv[0]**2 + g['d_22'][1]*ctrv[1]**2 + g['d_33'][1]*ctrv[2]**2))
    ret[2] = p[2] - pold[2] - dt*(qe/c*(ctrv[1]*dAth[2] + ctrv[2]*dAph[2]) - mu*dB[2] + m/2*(g['d_11'][2]*ctrv[0]**2 + g['d_22'][2]*ctrv[1]**2 + g['d_33'][2]*ctrv[2]**2))

    return ret



from time import time
tic = time()
for kt in range(nt):

    f.evaluate(z[0, kt], z[1, kt], z[2, kt])
    g = metric(z[:,kt])

    sol = root(implicit_p, p, method='hybr',tol=1e-12,args=(pold, f.Ath, f.Aph, f.dB, f.dAth, f.dAph, g))
    pold = p
    p = sol.x

    z[0,kt+1] = z[0,kt] + dt/m * g['^11']*(p[0])
    z[1,kt+1] = z[1,kt] + dt/m * g['^22']*(p[1] - qe/c*f.Ath)
    z[2,kt+1] = z[2,kt] + dt/m * g['^33']*(p[2] - qe/c*f.Aph)



plot_orbit(z)
plt.show()




#implicit





