from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import root
from scipy.interpolate import lagrange
import common
from common import (f,r0,th0,ph0,pph0,timesteps,get_val,get_der,mu,qe,m,c,)
from plotting import plot_orbit, plot_cost_function

z = np.zeros(3)
z[0] = r0
z[1] = th0
z[2] = ph0

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

g = metric(z[:])
f.evaluate(r0, th0, ph0)

x = np.sqrt(g['^22']*f.hth*f.hth + g['^33']*f.hph*f.hph)
y = np.sqrt(g['^33']*f.hph*f.hph)

print(x)
print(y)