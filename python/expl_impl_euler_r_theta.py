"""
Created:  2018-08-08
Modified: 2019-03-07
Author:   Christopher Albert <albert@alumni.tugraz.at>
"""

from numpy import array, zeros, arange
from scipy.optimize import root
from scipy.interpolate import lagrange
import common
from common import r0, th0, ph0, pph0, timesteps, get_der
from plotting import plot_orbit
  
dt, nt = timesteps(steps_per_bounce = 11, nbounce = 100)
nlag = 2

z = zeros([3, nt+1])
z[:,0] = [r0, th0, ph0]

[H, pth, vpar, dHdx, dHdpph, dpthdx, 
    dpthdpph, dvpardx, dvpardpph] = get_der(array([r0, th0, ph0, pph0]))

dtau = dt/dpthdx[0] # step in orbit parameter tau
t = zeros(nt+1)

def F(r, zold):
    """ Cost function for implicit step in r """
    
    [H, pth, vpar, dHdx, dHdpph, dpthdx, 
     dpthdpph, dvpardx, dvpardpph] = get_der(array([r[0],zold[1],zold[2],pph0]))
    
    return (r - zold[0]) + dtau*dHdx[1]

#%%
from time import time
tic = time()
for kt in range(nt):
    
    # Initialize via Lagrange extrapolation
    if(kt>=nlag):
        extrapr = lagrange(arange(-nlag, 1), z[0, kt-nlag:kt+1])
        r0 = extrapr(1)
    else:
        r0 = z[0,kt]
        
    # Implicit substep in r
    sol = root(F, r0, method='hybr', tol=1e-12, args=(z[:,kt]))
    z[0,kt+1] = sol.x
    [H, pth, vpar, dHdx, dHdpph, dpthdx, 
     dpthdpph, dvpardx, dvpardpph] = get_der(array([z[0,kt+1],z[1,kt],z[2,kt],pph0]))
    z[1,kt+1] = z[1,kt] + dtau*dHdx[0]
    t[kt+1] = t[kt] + dtau*dpthdx[0]
        
    
print('Field evaluations: {}'.format(common.neval))
print('Time taken: {}'.format(time()-tic))

plot_orbit(z)
