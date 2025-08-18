"""
Created:  2018-08-08
Modified: 2019-03-07
Author:   Christopher Albert <albert@alumni.tugraz.at>
"""

from numpy import array, zeros, arange
from scipy.optimize import root
from scipy.interpolate import lagrange
import common
from common import f, r0, th0, ph0, pph0, timesteps, get_val, get_der
from plotting import plot_orbit, plot_cost_function

dt, nt = timesteps(steps_per_bounce = 8, nbounce = 100)
nlag = 1 # order of Lagrange extrapolation

z = zeros([3,nt+1])
z[:,0] = [r0, th0, ph0]

def F(r, q, pthold):
    """ Cost function in r for axisymmetric field with pph=const. """
    [H, pth, vpar, dHdx, dHdpph, dpthdx, 
     dpthdpph, dvpardx, dvpardpph] = get_der(array([r[0], q[0], q[1], pph0]))
    
    return dpthdx[0]*(pth - pthold) + dt*(dHdx[1]*dpthdx[0]-dHdx[0]*dpthdx[1])

#%%
from time import time
tic = time()
[H, pth, vpar] = get_val(array([r0,th0,ph0,pph0]))
for kt in range(nt):
    pthold = pth
    
    # Initialize via Lagrange extrapolation
    if(kt>=nlag):
        extrapr = lagrange(arange(-nlag, 1), z[0, kt-nlag:kt+1])
        r0 = extrapr(1)
    else:
        r0 = z[0,kt]
        
    # Implicit substep in r
    sol = root(F, r0, method='hybr',tol=1e-12,args=(z[1:,kt], pthold))
    z[0,kt+1] = sol.x
    
    # Explicit substep in q = (th, ph)
    [H, pth, vpar, dHdx, dHdpph, dpthdx, dpthdpph, 
     dvpardx, dvpardpph] = get_der(array([sol.x[0], z[1,kt], z[2,kt], pph0]))
    z[1,kt+1] = z[1,kt] + dt*dHdx[0]/dpthdx[0]
    z[2,kt+1] = z[2,kt] + dt*vpar/f.hph
    
print('Field evaluations: {}'.format(common.neval))
print('Time taken: {}'.format(time()-tic))

plot_orbit(z)
plot_cost_function(F, z[:,-2], z[:,-1], pthold)
