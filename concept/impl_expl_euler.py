"""
Created:  2018-08-08
Modified: 2019-03-07
Author:   Christopher Albert <albert@alumni.tugraz.at>
"""

from numpy import array, zeros, arange
from scipy.optimize import root
from scipy.interpolate import lagrange
import common
from common import r0, th0, ph0, pph0, timesteps, get_val, get_der
from plotting import plot_orbit

dt, nt = timesteps(steps_per_bounce = 8, nbounce = 100)
nlag = 1 # order of Lagrange extrapolation

z = zeros([3,nt+1])
z[:,0] = [r0,th0,ph0]

def F(x, thold, pthold):
    """ Cost function in x=(r,th) ignoring toroidal angle ph in axisymmetry""" 
    [H, pth, vpar, dHdx, dHdpph, dpthdx, 
     dpthdpph, dvpardx, dvpardpph] = get_der(array([x[0],x[1],0.0,pph0]))
    
    ret = zeros(2)
    
    ret[0] = pth - pthold
    ret[1] = dpthdx[0]*(x[1] - thold) - dt*dHdx[0]
    
    return ret

#%%
from time import time
tic = time()
[H, pth, vpar] = get_val(array([r0,th0,ph0,pph0]))
for kt in range(nt):
    pthold = pth
    
    # Initialize via Lagrange extrapolation
    if(kt>=nlag):
        x0 = zeros(2)
        extrapr0 = lagrange(arange(-nlag, 1), z[0,kt-nlag:kt+1])
        extrapr1 = lagrange(arange(-nlag, 1), z[1,kt-nlag:kt+1])
        x0[0] = extrapr0(1.0)
        x0[1] = extrapr1(1.0)
    else:
        x0 = array([z[0,kt], z[1,kt]])
        
    sol = root(F, x0, method='hybr',tol=1e-12,args=(z[1,kt],pthold))
    z[:2,kt+1] = sol.x
    
    # Tokamak, no change in p_phi
    [H, pth, vpar, dHdx, dHdpph, dpthdx, 
     dpthdpph, dvpardx, dvpardpph] = get_der(array([sol.x[0],sol.x[1],0.0,pph0]))
    
    pth = pth - dt*(dHdx[1] - dHdx[0]*dpthdx[1]/dpthdx[0])
    
print('Field evaluations: {}'.format(common.neval))
print('Time taken: {}'.format(time()-tic))

plot_orbit(z)
