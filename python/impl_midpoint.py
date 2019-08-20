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

steps_per_bounce = 64
dt, nt = timesteps(steps_per_bounce, nbounce = 100)
nlag = 1 # order of Lagrange extrapolation

z = zeros([3,nt+1])
z[:,0] = [r0,th0,r0]
zold = z[:,0]

def F(x, xold, pthold):
    """ Cost function for implicit midpoint rule in axisymmetric field """   
    global pth
    ret = zeros(3)
    
     # evaluate at midpoint
    [H, pth, vpar, dHdx, dHdpph, dpthdx, 
     dpthdpph, dvpardx, dvpardpph] = get_der(
       array([x[2], 0.5*(x[1] + xold[1]), 0.0, pph0]))
    
    ret[0] = dpthdx[0]*(x[1] - xold[1]) - dt*dHdx[0]
    
    dpthdrmid = dpthdx[0]
    pthdotbar = dpthdrmid*dHdx[1]-dpthdx[1]*dHdx[0]
    ret[1] = dpthdrmid*(pth - pthold) + dt/2.0*pthdotbar
    
     # evaluate at final position
    [H, pth, vpar, dHdx, dHdpph, dpthdx, 
     dpthdpph, dvpardx, dvpardpph] = get_der(array([x[0],x[1],x[2],pph0]))
     
    ret[2] = dpthdrmid*(pth-pthold) + dt*pthdotbar
    
    return ret

#%%
from time import time
tic = time()
[H, pth, vpar] = get_val(array([r0,th0,ph0,pph0]))
for kt in range(nt):
    pthold = pth
    
    # Initialize via Lagrange extrapolation
    if(kt>=nlag):
        extrapr = lagrange(arange(-nlag,1), z[0,kt-nlag:kt+1])
        extrapth = lagrange(arange(-nlag,1), z[1,kt-nlag:kt+1])
        extraprmid = lagrange(arange(-nlag,1), z[2,kt-nlag:kt+1])
        z0 = array([extrapr(1.0),extrapth(1.0),extraprmid(1.0)])
    else:
        z0 = z[:,kt]
        
    sol = root(F, z0, method='hybr',tol=1e-12,args=(zold, pthold))
    z[:,kt+1] = sol.x
    zold = z[:,kt+1]
    pthold = pth
    
print('Field evaluations: {}'.format(common.neval))
print('Time taken: {}'.format(time()-tic))

plot_orbit(z)

