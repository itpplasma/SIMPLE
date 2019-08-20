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

dt, nt = timesteps(steps_per_bounce = 9, nbounce = 100)
nlag = 2 # order of Lagrange extrapolation

z = zeros([4,nt+1])
z[:,0] = [r0,th0,r0,th0]
zold = z[:,0]

def Fx(x,qold,pthold):
    """ Cost function in x=(r,th) ignoring toroidal angle ph in axisymmetry"""
    
    [H, pth, vpar, dHdx, dHdpph, dpthdx, 
     dpthdpph, dvpardx, dvpardpph] = get_der(array([x[0],x[1],0.0,pph0]))
    
    ret = zeros(2)
    ret[0] = pthold - pth
    ret[1] = qold[0] - x[1] + dt/2.0*dHdx[0]/dpthdx[0]
    
    return ret


def Fr(r,q,pthold,pdotold):
    """ Cost function in r for axisymmetric field with pph=const. """
    
    # Tokamak, no change in p_phi
    [H, pth, vpar, dHdx, dHdpph, dpthdx, 
     dpthdpph, dvpardx, dvpardpph] = get_der(array([r[0],q[0],0.0,pph0]))
    
    ret = pthold - pth + dt/2.0*pdotold - dt/2.0*(dHdx[1]-dHdx[0]/dpthdx[0]*dpthdx[1])
    
    return ret

#%%
from time import time
tic = time()
[H, pth, vpar] = get_val(array([r0,th0,ph0,pph0]))
for kt in range(nt):
    pthold = pth
    
    # Initialize via Lagrange extrapolation
    if(kt>=nlag):
         extrapri = lagrange(arange(-nlag, 1), z[2,kt-nlag:kt+1])
         extrapthi = lagrange(arange(-nlag, 1), z[3,kt-nlag:kt+1])
         x0i = array([extrapri(1.0), extrapthi(1.0)])
         extrapr = lagrange(arange(-nlag, 1), z[0,kt-nlag:kt+1])
         r0 = extrapr(1.0)
    else:
        x0i = z[:2,kt]
        r0 = z[0,kt]
        
    # Implicit step
    sol = root(Fx, x0i, method='hybr', tol=1e-12, args=(zold[1:2],pthold))
    xi = sol.x
    z[2:,kt+1] = xi
    
    [H, pth, vpar, dHdx, dHdpph, dpthdx, 
     dpthdpph, dvpardx, dvpardpph] = get_der(array([xi[0],xi[1],0.0,pph0]))
    pdotold = -(dHdx[1]-dHdx[0]/dpthdx[0]*dpthdx[1])
    pthold = pth
    
    sol = root(Fr, r0, method='hybr',tol=1e-12, args=(xi[1:],pthold,pdotold))
    z[0,kt+1] = sol.x
    
    [H, pth, vpar, dHdx, dHdpph, dpthdx, 
     dpthdpph, dvpardx, dvpardpph] = get_der(array([sol.x[0],xi[1],0.0,pph0]))
    z[1,kt+1] = xi[1] + dt/2.0*dHdx[0]/dpthdx[0]
    
    zold = z[:,kt+1]
    
print('Field evaluations: {}'.format(common.neval))
print('Time taken: {}'.format(time()-tic))

plot_orbit(z)

