"""
Created:  2018-08-08
Modified: 2019-03-07
Author:   Christopher Albert <albert@alumni.tugraz.at>
"""

from numpy import array, zeros, arange
from scipy.interpolate import lagrange
import common
from common import f, r0, th0, ph0, pph0, timesteps, get_val, get_der2, newton1
                      
from plotting import plot_orbit, plot_cost_function_jac

steps_per_bounce = 8
dt, nt = timesteps(steps_per_bounce, nbounce = 100)
nlag = 1 # order of Lagrange extrapolation

z = zeros([3,nt+1])
z[:,0] = [r0, th0, ph0]

Hplot = zeros(nt+1) # Hamiltonian for plotting
[Hplot[0], pth, vpar] = get_val(array([r0,th0,ph0,pph0]))

def F(r, q, pthold):
    global H, dHdx, dHdpph, pth, dpthdx, vpar, dvpardx, \
     d2pthdx2, d2pthdpphdz, d2Hdx2, d2Hdpphdz, \
     d2vpardx2, d2vpardpphdz
    
    [H, pth, vpar,  dHdx, dHdpph, dpthdx, dpthdpph, dvpardx, dvpardpph, d2pthdx2, 
     d2pthdpphdz, d2Hdx2, d2Hdpphdz, d2vpardx2, d2vpardpphdz] = get_der2(
         array([r[0],q[0],q[1],pph0]))
    
    ret = dpthdx[0]*(pth-pthold) - dt*(dHdx[0]*dpthdx[1]-dHdx[1]*dpthdx[0])
    
    jac = d2pthdx2[0]*(pth-pthold) + dpthdx[0]**2 - dt*(
            d2Hdx2[0]*dpthdx[1]-d2pthdx2[0]*dHdx[1]
            +dHdx[0]*d2pthdx2[3]-dpthdx[0]*d2Hdx2[3])
    
    return ret, [jac]


#%%
common.neval = 0

from time import time
tic = time()
nbounce = 0
for kt in range(nt):
    pthold = pth
    
    # Initialize via Lagrange extrapolation
    if(kt>=nlag):
        extrapr = lagrange(arange(-nlag,1), z[0,kt-nlag:kt+1])
        r0 = extrapr(1)
    else:
        r0 = z[0,kt]
        
    sol = newton1(F, r0, rtol=1e-7, atol=1e-15, args=(z[1:,kt], pthold))
    z[0,kt+1] = sol.x
    
    dHdx[0] = dHdx[0] + (sol.x[0] - sol.xold[0])*d2Hdx2[0]
    dpthdx[0] = dpthdx[0] + (sol.x[0] - sol.xold[0])*d2pthdx2[0]
    vpar = vpar + (sol.x[0] - sol.xold[0])*dvpardx[0]
    f.B = f.B + (sol.x[0] - sol.xold[0])*f.dB[0]
    f.hph = f.hph + (sol.x[0] - sol.xold[0])*f.dhph[0]
    pth = pth + (sol.x[0] - sol.xold[0])*dpthdx[0]
    H = H + (sol.x[0] - sol.xold[0])*dHdx[0]
    
    z[1,kt+1] = z[1,kt] + dt*dHdx[0]/dpthdx[0]
    z[2,kt+1] = z[2,kt] + dt*vpar/f.hph
    
    Hplot[kt+1] = H
    
print('Field evaluations: {}'.format(common.neval))
print('Time taken: {}'.format(time()-tic))

plot_orbit(z)
plot_cost_function_jac(F, z[:,-2], z[:,-1], pthold)
