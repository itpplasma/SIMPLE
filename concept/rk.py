"""
Created:  2018-08-08
Modified: 2019-03-07
Author:   Christopher Albert <albert@alumni.tugraz.at>
"""
from numpy import zeros, array, arange, append, hstack
from scipy.integrate import RK45
from common import r0, th0, ph0, pph0, timesteps, get_val, get_der
import common
from plotting import plot_orbit

dt, nt = timesteps(steps_per_bounce = 1, nbounce = 100)
z = zeros([2,1])
z[:,0] = [r0, th0]

Hplot = zeros(1)
[Hplot[0], pth, vpar] = get_val(array([r0,th0,ph0,pph0]))

def zdot(t, z):
    global Ath, dAth, Aph, dAph, Bmod, dBmod, hth, dhth, hph, dhph,\
           H, dH, pth, pph, dpth, vpar, dvpar    
    [H, pth, vpar, dHdx, dHdpph, dpthdx, 
     dpthdpph, dvpardx, dvpardpph] = get_der(array([z[0],z[1],ph0,pph0]))
    ret = zeros(2)
    ret[0] = -dHdx[1]/dpthdx[0]
    ret[1] = dHdx[0]/dpthdx[0]
    return ret
    
from time import time
tic = time()

tvec = zeros(1)
integ = RK45(zdot, 0.0, z[:,0], nt*dt, rtol=1e-6, atol=1e-12, max_step=dt)

common.neval = 0
maxnt = 10000000
for kt in arange(maxnt):
    if (integ.status == 'finished'):
        break
    integ.step()
    
    tvec = append(tvec, integ.t)
    z = hstack([z, array([integ.y]).T])
    Hplot = append(Hplot, H)
    
print('Field evaluations: {}'.format(common.neval))
print('Time taken: {}'.format(time()-tic))

plot_orbit(z)
