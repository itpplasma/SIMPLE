"""
Plot field standalone
"""

# %%
from numpy import *
from matplotlib.pyplot import *
from mpl_toolkits.mplot3d import Axes3D      # pylint: disable=unused-import
from netCDF4 import Dataset                  # pylint: disable=no-name-in-module

def fourier_cos(fmn, theta, zeta, xm, xn):
    f = zeros_like(theta)
    for k in np.arange(xm.size):
        f += fmn[k]*cos(xm[k]*theta - xn[k]*zeta)
    return f

# %%
data = Dataset("wout.nc", "r", format="NETCDF4")
ks = 16

s = data.variables['phi']/max(data.variables['phi'])
sh = (s[1:] + s[0:-1])/2.0

figure()
plot(s, data.variables['iotaf'])
plot(sh, data.variables['iotas'][1:])

rmnc = data.variables['rmnc']
zmns = data.variables['zmns']
lmns = data.variables['lmns']

# %%

#s0 = (sdat[ks+1] + sdat[ks])/2.0
#s0 = sdat[ks]
#print('s = {}'.format(s0))

nth = 72
nph = 36
th = linspace(-pi/2, 3*pi/2, nth, endpoint=False)
ph = linspace(0, pi, nph, endpoint=False)
[PH, TH] = np.meshgrid(ph, th)
TH = TH.flatten()
PH = PH.flatten()
TH2 = zeros_like(TH)

# %%

mn_mode_nyq = data.dimensions['mn_mode_nyq']
xm_nyq = data.variables['xm_nyq']
xn_nyq = data.variables['xn_nyq']
bmnc = data.variables['bmnc']

b = fourier_cos(bmnc[ks,:], TH, PH, xm_nyq, xn_nyq)

# %%

figure(figsize=(3.2, 3.2))
contour(PH.reshape(nth, nph), TH.reshape(nth, nph),
        b.reshape(nth, nph)/min(b), cmap='jet', levels=100)
xlabel(r'$\varphi$')
ylabel(r'$\vartheta$')

# %%
# b2 = np.sqrt(bsubu*bsupu + bsubv*bsupv)

# figure(figsize=(3.2, 3.2))
# contour(PH.reshape(nth, nph), TH.reshape(nth, nph),
#         b2.reshape(nth, nph)/min(b), cmap='jet', levels=100)
# xlabel(r'$\varphi$')
# ylabel(r'$\vartheta$')

# rootgrp.close()
