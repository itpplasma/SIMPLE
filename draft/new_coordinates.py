#%%
from pysimple import simple, params, orbit_symplectic
from pysimple import get_can_sub as coord
from pysimple import splint_vmec_data
from pysimple import magfie_sub

magfie_vmec = magfie_sub.magfie_vmec

import numpy as np
import matplotlib.pyplot as plt
#%%
tracy = params.Tracer()

simple.init_field(tracy, "wout.nc", 3, 3, 3, 1)
simple.init_params(tracy, 2, 4, 3.5e6, 256, 1, 1e-13)

#%%
ns = 32
nth = 32
nph = 1


s = np.linspace(1.0/ns, 1.0, ns)
th = np.arange(nth)*2*np.pi/nth
ph = np.arange(nph)*2*np.pi/nph

S, TH = np.meshgrid(s, th)

x = np.empty((ns*nth*nph, 3))

bder = np.empty(3)
hcovar = np.empty(3)
hctrvr = np.empty(3)
hcurl = np.empty(3)

bmods = np.empty(ns*nth*nph)

k = 0
for ks in np.arange(ns):
    for kth in np.arange(nth):
        for kph in np.arange(nph):
            x[k, :2] = coord.vmec_to_cyl(s[ks], th[kth], ph[kph])
            x[k, 2] = ph[kph]

            A_phi,A_theta,dA_phi_ds,dA_theta_ds,aiota, \
            R,Z,alam,dR_ds,dR_dt,dR_dp,dZ_ds,dZ_dt,dZ_dp, \
                dl_ds,dl_dt,dl_dp = splint_vmec_data(
                    s[ks], th[kth], ph[kph])

            bmod,sqrtg = magfie_vmec(
                [s[ks], th[kth], ph[kph]],bder,hcovar,hctrvr,hcurl)

            bmods[k] = bmod

            k = k+1
x = np.array(x)

plt.figure()
plt.scatter(x[:, 0], x[:, 1], c=bmods.flatten(), cmap='viridis')
plt.colorbar()

# %%
from scipy.interpolate import bisplrep, bisplev

rmin = 600.0
zmin = -400.0
rmax = 1400.0
zmax = 400.0

tx = np.linspace(rmin, rmax, 64)
ty = np.linspace(zmin, zmax, 64)

tck = bisplrep(x[:,0], x[:,1], bmods, kx=3, ky=3,
                  xb=rmin, xe=rmax, yb=zmin, ye=zmax, s=1e-4)

# %% evaluate and plot
nr = 64
nz = 64
rplot = np.linspace(rmin, rmax, nr)
zplot = np.linspace(zmin, zmax, nz)
RR, ZZ = np.meshgrid(rplot, zplot)
FF = bisplev(RR[0,:], ZZ[:,0], tck, dx=0, dy=0)  # Add dx and dy parameters
plt.figure()
plt.contourf(RR, ZZ, FF)
#plt.clim(bmods.min(), bmods.max())
plt.colorbar()


# %%
