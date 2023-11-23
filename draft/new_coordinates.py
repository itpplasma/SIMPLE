#%%
from pysimple import simple, params, orbit_symplectic
from pysimple import get_canonical_coordinates_sub as coord

import numpy as np
import matplotlib.pyplot as plt
#%%
tracy = params.Tracer()

simple.init_field(tracy, "wout.nc", 3, 3, 3, 1)
simple.init_params(tracy, 2, 4, 3.5e6, 256, 1, 1e-13)

#%%
ns = 32
nth = 32
nph = 32


s = np.linspace(0, 1, ns)
th = np.arange(nth)*2*np.pi/nth
ph = np.arange(nph)*2*np.pi/nph

S, TH = np.meshgrid(s, th)

x = np.zeros((ns*nth*nph, 3))
k = 0
for ks in np.arange(ns):
    for kth in np.arange(nth):
        for kph in np.arange(nph):
            x[k, :2] = coord.vmec_to_cyl(s[ks], th[kth], ph[kph])
            x[k, 2] = ph[kph]
            k = k+1
x = np.array(x)

plt.figure()
plt.plot(x[:, 0], x[:, 1], 'o')
# %%
import scipy as sp

fit = sp.interpolate.bisplrep(S.flatten(), TH.flatten(),
                              x[:, 0].T)

# %%