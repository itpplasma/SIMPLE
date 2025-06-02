
from pysimple import simple, simple_main, params, orbit_symplectic, params
from pysimple import get_can_sub as coord
from pysimple import field_can_mod

import numpy as np
import matplotlib.pyplot as plt
#%%
tracy = simple.Tracer()

simple_main.init_field(tracy, "wout.nc", 3, 3, 3, 1)
params.params_init()

z0_vmec = np.array([0.8, 1.0, 0.2, 1.0, 0.5])   # s, th, ph, v/v_th, v_par/v
z0_can = z0_vmec.copy()
z0_can[1:3] = coord.vmec_to_can(z0_vmec[0], z0_vmec[1], z0_vmec[2])
simple.init_integrator(tracy, z0_can)

nt = 10000

z_vmec = np.zeros([nt, 5])
x_cyl = np.zeros([nt, 3])
x_cyl = coord.vmec_to_cyl(z_vmec[0, 0], z_vmec[0, 1], z_vmec[0, 2])

for kt in range(nt-1):
    orbit_symplectic.orbit_timestep_sympl(tracy.si, tracy.f)
    z_integ[kt+1, :] = tracy.si.z
    z_vmec[kt+1, 0] = z_integ[kt+1, 0]
    z_vmec[kt+1, 1:3] = coord.can_to_vmec(
        z_integ[kt+1, 0], z_integ[kt+1, 1], z_integ[kt+1, 2])
    z_vmec[kt+1, 3] = np.sqrt(tracy.f.mu*tracy.f.bmod+0.5*tracy.f.vpar**2)
    z_vmec[kt+1, 4] = tracy.f.vpar/(z_vmec[kt+1, 3]*np.sqrt(2))
    z_cyl[kt+1, :2] = coord.vmec_to_cyl(z_vmec[kt+1, 0], z_vmec[kt+1, 1], z_vmec[kt+1, 2])
    z_cyl[kt+1, 2] = z_vmec[kt+1, 2]
