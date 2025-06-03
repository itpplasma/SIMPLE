#%%
from pysimple import simple, simple_main, params, orbit_symplectic, params
from pysimple import get_can_sub as coord
from pysimple import field_can_mod

import numpy as np
import matplotlib.pyplot as plt
#%%
tracy = simple.Tracer()

simple_main.init_field(tracy, "wout.nc", 3, 3, 3, 1)
params.params_init()
#%%

z0_vmec = np.array([0.8, 1.0, 0.2, 1.0, 0.5])   # s, th, ph, v/v_th, v_par/v
z0_can = z0_vmec.copy()
z0_can[1:3] = coord.vmec_to_can(z0_vmec[0], z0_vmec[1], z0_vmec[2])
simple.init_sympl(tracy.si, tracy.f, z0_can, params.dtaumin, params.dtaumin, 1e-13, tracy.integmode)

nt = 10000

z_can = np.zeros([nt, 4])  # s, th_c, ph_c, p_phi
z_vmec = np.zeros([nt, 5])
x_cyl = np.zeros([nt, 3])
R0, Z0 = coord.vmec_to_cyl(z_vmec[0, 0], z_vmec[0, 1], z_vmec[0, 2])
phi0 = z_vmec[0, 2]
x_cyl[0, 0] = R0
x_cyl[0, 1] = z_vmec[0, 2]
x_cyl[0, 2] = Z0

z_can[0,:] = tracy.si.z

for kt in range(nt-1):
    orbit_symplectic.orbit_timestep_sympl_expl_impl_euler(tracy.si, tracy.f)
    z_can[kt+1, :] = tracy.si.z
    z_vmec[kt+1, 0] = z_can[kt+1, 0]
    z_vmec[kt+1, 1:3] = coord.can_to_vmec(
        z_can[kt+1, 0], z_can[kt+1, 1], z_can[kt+1, 2])
    z_vmec[kt+1, 3] = np.sqrt(tracy.f.mu*tracy.f.bmod+0.5*tracy.f.vpar**2)
    z_vmec[kt+1, 4] = tracy.f.vpar/(z_vmec[kt+1, 3]*np.sqrt(2))
    R, Z = coord.vmec_to_cyl(z_vmec[kt+1, 0], z_vmec[kt+1, 1], z_vmec[kt+1, 2])
    phi = z_vmec[kt+1, 2]
    x_cyl[kt+1, 0] = R
    x_cyl[kt+1, 1] = phi
    x_cyl[kt+1, 2] = Z

plt.figure()
plt.plot(x_cyl[:, 0], x_cyl[:, 2])
plt.xlabel('R')
plt.ylabel('Z')

# %% now in 3D
from mpl_toolkits import mplot3d
fig = plt.figure()
ax = plt.axes(projection='3d')
ax.plot3D(x_cyl[:, 0]*np.cos(x_cyl[:, 1]), x_cyl[:, 0]*np.sin(x_cyl[:, 1]), x_cyl[:, 2])
ax.set_xlabel('R*cos(phi)')
ax.set_ylabel('R*sin(phi)')
ax.set_zlabel('Z')
plt.show()
