#%%
from pysimple import simple, params, orbit_symplectic, vmec_to_can,\
     can_to_vmec, vmec_to_cyl
import numpy as np
import matplotlib.pyplot as plt
#%%
tracy = params.Tracer()

simple.init_field(tracy, "wout.nc", 3, 3, 3, 1)
simple.init_params(tracy, 2, 4, 3.5e6, 256, 1, 1e-13)

#%% Initial conditions
z0_vmec = np.array([0.8, 1.0, 0.2, 1.0, 0.5])   # s, th, ph, v/v_th, v_par/v
z0_can = z0_vmec.copy()  # s, th_c, ph_c, v/v_th, v_par/v

z0_can[1:3] = vmec_to_can(z0_vmec[0], z0_vmec[1], z0_vmec[2])

simple.init_integrator(tracy, z0_can)

print(f'B = {tracy.f.bmod}')

nt = 10000
z_integ = np.zeros([nt, 4])  # s, th_c, ph_c, p_phi
z_vmec = np.zeros([nt, 5])  # s, th, ph, v/v_th, v_par/v
z_cyl = np.zeros([nt, 3])
z_integ[0,:] = tracy.si.z
z_vmec[0,:] = z0_vmec
z_cyl[0,:2] = vmec_to_cyl(z_vmec[0, 0], z_vmec[0, 1], z_vmec[0, 2])
z_cyl[0, 2] = z_vmec[0, 2]

for kt in range(nt-1):
    orbit_symplectic.orbit_timestep_sympl(tracy.si, tracy.f)
    z_integ[kt+1, :] = tracy.si.z
    z_vmec[kt+1, 0] = z_integ[kt+1, 0]
    z_vmec[kt+1, 1:3] = can_to_vmec(
        z_integ[kt+1, 0], z_integ[kt+1, 1], z_integ[kt+1, 2])
    z_vmec[kt+1, 3] = np.sqrt(tracy.f.mu*tracy.f.bmod+0.5*tracy.f.vpar**2)
    z_vmec[kt+1, 4] = tracy.f.vpar/(z_vmec[kt+1, 3]*np.sqrt(2))
    z_cyl[kt+1, :2] = vmec_to_cyl(z_vmec[kt+1, 0], z_vmec[kt+1, 1], z_vmec[kt+1, 2])
    z_cyl[kt+1, 2] = z_vmec[kt+1, 2]

plt.figure()
plt.plot(z_vmec[:, 0]*np.cos(z_vmec[:, 1]), z_vmec[:, 0]*np.sin(z_vmec[:, 1]))
plt.xlabel('s * cos(th)')
plt.ylabel('s * sin(th)')
plt.title('Poloidal orbit topology')


plt.figure()
plt.plot(z_vmec[:, 3])
plt.plot(z_vmec[:, 4])
plt.xlabel('Timestep')
plt.ylabel('Normalized velocity')
plt.legend(['v/v_0', 'v_par/v'])
plt.title('Velocities over time')

# 3D plot
from mpl_toolkits import mplot3d
# R = S0 + s*cos(th)
# Z = s*sin(th)
# X = R*cos(ph)
# Y = R*sin(ph)
# Z = Z

S0 = 5.0
plt.figure()
ax = plt.axes(projection='3d')
ax.plot3D((S0 + z_vmec[:, 0]*np.cos(z_vmec[:, 1]))*np.cos(z_vmec[:, 2]),
         (S0 + z_vmec[:, 0]*np.cos(z_vmec[:, 1]))*np.sin(z_vmec[:,2]),
          z_vmec[:, 0]*np.sin(z_vmec[:, 1]))
ax.set_xlabel('(3 + s * cos(th))*cos(ph)')
ax.set_ylabel('(3 + s * sin(th))*cos(ph)')
ax.set_zlabel('s*cos(ph)')
ax.set_title('3D orbit topology')
ax.set_xlim(-5,5)
ax.set_ylim(-5,5)
ax.set_zlim(-5,5)

# Poloidal orbit in RZ

plt.figure()
plt.plot(z_cyl[:, 0], z_cyl[:, 1])
plt.xlabel('R')
plt.ylabel('Z')
plt.title('Poloidal orbit projection')

# 3D orbit in RZ

plt.figure()
ax = plt.axes(projection='3d')
ax.plot(z_cyl[:, 0]*np.cos(z_cyl[:,2]),
    z_cyl[:, 0]*np.sin(z_cyl[:,2]),
    z_cyl[:,1])
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.set_title('Orbit in real space')
ax.set_xlim(-1400,1400)
ax.set_ylim(-1400,1400)
ax.set_zlim(-1400,1400)

# %%


plt.show()
