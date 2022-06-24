from pysimple import simple, orbit_symplectic, vmec_to_can, can_to_vmec
import numpy as np
import matplotlib.pyplot as plt

tracy = simple.Tracer()

simple.init_field(tracy, "wout.nc", 3, 3, 3, 1)
simple.init_params(tracy, 2, 4, 3.5e6, 256, 1, 1e-13)

# Initial conditions
z0_vmec = np.array([0.5, 0.3, 0.2, 1.0, 0.1])   # s, th, ph, v/v_th, v_par/v
z0_can = z0_vmec.copy()  # s, th_c, ph_c, v/v_th, v_par/v

z0_can[1:3] = vmec_to_can(z0_vmec[0], z0_vmec[1], z0_vmec[2])

simple.init_integrator(tracy, z0_can)

print(f'B = {tracy.f.bmod}')

nt = 10000
z_integ = np.zeros([nt, 4])  # s, th_c, ph_c, p_phi
z_vmec = np.zeros([nt, 5])  # s, th, ph, v/v_th, v_par/v
z_integ[0,:] = tracy.si.z
z_vmec[0,:] = z0_vmec

for kt in range(nt-1):
    orbit_symplectic.orbit_timestep_sympl(tracy.si, tracy.f)
    z_integ[kt+1, :] = tracy.si.z
    z_vmec[kt+1, 0] = z_integ[kt+1, 0]
    z_vmec[kt+1, 1:3] = can_to_vmec(
        z_integ[kt+1, 0], z_integ[kt+1, 1], z_integ[kt+1, 2])
    z_vmec[kt+1, 3] = np.sqrt(tracy.f.mu*tracy.f.bmod+0.5*tracy.f.vpar**2)
    z_vmec[kt+1, 4] = tracy.f.vpar/(z_vmec[kt+1, 3]*np.sqrt(2))

plt.figure()
plt.plot(z_vmec[:, 0]*np.cos(z_vmec[:, 1]), z_vmec[:, 0]*np.sin(z_vmec[:, 1]))
plt.xlabel('s * cos(th)')
plt.xlabel('s * sin(th)')
plt.title('Poloidal orbit topology')

plt.figure()
plt.plot(z_vmec[:, 3])
plt.plot(z_vmec[:, 4])
plt.xlabel('Timestep')
plt.ylabel('Normalized velocity')
plt.legend(['v/v_0', 'v_par/v'])
plt.title('Velocities over time')

plt.show()
