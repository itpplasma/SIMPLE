#%%
from numba import njit
from field import afield_simple, bfield_simple, bmod_simple
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np


@njit
def velo(t, z):
    B = np.empty(3)
    bfield_simple(z[0], z[1], z[2], B)

    zdot = np.empty(6)
    zdot[0] = z[3]
    zdot[1] = z[4]
    zdot[2] = z[5]
    zdot[3] = z[4] * B[2] - z[5] * B[1]
    zdot[4] = z[5] * B[0] - z[3] * B[2]
    zdot[5] = z[3] * B[1] - z[4] * B[0]

    return zdot


@njit
def velo_cpp(t, z, mu0ovm):
    B = np.empty(3)
    Bmod = np.empty(1)
    dBmod = np.empty(3)
    bfield_simple(z[0], z[1], z[2], B)
    bmod_simple(z[0], z[1], z[2], Bmod, dBmod)

    zdot = np.empty(6)
    zdot[0] = z[3]
    zdot[1] = z[4]
    zdot[2] = z[5]
    zdot[3] = z[4] * B[2] - z[5] * B[1] - mu0ovm * dBmod[0]
    zdot[4] = z[5] * B[0] - z[3] * B[2] - mu0ovm * dBmod[1]
    zdot[5] = z[3] * B[1] - z[4] * B[0] - mu0ovm * dBmod[2]
 
    return zdot


# %%
z0 = np.array([1.0, 2.0, 3.0, 0.1, 0.1, 0.1])

sol = solve_ivp(velo, [0, 500], z0, method="RK45", rtol=1e-9, atol=1e-9)

z0_cpp = np.array([1.04, 2.04, 2.94, 0.1, 0.1, 0.1])  # Guessed guiding center

# Set initial conditions for guiding center CPP
B = np.empty(3)
Bmod = np.empty(1)
dBmod = np.empty(3)
bfield_simple(z0_cpp[0], z0_cpp[1], z0_cpp[2], B)
bmod_simple(z0_cpp[0], z0_cpp[1], z0_cpp[2], Bmod, dBmod)

vpar = np.dot(z0_cpp[3:], B/Bmod)
vperp = np.sqrt(np.sum(z0_cpp[3:]**2) - vpar**2)
mu0ovm = 0.5 * vperp**2 / Bmod
print(mu0ovm)

z0_cpp[3:] = vpar*B/Bmod




sol_cpp = solve_ivp(
    velo_cpp, [0, 500], z0_cpp, method="RK45", rtol=1e-9, atol=1e-9, args=(mu0ovm)
)

# %%
plt.figure()
plt.plot(sol.y[0], sol.y[1])
plt.plot(z0_cpp[0], z0_cpp[1], "ro")
plt.plot(sol_cpp.y[0], sol_cpp.y[1])
plt.show()
# %%
plt.figure()
plt.plot(sol.y[1], sol.y[2])
plt.plot(z0_cpp[1], z0_cpp[2], "ro")
plt.plot(sol_cpp.y[1], sol_cpp.y[2])
plt.show()

print(len(sol_cpp.t))  # This will print the number of time steps the solver used

# %%