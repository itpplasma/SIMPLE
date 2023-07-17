#%%
import jax.numpy as np
from jax import grad, jit, vmap, value_and_grad
from jax import random
from jax.lax import fori_loop
from functools import partial
import matplotlib.pyplot as plt

#%%
@jit
def timestep(k, z):
    xnew = z[0] + 0.01*z[1]
    return np.array([
        xnew,
        z[1] + 0.01*np.sin(xnew)*(1.0 + z[3]*np.cos(z[2])),
        z[2] + 0.01,
        z[3]
    ])


def xend(x0, t0, eps, nt):
    z0 = np.hstack([x0, t0, eps])
    timesteps = fori_loop(
        0, nt, timestep, z0)
    return timesteps[0]
# %%
x0 = np.array([0.2, 0.3])
t0 = np.array(0.2)
eps = np.array(0.1)
#grad(xend, 2)(z0, t0, eps)
test = value_and_grad(xend, 2)
test(x0, t0, eps, 1000)
#dzdt(z0, t0, 0.1)
#%%
z0 = np.hstack([x0, t0, eps])
z = np.empty((4, 100))
z = z.at[:,0].set(z0)

for k in np.arange(100):
    z = z.at[:,k+1].set(timestep(k, z[k,:]))
# %%
fig, ax = plt.subplots()
ax.plot(z[0,:], z[1,:])

# %%

def f(x):
    if(x>0): return x
    return -x

fg = value_and_grad(f)

# %%
