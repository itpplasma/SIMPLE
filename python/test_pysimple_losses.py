#%%
from pysimple import simple, simple_main, params as p, new_vmec_stuff_mod as stuff
import numpy as np
import matplotlib.pyplot as plt

for item in dir(p):
    if item.startswith("__"): continue
    try:
        print(f'{item}: {getattr(p, item)}')
    except:
        print(f'{item}: NULL')

tracy = simple.Tracer()

print('wout.nc', stuff.ns_s, stuff.ns_tp, stuff.multharm, p.integmode)

stuff.multharm = 3     # Fast but inaccurate splines
p.ntestpart = 32
p.trace_time = 1e-3
p.contr_pp = -1e10     # Trace all passing passing
p.startmode = 2        # Read start.dat for initial conditions
simple_main.run()

#%%
print(p.times_lost)

#%% Time steps
t = np.linspace(p.dtau/p.v0, p.trace_time, p.ntimstep)

plt.figure()
plt.semilogx(t, p.confpart_pass + p.confpart_trap)
plt.xlim([1e-4, p.trace_time])
plt.xlabel('time')
plt.ylabel('confined fraction')

plt.figure()
condi = np.logical_and(p.times_lost > 0, p.times_lost < p.trace_time)
plt.semilogx(p.times_lost[condi], p.perp_inv[condi], 'x')
plt.xlim([1e-4, p.trace_time])
plt.xlabel('loss time')
plt.ylabel('perpendicular invariant')
plt.show()

# %%
