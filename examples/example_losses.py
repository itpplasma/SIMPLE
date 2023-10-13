#%%
from pysimple import simple, simple_main, params as p, new_vmec_stuff_mod as stuff
import numpy as np
import matplotlib.pyplot as plt

nth = 5
nph = 5

stuff.multharm = 5     # Fast but inaccurate splines
p.ntestpart = nth*nph
p.trace_time = 1e-2
p.contr_pp = -1e10     # Trace all passing passing
p.startmode = -1       # Manual start conditions
p.sbeg = [0.4]         # Starting flux surface

tracy = p.Tracer()

simple.init_field(tracy, 'wout.nc',
    stuff.ns_s, stuff.ns_tp, stuff.multharm, p.integmode)

print(stuff.ns_s, stuff.ns_tp, stuff.multharm, p.integmode)

p.params_init()

for item in dir(p):
    if item.startswith("__"): continue
    try:
        print(f'{item}: {getattr(p, item)}')
    except:
        print(f'{item}: NULL')
#%%
# s, th_vmec, ph_vmec, v/v0, v_par/v
p.zstart[0,:] = p.sbeg[0]
TH, PH = np.meshgrid(np.linspace(0,2*np.pi,nth+1)[:-1], np.linspace(0,2*np.pi,nph+1)[:-1])  # TODO: th and phi check order
p.zstart[2, :] = TH.flatten()
p.zstart[3,:] = PH.flatten() 
p.zstart[4, :] = 0.0
#%%
simple_main.run(tracy)

print(p.times_lost)

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
 