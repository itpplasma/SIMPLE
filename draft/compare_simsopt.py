#%%
import numpy as np
import matplotlib.pyplot as plt
from simsopt.mhd.vmec import Vmec
from simsopt.field.boozermagneticfield import BoozerRadialInterpolant
from simsopt.geo import SurfaceRZFourier

wout_file = "wout.nc"

vmec = Vmec(wout_file)
boozer = BoozerRadialInterpolant(vmec, order=3, mpol=16, ntor=16)
# %%
# surf = SurfaceRZFourier.from_wout(wout_file, s=0.5)
# surf.plot(close=True)
# %%
s = np.linspace(0, 1, 100)
th = np.ones_like(s)
ph = np.ones_like(s)

points = np.array([s, th, ph]).T.copy()

boozer.set_points(points)
B = boozer.modB()

plt.figure()
plt.plot(s, boozer.dmodBds(), label="dmodBds")
plt.plot(s, boozer.dmodBdtheta(), label="dmodBdtheta")
plt.plot(s, boozer.dmodBdzeta(), label="dmodBdzeta")
plt.legend()
# %%
import pysimple

tracer = pysimple.simple.Tracer()
pysimple.params.field_input = "wout.nc"
pysimple.simple_main.init_field(tracer, "wout.nc", 5, 5, 5, 1)
# %%
