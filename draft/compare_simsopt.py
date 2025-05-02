#%%
import numpy as np
import matplotlib.pyplot as plt

wout_file = "wout.nc"
# %%
import pysimple
BOOZER = pysimple.magfie_sub.boozer

tracer = pysimple.simple.Tracer()
pysimple.params.read_config("simple.in")
pysimple.simple_main.init_field(tracer, "wout.nc", 5, 5, 5, 1)
# %%
bder = np.empty(3)
hcovar = np.empty(3)
hctrvr = np.empty(3)
hcurl = np.empty(3)

s = np.linspace(0, 1, 100)
th = np.ones_like(s)
ph = np.ones_like(s)

modB = np.empty(len(s))
dmodBds = np.empty(len(s))
dmodBdtheta = np.empty(len(s))
dmodBdzeta = np.empty(len(s))
for i in range(len(s)):
    x = np.array([s[i], th[i], ph[i]])
    bmod, sqrtg = pysimple.magfie_sub.magfie_boozer(x, bder, hcovar, hctrvr, hcurl)
    modB[i] = bmod
    dmodBds[i] = bder[0]
    dmodBdtheta[i] = bder[1]
    dmodBdzeta[i] = bder[2]

plt.figure()
plt.plot(s, dmodBds, label="dmodBds")
plt.plot(s, dmodBdtheta, label="dmodBdtheta")
plt.plot(s, dmodBdzeta, label="dmodBdzeta")
plt.legend()


# %%
from simsopt.mhd.vmec import Vmec
from simsopt.field.boozermagneticfield import BoozerRadialInterpolant
from simsopt.geo import SurfaceRZFourier

vmec = Vmec(wout_file)
boozer = BoozerRadialInterpolant(vmec, order=3, mpol=16, ntor=16)
# %%
# surf = SurfaceRZFourier.from_wout(wout_file, s=0.5)
# surf.plot(close=True)
# %%

points = np.array([s, th, ph]).T.copy()

boozer.set_points(points)
B = boozer.modB()

plt.figure()
plt.plot(s, boozer.dmodBds(), label="dmodBds")
plt.plot(s, boozer.dmodBdtheta(), label="dmodBdtheta")
plt.plot(s, boozer.dmodBdzeta(), label="dmodBdzeta")
plt.legend()

# %%
