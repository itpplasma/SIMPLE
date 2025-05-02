#%%
import numpy as np
import matplotlib.pyplot as plt

GAUSS_IN_TESLA = 1e-4
CM_IN_M = 1e-2

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
sqrtg = np.empty(len(s))
dmodBds = np.empty(len(s))
dmodBdtheta = np.empty(len(s))
dmodBdzeta = np.empty(len(s))
Bscov = np.empty(len(s))
Bthetacov = np.empty(len(s))
Bzetacov = np.empty(len(s))
for i in range(len(s)):
    x = np.array([s[i], th[i], ph[i]])
    modB[i], sqrtg[i] = pysimple.magfie_sub.magfie_boozer(x, bder, hcovar, hctrvr, hcurl)
    modB[i] = modB[i]*GAUSS_IN_TESLA
    sqrtg[i] = sqrtg[i]*CM_IN_M**3
    dmodBds[i] = bder[0]*modB[i]
    dmodBdtheta[i] = bder[1]*modB[i]
    dmodBdzeta[i] = bder[2]*modB[i]
    Bscov[i] = hcovar[0]*modB[i]*CM_IN_M
    Bthetacov[i] = hcovar[1]*modB[i]*CM_IN_M
    Bzetacov[i] = hcovar[2]*modB[i]*CM_IN_M

# %%
from simsopt.mhd.vmec import Vmec
from simsopt.field.boozermagneticfield import BoozerRadialInterpolant

vmec = Vmec(wout_file)
boozer = BoozerRadialInterpolant(vmec, order=3, mpol=16, ntor=16)
# %%

points = np.array([s, th, ph]).T.copy()

boozer.set_points(points)
modB_simsopt = boozer.modB()
sqrtg_simsopt = (boozer.G()+boozer.iota()*boozer.I())/modB_simsopt**2

# %%
plt.figure()
plt.plot(s, boozer.dmodBds(), label="dmodBds SIMSOPT")
plt.plot(s, boozer.dmodBdtheta(), label="dmodBdtheta SIMSOPT")
plt.plot(s, boozer.dmodBdzeta(), label="dmodBdzeta SIMSOPT")
plt.plot(s, dmodBds, "--", label="dmodBds SIMPLE")
plt.plot(s, dmodBdtheta, "--", label="dmodBdtheta SIMPLE")
plt.plot(s, dmodBdzeta, "--", label="dmodBdzeta SIMPLE")
plt.gca().set_yscale("symlog", linthresh=1e-2)
plt.xlabel("s")
plt.ylabel("|B| derivatives [T]")
plt.legend()

# %%
plt.figure()
plt.plot(s, boozer.K(), label="Bscov SIMSOPT")
plt.plot(s, boozer.I(), label="Bthetacov SIMSOPT")
plt.plot(s, boozer.G(), label="Bzetacov SIMSOPT")
plt.plot(s, Bscov, "--", label="Bscov SIMPLE (set to zero)")
plt.plot(s, Bthetacov, "--", label="Bthetacov SIMPLE")
plt.plot(s, Bzetacov, "--", label="Bzetacov SIMPLE")
plt.gca().set_yscale("symlog", linthresh=1e-5)
plt.xlabel("s")
plt.ylabel("|B| derivatives [T]")
plt.legend()

# %%
plt.figure()
plt.plot(s, sqrtg_simsopt, label="sqrtg SIMSOPT")
plt.plot(s, -sqrtg/boozer.psi0, "--", label="-sqrtg SIMPLE")
#plt.gca().set_yscale("log")
plt.xlabel("s")
plt.ylabel("sqrtg [m^3]")
plt.legend()


# %%
