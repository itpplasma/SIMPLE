# %%
import numpy as np
import matplotlib.pyplot as plt
from simsopt.mhd.vmec import Vmec
from simsopt.field.boozermagneticfield import BoozerRadialInterpolant
import pysimple


GAUSS_IN_TESLA = 1e-4
CM_IN_M = 1e-2

wout_file = "wout.nc"

BOOZER = pysimple.magfie_sub.boozer

tracer = pysimple.simple.Tracer()
pysimple.params.read_config("simple.in")
pysimple.simple_main.init_field(tracer, "wout.nc", 5, 5, 5, 1)

# %% Check VMEC to real space
s = 0.5 * np.ones(100)
th = np.linspace(0, 2 * np.pi, 100)
ph = np.zeros(100)
R = np.empty_like(s)
Z = np.empty_like(s)
for i in range(len(s)):
    R[i], Z[i] = pysimple.get_can_sub.vmec_to_cyl(s[i], th[i], ph[i])
R = np.array(R) * CM_IN_M
Z = np.array(Z) * CM_IN_M

plt.figure()
plt.plot(R[0], Z[0], "o")
plt.quiver(
    R[:-1],
    Z[:-1],
    R[1:] - R[:-1],
    Z[1:] - Z[:-1],
    angles="xy",
    scale_units="xy",
    scale=1,
    color="r",
)

#%%
vmec = Vmec(wout_file)
boozer = BoozerRadialInterpolant(vmec, order=3, mpol=16, ntor=16)
points = np.array([s, th, ph]).T.copy()
boozer.set_points(points)
#%%
R = boozer.R()
Z = boozer.Z()

plt.figure()
plt.plot(R[0], Z[0], "o")
plt.quiver(
    R[:-1],
    Z[:-1],
    R[1:] - R[:-1],
    Z[1:] - Z[:-1],
    angles="xy",
    scale_units="xy",
    scale=1,
    color="r",
)

# %%
bder = np.empty(3)
hcovar = np.empty(3)
hctrvr = np.empty(3)
hcurl = np.empty(3)

s = np.linspace(0, 1, 100)
th = np.ones_like(s)
ph = np.ones_like(s)


def evaluate_magfie_simple(s, th, ph, magfie_routine):
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
        modB[i], sqrtg[i] = magfie_routine(x, bder, hcovar, hctrvr, hcurl)
        modB[i] = modB[i] * GAUSS_IN_TESLA
        sqrtg[i] = sqrtg[i] * CM_IN_M**3
        dmodBds[i] = bder[0] * modB[i]
        dmodBdtheta[i] = bder[1] * modB[i]
        dmodBdzeta[i] = bder[2] * modB[i]
        Bscov[i] = hcovar[0] * modB[i] * CM_IN_M
        Bthetacov[i] = hcovar[1] * modB[i] * CM_IN_M
        Bzetacov[i] = hcovar[2] * modB[i] * CM_IN_M

    return (
        modB,
        sqrtg,
        dmodBds,
        dmodBdtheta,
        dmodBdzeta,
        Bscov,
        Bthetacov,
        Bzetacov,
    )


print("Evaluating SIMPLE VMEC magnetic field")

modB, sqrtg, dmodBds, dmodBdtheta, dmodBdzeta, Bscov, Bthetacov, Bzetacov = (
    evaluate_magfie_simple(s, th, ph, pysimple.magfie_sub.magfie_vmec)
)

print("modB = ", modB[5])
print("sqrtg = ", sqrtg[5])
print("Bthetacov = ", Bthetacov[5])
print("Bzetacov = ", Bzetacov[5])

print("Evaluating SIMPLE Boozer magnetic field")

modB, sqrtg, dmodBds, dmodBdtheta, dmodBdzeta, Bscov, Bthetacov, Bzetacov = (
    evaluate_magfie_simple(s, th, ph, pysimple.magfie_sub.magfie_boozer)
)

print("modB = ", modB[5])
print("sqrtg = ", sqrtg[5])
print("Bthetacov = ", Bthetacov[5])
print("Bzetacov = ", Bzetacov[5])
# %%
# %%

points = np.array([s, th, ph]).T.copy()

boozer.set_points(points)
modB_simsopt = boozer.modB()
sqrtg_simsopt = (boozer.G() + boozer.iota() * boozer.I()) / modB_simsopt**2

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
plt.plot(s, -sqrtg_simsopt, label="-sqrtg SIMSOPT")
plt.plot(s, sqrtg / boozer.psi0, "--", label="sqrtg SIMPLE")
# plt.gca().set_yscale("log")
plt.xlabel("s")
plt.ylabel("sqrtg [m^3]")
plt.legend()
plt.show()

# %%
plt.figure()
plt.plot(s, boozer.iota(), label="iota SIMSOPT")

# %%
