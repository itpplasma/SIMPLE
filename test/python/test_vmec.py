import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

from pysimple import new_vmec_stuff_mod as v
from pysimple import spline_vmec_sub, vmecin_sub, vmec_alloc_sub, get_can_sub


v.netcdffile = 'wout.nc'
v.ns_s = 5
v.ns_tp = 5
v.multharm = 7

spline_vmec_sub.spline_vmec_data()

vmec_alloc_sub.new_allocate_vmec_stuff()

n_flux_surf = v.nsurfm
n_harmonics = v.nstrm
assert(v.rmnc.shape == (n_harmonics,n_flux_surf))

torflux = vmecin_sub.vmecin(v.rmnc,v.zmns,v.almns,v.rmns,v.zmnc,v.almnc,v.aiota,
    v.phi,v.sps,v.axm,v.axn,v.s,v.nsurfm,v.nstrm,v.kpar)
assert(np.allclose(torflux, v.phi[-1]))

#%%
def vmec_to_cyl(ks, th, ph):
    R, Z = get_can_sub.vmec_to_cyl(v.s[ks], th, ph)

    P = ph
    return R, P, Z

def vmec_to_cyl_ref(ks, th, ph):
    R = 0.0
    P = 0.0
    Z = 0.0

    for kmn in range(v.nstrm):
        m = int(v.axm[kmn])
        n = int(v.axn[kmn])
        cosphase = np.cos(m*th - n*ph)
        sinphase = np.sin(m*th - n*ph)
        R += v.rmnc[kmn,ks]*cosphase + v.rmns[kmn,ks]*sinphase
        P = ph
        Z += v.zmnc[kmn,ks]*cosphase + v.zmns[kmn,ks]*sinphase
    return R, P, Z

def cyl_to_cart(R, P, Z):
    x = R*np.cos(P)
    y = R*np.sin(P)
    z = Z
    return x, y, z


ks = 5

th = np.linspace(0,2*np.pi,10)
ph = np.linspace(0,2*np.pi,11)

TH, PH = np.meshgrid(th, ph)
TH = TH.flatten()
PH = PH.flatten()

fig, ax = plt.subplots(subplot_kw={'projection':'3d'})
for k, thk in enumerate(TH):
    phk = PH[k]

    R, P, Z = vmec_to_cyl(ks, thk, phk)
    x, y, z = cyl_to_cart(R, P, Z)

    ax.plot(x, y, z, '.')

    R, P, Z = vmec_to_cyl_ref(ks, thk, phk)
    x_ref, y_ref, z_ref = cyl_to_cart(R, P, Z)
    ax.plot(x_ref, y_ref, z_ref, 'x')

    assert(np.allclose([x, y, z], [x_ref, y_ref, z_ref]))

plt.axis('equal')
plt.savefig('test_vmec.png')
