"""
Plot field standalone
"""
##
using Plots
using NetCDF

function fourier_cos(fmn, theta, zeta, xm, xn)
    f = zero(theta)
    for k in eachindex(fmn)
        f .+= fmn[k]*cos.(xm[k]*theta - xn[k]*zeta)
    end
    return f
end
##
ncfile = "wout.nc"

println(ncinfo(ncfile))
iotaf = ncread(ncfile, "iotaf")
iotas = ncread(ncfile, "iotas")
psi_tor = ncread(ncfile, "phi") 
s = psi_tor/psi_tor[end]
sh = 0.5*(s[2:end] + s[1:end-1])

##

plt = plot(s, iotaf)
plot!(plt, sh, iotas[2:end])
display(plt)
##

rmnc = ncread(ncfile, "rmnc")
zmns = ncread(ncfile, "zmns")
lmns = ncread(ncfile, "lmns")

##
ks = 65
nth = 72
nph = 36

th = range(-pi/2, step=2*pi/nth, length=nth)
ph = range(0, step=pi/nph, length=nph)


##

[PH, TH] = np.meshgrid(ph, th)
TH = TH.flatten()
PH = PH.flatten()
TH2 = zeros_like(TH)


# %%

mn_mode_nyq = data.dimensions['mn_mode_nyq']
xm_nyq = data.variables['xm_nyq']
xn_nyq = data.variables['xn_nyq']
bmnc = data.variables['bmnc']
bsubumnc = data.variables['bsubumnc']
bsubvmnc = data.variables['bsubvmnc']
bsupumnc = data.variables['bsupumnc']
bsupvmnc = data.variables['bsupvmnc']


b = np.zeros(TH.shape)
bsubu = np.zeros(TH.shape)
bsubv = np.zeros(TH.shape)
bsupu = np.zeros(TH.shape)
bsupv = np.zeros(TH.shape)
R = np.zeros(TH.shape)
Z = np.zeros(TH.shape)

for k in np.arange(mn_mode_nyq.size):
    m = xm_nyq[k]
    n = xn_nyq[k]
    fac = np.cos(m*TH - n*PH)
    b += bmnc[ks, k]*fac
    bsubu += bsubumnc[ks, k]*fac
    bsubv += bsubvmnc[ks, k]*fac
    bsupu += bsupumnc[ks, k]*fac
    bsupv += bsupvmnc[ks, k]*fac

# %%

figure(figsize=(3.2, 3.2))
contour(PH.reshape(nth, nph), TH.reshape(nth, nph),
        b.reshape(nth, nph)/min(b), cmap='jet', levels=100)
xlabel(r'$\varphi$')
ylabel(r'$\vartheta$')

# %%
b2 = np.sqrt(bsubu*bsupu + bsubv*bsupv)

figure(figsize=(3.2, 3.2))
contour(PH.reshape(nth, nph), TH.reshape(nth, nph),
        b2.reshape(nth, nph)/min(b), cmap='jet', levels=100)
xlabel(r'$\varphi$')
ylabel(r'$\vartheta$')

# rootgrp.close()
