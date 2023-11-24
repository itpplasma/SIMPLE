# %%
using Plots
using PyCall
println(pwd())

py"""
import sys
sys.path.insert(0, "")
"""

# %%
@pyimport pysimple

# %%
tracy = pysimple.params.Tracer()
pysimple.simple.init_field(tracy, "wout.nc", 3, 3, 3, 1)
pysimple.simple.init_params(tracy, 2, 4, 3.5e6, 256, 1, 1e-13)

# %%
ns = 16
nth = 32
nph = 1

s = range(1.0/ns, stop=1.0, length=ns)
th = range(0.0, stop=2π - 2π/nth, length=nth)
ph = range(0.0, stop=2π - 2π/nph, length=nph)

x = zeros(ns*nth*nph, 3)

bder = zeros(3)
hcovar = zeros(3)
hctrvr = zeros(3)
hcurl = zeros(3)

bmod = 0.0
sqrtg = 0.0

bmods = zeros(ns*nth*nph)

k = 1
for ks in 1:ns
    for kth in 1:nth
        for kph in 1:nph
            x[k, 1:2] .= pysimple.get_canonical_coordinates_sub.vmec_to_cyl(s[ks], th[kth], ph[kph])
            x[k, 3] = ph[kph]

            A_phi, A_theta, dA_phi_ds, dA_theta_ds, aiota, R, Z, alam, dR_ds, dR_dt, dR_dp, dZ_ds, dZ_dt, dZ_dp, dl_ds, dl_dt, dl_dp = pysimple.splint_vmec_data(
                s[ks], th[kth], ph[kph])

            bmod, sqrtg = pysimple.magfie_sub.magfie_vmec([s[ks], th[kth], ph[kph]], bder, hcovar, hctrvr, hcurl)

            bmods[k] = bmod

            k += 1
        end
    end
end


scatter(x[:, 1], x[:, 2], st=:scatter, aspect_ratio=:equal,
    legend=false, marker_z=bmods, c=:viridis)

# %%
