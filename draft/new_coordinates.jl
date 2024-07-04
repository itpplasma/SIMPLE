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
            x[k, 1:2] .= pysimple.get_can_sub.vmec_to_cyl(s[ks], th[kth], ph[kph])
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

# %% https://discourse.julialang.org/t/interpolate-2d-function/96914/7

using BSplineKit
using LinearAlgebra

function fit_spline2D(fdata)
    ord = BSplineOrder(6)  # quintic splines

    # Create B-spline knots based on interpolation points (uses an internal function)
    ts_x = SplineInterpolations.make_knots(xs, ord, nothing)
    ts_y = SplineInterpolations.make_knots(ys, ord, nothing)

    # Create B-spline bases
    Bx = BSplineBasis(ord, ts_x; augment = Val(false))
    By = BSplineBasis(ord, ts_y; augment = Val(false))

    # Create and factorise interpolation matrices
    Cx = lu!(collocation_matrix(Bx, xs))
    Cy = lu!(collocation_matrix(By, ys))

    # 2D B-spline coefficients (output)
    coefs = similar(fdata)

    # Solve linear systems
    for j ∈ eachindex(ys)
        @views ldiv!(coefs[:, j], Cx, fdata[:, j])
    end
    for i ∈ eachindex(xs)
        @views ldiv!(Cy, coefs[i, :])
    end

    return coefs
end

# Evaluate 2D (tensor product) spline at point (x, y).
function eval_spline2D(coefs::AbstractMatrix, (Bx, By), (x, y))
    i, bx = Bx(x)  # evaluate B-spline basis at point x
    j, by = By(y)
    kx = order(Bx)
    ky = order(By)
    val = zero(eltype(coefs))  # spline evaluated at (x, y)
    for δj ∈ 1:ky, δi ∈ 1:kx
        coef = coefs[i - δi + 1, j - δj + 1]
        bi = bx[δi]
        bj = by[δj]
        val += coef * bi * bj
    end
    val
end


f(x, y) = 2.0 + sinpi(x) * cospi(y)

xs = 0:0.1:1
ys = 0:0.1:2
fdata = f.(xs, ys')

coefs = spline2D(fdata)

# Verification: evaluate 2D spline at data points
for m ∈ eachindex(ys), n ∈ eachindex(xs)
    x = xs[n]
    y = ys[m]
    val = eval_spline2D(coefs, (Bx, By), (x, y))
    fval = fdata[n, m]
    @assert val ≈ fval  # check that the value is correct
end

# Plot
using Plots
plot(fdata, st = :surface, color = :viridis, legend = false)
