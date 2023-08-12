#
# Julia. It supports UTF-8. (╯°□°）╯︵ ┻┻
#

using Plots: plot
using StaticArrays: SVector, SMatrix
using Printf: @printf
import Base.show

struct FArrayDims
    fstride::Cptrdiff_t
    lower_bound::Cptrdiff_t
    upper_bound::Cptrdiff_t
end

struct FDataType
    len::Csize_t
    ver::Cint
    rank::Cchar
    type::Cchar
    attribute::Cshort
end

struct FArray1D
    base_addr::Ptr{Cvoid}
    offset::Csize_t
    dtype::FDataType
    span::Cptrdiff_t
    dim::FArrayDims

    FArray1D(vector::Array{Float64,1}) = new(
                pointer(vector),
                0,
                FDataType(8, 0, 1, 3, 0),
                length(vector)*8,
                FArrayDims(1, 1, length(vector))
        )
end

struct FArray2D
    base_addr::Ptr{Cvoid}
    offset::Csize_t
    dtype::FDataType
    span::Cptrdiff_t
    dim1::FArrayDims
    dim2::FArrayDims

    FArray2D(arr::Array{Float64,2}) = new(
                pointer(arr),
                0,
                FDataType(8, 0, 2, 3, 0),
                size(arr)[0]*size(arr)[1]*8,
                FArrayDims(1, 1, size(arr)[0]),
                FArrayDims(size(arr)[0], 1, size(arr)[1])
        )
end


function show(io::IO, farr::FArray1D)
    data = unsafe_wrap(Array{Float64,1}, Ptr{Float64}(farr.base_addr),
        (farr.dim.upper_bound - farr.dim.lower_bound + 1,))
    print(io, data)
end


mutable struct Tracer
  fper::Float64
  dtau::Float64
  dtaumin::Float64
  v0::Float64
  n_e::Int32
  n_d::Int32
  integmode::Int32
  relerr::Float64

  #
  # The following structs cannot be stored in-line in a mutable struct.
  # This is why they are now inlined explicitly
  #

  ## f::FieldCan
  field_type::Int32  # -1: testing, 0: canonical, 2: Boozer

  Ath::Float64
  Aph::Float64
  hth::Float64
  hph::Float64
  Bmod::Float64

  dAth::SVector{3, Float64}
  dAph::SVector{3, Float64}
  dhth::SVector{3, Float64}
  dhph::SVector{3, Float64}
  dBmod::SVector{3, Float64}

  # second derivatives: drdr, drdhth, drdph, dthdth, dthdph, dphdph
  d2Ath::SVector{6, Float64}
  d2Aph::SVector{6, Float64}
  d2hth::SVector{6, Float64}
  d2hph::SVector{6, Float64}
  d2Bmod::SVector{6, Float64}

  H::Float64
  pth::Float64
  vpar::Float64

  dvpar::SVector{4, Float64}
  dH::SVector{4, Float64}
  dpth::SVector{4, Float64}

  # order of second derivatives:
  # d2dr2, d2drdth, d2drph, d2dth2, d2dthdph, d2dph2,
  # d2dpphdr, d2dpphdth, d2dpphdph, d2dpph2
  d2vpar::SVector{10, Float64}
  d2H::SVector{10, Float64}
  d2pth::SVector{10, Float64}

  mu::Float64
  ro0::Float64


  ## si::SymplecticIntegrator
  atol::Float64
  rtol::Float64

  z::SVector{4, Float64}
  pthold::Float64

  # Buffer for Lagrange polynomial interpolation
  kbuf::Int32
  kt::Int32
  k::Int32
  bufind::SVector{4, Int32}
  zbuf::SMatrix{4, 16, Float64, 64}
  coef::SMatrix{1, 4, Float64, 4}

  # Timestep and variables from z0
  ntau::Int32
  dt::Float64
  pabs::Float64

  # Integrator mode
  mode::Int32  # 1 = euler1, 2 = euler2, 3 = verlet


  ## mi::MultistageIntegrator
  s::Int32
  alpha::SVector{64, Float64}
  beta::SVector{64, Float64}
  stages::SVector{42496, Int8}  # Static 64*664 array of SymplecticIntegrators

  Tracer() = new()
end

# TODO
# mutable struct CutDetector
#   npl_half::Int32
#   alam_prev::Float64
#   par_inv::Float64
#   iper::Int32
#   itip::Int32
#   kper::Int32

#   orb_sten::SMatrix{6, 6, Float64, 36}
#   coef::SMatrix{1, 6, Float64, 6}
#   ipoi::SVector{6, Int32}

#   tracer::Tracer

# #   fper::Float64
# #   dtau::Float64
# #   dtaumax::Float64
# #   v0::Float64
# #   n_e::Int32
# #   n_d::Int32
# #   integmode::Int32
# #   relerr::Float64
# #   firstrun::Bool

#   CutDetector() = new()
# end


tracer = Tracer()
ccall((:__simple_MOD_init_field, "../build/libsimple"),
Cvoid, (Ref{Tracer}, Cstring, Ref{Int32}, Ref{Int32}, Ref{Int32}, Ref{Int32}),
tracer, "wout.nc", 3, 3, 3, 1)

Z_charge = 2
m_mass = 4
c=2.9979e10
e_charge=4.8032e-10
e_mass=9.1094e-28
p_mass=1.6726e-24
ev=1.6022e-12

bmod_ref = 5.0e4
Z_charge = 2
m_mass = 4
E_kin = 3.5e6
trace_time = 1.0e-2

v0 = sqrt(2.0*E_kin*ev/(m_mass*p_mass))    # alpha particle velocity, cm/s
rlarm=v0*m_mass*p_mass*c/(Z_charge*e_charge*bmod_ref) # reference gyroradius
tau=trace_time*v0

npoiper2 = 128

p_rmajor = cglobal((:__new_vmec_stuff_mod_MOD_rmajor, "../build/libsimple"), Float64)
rmajor = unsafe_load(p_rmajor)
println(rmajor)
rbig = rmajor*1.0e2
dtaumin = 2.0*pi*rbig/npoiper2
dtau = dtaumin

ccall((:__simple_MOD_init_params, "../build/libsimple"), Cvoid,
      (Ref{Tracer}, Ref{Int32}, Ref{Int32}, Ref{Float64}, Ref{Int32},
       Ref{Int32}, Ref{Float64}),
      tracer, Z_charge, m_mass, E_kin, npoiper2, 1, 1e-13)

s = [0.5]
th = [0.0]
ph = [0.314]
lam = [0.22]


th_c = [0.0]
ph_c = [0.314]


ccall((:vmec_to_can_, "../build/libsimple"), Cvoid,
      (Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}),
      s, th, ph, th_c, ph_c)

z0_can = [0.5, th_c[1], ph_c[1], 1.0, 0.22]

println("VMEC coordinates:")
println([s[1], th[1], ph[1], 1.0, lam[1]])
println("Canonical coordinates:")
println(z0_can)

z0_can_array = FArray1D(z0_can)

println("Before init_integrator:")
print("dAth = ")
println(tracer.dAth)
@printf "B = %e\n" tracer.Bmod

ccall((:__simple_MOD_init_integrator, "../build/libsimple"), Cvoid,
      (Ref{Tracer}, Ref{FArray1D}),
      tracer, z0_can_array)

ierr = 0


# This shows that Fortran can fiddle in the formally immutable SArray.
# Since Fortran changes the contents of the stack variables,
# this hack is required (MArray will not work).
println("After init_integrator:")
print("dAth = ")
println(tracer.dAth)
@printf "B = %e\n" tracer.Bmod
