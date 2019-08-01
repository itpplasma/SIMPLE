using StaticArrays
import Base.show

firstrun = true

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

a = FArray1D([1.0,2.0,3.0])

# mutable struct SymplecticIntegrator
#       atol::Float64
#       rtol::Float64
    
#       r::Float64
#       th::Float64
#       ph::Float64
#       pph::Float64
    
#       pph::Float64
#       ! Buffer for Lagrange polynomial interpolation
#       integer :: kbuf
#       integer :: kt
#       integer :: k
#       integer :: bufind(0:nlag)
#       double precision, dimension(4, nbuf) :: zbuf
#       double precision, dimension(0:0, nlag+1) :: coef
    
#       ! Timestep and variables from z0
#       integer :: ntau
#       double precision :: dt
#       double precision :: pabs
    
#       ! For initial conditions
#       double precision :: z0init(5)
    
#       ! Integrator mode
#       integer :: mode  ! 1 = euler1, 2 = euler2, 3 = verlet
    
#       ! Data for evaluating field, Hamiltonian and canonical poloidal momentum
#       type(FieldCan), pointer :: f
# end

mutable struct NeoOrb
  fper::Float64
  dtau::Float64
  dtaumax::Float64
  v0::Float64
  n_e::Int32
  n_d::Int32
  integmode::Int32
  relerr::Float64
  firstrun::Bool
#  si::Ptr{Cvoid}

  NeoOrb() = new()
end


mutable struct CutDetector
  npl_half::Int32
  alam_prev::Float64
  par_inv::Float64
  iper::Int32
  itip::Int32
  kper::Int32

  orb_sten::SMatrix{6, 6, Float64, 36}
  coef::SMatrix{1, 6, Float64, 6}
  ipoi::SVector{6, Int32}

  norb::NeoOrb

#   fper::Float64
#   dtau::Float64
#   dtaumax::Float64
#   v0::Float64
#   n_e::Int32
#   n_d::Int32
#   integmode::Int32
#   relerr::Float64
#   firstrun::Bool

  CutDetector() = new()
end

if(firstrun)
      norb = NeoOrb()

      ccall((:__neo_orb_MOD_init_field, "lib/libneo_orb"),
      Cvoid, (Ref{NeoOrb}, Ref{Int32}, Ref{Int32}, Ref{Int32}, Ref{Int32}), 
      norb, 5, 5, 3, 0)
end

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

npoiper2 = 64

p_rmajor = cglobal((:__new_vmec_stuff_mod_MOD_rmajor, "lib/libneo_orb"), Float64)
rmajor = unsafe_load(p_rmajor)
println(rmajor)
rbig = rmajor*1.0e2
dtaumax = 2.0*pi*rbig/npoiper2
dtau = dtaumax
ccall((:__neo_orb_MOD_init_params, "lib/libneo_orb"), Cvoid,
      (Ref{NeoOrb}, Ref{Int32}, Ref{Int32}, Ref{Float64}, Ref{Float64},
       Ref{Float64}, Ref{Float64}),
      norb, Z_charge, m_mass, E_kin, dtau, dtaumax, 1e-8)

s = [0.5]
th = [0.0]
ph = [0.314]
lam = [0.22]

z = [0.5, 0.0, 0.314, 1.0, 0.22]
za = FArray1D(z)

ierr = 0

ntimstep = 10000
@time begin
      for i = 1:ntimstep
            ccall((:__neo_orb_MOD_timestep_z, "lib/libneo_orb"), Cvoid,
                   (Ref{NeoOrb}, Ref{FArray1D}, Ref{Int32}), norb, za, ierr)
      end
end

println(z)

ncut = 1000
cutter = CutDetector()
var_tip = zeros(6)
var_tipa = FArray1D(var_tip)
cut_type = 0

println(cutter)

ccall((:__cut_detector_MOD_init, "lib/libneo_orb"), Cvoid,
       (Ref{CutDetector}, Ref{NeoOrb}, Ref{FArray1D}),
       cutter, norb, za)
         
#println(cutter)

@time begin
        for i = 1:ncut
                ccall((:__cut_detector_MOD_trace_to_cut, "lib/libneo_orb"), Cvoid,
                      (Ref{CutDetector}, Ref{FArray1D}, Ref{FArray1D}, Ref{Int32}, Ref{Int32}),
                      cutter, za, var_tipa, cut_type, ierr)
        end
end
