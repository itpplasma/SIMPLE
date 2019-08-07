using Plots
using StaticArrays

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
  
    NeoOrb() = new()
  end

mutable struct FieldCan
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

    FieldCan() = new()
end

norb = NeoOrb()
if !@isdefined(firstrun) | firstrun == true
    global firstrun, norb
    ccall((:__neo_orb_MOD_init_field, "lib/libneo_orb"),
    Cvoid, (Ref{NeoOrb}, Ref{Int32}, Ref{Int32}, Ref{Int32}, Ref{Int32}), 
    norb, 5, 5, 3, 0)
    firstrun = false
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
println("R0 = ", rmajor)
rbig = rmajor*1.0e2
dtaumin = 2.0*pi*rbig/npoiper2
dtau = dtaumin
ccall((:__neo_orb_MOD_init_params, "lib/libneo_orb"), Cvoid,
      (Ref{NeoOrb}, Ref{Int32}, Ref{Int32}, Ref{Float64}, Ref{Float64},
       Ref{Float64}, Ref{Float64}),
      norb, Z_charge, m_mass, E_kin, dtau, dtaumin, 1e-8)


s = [0.5]
th = [0.0]
ph = [0.314]
pabs = [1.0]
lam = [0.22]

z0 = [0.5, 0.0, 0.314, 1.0, 0.22]

rtol = 1e-12

kbuf = 0
kt = 0
k = 0

dt = dtaumin/sqrt(2e0) # factor 1/sqrt(2) due to velocity normalisation different from other modules

f = FieldCan()
ccall((:__field_can_mod_MOD_eval_field_can, "lib/libneo_orb"), Cvoid,
        (Ref{FieldCan}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Int32}),
        f, s, th, ph, 0)

mu = 0.5 * z0[4] * (1.0 - z0[5]^2)/f.Bmod*2.0 # mu by factor 2 from other modules

ro0_parmot = unsafe_load(
    cglobal((:__parmot_mod_MOD_ro0, "lib/libneo_orb"), Float64)
)

f.ro0 = ro0_parmot/sqrt(2.0) # ro0 = mc/e*v0, different by sqrt(2) from other modules
f.vpar = z0[4]*z0[5]*sqrt(2.0) # vpar_bar = vpar/sqrt(T/m), different by sqrt(2) from other modules

# ! initialize canonical variables
z = zeros(4)
z[1:3] = z0[1:3]  # r, th, ph
z[4] = f.vpar*f.hph + f.Aph/f.ro0
ccall((:__field_can_mod_MOD_get_val, "lib/libneo_orb"), Cvoid,
        (Ref{FieldCan}, Ref{Float64}),
        f, z[4])  # for pth

# call plag_coeff(nlag+1, 0, 1d0, 1d0*(/(k,k=-nlag,0)/), si%coef)

# ccall((:__field_can_mod_MOD_get_val, "lib/libneo_orb"), Cvoid,
#                    (Ref{FieldCan}, Ref{Float64}), f, pphi)

function timestep!(z)

end