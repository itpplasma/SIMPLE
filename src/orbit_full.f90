module orbit_full
  ! Full-orbit (gyro-resolved Lorentz) pushers for SIMPLE, alongside the
  ! existing symplectic guiding-center tracer. CGS Gaussian units throughout
  ! (see src/util.F90): the Lorentz force carries the explicit 1/c,
  !   m dv/dt = (q/c) v x B.
  !
  ! ORBIT_BORIS is the explicit gyro-resolved pusher (this module).
  ! ORBIT_PAULI is the variational CPP scheme; the seam (state slots, provider
  ! canonical-field method, non-fatal convergence error) is in place but the
  ! step is not implemented here yet.
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use orbit_full_provider, only: field_metric_provider_t, &
      FO_OK, FO_ERR_FIELD, FO_ERR_NO_CONVERGE, FO_ERR_OUT_OF_DOMAIN
  use orbit_full_mock_cart, only: cartesian_provider_t
  use neo_biotsavart, only: coils_t
  use util, only: c, p_mass, e_charge, twopi
  implicit none
  private

  ! orbit models (0 reserved for the existing symplectic guiding-center path)
  integer, parameter, public :: ORBIT_GC    = 0
  integer, parameter, public :: ORBIT_PAULI = 1   ! CPP variational, big dt, implicit
  integer, parameter, public :: ORBIT_BORIS = 2   ! gyro-resolved Lorentz, explicit

  ! coordinate kinds (3..5 reserved for the libneo PR: VMEC, Boozer, chartmap)
  integer, parameter, public :: COORD_CART = 1
  integer, parameter, public :: COORD_CYL  = 2

  ! error codes re-exported from the provider module
  public :: FO_OK, FO_ERR_FIELD, FO_ERR_NO_CONVERGE, FO_ERR_OUT_OF_DOMAIN

  ! Name kept as FullOrbitState (no _t) for API compatibility with
  ! test_full_orbit.f90.
  type, public :: FullOrbitState
    integer  :: orbit_model = ORBIT_BORIS
    integer  :: ncoord      = COORD_CART
    ! Boris/Cartesian:   z = (x1,x2,x3, v1,v2,v3) physical position and velocity.
    ! Boris/curvilinear: z = (u1,u2,u3, v1,v2,v3) with v contravariant.
    ! Pauli (CPP):       z = (u1,u2,u3, vpar,0,0); z(1:4) carry state.
    real(dp) :: z(6)   = 0.0_dp
    real(dp) :: mu     = 0.0_dp   ! m*vperp^2/(2 B), adiabatic invariant
    real(dp) :: dt     = 0.0_dp
    real(dp) :: qm     = 0.0_dp   ! charge/mass (esu/g)
    real(dp) :: mass   = 0.0_dp   ! particle mass (g)
    real(dp) :: charge = 0.0_dp   ! particle charge (esu)
    class(field_metric_provider_t), pointer :: prov => null()
  end type FullOrbitState

  public :: init_full_orbit_state, timestep_full_orbit, &
            convert_full_to_gc, compute_energy

  interface init_full_orbit_state
    module procedure init_full_orbit_state_prov
    module procedure init_full_orbit_state_coils
  end interface init_full_orbit_state

contains

  subroutine init_full_orbit_state_prov(state, x0, v0, orbit_model, ncoord, &
                                        mass, charge, dt, prov)
    type(FullOrbitState), intent(out) :: state
    real(dp), intent(in) :: x0(3), v0(3)
    integer,  intent(in) :: orbit_model, ncoord
    real(dp), intent(in) :: mass, charge, dt
    class(field_metric_provider_t), intent(in), target :: prov

    state%orbit_model = orbit_model
    state%ncoord      = ncoord
    state%mass        = mass
    state%charge      = charge
    state%qm          = charge / mass
    state%dt          = dt
    state%prov        => prov
    state%z(1:3)      = x0
    state%z(4:6)      = v0
    call set_mu_from_state(state)
  end subroutine init_full_orbit_state_prov

  ! Convenience init matching test_full_orbit.f90: flux-like IC (s,theta,phi)
  ! plus a coil set. Builds a Cartesian Biot-Savart provider, places the
  ! particle on a circle of radius proportional to s, and sets the velocity
  ! from pitch lambda=vpar/v with speed v. mass/charge are atomic units here
  ! (e.g. 4.0 amu, 2.0 e) and converted to CGS internally.
  subroutine init_full_orbit_state_coils(state, s, theta, phi, lambda, v, &
                                         orbit_model, mass_amu, charge_e, dt, coils)
    type(FullOrbitState), intent(out) :: state
    real(dp), intent(in) :: s, theta, phi, lambda, v
    integer,  intent(in) :: orbit_model
    real(dp), intent(in) :: mass_amu, charge_e, dt
    type(coils_t), intent(in), target :: coils
    type(cartesian_provider_t), allocatable, save :: cart_prov
    real(dp) :: mass_cgs, charge_cgs, vpar, vperp, R, x0(3), v0(3)
    real(dp), parameter :: R_AXIS = 1200.0_dp  ! cm, reactor-scale placeholder

    if (allocated(cart_prov)) deallocate(cart_prov)
    allocate(cart_prov)
    cart_prov%field_kind = 3                     ! FIELD_COILS
    cart_prov%coils => coils

    mass_cgs   = mass_amu * p_mass
    charge_cgs = charge_e * e_charge

    R = R_AXIS * (1.0_dp + 0.1_dp * (s - 0.5_dp))
    x0 = [R * cos(phi), R * sin(phi), 0.0_dp]
    vpar  = lambda * v
    vperp = sqrt(max(v*v - vpar*vpar, 0.0_dp))
    ! Velocity: vpar along toroidal e_phi, vperp along e_R for a generic start.
    v0 = [ vperp * cos(phi) - vpar * sin(phi) * sin(theta), &
           vperp * sin(phi) + vpar * cos(phi) * sin(theta), &
           vpar  * cos(theta) ]

    call init_full_orbit_state_prov(state, x0, v0, orbit_model, COORD_CART, &
                                    mass_cgs, charge_cgs, dt, cart_prov)
  end subroutine init_full_orbit_state_coils

  subroutine set_mu_from_state(state)
    type(FullOrbitState), intent(inout) :: state
    real(dp) :: Bvec(3), Bmod, hcov(3), vperp2, vpar
    integer :: ierr

    call state%prov%eval_field(state%z(1:3), Bvec, Bmod, hcov, ierr)
    if (ierr /= FO_OK .or. Bmod <= 0.0_dp) then
      state%mu = 0.0_dp
      return
    end if
    select case (state%ncoord)
    case (COORD_CYL)
      vperp2 = vperp_sq_cyl(state%z(1), state%z(4:6), hcov)
    case default
      vpar   = dot_product(state%z(4:6), hcov)
      vperp2 = max(dot_product(state%z(4:6), state%z(4:6)) - vpar*vpar, 0.0_dp)
    end select
    state%mu = state%mass * vperp2 / (2.0_dp * Bmod)
  end subroutine set_mu_from_state

  subroutine timestep_full_orbit(state, ierr)
    type(FullOrbitState), intent(inout) :: state
    integer, intent(out) :: ierr

    select case (state%orbit_model)
    case (ORBIT_BORIS)
      select case (state%ncoord)
      case (COORD_CART)
        call boris_step_cart(state, ierr)
      case (COORD_CYL)
        call boris_step_cyl(state, ierr)
      case default
        ierr = FO_ERR_OUT_OF_DOMAIN
      end select
    case (ORBIT_PAULI)
      ! CPP variational step not implemented yet. The seam (canonical-field
      ! provider method, vpar in z(4), mu, FO_ERR_NO_CONVERGE channel) is in
      ! place; flag not-implemented without touching the Boris branch.
      ierr = FO_ERR_NO_CONVERGE
    case default
      ierr = FO_ERR_OUT_OF_DOMAIN
    end select
  end subroutine timestep_full_orbit

  ! Cartesian Boris, drift-kick-drift (leapfrog), CGS:
  !   m dv/dt = (q/c) v x B,  Omega = q B/(m c).
  subroutine boris_step_cart(state, ierr)
    type(FullOrbitState), intent(inout) :: state
    integer, intent(out) :: ierr
    real(dp) :: x(3), v(3), Bvec(3), Bmod, hcov(3)
    real(dp) :: t(3), s(3), vprime(3), tmag2
    real(dp) :: Phi, dPhi(3), Efield(3)
    real(dp) :: dt, qmc

    dt  = state%dt
    qmc = state%qm / c
    x = state%z(1:3)
    v = state%z(4:6)

    ! half position drift
    x = x + 0.5_dp * dt * v

    call state%prov%eval_field(x, Bvec, Bmod, hcov, ierr)
    if (ierr /= FO_OK) return

    ! half electric kick (E = -grad Phi); mocks return Phi=0
    call state%prov%eval_potential(x, Phi, dPhi)
    Efield = -dPhi
    v = v + 0.5_dp * dt * state%qm * Efield

    ! magnetic rotation (exact for constant B over the step)
    t = qmc * Bvec * 0.5_dp * dt
    tmag2 = dot_product(t, t)
    s = 2.0_dp * t / (1.0_dp + tmag2)
    vprime = v + cross(v, t)
    v = v + cross(vprime, s)

    ! second half electric kick
    v = v + 0.5_dp * dt * state%qm * Efield

    ! half position drift
    x = x + 0.5_dp * dt * v

    state%z(1:3) = x
    state%z(4:6) = v
    ierr = FO_OK
  end subroutine boris_step_cart

  ! Cylindrical Boris in coordinates (R,phi,Z), contravariant velocity in
  ! z(4:6). Symmetric splitting: half geodesic kick / half drift / orthonormal
  ! magnetic rotation / half drift / half geodesic kick. The geodesic kick
  ! applies -Gamma^i_{mn} v^m v^n from the provider's Christoffel symbols.
  subroutine boris_step_cyl(state, ierr)
    type(FullOrbitState), intent(inout) :: state
    integer, intent(out) :: ierr
    real(dp) :: u(3), v(3), Bvec(3), Bmod, hcov(3)
    real(dp) :: dt, qmc, R
    real(dp) :: vorth(3), t(3), sv(3), vprime(3), tmag2

    dt  = state%dt
    qmc = state%qm / c
    u = state%z(1:3)
    v = state%z(4:6)

    call geodesic_half_kick(state, u, v, 0.5_dp*dt)

    u = u + 0.5_dp * dt * v
    R = u(1)
    if (R <= 0.0_dp) then
      ierr = FO_ERR_OUT_OF_DOMAIN
      return
    end if

    call state%prov%eval_field(u, Bvec, Bmod, hcov, ierr)
    if (ierr /= FO_OK) return

    ! magnetic rotation in the local orthonormal frame
    vorth = contravar_to_orth(v, R)
    t = qmc * Bvec * 0.5_dp * dt
    tmag2 = dot_product(t, t)
    sv = 2.0_dp * t / (1.0_dp + tmag2)
    vprime = vorth + cross(vorth, t)
    vorth = vorth + cross(vprime, sv)
    v = orth_to_contravar(vorth, R)

    u = u + 0.5_dp * dt * v
    if (u(1) <= 0.0_dp) then
      ierr = FO_ERR_OUT_OF_DOMAIN
      return
    end if

    call geodesic_half_kick(state, u, v, 0.5_dp*dt)

    state%z(1:3) = u
    state%z(4:6) = v
    ierr = FO_OK
  end subroutine boris_step_cyl

  subroutine geodesic_half_kick(state, u, v, h)
    type(FullOrbitState), intent(in) :: state
    real(dp), intent(in)    :: u(3), h
    real(dp), intent(inout) :: v(3)
    real(dp) :: Gamma(3,3,3), acc(3)
    integer :: i, m, n

    call state%prov%christoffel(u, Gamma)
    acc = 0.0_dp
    do i = 1, 3
      do m = 1, 3
        do n = 1, 3
          acc(i) = acc(i) - Gamma(i,m,n) * v(m) * v(n)
        end do
      end do
    end do
    v = v + h * acc
  end subroutine geodesic_half_kick

  ! Orthonormal cylindrical components (vR,vphi_phys,vZ) from contravariant
  ! (v^R, v^phi, v^Z): vphi_phys = R v^phi.
  pure function contravar_to_orth(v, R) result(vorth)
    real(dp), intent(in) :: v(3), R
    real(dp) :: vorth(3)
    vorth = [v(1), R * v(2), v(3)]
  end function contravar_to_orth

  pure function orth_to_contravar(vorth, R) result(v)
    real(dp), intent(in) :: vorth(3), R
    real(dp) :: v(3)
    v = [vorth(1), vorth(2) / R, vorth(3)]
  end function orth_to_contravar

  pure function vperp_sq_cyl(R, v, hcov) result(vperp2)
    real(dp), intent(in) :: R, v(3), hcov(3)
    real(dp) :: vperp2, vorth(3), hunit(3), vpar, vsq
    vorth = contravar_to_orth(v, R)
    ! hcov are covariant unit-field comps; orthonormal unit field is hcov/g_ii.
    ! For the toroidal mock h_orth = (0,1,0).
    hunit = [hcov(1), hcov(2) / R, hcov(3)]
    vpar = dot_product(vorth, hunit)
    vsq  = dot_product(vorth, vorth)
    vperp2 = max(vsq - vpar*vpar, 0.0_dp)
  end function vperp_sq_cyl

  pure function cross(a, b) result(c)
    real(dp), intent(in) :: a(3), b(3)
    real(dp) :: c(3)
    c(1) = a(2)*b(3) - a(3)*b(2)
    c(2) = a(3)*b(1) - a(1)*b(3)
    c(3) = a(1)*b(2) - a(2)*b(1)
  end function cross

  ! Gyro-average to guiding-center scalars. For the analytic mocks this returns
  ! position-coordinate scalars and the instantaneous pitch/speed; full
  ! gyro-averaging over coordinate maps lands with the libneo PR.
  subroutine convert_full_to_gc(state, s, theta, phi, lambda, v)
    type(FullOrbitState), intent(in)  :: state
    real(dp), intent(out) :: s, theta, phi, lambda, v
    real(dp) :: Bvec(3), Bmod, hcov(3), vpar, speed
    real(dp) :: vorth(3), hunit(3)
    integer :: ierr

    call state%prov%eval_field(state%z(1:3), Bvec, Bmod, hcov, ierr)
    select case (state%ncoord)
    case (COORD_CYL)
      s = state%z(1); theta = state%z(3); phi = state%z(2)
      vorth = contravar_to_orth(state%z(4:6), state%z(1))
      hunit = [hcov(1), hcov(2)/state%z(1), hcov(3)]
      speed = sqrt(dot_product(vorth, vorth))
      vpar  = dot_product(vorth, hunit)
    case default
      s = state%z(1); theta = state%z(2); phi = state%z(3)
      speed = sqrt(dot_product(state%z(4:6), state%z(4:6)))
      vpar  = dot_product(state%z(4:6), hcov)
    end select
    v = speed
    if (speed > 0.0_dp) then
      lambda = vpar / speed
    else
      lambda = 0.0_dp
    end if
  end subroutine convert_full_to_gc

  ! 0.5*m*|v|^2 (+ q*Phi). Velocity magnitude uses the metric in curvilinear
  ! coordinates so it is the physical speed.
  function compute_energy(state) result(energy)
    type(FullOrbitState), intent(in) :: state
    real(dp) :: energy
    real(dp) :: speed2, Phi, dPhi(3), vorth(3)

    select case (state%ncoord)
    case (COORD_CYL)
      vorth = contravar_to_orth(state%z(4:6), state%z(1))
      speed2 = dot_product(vorth, vorth)
    case default
      speed2 = dot_product(state%z(4:6), state%z(4:6))
    end select
    call state%prov%eval_potential(state%z(1:3), Phi, dPhi)
    energy = 0.5_dp * state%mass * speed2 + state%charge * Phi
  end function compute_energy

end module orbit_full
