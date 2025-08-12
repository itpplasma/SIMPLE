module simple
  use util, only: c, e_charge, p_mass, ev, twopi
  use new_vmec_stuff_mod, only : netcdffile, multharm, ns_s, ns_tp, &
                                 vmec_B_scale, vmec_RZ_scale

  use parmot_mod, only : rmu, ro0
  use velo_mod,   only : isw_field_type
  use field_can_mod, only : FieldCan
  use orbit_symplectic, only : SymplecticIntegrator, MultistageIntegrator, &
    orbit_sympl_init, orbit_timestep_sympl
  use field, only : VmecField
  use field_can_mod, only : eval_field => evaluate, init_field_can, FieldCan
  use diag_mod, only : icounter
  use chamb_sub, only : chamb_can

  implicit none

  ! Define real(dp) kind parameter
  integer, parameter :: dp = kind(1.0d0)
save

public

  type :: Tracer
    ! Physical parameters
    real(dp) :: fper              ! Field period
    real(dp) :: v0                ! Reference velocity
    integer  :: n_e, n_d          ! Charge and mass numbers
    
    ! Integration parameters
    real(dp) :: dtau, dtaumin     ! Time steps
    real(dp) :: relerr            ! Relative error tolerance
    integer  :: integmode = 0     ! 0 = RK, 1 = Euler1, 2 = Euler2, 3 = Verlet
    
    ! Field and integrators
    type(FieldCan) :: f
    type(SymplecticIntegrator) :: si
    type(MultistageIntegrator) :: mi
  contains
    ! Improved type-bound procedures for better encapsulation
    procedure :: init_parameters => tracer_init_parameters
    procedure :: validate_timestep => tracer_validate_timestep
    procedure :: normalize_velocity => tracer_normalize_velocity
    procedure :: compute_larmor_radius => tracer_compute_larmor_radius
  end type Tracer

  interface tstep
      module procedure timestep
      module procedure timestep_z
      module procedure timestep_sympl_z
   end interface tstep

contains


  subroutine init_vmec(vmec_file, ans_s, ans_tp, amultharm, fper)
    use spline_vmec_sub, only : spline_vmec_data

    character(*), intent(in) :: vmec_file
    integer, intent(in) :: ans_s, ans_tp, amultharm
    real(dp), intent(out) :: fper

    integer :: L1i
    real(dp) :: RT0, R0i, cbfi, bz0i, bf0

    ! Initialize VMEC parameters
    call init_vmec_parameters(vmec_file, ans_s, ans_tp, amultharm)

    call spline_vmec_data ! initialize splines for VMEC field
    
    ! Compute field geometry
    call compute_field_geometry(RT0, R0i, L1i, cbfi, bz0i, bf0, fper)
  end subroutine init_vmec


  subroutine init_params(self, Z_charge, m_mass, E_kin, npoints, store_step, relerr)
    ! Initializes normalization for velocity and Larmor radius based on kinetic energy
    ! of plasma particles (= temperature for thermal particles).

    type(Tracer) :: self
    integer, intent(in), optional :: Z_charge, m_mass
    real(dp), intent(in), optional :: E_kin
    integer, intent(in), optional :: npoints  ! Integrator resolution. Number of
    ! integration points for strongly passing particles around the torus

    integer, intent(in), optional :: store_step       ! Store every X timesteps
    real(dp), intent(in), optional :: relerr  ! Relative error

    ! Use the new type-bound procedure for better encapsulation
    call self%init_parameters(Z_charge, m_mass, E_kin, npoints, store_step, relerr)
  end subroutine init_params

  subroutine init_sympl(si, f, z0, dtau, dtaumin, rtol_init, mode_init)
    !
    type(SymplecticIntegrator), intent(inout) :: si
    type(FieldCan), intent(inout) :: f
    real(dp), intent(in) :: z0(:)
    real(dp), intent(in) :: dtau, dtaumin
    real(dp), intent(in) :: rtol_init
    integer, intent(in) :: mode_init ! 1 = expl.-impl. Euler, 2 = impl.-expl. Euler, 3 = Midpoint

    real(dp) :: z(4)

    ! Validate timestep parameters
    if(min(dabs(mod(dtau, dtaumin)), dabs(mod(dtau, dtaumin)-dtaumin)) > 1d-9*dtaumin) then
      print *, 'dtau = ', dtau, ' dtaumin = ', dtaumin
      error stop 'orbit_sympl_init - error: dtau/dtaumin not integer'
    endif

    ! Initialize symplectic field configuration
    call init_symplectic_field(f, z0, si)

    z(1:3) = z0(1:3)  ! s, th, ph
    z(4) = f%vpar*f%hph + f%Aph/f%ro0 ! pphi

    ! factor 1/sqrt(2) due to velocity normalisation different from other modules
    call orbit_sympl_init(si, f, z, dtaumin/dsqrt(2d0), nint(dtau/dtaumin), &
                          rtol_init, mode_init)
  end subroutine init_sympl

  subroutine timestep(self, s, th, ph, lam, ierr)
    type(Tracer), intent(inout) :: self
    real(dp), intent(inout) :: s, th, ph, lam
    integer, intent(out) :: ierr

    real(dp), dimension(5) :: z

    z(1) = s
    z(2) = th
    z(3) = ph
    z(4) = 1d0
    z(5) = lam

    call timestep_z(self, z, ierr)

    s = z(1)
    th = z(2)
    ph = z(3)
    lam = z(5)
  end subroutine timestep

  subroutine timestep_z(self, z, ierr)
    use alpha_lifetime_sub, only : orbit_timestep_axis

    type(Tracer), intent(inout) :: self
    real(dp), intent(inout) :: z(:)
    integer, intent(out) :: ierr

    call orbit_timestep_axis(z, self%dtau, self%dtau, self%relerr, ierr)
  end subroutine timestep_z

  subroutine timestep_sympl_z(si, f, z, ierr)
    type(SymplecticIntegrator), intent(inout) :: si
    type(FieldCan), intent(inout) :: f
    real(dp), intent(inout) :: z(:)
    integer, intent(out) :: ierr

    if (z(1) < 0.0 .or. z(1) > 1.0) then
      ierr = 1
      return
    end if

    call orbit_timestep_sympl(si, f, ierr)

    z(1:3) = si%z(1:3)
    z(5) = f%vpar/(si%pabs*dsqrt(2d0))  ! alambda

  end subroutine timestep_sympl_z

  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! Helper: Initialize Tracer parameters with validation
  subroutine tracer_init_parameters(self, Z_charge, m_mass, E_kin, npoints, store_step, relerr)
    use new_vmec_stuff_mod, only : rmajor
    
    class(Tracer), intent(inout) :: self
    integer, intent(in), optional :: Z_charge, m_mass, npoints, store_step
    real(dp), intent(in), optional :: E_kin, relerr
    
    ! Set parameters with defaults
    if (present(Z_charge)) self%n_e = Z_charge
    if (present(m_mass)) self%n_d = m_mass
    if (present(relerr)) self%relerr = relerr
    
    ! Compute velocity and time stepping parameters
    call self%normalize_velocity(E_kin)
    call self%compute_larmor_radius()
    call compute_timestep_parameters(self, npoints, store_step)
  end subroutine tracer_init_parameters
  
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! Helper: Validate timestep parameters
  logical function tracer_validate_timestep(self, dtau, dtaumin) result(is_valid)
    class(Tracer), intent(in) :: self
    real(dp), intent(in) :: dtau, dtaumin
    
    real(dp) :: tolerance, r
    
    tolerance = 1.0d-9
    if (dtaumin <= 0.0d0) then
      is_valid = .false.
      return
    end if
    
    r = mod(dtau, dtaumin)
    if (r < 0.0d0) r = r + dtaumin
    is_valid = min(dabs(r), dabs(r - dtaumin)) <= tolerance*dtaumin
  end function tracer_validate_timestep
  
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! Helper: Normalize velocity based on kinetic energy
  subroutine tracer_normalize_velocity(self, E_kin)
    class(Tracer), intent(inout) :: self
    real(dp), intent(in), optional :: E_kin
    
    real(dp) :: energy
    
    ! Set energy with default for thermal particles
    energy = 3.5d6*ev  ! Default thermal energy
    if (present(E_kin)) energy = E_kin
    
    ! Neglect relativistic effects by large inverse relativistic temperature
    rmu = 1d8
    
    ! Reference velocity
    self%v0 = sqrt(2.d0*energy/(self%n_d*p_mass))
  end subroutine tracer_normalize_velocity
  
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! Helper: Compute Larmor radius
  subroutine tracer_compute_larmor_radius(self)
    class(Tracer), intent(in) :: self
    
    ! Larmor radius times magnetic field:
    ro0 = self%v0*self%n_d*p_mass*c/(self%n_e*e_charge)
  end subroutine tracer_compute_larmor_radius
  
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! Helper: Compute timestep parameters for integration
  subroutine compute_timestep_parameters(self, npoints, store_step)
    use new_vmec_stuff_mod, only : rmajor
    
    type(Tracer), intent(inout) :: self
    integer, intent(in), optional :: npoints, store_step
    
    integer :: n_integration_points = 256
    
    ! Set integration resolution
    if (present(npoints)) n_integration_points = npoints
    
    self%dtaumin = twopi*rmajor*1.0d2/n_integration_points
    
    ! Set output timestep
    if (present(store_step)) then
      self%dtau = store_step*self%dtaumin
    else
      self%dtau = self%dtaumin
    end if
  end subroutine compute_timestep_parameters
  
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! Helper: Initialize VMEC field parameters
  subroutine init_vmec_parameters(vmec_file, ans_s, ans_tp, amultharm)
    character(*), intent(in) :: vmec_file
    integer, intent(in) :: ans_s, ans_tp, amultharm
    
    ! TODO: Remove side effects - these should be passed as parameters
    netcdffile = vmec_file
    ns_s = ans_s
    ns_tp = ans_tp
    multharm = amultharm
  end subroutine init_vmec_parameters
  
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! Helper: Compute field geometry parameters
  subroutine compute_field_geometry(RT0, R0i, L1i, cbfi, bz0i, bf0, fper)
    use spline_vmec_sub, only : volume_and_B00
    use vmecin_sub, only : stevvo
    
    real(dp), intent(out) :: RT0, R0i, cbfi, bz0i, bf0, fper
    integer, intent(out) :: L1i
    real(dp) :: volume, B00
    
    call stevvo(RT0, R0i, L1i, cbfi, bz0i, bf0) ! initialize periods and major radius
    fper = twopi/dble(L1i)   !<= field period
    print *, 'R0 = ', RT0, ' cm, fper = ', fper
    
    call volume_and_B00(volume, B00)
    print *, 'volume = ', volume, ' cm^3,  B_00 = ', B00, ' G'
  end subroutine compute_field_geometry
  
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! Helper: Initialize symplectic field configuration
  subroutine init_symplectic_field(f, z0, si)
    type(FieldCan), intent(inout) :: f
    real(dp), intent(in) :: z0(:)
    type(SymplecticIntegrator), intent(inout) :: si
    
    ! Initialize symplectic integrator
    call eval_field(f, z0(1), z0(2), z0(3), 0)
    
    si%pabs = z0(4)
    
    f%mu = .5d0*z0(4)**2*(1.d0-z0(5)**2)/f%Bmod*2d0 ! mu by factor 2 from other modules
    f%ro0 = ro0/dsqrt(2d0) ! ro0 = mc/e*v0, different by sqrt(2) from other modules
    f%vpar = z0(4)*z0(5)*dsqrt(2d0) ! vpar_bar = vpar/sqrt(T/m), different by sqrt(2) from other modules
  end subroutine init_symplectic_field

end module simple
