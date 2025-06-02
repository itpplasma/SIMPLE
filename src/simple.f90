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
save

public

  type :: Tracer
    double precision :: fper
    double precision :: dtau, dtaumin, v0
    integer          :: n_e, n_d

    integer :: integmode = 0 ! 0 = RK, 1 = Euler1, 2 = Euler2, 3 = Verlet
    double precision :: relerr

    type(FieldCan) :: f
    type(SymplecticIntegrator) :: si
    type(MultistageIntegrator) :: mi
  end type Tracer

  interface tstep
      module procedure timestep
      module procedure timestep_z
      module procedure timestep_sympl_z
   end interface tstep

contains


  subroutine init_vmec(vmec_file, ans_s, ans_tp, amultharm, fper)
    use spline_vmec_sub, only : spline_vmec_data, volume_and_B00
    use vmecin_sub, only : stevvo

    character(*), intent(in) :: vmec_file
    integer, intent(in) :: ans_s, ans_tp, amultharm
    double precision, intent(out) :: fper

    integer             :: L1i
    double precision    :: RT0, R0i, cbfi, bz0i, bf0, volume, B00

    ! TODO: Remove side effects
    netcdffile = vmec_file
    ns_s = ans_s
    ns_tp = ans_tp
    multharm = amultharm

    call spline_vmec_data ! initialize splines for VMEC field
    call stevvo(RT0, R0i, L1i, cbfi, bz0i, bf0) ! initialize periods and major radius
    fper = twopi/dble(L1i)   !<= field period
    print *, 'R0 = ', RT0, ' cm, fper = ', fper
    call volume_and_B00(volume,B00)
    print *,'volume = ',volume,' cm^3,  B_00 = ',B00,' G'
  end subroutine init_vmec


  subroutine init_params(self, Z_charge, m_mass, E_kin, npoints, store_step, relerr)
    ! Initializes normalization for velocity and Larmor radius based on kinetic energy
    ! of plasma particles (= temperature for thermal particles).
    use new_vmec_stuff_mod, only : rmajor

    type(Tracer) :: self
    integer, intent(in), optional :: Z_charge, m_mass
    double precision, intent(in), optional :: E_kin
    integer, intent(in), optional :: npoints  ! Integrator resolution. Number of
    ! integration points for strongly passing particles around the torus

    integer, intent(in), optional :: store_step       ! Store every X timesteps
    double precision, intent(in), optional :: relerr  ! Relative error

    if (present(Z_charge)) self%n_e = Z_charge
    if (present(m_mass)) self%n_d = m_mass
    if (present(relerr)) self%relerr = relerr

    ! Neglect relativistic effects by large inverse relativistic temperature
    rmu=1d8

    ! Reference velocity and normalized Larmor radius
    if (present(E_kin)) then
      self%v0 = sqrt(2.d0*E_kin*ev/(self%n_d*p_mass))
    else
      self%v0 = sqrt(2.d0*3.5d6*ev/(self%n_d*p_mass))
    end if

    ! Larmor radius times magnetic field:
    ro0=self%v0*self%n_d*p_mass*c/(self%n_e*e_charge)

    if (present(npoints)) then
      self%dtaumin=twopi*rmajor*1.0d2/npoints
    else
      self%dtaumin=twopi*rmajor*1.0d2/256.0d0
    end if

   ! timestep where to get results
    if (present(store_step)) then
      self%dtau = store_step*self%dtaumin
    else
      self%dtau = self%dtaumin
    end if

  end subroutine init_params

  subroutine init_sympl(si, f, z0, dtau, dtaumin, rtol_init, mode_init)
    !
    type(SymplecticIntegrator), intent(inout) :: si
    type(FieldCan), intent(inout) :: f
    double precision, intent(in) :: z0(:)
    double precision, intent(in) :: dtau, dtaumin
    double precision, intent(in) :: rtol_init
    integer, intent(in) :: mode_init ! 1 = expl.-impl. Euler, 2 = impl.-expl. Euler, 3 = Midpoint

    double precision :: z(4)

    if(min(dabs(mod(dtau, dtaumin)), dabs(mod(dtau, dtaumin)-dtaumin)) > 1d-9*dtaumin) then
      print *, 'dtau = ', dtau, ' dtaumin = ', dtaumin
      error stop 'orbit_sympl_init - error: dtau/dtaumin not integer'
    endif

    ! Initialize symplectic integrator
    call eval_field(f, z0(1), z0(2), z0(3), 0)

    si%pabs = z0(4)

    f%mu = .5d0*z0(4)**2*(1.d0-z0(5)**2)/f%Bmod*2d0 ! mu by factor 2 from other modules
    f%ro0 = ro0/dsqrt(2d0) ! ro0 = mc/e*v0, different by sqrt(2) from other modules
    f%vpar = z0(4)*z0(5)*dsqrt(2d0) ! vpar_bar = vpar/sqrt(T/m), different by sqrt(2) from other modules

    z(1:3) = z0(1:3)  ! s, th, ph
    z(4) = f%vpar*f%hph + f%Aph/f%ro0 ! pphi

    ! factor 1/sqrt(2) due to velocity normalisation different from other modules
    call orbit_sympl_init(si, f, z, dtaumin/dsqrt(2d0), nint(dtau/dtaumin), &
                          rtol_init, mode_init)
  end subroutine init_sympl

  subroutine init_integrator(self, z0)
    type(Tracer), intent(inout) :: self
    double precision, intent(in) :: z0(:)

    call init_sympl(self%si, self%f, z0, self%dtau, self%dtaumin, &
      self%relerr, self%integmode)
  end subroutine init_integrator

  subroutine timestep(self, s, th, ph, lam, ierr)
    type(Tracer), intent(inout) :: self
    double precision, intent(inout) :: s, th, ph, lam
    integer, intent(out) :: ierr

    double precision, dimension(5) :: z

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
    double precision, intent(inout) :: z(:)
    integer, intent(out) :: ierr

    call orbit_timestep_axis(z, self%dtau, self%dtau, self%relerr, ierr)
  end subroutine timestep_z

  subroutine timestep_sympl_z(si, f, z, ierr)
    type(SymplecticIntegrator), intent(inout) :: si
    type(FieldCan), intent(inout) :: f
    double precision, intent(inout) :: z(:)
    integer, intent(out) :: ierr

    if (z(1) < 0.0 .or. z(1) > 1.0) then
      ierr = 1
      return
    end if

    call orbit_timestep_sympl(si, f, ierr)

    z(1:3) = si%z(1:3)
    z(5) = f%vpar/(si%pabs*dsqrt(2d0))  ! alambda

  end subroutine timestep_sympl_z

end module simple
