module simple
  use util, only: c, e_charge, p_mass, ev, twopi
  use new_vmec_stuff_mod, only : netcdffile, multharm, ns_s, ns_tp, &
                                 vmec_B_scale, vmec_RZ_scale

  use parmot_mod, only : rmu, ro0
  use velo_mod,   only : isw_field_type
  use field_can_mod, only : field_can_t
  use orbit_symplectic, only : symplectic_integrator_t, multistage_integrator_t, &
    orbit_sympl_init, orbit_timestep_sympl
  use field, only : vmec_field_t
  use field_can_mod, only : eval_field => evaluate, init_field_can, field_can_t
  use orbit_cpp_canonical, only : cpp_canon_state_t, cpp_canon_init, &
    cpp_canon_step, cpp_canon_to_gc, MODEL_CPP_SYM, COORD_CHARTMAP
  use diag_mod, only : icounter
  use chamb_sub, only : chamb_can

  implicit none

  ! Define real(dp) kind parameter
  integer, parameter :: dp = kind(1.0d0)
save

public

  type :: tracer_t
    real(dp) :: fper
    real(dp) :: dtau, dtaumin, v0
    integer          :: n_e, n_d

    integer :: integmode = 0 ! 0 = RK, 1 = Euler1, 2 = Euler2, 3 = Midpoint, 4-7 = Gauss1-4, 15 = Lobatto3
    real(dp) :: relerr

    type(field_can_t) :: f
    type(symplectic_integrator_t) :: si
    type(multistage_integrator_t) :: mi
    type(cpp_canon_state_t) :: cpp  ! genuine 6D CPP state (orbit_model=ORBIT_CPP6D)
  end type tracer_t

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
    real(dp), intent(out) :: fper

    integer             :: L1i
    real(dp)    :: RT0, R0i, cbfi, bz0i, bf0, volume, B00

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

    type(tracer_t) :: self
    integer, intent(in), optional :: Z_charge, m_mass
    real(dp), intent(in), optional :: E_kin
    integer, intent(in), optional :: npoints  ! Integrator resolution. Number of
    ! integration points for strongly passing particles around the torus

    integer, intent(in), optional :: store_step       ! Store every X timesteps
    real(dp), intent(in), optional :: relerr  ! Relative error

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
    type(symplectic_integrator_t), intent(inout) :: si
    type(field_can_t), intent(inout) :: f
    real(dp), intent(in) :: z0(:)
    real(dp), intent(in) :: dtau, dtaumin
    real(dp), intent(in) :: rtol_init
    integer, intent(in) :: mode_init ! 1 = expl.-impl. Euler, 2 = impl.-expl. Euler, 3 = Midpoint

    real(dp) :: z(4)

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

  subroutine init_cpp(cpp, f, z0, dtaumin)
    ! Initialize the genuine 6D canonical CPP state (orbit_model=ORBIT_CPP6D) from
    ! the SAME (s,theta,phi,vpar,mu) GC start as init_sympl, in NORMALIZED TIME.
    ! Replicates the GC sqrt(2) convention verbatim (init_sympl lines above), then
    ! maps onto the dimensionless 6D Hamiltonian on the production Boozer/chartmap
    ! chart: the 6D state runs in u=(rho,theta_B,phi_B) with rho=sqrt(s) so the
    ! libneo chartmap metric is native; field_can supplies A_i,|B|,h_i in s=rho^2.
    ! The magnetic coupling qc=1/ro0_bar=sqrt(2)/ro0 is threaded via st%ro0=ro0_bar,
    ! so the canonical momentum p_i=vpar*h_i+A_i/ro0_bar matches the GC pphi seed.
    type(cpp_canon_state_t), intent(out) :: cpp
    type(field_can_t), intent(inout) :: f
    real(dp), intent(in) :: z0(:)
    real(dp), intent(in) :: dtaumin

    real(dp) :: ro0_bar, x0(3)

    call eval_field(f, z0(1), z0(2), z0(3), 0)

    f%mu = .5d0*z0(4)**2*(1.d0-z0(5)**2)/f%Bmod*2d0 ! mu by factor 2 (GC convention)
    ro0_bar = ro0/dsqrt(2d0)                          ! ro0 smaller by sqrt(2)
    f%vpar = z0(4)*z0(5)*dsqrt(2d0)                   ! vpar_bar = vpar/sqrt(T/m)

    ! 6D state in the metric chart: u=(rho,theta_B,phi_B), rho=sqrt(s).
    x0(1) = dsqrt(max(z0(1), 0d0))
    x0(2) = z0(2)
    x0(3) = z0(3)

    ! mass=charge=1 (thesis e=m=1); dt=dtaumin/sqrt(2) (SAME as GC).
    ! st%ro0=ro0_bar gives qc=1/ro0_bar so p_i seeds match the GC pphi convention;
    ! p_s carries only the O(rho*) g_si v^i metric term (the genuine 6D start).
    call cpp_canon_init(cpp, MODEL_CPP_SYM, COORD_CHARTMAP, x0, vpar0=f%vpar, &
      vperp0=0d0, mu_in=f%mu, mass=1d0, charge=1d0, dt=dtaumin/dsqrt(2d0), &
      ro0_in=ro0_bar)
    cpp%pabs = z0(4)   ! normalized speed; z(4) on write-back, conserved
  end subroutine init_cpp

  subroutine orbit_timestep_cpp_canonical(cpp, f, z, ierr)
    ! Advance the genuine 6D CPP one normalized step (dtaumin/sqrt(2)) and write
    ! back the standard SIMPLE z(1:5) so times_lost/confined_fraction/output read
    ! it identically to the GC path. The wrapper does NOT call
    ! to_standard_z_coordinates (that reads the sympl path); it builds z itself.
    type(cpp_canon_state_t), intent(inout) :: cpp
    type(field_can_t), intent(inout) :: f
    real(dp), intent(inout) :: z(:)
    integer, intent(out) :: ierr

    real(dp) :: r, th, ph, vpar

    if (z(1) < 0.0d0 .or. z(1) > 1.0d0) then
      ierr = 1
      return
    end if

    call cpp_canon_step(cpp, ierr)
    ! cpp ierr: 2 = rho>=1 (s>=1 loss), 1 = LU fail, 3 = non-converge. All map to
    ! a nonzero orbit error consistent with the sympl loss/abort semantics.
    if (ierr /= 0) return

    ! Write back z. State runs in rho; output uses s=rho^2 (loss test, classifier).
    call cpp_canon_to_gc(cpp, r, th, ph, vpar)
    z(1) = cpp%z(1)**2   ! s = rho^2
    z(2) = cpp%z(2)
    z(3) = cpp%z(3)
    ! z(4)=pabs is the normalized speed (conserved); z(5)=vpar/(pabs*sqrt2) matches
    ! to_standard_z_coordinates so classification/output read z(4:5) unchanged.
    z(4) = cpp%pabs
    z(5) = vpar/(z(4)*dsqrt(2d0))
  end subroutine orbit_timestep_cpp_canonical

  subroutine timestep(self, s, th, ph, lam, ierr)
    type(tracer_t), intent(inout) :: self
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

    type(tracer_t), intent(inout) :: self
    real(dp), intent(inout) :: z(:)
    integer, intent(out) :: ierr

    call orbit_timestep_axis(z, self%dtau, self%dtau, self%relerr, ierr)
  end subroutine timestep_z

  subroutine timestep_sympl_z(si, f, z, ierr)
    type(symplectic_integrator_t), intent(inout) :: si
    type(field_can_t), intent(inout) :: f
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

end module simple
