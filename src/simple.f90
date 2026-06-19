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
    cpp_canon_step, cpp_canon_to_gc, MODEL_CPP_SYM, COORD_CHARTMAP, COORD_VMEC
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
    ! the SAME (s,theta,phi,v/v0,lambda) GC start as init_sympl.
    !
    ! Coordinate route: REAL VMEC flux coordinates (COORD_VMEC). The diagnosis on
    ! the Cartesian-storage Boozer chartmap (DOC/coordinates-and-fields.md, "6D
    ! canonical CPP") found its libneo periodic-Cartesian spline destroys the
    ! secular toroidal rotation for nfp>1, so the spline metric is inconsistent
    ! with the Boozer covariant field (h_i g^ij h_j ~ nfp^2, not 1). The VMEC
    ! flux metric from libneo is consistent (test_cpp_vmec: |g g^-1 - I| < 1e-10,
    ! h_i g^ij h_j ~ 1 to FD accuracy), so the production loss path runs there. The
    ! 6D state runs natively in u=(s,vartheta,varphi); s is the chart-independent
    ! flux label, so the s>=1 loss test and the s-binned confined fraction carry
    ! over even though Boozer and VMEC angles differ.
    !
    ! Units: the SIMPLE GC normalization (same as init_sympl). With the CONSISTENT
    ! VMEC metric the covariant unit field obeys h_i g^ij h_j = |h|^2 = 1, so the
    ! 6D Hamiltonian H = (1/2m)(p-qcA)g^ij(p-qcA) + mu|B| reduces to the GC
    ! H = vpar_bar^2/2 + mu_bar|B| with mass=1 and the seed p_i = vpar_bar h_i +
    ! A_i/ro0_bar: along the field (p-qcA) = vpar_bar h, so the kinetic term is
    ! (vpar_bar^2/2m)|h|^2 = vpar_bar^2/2. (This identity FAILED on the chartmap,
    ! whose |h|^2 was ~nfp^2.) Keeping mass=1 also keeps the velocities O(vpar_bar)
    ! ~ O(1), so the canonical-midpoint Newton stays well conditioned -- physical
    ! CGS mass ~ 1e-24 would blow up v^i = g^ij(...)/m and wreck the solve.
    ! qc = 1/ro0_bar = sqrt(2)/ro0, dt = dtaumin/sqrt(2): both identical to GC.
    use orbit_cpp_vmec_metric, only: vmec_metric_attach, vmec_metric_ready, &
      vmec_eval_field
    type(cpp_canon_state_t), intent(out) :: cpp
    type(field_can_t), intent(inout) :: f
    real(dp), intent(in) :: z0(:)
    real(dp), intent(in) :: dtaumin

    real(dp) :: ro0_bar, x0(3), Acov(3), Bmod, dBmod(3), hcov(3), mu, vpar_bar

    if (.not. vmec_metric_ready()) call vmec_metric_attach()

    ! 6D state in the VMEC flux chart: u=(s,vartheta,varphi), s direct (no rho).
    x0(1) = min(max(z0(1), 0d0), 1d0)
    x0(2) = z0(2)
    x0(3) = z0(3)

    call vmec_eval_field(x0, Acov, Bmod, dBmod, hcov)

    mu = .5d0*z0(4)**2*(1.d0-z0(5)**2)/Bmod*2d0      ! mu by factor 2 (GC convention)
    ro0_bar = ro0/dsqrt(2d0)                          ! ro0 smaller by sqrt(2)
    vpar_bar = z0(4)*z0(5)*dsqrt(2d0)                 ! vpar_bar = vpar/sqrt(T/m)

    ! mass=1 (see header): the consistent |h|^2=1 metric makes the GC reduction
    ! exact; st%ro0=ro0_bar gives qc=1/ro0_bar so p_i seeds match the GC pphi.
    call cpp_canon_init(cpp, MODEL_CPP_SYM, COORD_VMEC, x0, vpar0=vpar_bar, &
      vperp0=0d0, mu_in=mu, mass=1d0, charge=1d0, dt=dtaumin/dsqrt(2d0), &
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
    ! cpp ierr: 2 = z(1)>=1 (s>=1 loss), 1 = LU fail, 3 = non-converge. All map to
    ! a nonzero orbit error consistent with the sympl loss/abort semantics.
    if (ierr /= 0) return

    ! Write back z. COORD_VMEC runs in s directly; COORD_CHARTMAP in rho (s=rho^2).
    ! z(4)=pabs is the conserved normalized speed; z(5)=lambda (vpar is the
    ! normalized vpar_bar in both wires) so classification/output read z(4:5) the
    ! same as to_standard_z_coordinates.
    call cpp_canon_to_gc(cpp, r, th, ph, vpar)
    z(4) = cpp%pabs
    z(2) = cpp%z(2)
    z(3) = cpp%z(3)
    z(5) = vpar/(z(4)*dsqrt(2d0))
    if (cpp%coord == COORD_CHARTMAP) then
      z(1) = cpp%z(1)**2   ! s = rho^2 (chartmap chart)
    else
      z(1) = cpp%z(1)      ! s direct (VMEC flux chart)
    end if
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
