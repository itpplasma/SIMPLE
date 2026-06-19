program test_fo_symplectic
  ! Behavioral tests for the implicit-midpoint curvilinear full-orbit pusher
  ! (ORBIT_FOSYMPL, foimpl_step) on the analytic mock providers. Oracles:
  !   1. uniform-B: |v| and energy conserved to solver tolerance; closed circle.
  !   2. cylindrical 1/R curvature drift vs analytic v_d (geodesic + Lorentz).
  !   3. long-run energy: no secular drift (structure-preserving midpoint), the
  !      property an explicit RK does not have.
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use util, only: c, twopi, p_mass, e_charge
  use orbit_full, only: FullOrbitState, init_full_orbit_state, &
      timestep_full_orbit, compute_energy, ORBIT_FOSYMPL, COORD_CART, COORD_CYL, &
      FO_OK
  use orbit_full_mock_cart, only: cartesian_provider_t, FIELD_UNIFORM
  use orbit_full_mock_cyl, only: cylindrical_provider_t
  implicit none

  integer :: nfail
  nfail = 0

  call test_uniform_gyration(nfail)
  call test_cyl_curvature_drift(nfail)
  call test_energy_no_secular_drift(nfail)

  if (nfail == 0) then
    print *, 'ALL FO-SYMPLECTIC TESTS PASSED'
  else
    print *, 'FO-SYMPLECTIC TESTS FAILED: ', nfail
    error stop 1
  end if

contains

  subroutine check(name, ok, nfail)
    character(*), intent(in) :: name
    logical, intent(in) :: ok
    integer, intent(inout) :: nfail
    if (ok) then
      print '(A,A)', 'PASS  ', name
    else
      print '(A,A)', 'FAIL  ', name
      nfail = nfail + 1
    end if
  end subroutine check

  subroutine test_uniform_gyration(nfail)
    integer, intent(inout) :: nfail
    type(cartesian_provider_t), target :: prov
    type(FullOrbitState) :: st
    real(dp) :: mass, charge, B0, vperp, vpar, speed0
    real(dp) :: Omega, period, dt, rL
    real(dp) :: x0(3), v0(3), xstart(3)
    real(dp) :: e0, e1, errpos, errv
    integer :: i, nstep, ierr

    mass   = 4.0_dp * p_mass
    charge = 2.0_dp * e_charge
    B0     = 1.0d4
    prov%field_kind = FIELD_UNIFORM
    prov%B0 = [0.0_dp, 0.0_dp, B0]

    vperp = 1.0d7
    vpar  = 3.0d6
    speed0 = sqrt(vperp**2 + vpar**2)
    Omega = charge * B0 / (mass * c)
    period = twopi / Omega
    rL = mass * c * vperp / (charge * B0)
    nstep = 400
    dt = period / nstep

    x0 = [0.0_dp, 0.0_dp, 0.0_dp]
    v0 = [vperp, 0.0_dp, vpar]
    call init_full_orbit_state(st, x0, v0, ORBIT_FOSYMPL, COORD_CART, &
                               mass, charge, dt, prov)
    xstart = st%z(1:3)
    e0 = compute_energy(st)

    do i = 1, nstep
      call timestep_full_orbit(st, ierr)
      if (ierr /= FO_OK) then
        call check('uniform: timestep ierr', .false., nfail)
        return
      end if
    end do

    e1 = compute_energy(st)
    errv = abs(sqrt(dot_product(st%z(4:6), st%z(4:6))) - speed0) / speed0
    errpos = sqrt((st%z(1)-xstart(1))**2 + (st%z(2)-xstart(2))**2)

    print '(A,ES12.4,A,ES12.4)', '  uniform: |v| relerr=', errv, &
        ' return-pos err=', errpos
    print '(A,ES12.4)', '  uniform: dE/E=', abs(e1-e0)/e0

    call check('uniform: |v| constant', errv < 1d-9, nfail)
    call check('uniform: energy constant', abs(e1-e0)/e0 < 1d-9, nfail)
    call check('uniform: return to start', errpos < 1d-2*rL, nfail)
  end subroutine test_uniform_gyration

  subroutine test_cyl_curvature_drift(nfail)
    integer, intent(inout) :: nfail
    type(cylindrical_provider_t), target :: prov
    type(FullOrbitState) :: st
    real(dp) :: mass, charge, B0, R0, R, Bloc
    real(dp) :: Omega, period, dt, vd_exact, vd_meas
    real(dp) :: u0(3), w0(3), z0, t_total, vpar_in, vperp_in
    integer :: i, nstep_per, nper_run, nstep, ierr

    mass   = 4.0_dp * p_mass
    charge = 2.0_dp * e_charge
    B0     = 1.0d4
    R0     = 200.0_dp
    R      = 200.0_dp
    prov%B0 = B0
    prov%R0 = R0

    vpar_in  = 1.0d7
    vperp_in = 1.0d5
    Bloc = B0 * R0 / R
    Omega = charge * Bloc / (mass * c)
    period = twopi / Omega
    nstep_per = 300
    nper_run = 40
    nstep = nstep_per * nper_run
    dt = period / nstep_per

    vd_exact = (mass * c) / (charge * Bloc * R) * &
               (vpar_in**2 + 0.5_dp * vperp_in**2)

    u0 = [R, 0.0_dp, 0.0_dp]
    w0 = [vperp_in, vpar_in / R, 0.0_dp]
    call init_full_orbit_state(st, u0, w0, ORBIT_FOSYMPL, COORD_CYL, &
                               mass, charge, dt, prov)
    z0 = st%z(3)
    do i = 1, nstep
      call timestep_full_orbit(st, ierr)
      if (ierr /= FO_OK) then
        call check('cyl curvature: timestep ierr', .false., nfail)
        return
      end if
    end do
    t_total = nstep * dt
    vd_meas = (st%z(3) - z0) / t_total

    print '(A,ES12.4,A,ES12.4)', '  cyl curvature: vd_exact=', vd_exact, &
        ' vd_meas=', vd_meas
    call check('cyl curvature: vertical drift', &
        abs(vd_meas - vd_exact) < 0.10_dp*abs(vd_exact), nfail)
  end subroutine test_cyl_curvature_drift

  subroutine test_energy_no_secular_drift(nfail)
    ! Long run in uniform B: the midpoint energy stays bounded with no secular
    ! growth. Compare mean energy over first vs last quarter of the run.
    integer, intent(inout) :: nfail
    type(cartesian_provider_t), target :: prov
    type(FullOrbitState) :: st
    real(dp) :: mass, charge, B0, vperp, vpar
    real(dp) :: Omega, period, dt, e0, e, emin, emax
    real(dp) :: efirst, elast, x0(3), v0(3)
    integer :: i, nstep, nq, ierr

    mass   = 4.0_dp * p_mass
    charge = 2.0_dp * e_charge
    B0     = 1.0d4
    prov%field_kind = FIELD_UNIFORM
    prov%B0 = [0.0_dp, 0.0_dp, B0]

    vperp = 1.0d7
    vpar  = 2.0d6
    Omega = charge * B0 / (mass * c)
    period = twopi / Omega
    nstep = 200 * 200
    dt = period / 200

    x0 = [0.0_dp, 0.0_dp, 0.0_dp]
    v0 = [vperp, 0.0_dp, vpar]
    call init_full_orbit_state(st, x0, v0, ORBIT_FOSYMPL, COORD_CART, &
                               mass, charge, dt, prov)
    e0 = compute_energy(st)
    emin = e0; emax = e0
    nq = nstep/4
    efirst = 0.0_dp; elast = 0.0_dp

    do i = 1, nstep
      call timestep_full_orbit(st, ierr)
      if (ierr /= FO_OK) then
        call check('no-drift: timestep ierr', .false., nfail)
        return
      end if
      e = compute_energy(st)
      emin = min(emin, e); emax = max(emax, e)
      if (i <= nq) efirst = efirst + e
      if (i > nstep - nq) elast = elast + e
    end do
    efirst = efirst / nq
    elast  = elast / nq

    print '(A,ES12.4)', '  no-drift: energy bound (emax-emin)/e0=', (emax-emin)/e0
    print '(A,ES12.4)', '  no-drift: secular |<e>_last-<e>_first|/e0=', &
        abs(elast - efirst)/e0

    call check('fosympl: energy bounded long run', (emax-emin)/e0 < 1d-8, nfail)
    call check('fosympl: no secular energy drift', &
        abs(elast - efirst)/e0 < 1d-10, nfail)
  end subroutine test_energy_no_secular_drift

end program test_fo_symplectic
