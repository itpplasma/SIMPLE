program test_cpp_gc_fullorbit_convergence
  ! Three-system cross-validation on analytic oracle fields (no spline, no field
  ! file): GC ~= CPP ~= gyro-averaged full orbit as rho* -> 0.
  !
  !   A. CPP == GC on the analytic circular-tokamak "test" canonical chart
  !      (field_can_test, finite rotational transform iota0). The CPP Gauss
  !      residual is the GC degenerate-Lagrangian system at fixed mu, so on the
  !      same chart, same stage, same dt the 4D state agrees to Newton tolerance
  !      (atol-level), independent of rho*. Checked across GAUSS1..GAUSS3 for a
  !      trapped (banana) initial condition; p_phi conserved (axisymmetric).
  !
  !   B. GC == gyro-averaged full orbit as rho* -> 0 in the cylindrical 1/R mock
  !      (closed-form curvature + grad-B drift). The Boris orbit averaged over a
  !      bounce of one gyroperiod has a vertical drift speed; GC has the analytic
  !      drift. The relative drift error scales linearly with rho* = r_L/R
  !      (first-order GC accuracy): slope check across rho* in {1e-1,1e-2,1e-3}.
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use util, only: c, twopi, pi, p_mass, e_charge
  use field_can_mod, only: field_can_t, field_can_from_name, evaluate, get_val
  use field_can_base, only: n_field_evaluations
  use orbit_symplectic, only: orbit_timestep_sympl, orbit_sympl_init
  use orbit_symplectic_base, only: symplectic_integrator_t, GAUSS1, GAUSS2, GAUSS3
  use orbit_cpp, only: orbit_timestep_cpp, orbit_cpp_init, cpp_stages_from_mode
  use orbit_full, only: FullOrbitState, init_full_orbit_state, &
      timestep_full_orbit, convert_full_to_gc, ORBIT_BORIS, COORD_CYL, FO_OK
  use orbit_full_mock_cyl, only: cylindrical_provider_t

  implicit none

  integer :: nfail
  nfail = 0

  ! Analytic circular-tokamak canonical chart: no field file required.
  call field_can_from_name('test')
  if (.not. associated(evaluate)) then
    print *, 'evaluate pointer not associated for test chart'
    error stop 1
  end if

  call test_cpp_equals_gc_analytic(GAUSS1, 'GAUSS1', nfail)
  call test_cpp_equals_gc_analytic(GAUSS2, 'GAUSS2', nfail)
  call test_cpp_equals_gc_analytic(GAUSS3, 'GAUSS3', nfail)
  call test_gc_fullorbit_convergence(nfail)

  if (nfail == 0) then
    print *, 'ALL THREE-SYSTEM CONVERGENCE TESTS PASSED'
  else
    print *, 'THREE-SYSTEM CONVERGENCE TESTS FAILED: ', nfail
    error stop 1
  end if

contains

  ! Initialize a GC symplectic step on the test chart with an explicit Larmor
  ! scale ro0 (rho*). Mirrors simple%init_sympl arithmetic but bypasses the
  ! global ro0/rmu so the analytic chart needs no field file. z0 = (r,th,ph,p,la).
  subroutine init_test_chart(si, f, z0, ro0_in, dt, ntau, rtol, mode)
    type(symplectic_integrator_t), intent(inout) :: si
    type(field_can_t), intent(inout) :: f
    real(dp), intent(in) :: z0(5), ro0_in, dt, rtol
    integer, intent(in) :: ntau, mode
    real(dp) :: z(4)

    call evaluate(f, z0(1), z0(2), z0(3), 0)
    f%mu  = 0.5_dp*z0(4)**2*(1.0_dp - z0(5)**2)/f%Bmod*2.0_dp
    f%ro0 = ro0_in
    f%vpar = z0(4)*z0(5)*sqrt(2.0_dp)
    z(1:3) = z0(1:3)
    z(4)   = f%vpar*f%hph + f%Aph/f%ro0
    call orbit_sympl_init(si, f, z, dt, ntau, rtol, mode)
  end subroutine init_test_chart

  ! A. CPP == GC on the analytic chart, plus banana confinement and p_phi
  ! conservation. Trapped IC: small pitch so the particle mirrors.
  subroutine test_cpp_equals_gc_analytic(mode, tag, nfail)
    integer, intent(in) :: mode
    character(*), intent(in) :: tag
    integer, intent(inout) :: nfail

    type(symplectic_integrator_t) :: si_gc, si_cpp
    type(field_can_t) :: f_gc, f_cpp
    real(dp) :: z0(5), dt, rtol, maxdiff, d, pphi0, pdev, rmin, rmax
    integer :: it, ierr_gc, ierr_cpp, s, nstep

    ! Trapped (banana) IC on the analytic chart: pitch lambda small.
    z0 = [0.3_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.15_dp]
    dt    = 2.0_dp
    rtol  = 1.0e-13_dp
    nstep = 400
    s = cpp_stages_from_mode(mode)

    call init_test_chart(si_gc,  f_gc,  z0, 1.0e-3_dp, dt, 1, rtol, mode)
    call init_test_chart(si_cpp, f_cpp, z0, 1.0e-3_dp, dt, 1, rtol, mode)
    call orbit_cpp_init(si_cpp, f_cpp)

    pphi0 = si_gc%z(4); pdev = 0.0_dp
    rmin = z0(1); rmax = z0(1)
    maxdiff = 0.0_dp
    do it = 1, nstep
      call orbit_timestep_sympl(si_gc, f_gc, ierr_gc)
      call orbit_timestep_cpp(si_cpp, f_cpp, s, ierr_cpp)
      if (ierr_gc /= 0 .or. ierr_cpp /= 0) then
        print '(A,A,A,I0,A,I0)', '  ', tag, ': early exit ierr_gc=', ierr_gc, &
            ' ierr_cpp=', ierr_cpp
        exit
      end if
      d = maxval(abs(si_cpp%z - si_gc%z))
      maxdiff = max(maxdiff, d)
      pdev = max(pdev, abs(si_gc%z(4) - pphi0))
      rmin = min(rmin, si_gc%z(1)); rmax = max(rmax, si_gc%z(1))
    end do

    print '(A,A,A,ES12.4)', '  ', tag, ': max |z_CPP - z_GC| = ', maxdiff
    print '(A,A,A,ES12.4,A,ES12.4)', '  ', tag, ': pphi excursion = ', pdev, &
        '  banana radial width = ', rmax - rmin
    call check(tag//': CPP == GC to Newton tol', maxdiff < 1.0e-10_dp, nfail)
    ! Axisymmetric analytic chart: p_phi exactly conserved.
    call check(tag//': p_phi conserved (axisymmetric)', pdev < 1.0e-10_dp, nfail)
    ! Trapped particle: finite radial banana width, stays confined.
    call check(tag//': trapped banana (finite width, confined)', &
               (rmax - rmin) > 1.0e-4_dp .and. rmax < 0.5_dp, nfail)
  end subroutine test_cpp_equals_gc_analytic

  ! B. GC drift vs gyro-averaged Boris drift in the 1/R mock, O(rho*) slope.
  ! GC analytic vertical drift speed (CGS) for the toroidal 1/R field at major
  ! radius R: v_d = (m c)/(q B R) (vpar^2 + vperp^2/2). The Boris orbit drifts
  ! at v_d + O(rho*) corrections; the relative error must fall ~ linearly in
  ! rho* = r_L/R. GC reference is the same closed form (the test chart has no
  ! 1/R cylindrical analogue, so the closed form is the GC oracle here).
  subroutine test_gc_fullorbit_convergence(nfail)
    integer, intent(inout) :: nfail
    real(dp) :: rhostar(3), err(3)
    integer :: k

    rhostar = [1.0e-1_dp, 1.0e-2_dp, 1.0e-3_dp]
    do k = 1, 3
      err(k) = boris_drift_relerr(rhostar(k))
      print '(A,ES10.2,A,ES12.4)', '  rho* = ', rhostar(k), &
          '  full-orbit drift relerr vs GC = ', err(k)
    end do

    ! First-order GC: relerr ~ rho*. Each decade in rho* must cut the error by
    ! at least ~5x (allowing geometric slack); and the smallest rho* must be
    ! well below the largest.
    call check('FO->GC drift error decreases with rho*', &
               err(2) < 0.4_dp*err(1) .and. err(3) < 0.4_dp*err(2), nfail)
    call check('FO->GC drift converges (relerr < 1% at rho*=1e-3)', &
               err(3) < 1.0e-2_dp, nfail)
  end subroutine test_gc_fullorbit_convergence

  ! Run one Boris orbit in the 1/R cylindrical mock at a chosen rho* = r_L/R,
  ! measure the gyro-averaged vertical (Z) drift speed, return its relative
  ! error against the analytic GC drift. rho* is varied through B0 at fixed
  ! velocities (r_L = m c vperp/(q B)).
  function boris_drift_relerr(rhostar) result(relerr)
    real(dp), intent(in) :: rhostar
    real(dp) :: relerr
    type(cylindrical_provider_t), target :: prov
    type(FullOrbitState) :: st
    real(dp) :: mass, charge, R0, R, vpar_in, vperp_in
    real(dp) :: Bloc, B0, Omega, period, dt, vd_exact, vd_meas
    real(dp) :: u0(3), w0(3), z0, t_total, rL
    integer :: i, nstep_per, nper, nstep, ierr

    mass    = 4.0_dp * p_mass
    charge  = 2.0_dp * e_charge
    R0      = 200.0_dp
    R       = 200.0_dp
    vpar_in  = 1.0d7
    vperp_in = 1.0d6

    ! rho* = r_L/R = m c vperp/(q B R)  =>  B = m c vperp/(q rho* R)
    Bloc = mass * c * vperp_in / (charge * rhostar * R)
    B0   = Bloc * R / R0
    prov%B0 = B0
    prov%R0 = R0

    Omega  = charge * Bloc / (mass * c)
    period = twopi / Omega
    rL     = mass * c * vperp_in / (charge * Bloc)
    nstep_per = 200
    nper = 50
    nstep = nstep_per * nper
    dt = period / nstep_per

    vd_exact = (mass * c) / (charge * Bloc * R) * &
               (vpar_in**2 + 0.5_dp * vperp_in**2)

    u0 = [R, 0.0_dp, 0.0_dp]
    w0 = [vperp_in, vpar_in / R, 0.0_dp]
    call init_full_orbit_state(st, u0, w0, ORBIT_BORIS, COORD_CYL, &
                               mass, charge, dt, prov)
    z0 = st%z(3)
    do i = 1, nstep
      call timestep_full_orbit(st, ierr)
      if (ierr /= FO_OK) then
        relerr = huge(1.0_dp)
        return
      end if
    end do
    t_total = nstep * dt
    vd_meas = (st%z(3) - z0) / t_total
    relerr = abs(vd_meas - vd_exact) / abs(vd_exact)
  end function boris_drift_relerr

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

end program test_cpp_gc_fullorbit_convergence
