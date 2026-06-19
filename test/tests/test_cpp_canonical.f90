program test_cpp_canonical
  ! Behavioral validation of the 6D canonical-midpoint port against the python
  ! reference oracle (DVI_python). The oracle was regenerated with the CORRECTED
  ! field and metric (numpy 2.4.6 / scipy 1.17.1, scipy.optimize.root hybr,
  ! tol=1e-12): the metric d_33 theta-derivative carries the factor r, and the
  ! field d|B| is the exact closed-form gradient of |B|=sqrt(W) (the previous
  ! python listing dropped the A_th,rth chain-rule term in d|B|/dtheta). The
  ! Fortran uses the same corrected field/metric and reproduces z(t):
  !   CP   dt=1   : per-step z to 1e-10 (CP has no mu|B| term, unaffected by d|B|)
  !   CPP-sym dt=80: per-step z to 1e-10; small, dt^2-shrinking energy oscillation
  !   CPP-var dt=800, ph0=1.0: per-step z to 1e-7 over 2000 steps
  ! With the corrected field the CPP-sym energy oscillation converges as dt^2
  ! (dt=80,40,20,10 -> ~4.2e-5, 1.0e-5, 2.2e-6, 4.8e-7, each halving ~/4), the
  ! true symplectic signature. The earlier buggy d|B| produced a flat ~1e-3
  ! plateau that did NOT converge; the dt^2 test below asserts the correct one.
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use orbit_cpp_canonical, only: cpp_canon_state_t, cpp_canon_init, cpp_canon_step, &
       cpp_canon_energy, residual, jacobian, &
       MODEL_CP, MODEL_CPP_SYM, MODEL_CPP_VAR, COORD_TOK
  implicit none

  integer :: nfail
  real(dp), parameter :: mu = 1.0e-5_dp, mass = 1.0_dp, charge = 1.0_dp
  real(dp), parameter :: x0(3) = [0.1_dp, 1.5_dp, 0.0_dp]

  nfail = 0

  call test_cp(nfail)
  call test_cpp_sym(nfail)
  call test_cpp_var(nfail)
  call test_cpp_sym_convergence(nfail)
  call test_cpp_banana(nfail)
  call test_jacobian_fd(nfail)

  if (nfail == 0) then
    print *, 'ALL CANONICAL 6D PORT TESTS PASSED'
  else
    print *, 'CANONICAL 6D PORT TESTS FAILED: ', nfail
    error stop 1
  end if

contains

  subroutine test_cp(nfail)
    integer, intent(inout) :: nfail
    type(cpp_canon_state_t) :: st
    real(dp) :: vperp0, B0ic, Emax, E0, dE
    integer :: it, ierr
    ! vperp0 chosen so mu = m vperp^2/(2|B|) = 1e-5 at the IC (|B|=0.99749...).
    B0ic = 0.99749164651988482_dp
    vperp0 = sqrt(2.0_dp*mu*B0ic/mass)

    call cpp_canon_init(st, MODEL_CP, COORD_TOK, x0, 0.0_dp, vperp0, mu, &
                        mass, charge, 1.0_dp)
    E0 = cpp_canon_energy(st); Emax = 0.0_dp

    call step_check(st, 1, [1.035959960987658e-01_dp, 1.482618475162977e+00_dp, &
         1.712482822602558e-04_dp], 'CP step1', nfail, ierr)
    call step_check(st, 1, [1.044428104384642e-01_dp, 1.445034902131783e+00_dp, &
         5.556086308662086e-04_dp], 'CP step2', nfail, ierr)
    call step_check(st, 1, [1.019119149294356e-01_dp, 1.414921956768324e+00_dp, &
         8.562929696301598e-04_dp], 'CP step3', nfail, ierr)
    call step_check(st, 2, [9.559580783940201e-02_dp, 1.448795795139552e+00_dp, &
         5.490863655870221e-04_dp], 'CP step5', nfail, ierr)

    ! Energy bound over 1000 steps (re-init for a clean run).
    call cpp_canon_init(st, MODEL_CP, COORD_TOK, x0, 0.0_dp, vperp0, mu, &
                        mass, charge, 1.0_dp)
    E0 = cpp_canon_energy(st)
    do it = 1, 1000
      call cpp_canon_step(st, ierr)
      if (ierr /= 0) then
        call check('CP 1000-step ierr==0', .false., nfail); return
      end if
      dE = abs((cpp_canon_energy(st) - E0)/E0)
      Emax = max(Emax, dE)
    end do
    print '(A,ES12.4)', '  CP max|dE/E0| (1000 steps) = ', Emax
    call check('CP energy bounded (<5e-2)', Emax < 5.0e-2_dp, nfail)
  end subroutine test_cp

  subroutine test_cpp_sym(nfail)
    integer, intent(inout) :: nfail
    type(cpp_canon_state_t) :: st
    real(dp) :: Emax, Emin, E0, E, drift, Eend
    integer :: it, ierr

    call cpp_canon_init(st, MODEL_CPP_SYM, COORD_TOK, x0, 0.0_dp, 0.0_dp, mu, &
                        mass, charge, 80.0_dp)
    call step_check(st, 1, [9.920194206034304e-02_dp, 1.496968003781973e+00_dp, &
         -3.015381629691890e-03_dp], 'CPPsym step1', nfail, ierr)
    call step_check(st, 1, [9.839799891542926e-02_dp, 1.488463997662909e+00_dp, &
         -1.204121153288971e-02_dp], 'CPPsym step2', nfail, ierr)
    call step_check(st, 3, [9.598548498685663e-02_dp, 1.426793705958419e+00_dp, &
         -7.446170518132635e-02_dp], 'CPPsym step5', nfail, ierr)
    call step_check(st, 5, [9.190808176109883e-02_dp, 1.206541898048687e+00_dp, &
         -2.899409392554551e-01_dp], 'CPPsym step10', nfail, ierr)

    ! Symplectic energy bound and end-drift over 1000 steps.
    call cpp_canon_init(st, MODEL_CPP_SYM, COORD_TOK, x0, 0.0_dp, 0.0_dp, mu, &
                        mass, charge, 80.0_dp)
    E0 = cpp_canon_energy(st); Emin = E0; Emax = E0; Eend = E0
    do it = 1, 1000
      call cpp_canon_step(st, ierr)
      if (ierr /= 0) then
        call check('CPPsym 1000-step ierr==0', .false., nfail); return
      end if
      E = cpp_canon_energy(st)
      Emin = min(Emin, E); Emax = max(Emax, E); Eend = E
    end do
    drift = (Eend - E0)/E0
    print '(A,ES12.4,A,ES12.4)', '  CPPsym max|dE/E0| = ', (Emax - Emin)/abs(E0), &
         '  end-drift = ', drift
    ! Corrected-field oracle (dt=80, 1000 steps): max ~4.17e-5, end-drift ~7.3e-7.
    ! No secular growth; far below the buggy-field ~1e-3 plateau.
    call check('CPPsym energy bound ~4e-5', (Emax - Emin)/abs(E0) < 1.0e-4_dp, nfail)
    call check('CPPsym end-drift tiny (<1e-5)', abs(drift) < 1.0e-5_dp, nfail)
  end subroutine test_cpp_sym

  subroutine test_cpp_var(nfail)
    integer, intent(inout) :: nfail
    type(cpp_canon_state_t) :: st
    real(dp) :: xv0(3)
    integer :: it, ierr
    xv0 = [0.1_dp, 1.5_dp, 1.0_dp]   ! cpp_var overrides ph0 = 1.0

    call cpp_canon_init(st, MODEL_CPP_VAR, COORD_TOK, xv0, 0.0_dp, 0.0_dp, mu, &
                        mass, charge, 800.0_dp)
    call step_check_tol(st, 1, [8.417248911051033e-02_dp, 9.254752624927132e-01_dp, &
         4.504656792084775e-01_dp], 'CPPvar step1', nfail, ierr, 1.0e-9_dp)
    call step_check_tol(st, 1, [8.807902711985266e-02_dp, -7.485694411011361e-02_dp, &
         -4.067197600076082e-01_dp], 'CPPvar step2', nfail, ierr, 1.0e-9_dp)
    ! Long-run reference at step 2000.
    do it = 3, 2000
      call cpp_canon_step(st, ierr)
      if (ierr /= 0) then
        call check('CPPvar 2000-step ierr==0', .false., nfail); return
      end if
    end do
    call assert_vec(st%z(1:3), [8.462270505499311e-02_dp, -2.362182109854340e-01_dp, &
         8.265908407137218e+01_dp], 'CPPvar step2000', nfail, 1.0e-7_dp)
  end subroutine test_cpp_var

  subroutine test_cpp_sym_convergence(nfail)
    ! A symplectic midpoint integrator has a bounded energy oscillation that
    ! shrinks as dt^2 with no secular drift. With the corrected field the
    ! oscillation over a fixed time window (80000) converges:
    !   dt=80 -> ~4.2e-5, dt=40 -> ~1.0e-5, dt=20 -> ~2.2e-6, dt=10 -> ~4.8e-7,
    ! each halving of dt reducing the bound by ~4. The earlier buggy d|B| gave a
    ! flat ~1e-3 plateau that did NOT converge; assert the correct dt^2 instead.
    integer, intent(inout) :: nfail
    real(dp) :: dts(4), osc(4), ratio
    integer :: i, n, it, ierr
    type(cpp_canon_state_t) :: st
    real(dp) :: E0, E, Emin, Emax
    logical :: dt2_ok
    dts = [80.0_dp, 40.0_dp, 20.0_dp, 10.0_dp]

    do i = 1, 4
      call cpp_canon_init(st, MODEL_CPP_SYM, COORD_TOK, x0, 0.0_dp, 0.0_dp, mu, &
                          mass, charge, dts(i))
      E0 = cpp_canon_energy(st); Emin = E0; Emax = E0
      n = nint(80000.0_dp/dts(i))
      do it = 1, n
        call cpp_canon_step(st, ierr)
        if (ierr /= 0) exit
        E = cpp_canon_energy(st); Emin = min(Emin, E); Emax = max(Emax, E)
      end do
      osc(i) = (Emax - Emin)/abs(E0)
      print '(A,F6.1,A,ES12.4)', '  CPPsym dt=', dts(i), '  max|dE/E0|=', osc(i)
    end do

    ! Each dt-halving must reduce the oscillation by a factor in [3, 5] (dt^2 ~ 4),
    ! and the finest dt must be well below the buggy-field 1e-3 plateau.
    dt2_ok = .true.
    do i = 1, 3
      ratio = osc(i)/osc(i+1)
      print '(A,I0,A,I0,A,F6.3)', '  CPPsym ratio dt', nint(dts(i)), '/dt', &
           nint(dts(i+1)), ' = ', ratio
      if (ratio < 3.0_dp .or. ratio > 5.0_dp) dt2_ok = .false.
    end do
    call check('CPPsym energy oscillation converges as dt^2', dt2_ok, nfail)
    call check('CPPsym finest-dt osc below buggy plateau (<1e-4)', &
               osc(4) < 1.0e-4_dp, nfail)
  end subroutine test_cpp_sym_convergence

  subroutine test_cpp_banana(nfail)
    ! At a guiding-center-sized dt (=80) the CPP-sym orbit stays confined on a
    ! bounded radial band (the GC banana, not lost to the wall) and conserves the
    ! canonical toroidal momentum p_phi = z(6) exactly: the analytic tokamak is
    ! axisymmetric (A, |B| and the metric have no phi dependence), so p_phi is an
    ! exact invariant of the canonical map, to machine precision. This is the GC
    ! banana signature the big-dt CPP must reproduce.
    integer, intent(inout) :: nfail
    type(cpp_canon_state_t) :: st
    real(dp) :: rmin, rmax, pph0, pphdev
    integer :: it, ierr

    call cpp_canon_init(st, MODEL_CPP_SYM, COORD_TOK, x0, 0.0_dp, 0.0_dp, mu, &
                        mass, charge, 80.0_dp)
    rmin = st%z(1); rmax = st%z(1); pph0 = st%z(6); pphdev = 0.0_dp
    do it = 1, 1000
      call cpp_canon_step(st, ierr)
      if (ierr /= 0) then
        call check('banana ierr==0', .false., nfail); return
      end if
      rmin = min(rmin, st%z(1)); rmax = max(rmax, st%z(1))
      pphdev = max(pphdev, abs(st%z(6) - pph0))
    end do
    print '(A,ES12.4,A,ES12.4)', '  banana r band = ', rmax - rmin, &
         '  p_phi drift = ', pphdev
    ! Bounded banana band, well inside the wall; p_phi conserved to ~machine eps.
    call check('banana r confined in (0,1)', rmin > 0.05_dp .and. rmax < 0.2_dp, nfail)
    call check('banana r oscillates (band > 1e-3)', rmax - rmin > 1.0e-3_dp, nfail)
    call check('banana p_phi invariant (<1e-12)', pphdev < 1.0e-12_dp, nfail)
  end subroutine test_cpp_banana

  subroutine test_jacobian_fd(nfail)
    ! The analytic 6x6 Jacobian must match a central finite difference of the
    ! residual for every model -- the GPU-portable path has no FD fallback, so a
    ! wrong analytic Jacobian would silently degrade Newton convergence. Check at
    ! a generic displaced point for all three models.
    integer, intent(inout) :: nfail
    integer :: models(3), im, j
    real(dp) :: dts(3)
    type(cpp_canon_state_t) :: st
    real(dp) :: zold(6), z(6), jan(6,6), jfd(6,6), rp(6), rm(6), zp(6), zm(6), h, err
    character(12) :: names(3)
    models = [MODEL_CP, MODEL_CPP_SYM, MODEL_CPP_VAR]
    dts = [1.0_dp, 80.0_dp, 800.0_dp]
    names = ['CP-jac      ', 'CPPsym-jac  ', 'CPPvar-jac  ']

    do im = 1, 3
      call cpp_canon_init(st, models(im), COORD_TOK, x0, 0.0_dp, 1.0e-3_dp, mu, &
                          mass, charge, dts(im))
      zold = st%z
      z = zold + [0.003_dp, -0.008_dp, 0.0004_dp, 0.0_dp, 0.0_dp, 0.0_dp]
      call jacobian(st, zold, z, jan)
      h = 1.0e-7_dp
      do j = 1, 6
        zp = z; zm = z; zp(j) = zp(j) + h; zm(j) = zm(j) - h
        call residual(st, zold, zp, rp)
        call residual(st, zold, zm, rm)
        jfd(:,j) = (rp - rm)/(2.0_dp*h)
      end do
      err = maxval(abs(jan - jfd))
      print '(A,A,A,ES10.2)', '  ', trim(names(im)), ' max|Jan-Jfd| = ', err
      call check(trim(names(im))//' analytic==FD', err < 1.0e-6_dp, nfail)
    end do
  end subroutine test_jacobian_fd

  subroutine step_check(st, nstep, ref, name, nfail, ierr)
    type(cpp_canon_state_t), intent(inout) :: st
    integer, intent(in) :: nstep
    real(dp), intent(in) :: ref(3)
    character(*), intent(in) :: name
    integer, intent(inout) :: nfail
    integer, intent(out) :: ierr
    integer :: i
    do i = 1, nstep
      call cpp_canon_step(st, ierr)
      if (ierr /= 0) then
        call check(name//' ierr==0', .false., nfail); return
      end if
    end do
    call assert_vec(st%z(1:3), ref, name, nfail, 1.0e-10_dp)
  end subroutine step_check

  subroutine step_check_tol(st, nstep, ref, name, nfail, ierr, tol)
    type(cpp_canon_state_t), intent(inout) :: st
    integer, intent(in) :: nstep
    real(dp), intent(in) :: ref(3), tol
    character(*), intent(in) :: name
    integer, intent(inout) :: nfail
    integer, intent(out) :: ierr
    integer :: i
    do i = 1, nstep
      call cpp_canon_step(st, ierr)
      if (ierr /= 0) then
        call check(name//' ierr==0', .false., nfail); return
      end if
    end do
    call assert_vec(st%z(1:3), ref, name, nfail, tol)
  end subroutine step_check_tol

  subroutine assert_vec(got, ref, name, nfail, tol)
    real(dp), intent(in) :: got(3), ref(3), tol
    character(*), intent(in) :: name
    integer, intent(inout) :: nfail
    real(dp) :: err
    err = maxval(abs(got - ref))
    if (err <= tol) then
      print '(A,A,A,ES10.2)', 'PASS  ', name, '  maxerr=', err
    else
      print '(A,A,A,ES10.2,A,ES10.2)', 'FAIL  ', name, '  maxerr=', err, ' tol=', tol
      print '(A,3ES23.15)', '   got = ', got
      print '(A,3ES23.15)', '   ref = ', ref
      nfail = nfail + 1
    end if
  end subroutine assert_vec

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

end program test_cpp_canonical
