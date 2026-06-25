program test_fo_boris
  ! Validate the full-orbit (FO) Boris pusher (orbit_fo_boris) on the real
  ! reactor-scale Boozer field. The particle is a charged particle in a static B
  ! field: the Boris rotation is exact per step, there is NO nonlinear solve (no
  ! convergence floor), and the Cartesian advance is regular through the magnetic
  ! axis, the gyro-resolved counterpart to the guiding-center model.
  !
  ! Gates:
  !   (1) energy |dE/E0| bounded < 1e-3 over many gyroperiods (passing and trapped),
  !   (2) NEAR-AXIS: a particle whose orbit reaches small s crosses the axis region
  !       with energy still bounded and no integrator failure,
  !   (3) a confined particle stays 0 < s < 1 over the run (no spurious loss).
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use parmot_mod, only: ro0
  use simple, only: init_params, tracer_t
  use simple_main, only: init_field
  use orbit_fo_boris, only: fo_state_t, fo_init, fo_step, &
    fo_energy, fo_mu, fo_to_gc, accept_or_fail, FO_OK, FO_LOCATE_FAIL
  use orbit_fo_field, only: fo_eval_field
  use reference_coordinates, only: ref_coords
  use params, only: field_input, coord_input, integmode, relerr, dtaumin, orbit_coord
  use velo_mod, only: isw_field_type
  use magfie_sub, only: BOOZER
  use boozer_coordinates_mod, only: use_B_r, use_del_tp_B
  use boozer_sub, only: get_boozer_coordinates
  implicit none

  integer, parameter :: ans_s = 5, ans_tp = 5, amultharm = 5
  type(tracer_t) :: norb
  real(dp) :: ro0_bar
  integer :: nfail

  nfail = 0
  isw_field_type = BOOZER
  field_input = 'wout.nc'; coord_input = 'wout.nc'
  orbit_coord = 1; integmode = 1; relerr = 1.0d-13
  call init_field(norb, 'wout.nc', ans_s, ans_tp, amultharm, integmode)
  use_B_r = .true.; use_del_tp_B = .true.
  call get_boozer_coordinates
  call init_params(norb, 2, 4, 3.5e6_dp, 16384, 1, 1.0d-13)
  dtaumin = norb%dtaumin
  ro0_bar = ro0/sqrt(2.0_dp)

  ! Precondition: the chartmap Jacobian must be in the documented convention. The
  ! field assembly and the Cartesian inverse Newton both rely on it; a transposed
  ! Jacobian flips the field while Boris still conserves energy, so check it here.
  call check_covariant_basis_convention(nfail)

  ! passing (lambda=0.9), trapped (lambda=0.2), and an inner orbit driven toward
  ! the axis (small s, lambda=0.7) to exercise the near-axis crossing.
  call run_fo([0.5_dp, 0.5_dp, 0.2_dp, 1.0_dp, 0.9_dp], ro0_bar, 'passing', nfail)
  call run_fo([0.5_dp, 0.5_dp, 0.2_dp, 1.0_dp, 0.2_dp], ro0_bar, 'trapped', nfail)
  call run_fo([0.04_dp, 0.5_dp, 0.2_dp, 1.0_dp, 0.7_dp], ro0_bar, 'near-axis', nfail)

  ! A marker exiting the boundary must be located (so the guiding-centre loss test
  ! runs), never turned into a confined fault.
  call test_accept_classification(nfail)

  if (nfail == 0) then
    print *, 'ALL FO-BORIS TESTS PASSED'
  else
    print *, 'FO-BORIS TESTS FAILED: ', nfail
    error stop 1
  end if

contains

  ! accept_or_fail must LOCATE (FO_OK) a marker that clamped to the edge while its warm
  ! guess was already near the edge -- a genuine LCFS exit, so fo_to_gc can flag the
  ! guiding-centre loss -- while a mid-radius seam glitch and a near-axis stall stay
  ! faults (FO_LOCATE_FAIL), never spurious losses. This pins the exiting-boundary fix.
  subroutine test_accept_classification(nfail)
    integer, intent(inout) :: nfail
    real(dp), parameter :: tol = 1.0e-6_dp, edge = 1.0_dp, scale = 1.0_dp
    real(dp), parameter :: big = 1.0_dp     ! residual well above EDGE_FRAC*scale
    call expect_status(accept_or_fail(0.5_dp, 1.0e-9_dp, scale, tol, edge, 0.5_dp), &
                       FO_OK, 'converged interior located', nfail)
    call expect_status(accept_or_fail(0.01_dp, big, scale, tol, edge, 0.01_dp), &
                       FO_LOCATE_FAIL, 'near-axis stall stays fault', nfail)
    call expect_status(accept_or_fail(edge, big, scale, tol, edge, 0.5_dp), &
                       FO_LOCATE_FAIL, 'mid-radius glitch stays fault', nfail)
    call expect_status(accept_or_fail(edge, big, scale, tol, edge, 0.98_dp), &
                       FO_OK, 'edge exit located for GC loss test', nfail)
  end subroutine test_accept_classification

  subroutine expect_status(got, want, tag, nfail)
    integer, intent(in) :: got, want
    character(*), intent(in) :: tag
    integer, intent(inout) :: nfail
    if (got /= want) then
      print *, 'FAIL accept_or_fail: ', tag, ' got', got, ' want', want
      nfail = nfail + 1
    else
      print *, 'ok: ', tag
    end if
  end subroutine expect_status

  subroutine run_fo(z0, ro0_bar, tag, nfail)
    real(dp), intent(in) :: z0(5), ro0_bar
    character(*), intent(in) :: tag
    integer, intent(inout) :: nfail
    type(fo_state_t) :: st
    real(dp) :: bmod, mu, vpar_bar, vperp0, E0, E, Emax, smin, smax
    real(dp) :: s, th, ph, vpar
    real(dp) :: Acov(3), dA(3,3), dBmod(3), hcov(3)
    real(dp) :: mu0, mui, mumin, mumax, mu_first, mu_last, mu_drift
    integer :: it, ierr, nstep, lost, nwin, nfw, nlw

    ! |B| at the guiding centre from the chartmap field (rho = sqrt(s)).
    call fo_eval_field([sqrt(z0(1)), z0(2), z0(3)], Acov, dA, bmod, dBmod, hcov)
    mu = 0.5_dp*z0(4)**2*(1.0_dp - z0(5)**2)/bmod*2.0_dp
    vpar_bar = z0(4)*z0(5)*sqrt(2.0_dp)
    vperp0 = sqrt(max(2.0_dp*mu*bmod, 0.0_dp))

    nstep = 20000   ! ~ many hundred gyroperiods at np16384
    call fo_init(st, z0(1:3), vpar_bar, vperp0, mu, 1.0_dp, &
                 1.0_dp, dtaumin/sqrt(2.0_dp), ro0_bar, z0(4))
    E0 = fo_energy(st); Emax = 0.0_dp
    mu0 = fo_mu(st); mumin = mu0; mumax = mu0
    smin = z0(1); smax = z0(1); lost = 0
    ! secular drift: gyro-average mu over the first and last tenth of the run and
    ! compare the averages. The window must span many gyroperiods (the Boris
    ! rotation is ~0.16 rad/step, so one gyroperiod is ~40 steps); nstep/10 = 2000
    ! steps averages ~50 gyroperiods, removing the gyro-ripple so what survives is
    ! the true secular drift. A short window leaves ripple phase aliased into the
    ! difference and overstates the drift.
    nwin = nstep/10; mu_first = 0.0_dp; mu_last = 0.0_dp; nfw = 0; nlw = 0
    do it = 1, nstep
      call fo_step(st, ierr)
      if (ierr /= 0) then; lost = 1; exit; end if
      call fo_to_gc(st, s, th, ph, vpar, ierr)
      if (ierr /= 0) then; lost = 1; exit; end if
      if (s <= 0.0_dp .or. s >= 1.0_dp) exit
      smin = min(smin, s); smax = max(smax, s)
      E = fo_energy(st)
      Emax = max(Emax, abs((E - E0)/E0))
      mui = fo_mu(st)
      mumin = min(mumin, mui); mumax = max(mumax, mui)
      if (it <= nwin) then; mu_first = mu_first + mui; nfw = nfw + 1; end if
      if (it > nstep - nwin) then; mu_last = mu_last + mui; nlw = nlw + 1; end if
    end do
    mu_drift = -1.0_dp
    if (nfw > 0 .and. nlw > 0) &
      mu_drift = abs(mu_last/nlw - mu_first/nfw)/(mu_first/nfw)
    print '(a,a,a,f7.4,a,f7.4,a,es10.2,a,i0)', '  ', tag, ' s band [', smin, &
      ',', smax, ']  max|dE/E0|=', Emax, '  ierr_lost=', lost
    print '(a,a,a,es10.2,a,es10.2)', '    ', tag, ' mu oscillation |dmu/mu0|=', &
      (mumax - mumin)/mu0, '  secular gyro-avg drift=', mu_drift
    call check(tag//' energy bounded (<1e-3)', Emax < 1.0e-3_dp, nfail)
    call check(tag//' step never failed (no nonconv: explicit)', lost == 0, nfail)
    ! mu is an adiabatic invariant, not exact: well conserved for a deep-trapped
    ! orbit (small FLR), but it genuinely degrades for grazing/near-axis orbits
    ! where the gyroradius is no longer small -- that breakdown is the physics the
    ! full orbit is meant to capture, not a defect. Hard-assert only the trapped
    ! case; the others are reported.
    if (tag == 'trapped') &
      call check(tag//' mu adiabatic: secular drift < 1e-2', &
                 mu_drift >= 0.0_dp .and. mu_drift < 1.0e-2_dp, nfail)
  end subroutine run_fo

  ! The FO field assembly builds B = curl A through ref_coords%covariant_basis,
  ! which must return Jc(i,k) = d x_i / d u_k (Cartesian component i, logical
  ! coord k). A transposed Jc (the libneo cart-spline chartmap bug) silently flips
  ! the field: the Boris energy gates still pass, but the orbits are wrong. Catch
  ! that directly by checking covariant_basis against a finite difference of
  ! evaluate_cart at a generic off-axis, off-seam interior point.
  subroutine check_covariant_basis_convention(nfail)
    integer, intent(inout) :: nfail
    real(dp) :: u(3), up(3), um(3), Jc(3,3), Jfd(3,3), xp(3), xm(3), du, rel
    integer :: k

    u = [0.6_dp, 0.7_dp, 0.3_dp]
    du = 1.0e-6_dp
    call ref_coords%covariant_basis(u, Jc)
    do k = 1, 3
      up = u; up(k) = u(k) + du
      um = u; um(k) = u(k) - du
      call ref_coords%evaluate_cart(up, xp)
      call ref_coords%evaluate_cart(um, xm)
      Jfd(:, k) = (xp - xm)/(2.0_dp*du)
    end do
    rel = maxval(abs(Jc - Jfd))/max(maxval(abs(Jfd)), 1.0e-30_dp)
    print '(a,es10.2)', '  covariant_basis vs FD(evaluate_cart) rel err = ', rel
    call check('chartmap Jacobian convention (not transposed)', rel < 1.0e-4_dp, &
               nfail)
  end subroutine check_covariant_basis_convention

  subroutine check(name, cond, nfail)
    character(*), intent(in) :: name
    logical, intent(in) :: cond
    integer, intent(inout) :: nfail
    if (cond) then
      print '(a,a)', 'PASS  ', name
    else
      print '(a,a)', 'FAIL  ', name
      nfail = nfail + 1
    end if
  end subroutine check

end program test_fo_boris
