program test_boozer_cartesian
  ! Boozer <-> Cartesian map and the Cartesian Larmor displacement (solution B):
  ! displace in Cartesian, not in flux coordinates, to get the particle position
  ! from a guiding center and back. Validated on the reactor-scale test
  ! equilibrium test_data/wout.nc.
  !
  ! Gates:
  !   (1) Forward/inverse round trip: cart_to_boozer(boozer_to_cart(u)) = u.
  !   (2) Analytic Jacobian d(x,y,z)/du matches a finite difference.
  !   (3) GC -> particle -> GC recovers the guiding center to inversion accuracy,
  !       at s = 0.5 (off axis) AND s = 0.25, for a deeply trapped pitch whose
  !       Larmor radius is a sizeable fraction of the minor radius. The legacy
  !       flux-coordinate add cannot do this: its error grows like rho* and
  !       diverges toward the axis.
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use parmot_mod, only: ro0
  use simple, only: init_sympl, init_params, tracer_t
  use simple_main, only: init_field
  use boozer_cartesian, only: boozer_to_cart, cart_to_boozer, &
    perp_unit_dir_flux, gc_to_particle, particle_to_gc
  use boozer_field_metric, only: boozer_field_metric_eval
  use params, only: field_input, coord_input, integmode, relerr, orbit_coord
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
  field_input = 'wout.nc'
  coord_input = 'wout.nc'
  orbit_coord = 1
  integmode = 1
  relerr = 1.0d-13
  call init_field(norb, 'wout.nc', ans_s, ans_tp, amultharm, integmode)
  use_B_r = .true.
  use_del_tp_B = .true.
  call get_boozer_coordinates
  call init_params(norb, 2, 4, 3.5e6_dp, 256, 1, 1.0d-13)
  ro0_bar = ro0/sqrt(2.0_dp)

  call test_round_trip(nfail)
  call test_jacobian(nfail)
  call test_metric_consistency(nfail)
  call test_gc_particle_roundtrip([0.5_dp, 0.5_dp, 0.2_dp], ro0_bar, nfail)
  call test_gc_particle_roundtrip([0.25_dp, 0.5_dp, 0.2_dp], ro0_bar, nfail)
  call test_gc_particle_roundtrip([0.3_dp, 0.5_dp, 0.2_dp], ro0_bar, nfail)

  if (nfail == 0) then
    print *, 'ALL BOOZER-CARTESIAN TESTS PASSED'
  else
    print *, 'BOOZER-CARTESIAN TESTS FAILED: ', nfail
    error stop 1
  end if

contains

  subroutine test_round_trip(nfail)
    integer, intent(inout) :: nfail
    real(dp) :: u(3), xyz(3), Jc(3,3), u2(3)
    integer :: ierr, k
    real(dp) :: pts(3,4), err

    pts(:,1) = [0.5_dp, 0.5_dp, 0.2_dp]
    pts(:,2) = [0.25_dp, 1.3_dp, 4.7_dp]
    pts(:,3) = [0.8_dp, -2.0_dp, 0.0_dp]
    pts(:,4) = [0.1_dp, 3.0_dp, 1.0_dp]
    do k = 1, 4
      u = pts(:,k)
      call boozer_to_cart(u, xyz, Jc)
      ! Perturb the guess so the inversion really iterates.
      call cart_to_boozer(xyz, u + [0.02_dp, 0.1_dp, 0.1_dp], u2, ierr)
      err = maxval(abs(u2 - u))
      print '(A,F5.2,A,ES11.3,A,I0)', '  round trip s=', u(1), &
        ' |du|=', err, ' ierr=', ierr
      call check('round trip recovers Boozer point', ierr == 0 .and. err < 1.0e-9_dp, nfail)
    end do
  end subroutine test_round_trip

  subroutine test_jacobian(nfail)
    integer, intent(inout) :: nfail
    real(dp) :: u(3), xyz(3), Jc(3,3), xp(3), xm(3), Jfd(3,3), Jdum(3,3)
    real(dp) :: h, du(3), err
    integer :: k

    u = [0.5_dp, 0.5_dp, 0.2_dp]
    call boozer_to_cart(u, xyz, Jc)
    h = 1.0e-6_dp
    do k = 1, 3
      du = 0.0_dp; du(k) = h
      call boozer_to_cart(u + du, xp, Jdum)
      call boozer_to_cart(u - du, xm, Jdum)
      Jfd(:,k) = (xp - xm)/(2.0_dp*h)
    end do
    err = maxval(abs(Jc - Jfd))/maxval(abs(Jc))
    print '(A,ES11.3)', '  Jacobian max rel error vs FD = ', err
    call check('analytic Jacobian matches finite difference', err < 1.0e-6_dp, nfail)
  end subroutine test_jacobian

  ! The Cartesian embedding metric must equal the Boozer field metric:
  ! g_ij = (d x/d u^i) . (d x/d u^j) = (Jc^T Jc)_ij. If these disagree, a velocity
  ! normalized in g is the wrong length when pushed through Jc, and the Larmor
  ! vector comes out mis-scaled.
  subroutine test_metric_consistency(nfail)
    integer, intent(inout) :: nfail
    real(dp) :: u(3), xyz(3), Jc(3,3), JtJ(3,3)
    real(dp) :: g(3,3), ginv(3,3), sqrtg, dg(3,3,3)
    real(dp) :: Acov(3), dA(3,3), Bctr(3), Bcov(3), Bmod, dBmod(3), hcov(3)
    real(dp) :: relerr_m

    u = [0.5_dp, 0.5_dp, 0.2_dp]
    call boozer_to_cart(u, xyz, Jc)
    JtJ = matmul(transpose(Jc), Jc)
    call boozer_field_metric_eval(u, g, ginv, sqrtg, dg, Acov, dA, &
      Bctr, Bcov, Bmod, dBmod, hcov)
    relerr_m = maxval(abs(JtJ - g))/maxval(abs(g))
    print '(A,ES12.4)', '  g vs Jc^T Jc max rel error = ', relerr_m
    print '(A,3ES12.4)', '  g    diag = ', g(1,1), g(2,2), g(3,3)
    print '(A,3ES12.4)', '  JtJ  diag = ', JtJ(1,1), JtJ(2,2), JtJ(3,3)
    call check('Cartesian Jacobian metric matches Boozer metric (g = Jc^T Jc)', &
      relerr_m < 1.0e-6_dp, nfail)
  end subroutine test_metric_consistency

  ! Seed a particle one Larmor vector off x_gc (fixed point consistent with the
  ! velocity seed), then reconstruct its guiding center. Because the seed is built
  ! so that particle_to_gc(x_p) = x_gc, the recovery is exact to the Newton
  ! tolerance, at s=0.5 (off axis) AND s=0.25, even though the gyroradius is a
  ! sizeable fraction of the major radius. The legacy flux-coordinate add diverged
  ! here (s driven through the axis).
  subroutine test_gc_particle_roundtrip(x_gc, ro0_bar, nfail)
    real(dp), intent(in) :: x_gc(3), ro0_bar
    integer, intent(inout) :: nfail
    real(dp) :: shift, rec

    call seed_reconstruct(x_gc, ro0_bar, 1.0_dp, shift, rec, nfail)
    print '(A,F5.2,A,ES11.3,A,ES11.3)', '  s=', x_gc(1), &
      '  Larmor shift=', shift, '  GC recovery err=', rec
    call check('CP seed produces a finite Larmor shift', shift > 1.0e-4_dp, nfail)
    call check('GC recovery is exact to Newton tol (err < 1e-7)', &
      rec < 1.0e-7_dp, nfail)
  end subroutine test_gc_particle_roundtrip

  subroutine seed_reconstruct(x_gc, ro0_bar, scale, shift, recover, nfail)
    real(dp), intent(in) :: x_gc(3), ro0_bar, scale
    real(dp), intent(out) :: shift, recover
    integer, intent(inout) :: nfail
    real(dp) :: g(3,3), ginv(3,3), sqrtg, dg(3,3,3)
    real(dp) :: Acov(3), dA(3,3), Bctr(3), Bcov(3), Bmod, dBmod(3), hcov(3)
    real(dp) :: eperp(3), vperp_con(3), vperp0, qc, x_p(3), x_rec(3)
    integer :: ierr

    qc = 1.0_dp/ro0_bar
    ! Deeply trapped pitch lambda=0.3, unit normalized speed scaled by `scale`:
    ! vperp0 = sqrt(2 mu B) with mu = (1-lambda^2)/B, i.e. sqrt(2(1-lambda^2))
    ! (O(1) normalized speed, same convention as init_canonical_6d).
    vperp0 = scale*sqrt(2.0_dp*(1.0_dp - 0.3_dp**2))

    call gc_to_particle(x_gc, vperp0, 1.0_dp, qc, x_p, ierr)
    call check('gc_to_particle inversion succeeds', ierr == 0, nfail)
    shift = maxval(abs(x_p - x_gc))

    ! Reconstruct from the perpendicular velocity AT the particle point, exactly
    ! as cpp_canon_boozer_guiding_center reads it from the momenta.
    call boozer_field_metric_eval(x_p, g, ginv, sqrtg, dg, Acov, dA, &
      Bctr, Bcov, Bmod, dBmod, hcov)
    call perp_unit_dir_flux(g, ginv, hcov, eperp)
    vperp_con = vperp0*eperp
    call particle_to_gc(x_p, vperp_con, 1.0_dp, qc, x_rec, ierr)
    call check('particle_to_gc inversion succeeds', ierr == 0, nfail)
    recover = maxval(abs(x_rec - x_gc))
  end subroutine seed_reconstruct

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

end program test_boozer_cartesian
