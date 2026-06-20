program test_boozer_field_metric
  ! GATE for the single-source Boozer boozer_field_metric (Boozer analogue of
  ! test_vmec_field_metric). On the real QA equilibrium test_data/wout.nc, at
  ! several interior points in Boozer coordinates u = (s, vartheta_B, varphi_B):
  !
  !   (a) h_i g^ij h_j = 1. The metric is pulled back from the VMEC chart through
  !       the angle Jacobian; the field (Bcov, |B|) is the production Boozer
  !       spline. The two are INDEPENDENT spline families, unlike the VMEC case
  !       where |B| is DEFINED from the metric (and the identity holds to 1e-13).
  !       Here the identity is limited by the production Boozer field's own
  !       internal consistency: the angular covariant components I(s), g(s) are
  !       flux functions, but the radial covariant B_s is reconstructed on a grid
  !       (compute_br_from_symflux), so |B|^2 = g^ij B_i B_j matches the splined
  !       |B|^2 only to ~1e-4 relative. That is a field-construction floor, NOT a
  !       metric error: check (e) shows the metric itself is geometrically exact.
  !   (b) |B| from boozer_field_metric matches the production Boozer |B| from
  !       splint_boozer_coord to 1e-12 (taken from the same spline).
  !   (c) analytic dg vs 4th-order central difference of g to 1e-6 (machine
  !       precision: the geometry point and the Jacobian use one angle map).
  !   (d) g . ginv = I to 1e-10.
  !   (e) geometric |B| = sqrt(g^ij B_i B_j) from the pulled-back metric vs VMEC
  !       |B| at the mapped angles to 1e-4. VMEC |B| is geometrically exact, so
  !       this is the strict proof the metric transform is right and isolates the
  !       residual in (a) to the production field spline.
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use simple, only: tracer_t
  use simple_main, only: init_field
  use boozer_sub, only: get_boozer_coordinates, splint_boozer_coord, boozer_to_vmec
  use boozer_coordinates_mod, only: use_B_r, use_del_tp_B
  use boozer_field_metric, only: boozer_field_metric_eval
  use vmec_field_metric, only: vmec_field_metric_eval
  use params, only: coord_input
  implicit none

  integer, parameter :: npts = 5
  real(dp), parameter :: pts(3, npts) = reshape([ &
       0.15_dp, 0.6_dp, 0.2_dp, &
       0.30_dp, 1.7_dp, 0.9_dp, &
       0.50_dp, 3.1_dp, 2.4_dp, &
       0.70_dp, 4.8_dp, 1.1_dp, &
       0.90_dp, 5.9_dp, 0.4_dp], [3, npts])

  type(tracer_t) :: norb
  integer :: nfail, ip
  real(dp) :: worst_hgh, worst_dg, worst_bmod, worst_ggi, worst_geo

  coord_input = 'wout.nc'
  call init_field(norb, 'wout.nc', 5, 5, 5, -1)
  ! The metric transform needs the Boozer-side angle map (use_del_tp_B) and the
  ! radial covariant field component B_s (use_B_r). Set before building data.
  use_B_r = .true.
  use_del_tp_B = .true.
  call get_boozer_coordinates
  print *, 'Built Boozer coordinates from wout.nc'

  nfail = 0
  worst_hgh = 0.0_dp
  worst_dg = 0.0_dp
  worst_bmod = 0.0_dp
  worst_ggi = 0.0_dp
  worst_geo = 0.0_dp

  print '(A)', '  point (s, vth_B, vph_B)        h_i g^ij h_j        |h.g.h - 1|'
  do ip = 1, npts
    call check_point(pts(:, ip))
  end do

  print '(A,ES12.4)', '  worst |h_i g^ij h_j - 1| (field-spline floor) = ', worst_hgh
  print '(A,ES12.4)', '  worst ||B|_metric - |B|_spline| / |B|         = ', worst_bmod
  print '(A,ES12.4)', '  worst |dg analytic - dg FD| (relative)        = ', worst_dg
  print '(A,ES12.4)', '  worst |g.ginv - I|                            = ', worst_ggi
  print '(A,ES12.4)', '  worst |sqrt(ginv B B) - |B|_VMEC| / |B| (geo) = ', worst_geo

  ! Strict metric/derivative checks at machine precision.
  call check('|B| metric vs spline to 1e-12', worst_bmod < 1.0e-12_dp, nfail)
  call check('dg analytic vs central FD to 1e-6', worst_dg < 1.0e-6_dp, nfail)
  call check('g.ginv = I to 1e-10', worst_ggi < 1.0e-10_dp, nfail)
  ! Strict geometric proof that the pulled-back metric is correct.
  call check('geometric |B| vs VMEC |B| to 1e-4', worst_geo < 1.0e-4_dp, nfail)
  ! h.g.h is bounded by the production Boozer field's covariant-component
  ! accuracy (B_s grid reconstruction), not by the metric: assert that floor.
  call check('h_i g^ij h_j = 1 to 2e-4 (field floor)', worst_hgh < 2.0e-4_dp, nfail)

  if (nfail == 0) then
    print *, 'ALL boozer_field_metric GATE TESTS PASSED'
  else
    print *, 'boozer_field_metric GATE TESTS FAILED: ', nfail
    error stop 1
  end if

contains

  subroutine check_point(u)
    real(dp), intent(in) :: u(3)

    real(dp) :: g(3,3), ginv(3,3), sqrtg, dg(3,3,3)
    real(dp) :: Acov(3), dA(3,3), Bctr(3), Bcov(3), Bmod, dBmod(3), hcov(3)
    real(dp) :: hgh, rel_dg, ggi, bmod_ref, bmod_err, bgeo, geo_err
    integer :: i, j

    ! Production Boozer |B| reference.
    bmod_ref = booz_bmod(u)

    call boozer_field_metric_eval(u, g, ginv, sqrtg, dg, Acov, dA, &
                                  Bctr, Bcov, Bmod, dBmod, hcov)

    ! (a) h_i g^ij h_j with covariant hcov contracted with the inverse metric.
    hgh = 0.0_dp
    do i = 1, 3
      do j = 1, 3
        hgh = hgh + hcov(i)*ginv(i,j)*hcov(j)
      end do
    end do
    worst_hgh = max(worst_hgh, abs(hgh - 1.0_dp))
    print '(A,3F7.3,A,F18.15,A,ES10.2)', '  (', u, ')  ', hgh, '   ', abs(hgh - 1.0_dp)

    ! (b) |B| match.
    bmod_err = abs(Bmod - bmod_ref)/abs(bmod_ref)
    worst_bmod = max(worst_bmod, bmod_err)

    ! (c) dg vs central FD.
    rel_dg = max_rel_dg_error(u, dg)
    worst_dg = max(worst_dg, rel_dg)

    ! (d) g . ginv = I.
    ggi = 0.0_dp
    do i = 1, 3
      do j = 1, 3
        ggi = max(ggi, abs(dot3(g(i,:), ginv(:,j)) - id(i,j)))
      end do
    end do
    worst_ggi = max(worst_ggi, ggi)

    ! (e) geometric |B| from metric+covariant field vs VMEC |B| at mapped angles.
    bgeo = 0.0_dp
    do i = 1, 3
      do j = 1, 3
        bgeo = bgeo + ginv(i,j)*Bcov(i)*Bcov(j)
      end do
    end do
    bgeo = sqrt(bgeo)
    geo_err = abs(bgeo - vmec_bmod_at_mapped(u))/abs(bmod_ref)
    worst_geo = max(worst_geo, geo_err)
  end subroutine check_point

  real(dp) function dot3(a, b) result(s)
    real(dp), intent(in) :: a(3), b(3)
    s = a(1)*b(1) + a(2)*b(2) + a(3)*b(3)
  end function dot3

  real(dp) function id(i, j) result(v)
    integer, intent(in) :: i, j
    v = 0.0_dp
    if (i == j) v = 1.0_dp
  end function id

  ! Production Boozer |B| from splint_boozer_coord at u = (s, vartheta_B, varphi_B).
  real(dp) function booz_bmod(u) result(bm)
    real(dp), intent(in) :: u(3)
    real(dp) :: A_theta, A_phi, dA_theta_dr, dA_phi_dr, d2A_phi_dr2, d3A_phi_dr3
    real(dp) :: B_vth, dB_vth, d2B_vth, B_vph, dB_vph, d2B_vph
    real(dp) :: dBmod_B(3), d2Bmod_B(6), B_r, dB_r(3), d2B_r(6)
    call splint_boozer_coord(u(1), u(2), u(3), 0, &
                             A_theta, A_phi, dA_theta_dr, dA_phi_dr, &
                             d2A_phi_dr2, d3A_phi_dr3, &
                             B_vth, dB_vth, d2B_vth, B_vph, dB_vph, d2B_vph, &
                             bm, dBmod_B, d2Bmod_B, B_r, dB_r, d2B_r)
  end function booz_bmod

  ! Geometrically exact VMEC |B| at the VMEC angles that map to the Boozer point.
  real(dp) function vmec_bmod_at_mapped(u) result(bm)
    real(dp), intent(in) :: u(3)
    real(dp) :: theta_V, varphi_V, uV(3)
    real(dp) :: g(3,3), ginv(3,3), sqrtg, dg(3,3,3)
    real(dp) :: Acov(3), dA(3,3), Bctr(3), Bcov(3), dBmod(3), hcov(3)
    call boozer_to_vmec(u(1), u(2), u(3), theta_V, varphi_V)
    uV = [u(1), theta_V, varphi_V]
    call vmec_field_metric_eval(uV, g, ginv, sqrtg, dg, Acov, dA, &
                                Bctr, Bcov, bm, dBmod, hcov)
  end function vmec_bmod_at_mapped

  ! Relative max error of analytic dg against a 4th-order central difference of g.
  real(dp) function max_rel_dg_error(u, dg) result(maxerr)
    real(dp), intent(in) :: u(3), dg(3,3,3)
    real(dp) :: gscale, dgfd, hstep(3)
    real(dp) :: gp1(3,3), gm1(3,3), gp2(3,3), gm2(3,3)
    integer :: i, j, k

    hstep = [1.0e-4_dp, 1.0e-4_dp, 1.0e-4_dp]
    gscale = 0.0_dp
    maxerr = 0.0_dp
    do k = 1, 3
      call eval_g(shift(u, k, hstep(k)), gp1)
      call eval_g(shift(u, k, -hstep(k)), gm1)
      call eval_g(shift(u, k, 2.0_dp*hstep(k)), gp2)
      call eval_g(shift(u, k, -2.0_dp*hstep(k)), gm2)
      do j = 1, 3
        do i = 1, 3
          gscale = max(gscale, abs(gp1(i,j)))
        end do
      end do
      do j = 1, 3
        do i = 1, 3
          dgfd = (-gp2(i,j) + 8.0_dp*gp1(i,j) - 8.0_dp*gm1(i,j) + gm2(i,j)) &
                 / (12.0_dp*hstep(k))
          maxerr = max(maxerr, abs(dg(i,j,k) - dgfd))
        end do
      end do
    end do
    maxerr = maxerr/max(gscale, 1.0_dp)
  end function max_rel_dg_error

  function shift(u, k, d) result(uu)
    real(dp), intent(in) :: u(3), d
    integer, intent(in) :: k
    real(dp) :: uu(3)
    uu = u
    uu(k) = uu(k) + d
  end function shift

  subroutine eval_g(u, g)
    real(dp), intent(in) :: u(3)
    real(dp), intent(out) :: g(3,3)
    real(dp) :: ginv(3,3), sqrtg, dg(3,3,3)
    real(dp) :: Acov(3), dA(3,3), Bctr(3), Bcov(3), Bmod, dBmod(3), hcov(3)
    call boozer_field_metric_eval(u, g, ginv, sqrtg, dg, Acov, dA, &
                                  Bctr, Bcov, Bmod, dBmod, hcov)
  end subroutine eval_g

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

end program test_boozer_field_metric
