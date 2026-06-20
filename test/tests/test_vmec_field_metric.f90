program test_vmec_field_metric
  ! GATE for the single-source plain vmec_field_metric (FACTS design 1).
  !
  ! 1. h_i g^ij h_j = 1 to ~1e-13 at several interior points. This is the
  !    consistency that the dual-source path fails (it gave 1.009): because
  !    h_i = g_ij B^j / |B| and |B| = sqrt(g_ij B^i B^j) come from the SAME g,
  !    the identity must hold to round-off.
  ! 2. dg_ij,k analytic vs central finite difference of g_ij to ~1e-8.
  !
  ! Runs on the real QA equilibrium test_data/wout.nc (symlinked into the test
  ! binary dir) at interior points away from the axis (s -> 0 is metric-singular).
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use new_vmec_stuff_mod, only: netcdffile, multharm, ns_s, ns_tp
  use spline_vmec_sub, only: spline_vmec_data
  use vmec_field_metric, only: vmec_field_metric_eval
  implicit none

  integer, parameter :: npts = 5
  real(dp), parameter :: pts(3, npts) = reshape([ &
       0.15_dp, 0.6_dp, 0.2_dp, &
       0.30_dp, 1.7_dp, 0.9_dp, &
       0.50_dp, 3.1_dp, 2.4_dp, &
       0.70_dp, 4.8_dp, 1.1_dp, &
       0.90_dp, 5.9_dp, 0.4_dp], [3, npts])

  integer :: nfail, ip
  real(dp) :: worst_hgh, worst_dg

  netcdffile = 'wout.nc'
  ns_s = 5
  ns_tp = 5
  multharm = 3
  call spline_vmec_data
  print *, 'Splined VMEC data from wout.nc'

  nfail = 0
  worst_hgh = 0.0_dp
  worst_dg = 0.0_dp

  print '(A)', '  point (s, theta, phi)         h_i g^ij h_j        |h.g.h - 1|'
  do ip = 1, npts
    call check_point(pts(:, ip), nfail, worst_hgh, worst_dg)
  end do

  print '(A,ES12.4)', '  worst |h_i g^ij h_j - 1| over all points = ', worst_hgh
  print '(A,ES12.4)', '  worst |dg analytic - dg FD| (relative)    = ', worst_dg

  call check('h_i g^ij h_j = 1 to 1e-13', worst_hgh < 1.0e-13_dp, nfail)
  call check('dg analytic vs central FD to 1e-8', worst_dg < 1.0e-8_dp, nfail)

  if (nfail == 0) then
    print *, 'ALL vmec_field_metric GATE TESTS PASSED'
  else
    print *, 'vmec_field_metric GATE TESTS FAILED: ', nfail
    error stop 1
  end if

contains

  subroutine check_point(u, nfail, worst_hgh, worst_dg)
    real(dp), intent(in) :: u(3)
    integer, intent(inout) :: nfail
    real(dp), intent(inout) :: worst_hgh, worst_dg

    real(dp) :: g(3,3), ginv(3,3), sqrtg, dg(3,3,3)
    real(dp) :: Acov(3), dA(3,3), Bctr(3), Bcov(3), Bmod, dBmod(3), hcov(3)
    real(dp) :: hgh, rel_dg
    integer :: i, j

    call vmec_field_metric_eval(u, g, ginv, sqrtg, dg, Acov, dA, &
                                Bctr, Bcov, Bmod, dBmod, hcov)

    ! h_i g^ij h_j with hcov covariant -> contract with the inverse metric.
    hgh = 0.0_dp
    do i = 1, 3
      do j = 1, 3
        hgh = hgh + hcov(i)*ginv(i,j)*hcov(j)
      end do
    end do
    worst_hgh = max(worst_hgh, abs(hgh - 1.0_dp))
    print '(A,3F7.3,A,F18.15,A,ES10.2)', '  (', u, ')  ', hgh, '   ', abs(hgh - 1.0_dp)

    rel_dg = max_rel_dg_error(u, dg)
    worst_dg = max(worst_dg, rel_dg)
    if (.not. (abs(hgh - 1.0_dp) < 1.0e-13_dp)) nfail = nfail + 1
  end subroutine check_point

  ! Relative max error of the analytic dg against a 4th-order central difference
  ! of g_ij in each direction. Step sized for ~1e-8 relative truncation; the
  ! denominator is the local metric scale so cm^2-sized entries are compared fairly.
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
          ! 4th-order central difference: (-f2 + 8 f1 - 8 f-1 + f-2)/(12 h).
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
    call vmec_field_metric_eval(u, g, ginv, sqrtg, dg, Acov, dA, &
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

end program test_vmec_field_metric
