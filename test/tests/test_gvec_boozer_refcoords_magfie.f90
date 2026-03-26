program test_gvec_boozer_refcoords_magfie
   use, intrinsic :: iso_fortran_env, only: dp => real64
   use field_gvec, only: gvec_field_t, create_gvec_field
   use magfie_sub, only: REFCOORDS, init_magfie, magfie, set_magfie_refcoords_field

   implicit none

   class(gvec_field_t), allocatable :: raw_field
   integer, parameter :: n_r = 4
   integer, parameter :: n_th = 4
   integer, parameter :: n_ph = 3
   real(dp), parameter :: hs = 1.0e-4_dp
   real(dp), parameter :: ht = 2.5e-4_dp
   real(dp), parameter :: hp = 2.5e-4_dp
   real(dp), parameter :: tol_bmod = 7.0e-2_dp
   real(dp), parameter :: tol_sqrtg = 1.0e-1_dp
   real(dp), parameter :: tol_hcov = 1.2e-1_dp
   real(dp), parameter :: tol_hctr = 1.2e-1_dp
   real(dp), parameter :: tol_bder = 2.0e-1_dp
   real(dp), parameter :: tol_hcurl = 3.0e-1_dp

   real(dp) :: rho
   real(dp) :: theta
   real(dp) :: phi
   real(dp) :: phi_period
   real(dp) :: x_r(3)
   real(dp) :: bmod_ref
   real(dp) :: sqrtg_ref
   real(dp) :: bder_ref(3)
   real(dp) :: hcov_ref(3)
   real(dp) :: hctr_ref(3)
   real(dp) :: hcurl_ref(3)
   real(dp) :: bmod_raw
   real(dp) :: sqrtg_raw
   real(dp) :: bder_raw(3)
   real(dp) :: hcov_raw(3)
   real(dp) :: hctr_raw(3)
   real(dp) :: hcurl_raw(3)
   real(dp) :: max_re_bmod
   real(dp) :: max_re_sqrtg
   real(dp) :: max_re_bder(3)
   real(dp) :: max_re_hcov(3)
   real(dp) :: max_re_hctr(3)
   real(dp) :: max_re_hcurl(3)
   integer :: i_r
   integer :: i_th
   integer :: i_ph
   integer :: n_failed

   call create_gvec_field('wout.gvec_boozer_export.nc', raw_field)
   call set_magfie_refcoords_field(raw_field)
   call init_magfie(REFCOORDS)

   phi_period = raw_field%data%phi_period()
   max_re_bmod = 0.0_dp
   max_re_sqrtg = 0.0_dp
   max_re_bder = 0.0_dp
   max_re_hcov = 0.0_dp
   max_re_hctr = 0.0_dp
   max_re_hcurl = 0.0_dp
   n_failed = 0

   do i_r = 1, n_r
      rho = 0.25_dp + 0.5_dp*real(i_r - 1, dp)/real(n_r - 1, dp)
      do i_th = 1, n_th
         theta = 0.2_dp + (2.0_dp*acos(-1.0_dp) - 0.4_dp)*real(i_th - 1, dp)/ &
                 real(n_th - 1, dp)
         do i_ph = 1, n_ph
            phi = 0.1_dp + (phi_period - 0.2_dp)*real(i_ph - 1, dp)/ &
                  real(n_ph - 1, dp)

            x_r = [rho, theta, phi]
            call magfie(x_r, bmod_ref, sqrtg_ref, bder_ref, hcov_ref, hctr_ref, hcurl_ref)
            call raw_magfie_oracle(raw_field, x_r, hs, ht, hp, phi_period, bmod_raw, &
                                   sqrtg_raw, bder_raw, hcov_raw, hctr_raw, hcurl_raw)

            call update_max_scalar(bmod_ref, bmod_raw, max_re_bmod)
            call update_max_scalar(sqrtg_ref, sqrtg_raw, max_re_sqrtg)
            call update_max_vec(bder_ref, bder_raw, max_re_bder)
            call update_max_vec(hcov_ref, hcov_raw, max_re_hcov)
            call update_max_vec(hctr_ref, hctr_raw, max_re_hctr)
            call update_max_vec(hcurl_ref, hcurl_raw, max_re_hcurl)

            if (.not. approx_rel_or_abs(bmod_ref, bmod_raw, tol_bmod, 1.0e-10_dp)) then
               call report_fail('bmod', x_r, bmod_ref, bmod_raw, tol_bmod, n_failed)
            end if
            if (.not. approx_rel_or_abs(sqrtg_ref, sqrtg_raw, tol_sqrtg, 1.0e-10_dp)) then
               call report_fail('sqrtg', x_r, sqrtg_ref, sqrtg_raw, tol_sqrtg, n_failed)
            end if
            call check_vec('bder', x_r, bder_ref, bder_raw, tol_bder, 1.0e-10_dp, n_failed)
            call check_vec('hcov', x_r, hcov_ref, hcov_raw, tol_hcov, 1.0e-10_dp, n_failed)
            call check_vec('hctr', x_r, hctr_ref, hctr_raw, tol_hctr, 1.0e-10_dp, n_failed)
            call check_vec('hcurl', x_r, hcurl_ref, hcurl_raw, tol_hcurl, 1.0e-10_dp, n_failed)
         end do
      end do
   end do

   call print_summary(max_re_bmod, max_re_sqrtg, max_re_bder, max_re_hcov, &
                      max_re_hctr, max_re_hcurl)
   if (n_failed > 0) then
      error stop 'test_gvec_boozer_refcoords_magfie: refcoords mismatch'
   end if

contains

subroutine raw_magfie_oracle(field, x_r, hs, ht, hp, phi_period, bmod, sqrtg, bder, hcov, &
                             hctr, hcurl)
   use, intrinsic :: iso_fortran_env, only: dp => real64
   use field_gvec, only: gvec_field_t

   class(gvec_field_t), intent(in) :: field
   real(dp), intent(in) :: x_r(3)
   real(dp), intent(in) :: hs
   real(dp), intent(in) :: ht
   real(dp), intent(in) :: hp
   real(dp), intent(in) :: phi_period
   real(dp), intent(out) :: bmod
   real(dp), intent(out) :: sqrtg
   real(dp), intent(out) :: bder(3)
   real(dp), intent(out) :: hcov(3)
   real(dp), intent(out) :: hctr(3)
   real(dp), intent(out) :: hcurl(3)

   real(dp) :: Acov_tmp(3)
   real(dp) :: hcov_s(3)
   real(dp) :: g_s(3, 3)
   real(dp) :: ginv_s(3, 3)
   real(dp) :: ginv_r(3, 3)
   real(dp) :: e_cov_s(3, 3)
   real(dp) :: sqrtg_abs
   real(dp) :: dh(3, 3)
   real(dp) :: dBdx(3)
   real(dp) :: J
   real(dp) :: x_plus(3)
   real(dp) :: x_minus(3)
   real(dp) :: Bp
   real(dp) :: Bm
   real(dp) :: hpv(3)
   real(dp) :: hmv(3)
   real(dp) :: sqgBctr(3)
   integer :: i

   J = 2.0_dp*x_r(1)

   call field%evaluate([x_r(1)**2, x_r(2), x_r(3)], Acov_tmp, hcov_s, bmod)
   call field%coords%metric_tensor([x_r(1)**2, x_r(2), x_r(3)], g_s, ginv_s, sqrtg_abs)
   call field%coords%covariant_basis([x_r(1)**2, x_r(2), x_r(3)], e_cov_s)
   call metric_inverse_scaled(J, ginv_s, ginv_r)

   hcov(1) = hcov_s(1)*J
   hcov(2:3) = hcov_s(2:3)
   sqrtg = signed_jacobian(e_cov_s)*J
   call field%evaluate([x_r(1)**2, x_r(2), x_r(3)], Acov_tmp, hcov_s, Bp, sqgBctr)
   hctr(1) = sqgBctr(1)/(sqrtg*bmod)
   hctr(2) = J*sqgBctr(2)/(sqrtg*bmod)
   hctr(3) = J*sqgBctr(3)/(sqrtg*bmod)

   do i = 1, 3
      x_plus = x_r
      x_minus = x_r
      select case (i)
      case (1)
         x_plus(1) = x_r(1) + hs
         x_minus(1) = max(1.0e-6_dp, x_r(1) - hs)
      case (2)
         x_plus(2) = wrap_angle(x_r(2) + ht, 2.0_dp*acos(-1.0_dp))
         x_minus(2) = wrap_angle(x_r(2) - ht, 2.0_dp*acos(-1.0_dp))
      case (3)
         x_plus(3) = wrap_angle(x_r(3) + hp, phi_period)
         x_minus(3) = wrap_angle(x_r(3) - hp, phi_period)
      end select

      call field%evaluate([x_plus(1)**2, x_plus(2), x_plus(3)], Acov_tmp, hpv, Bp)
      call field%evaluate([x_minus(1)**2, x_minus(2), x_minus(3)], Acov_tmp, hmv, Bm)
      hpv(1) = hpv(1)*2.0_dp*x_plus(1)
      hmv(1) = hmv(1)*2.0_dp*x_minus(1)

      select case (i)
      case (1)
         dBdx(i) = (Bp - Bm)/(x_plus(1) - x_minus(1))
         dh(i, :) = (hpv - hmv)/(x_plus(1) - x_minus(1))
      case (2)
         dBdx(i) = (Bp - Bm)/(2.0_dp*ht)
         dh(i, :) = (hpv - hmv)/(2.0_dp*ht)
      case (3)
         dBdx(i) = (Bp - Bm)/(2.0_dp*hp)
         dh(i, :) = (hpv - hmv)/(2.0_dp*hp)
      end select
   end do

   bder = dBdx/max(bmod, 1.0e-30_dp)
   call compute_hcurl(sqrtg, dh, hcurl)
end subroutine raw_magfie_oracle

subroutine metric_inverse_scaled(J, ginv_s, ginv_r)
   use, intrinsic :: iso_fortran_env, only: dp => real64

   real(dp), intent(in) :: J
   real(dp), intent(in) :: ginv_s(3, 3)
   real(dp), intent(out) :: ginv_r(3, 3)

   ginv_r = ginv_s
   ginv_r(1, 1) = ginv_s(1, 1)/(J*J)
   ginv_r(1, 2) = ginv_s(1, 2)/J
   ginv_r(1, 3) = ginv_s(1, 3)/J
   ginv_r(2, 1) = ginv_s(2, 1)/J
   ginv_r(3, 1) = ginv_s(3, 1)/J
end subroutine metric_inverse_scaled

subroutine compute_hcurl(sqrtg, dh, hcurl)
   use, intrinsic :: iso_fortran_env, only: dp => real64

   real(dp), intent(in) :: sqrtg
   real(dp), intent(in) :: dh(3, 3)
   real(dp), intent(out) :: hcurl(3)

   hcurl(1) = (dh(2, 3) - dh(3, 2))/sqrtg
   hcurl(2) = (dh(3, 1) - dh(1, 3))/sqrtg
   hcurl(3) = (dh(1, 2) - dh(2, 1))/sqrtg
end subroutine compute_hcurl

pure real(dp) function wrap_angle(angle, period) result(wrapped)
   use, intrinsic :: iso_fortran_env, only: dp => real64

   real(dp), intent(in) :: angle
   real(dp), intent(in) :: period

   wrapped = modulo(angle, period)
   if (wrapped < 0.0_dp) wrapped = wrapped + period
end function wrap_angle

pure real(dp) function signed_jacobian(e_cov) result(jac)
   use, intrinsic :: iso_fortran_env, only: dp => real64

   real(dp), intent(in) :: e_cov(3, 3)
   real(dp) :: c(3)

   c = cross_product(e_cov(:, 2), e_cov(:, 3))
   jac = dot_product(e_cov(:, 1), c)
end function signed_jacobian

pure function cross_product(a, b) result(c)
   use, intrinsic :: iso_fortran_env, only: dp => real64

   real(dp), intent(in) :: a(3)
   real(dp), intent(in) :: b(3)
   real(dp) :: c(3)

   c(1) = a(2)*b(3) - a(3)*b(2)
   c(2) = a(3)*b(1) - a(1)*b(3)
   c(3) = a(1)*b(2) - a(2)*b(1)
end function cross_product

pure logical function approx_rel_or_abs(a, b, tol_rel, tol_abs) result(ok)
   use, intrinsic :: iso_fortran_env, only: dp => real64

   real(dp), intent(in) :: a
   real(dp), intent(in) :: b
   real(dp), intent(in) :: tol_rel
   real(dp), intent(in) :: tol_abs

   if (abs(b) <= tol_abs) then
      ok = abs(a - b) <= tol_abs
   else
      ok = abs(a - b)/abs(b) <= tol_rel
   end if
end function approx_rel_or_abs

subroutine check_vec(name, x, a, b, tol_rel, tol_abs, n_failed)
   use, intrinsic :: iso_fortran_env, only: dp => real64

   character(*), intent(in) :: name
   real(dp), intent(in) :: x(3)
   real(dp), intent(in) :: a(3)
   real(dp), intent(in) :: b(3)
   real(dp), intent(in) :: tol_rel
   real(dp), intent(in) :: tol_abs
   integer, intent(inout) :: n_failed
   integer :: i

   do i = 1, 3
      if (.not. approx_rel_or_abs(a(i), b(i), tol_rel, tol_abs)) then
         call report_fail(trim(name)//'('//trim(int_to_str(i))//')', x, a(i), b(i), &
                          tol_rel, n_failed)
      end if
   end do
end subroutine check_vec

subroutine report_fail(name, x, a, b, tol, n_failed)
   use, intrinsic :: iso_fortran_env, only: dp => real64

   character(*), intent(in) :: name
   real(dp), intent(in) :: x(3)
   real(dp), intent(in) :: a
   real(dp), intent(in) :: b
   real(dp), intent(in) :: tol
   integer, intent(inout) :: n_failed
   real(dp) :: rel

   rel = abs(a - b)/max(abs(b), 1.0e-30_dp)
   print *, 'FAIL: ', trim(name), ' at x=(rho,th,ph)=', x
   print *, '  got=', a, ' expected=', b, ' rel=', rel, ' tol=', tol
   n_failed = n_failed + 1
end subroutine report_fail

pure function int_to_str(i) result(s)
   integer, intent(in) :: i
   character(len=16) :: s

   write (s, '(i0)') i
end function int_to_str

subroutine update_max_scalar(a, b, max_rel)
   use, intrinsic :: iso_fortran_env, only: dp => real64

   real(dp), intent(in) :: a
   real(dp), intent(in) :: b
   real(dp), intent(inout) :: max_rel
   real(dp) :: rel

   rel = abs(a - b)/max(abs(b), 1.0e-14_dp)
   max_rel = max(max_rel, rel)
end subroutine update_max_scalar

subroutine update_max_vec(a, b, max_rel)
   use, intrinsic :: iso_fortran_env, only: dp => real64

   real(dp), intent(in) :: a(3)
   real(dp), intent(in) :: b(3)
   real(dp), intent(inout) :: max_rel(3)
   integer :: i

   do i = 1, 3
      call update_max_scalar(a(i), b(i), max_rel(i))
   end do
end subroutine update_max_vec

subroutine print_summary(max_re_bmod, max_re_sqrtg, max_re_bder, max_re_hcov, max_re_hctr, &
                         max_re_hcurl)
   use, intrinsic :: iso_fortran_env, only: dp => real64

   real(dp), intent(in) :: max_re_bmod
   real(dp), intent(in) :: max_re_sqrtg
   real(dp), intent(in) :: max_re_bder(3)
   real(dp), intent(in) :: max_re_hcov(3)
   real(dp), intent(in) :: max_re_hctr(3)
   real(dp), intent(in) :: max_re_hcurl(3)

   print *, 'Max relative errors for GVEC Boozer refcoords magfie:'
   print '(a,1x,es12.4)', '  bmod  :', max_re_bmod
   print '(a,1x,es12.4)', '  sqrtg :', max_re_sqrtg
   print '(a,3(1x,es12.4))', '  bder  :', max_re_bder
   print '(a,3(1x,es12.4))', '  hcov  :', max_re_hcov
   print '(a,3(1x,es12.4))', '  hctr  :', max_re_hctr
   print '(a,3(1x,es12.4))', '  hcurl :', max_re_hcurl
end subroutine print_summary

end program test_gvec_boozer_refcoords_magfie
