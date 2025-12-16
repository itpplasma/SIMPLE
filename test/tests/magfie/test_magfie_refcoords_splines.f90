program test_magfie_refcoords_splines
    !> Verify magfie_refcoords assembly against finite-difference oracles.
    !>
    !> Oracle sources:
    !> - Bmod, hcov: splined_field_t%evaluate
    !> - bder, dhcov: 5-point central finite differences of evaluate()
    !> - sqrtg, ginv: finite differences of coordinate_system_t%evaluate_cart
    !>
    !> Target:
    !> - magfie_refcoords (via init_magfie(REFCOORDS))

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use simple, only: init_vmec
    use util, only: twopi
    use new_vmec_stuff_mod, only: nper
    use field_coils, only: coils_field_t, create_coils_field
    use field_splined, only: splined_field_t, create_splined_field
    use reference_coordinates, only: init_reference_coordinates, ref_coords
    use magfie_sub, only: init_magfie, magfie, REFCOORDS, set_magfie_refcoords_field
    use libneo_coordinates, only: coordinate_system_t

    implicit none

    type(coils_field_t) :: raw_coils
    type(splined_field_t) :: f
    real(dp) :: dummy

    integer, parameter :: n_points = 240
    real(dp), parameter :: h_rho = 1.0e-4_dp
    real(dp), parameter :: h_ang = 1.0e-5_dp

    real(dp), parameter :: tol_bmod_rel = 1.0e-12_dp
    real(dp), parameter :: tol_bmod_abs = 1.0e-12_dp
    real(dp), parameter :: tol_sqrtg_rel = 5.0e-7_dp
    real(dp), parameter :: tol_sqrtg_abs = 1.0e-10_dp
    real(dp), parameter :: tol_bder_rel = 5.0e-8_dp
    real(dp), parameter :: tol_bder_abs = 1.0e-10_dp
    real(dp), parameter :: tol_hcov_rel = 1.0e-12_dp
    real(dp), parameter :: tol_hcov_abs = 1.0e-12_dp
    real(dp), parameter :: tol_hctr_rel = 1.0e-8_dp
    real(dp), parameter :: tol_hctr_abs = 5.0e-11_dp
    real(dp), parameter :: tol_hcurl_rel = 1.0e-6_dp
    real(dp), parameter :: tol_hcurl_abs = 1.0e-12_dp

    real(dp) :: phi_period
    integer :: n_failed, i

    n_failed = 0

    call init_vmec('wout_ncsx.nc', 5, 5, 5, dummy)
    call init_reference_coordinates('wout_ncsx.nc')

    call create_coils_field('coils.simple', raw_coils)
    call create_splined_field(raw_coils, ref_coords, f)

    call set_magfie_refcoords_field(f)
    call init_magfie(REFCOORDS)

    phi_period = twopi/real(nper, dp)

    do i = 1, n_points
        call check_point(i, phi_period, f, n_failed)
    end do

    if (n_failed > 0) then
        print *, '================================'
        print *, 'FAILED: ', n_failed, ' mismatch(es) in magfie_refcoords test'
        print *, '================================'
        error stop 1
    end if

    print *, '================================'
    print *, 'PASSED: magfie_refcoords matches FD oracles'
    print *, '================================'

contains

    subroutine check_point(i, phi_period, f, n_failed)
        integer, intent(in) :: i
        real(dp), intent(in) :: phi_period
        type(splined_field_t), intent(in) :: f
        integer, intent(inout) :: n_failed

        real(dp) :: x(3)
        real(dp) :: bmod_t, sqrtg_t, bder_t(3), hcov_t(3), hctr_t(3), hcurl_t(3)
        real(dp) :: B0, h0(3)
        real(dp) :: bder_fd(3), hcov_fd(3), hctr_fd(3), hcurl_fd(3)
        real(dp) :: dhcov_fd(3, 3)
        real(dp) :: ginv_fd(3, 3), sqrtg_fd

        x(1) = 0.25_dp + 0.55_dp*real(mod(13*i, 997), dp)/997.0_dp
        x(2) = 0.1_dp + (twopi - 0.2_dp)*real(mod(29*i, 991), dp)/991.0_dp
        x(3) = 0.1_dp + (phi_period - 0.2_dp)*real(mod(31*i, 983), dp)/983.0_dp

        call oracle_field_fd(f, x, phi_period, B0, h0, bder_fd, dhcov_fd)
        hcov_fd = h0

        call metric_from_cart_fd(ref_coords, x, phi_period, ginv_fd, sqrtg_fd)
        hctr_fd = matmul(ginv_fd, hcov_fd)
        call compute_hcurl(sqrtg_fd, dhcov_fd, hcurl_fd)

        call magfie(x, bmod_t, sqrtg_t, bder_t, hcov_t, hctr_t, hcurl_t)

        call check_scalar('bmod', x, bmod_t, B0, tol_bmod_abs, tol_bmod_rel, n_failed)
        call check_scalar('sqrtg', x, sqrtg_t, sqrtg_fd, tol_sqrtg_abs, &
                          tol_sqrtg_rel, n_failed)
        call check_vec('bder', x, bder_t, bder_fd, tol_bder_abs, tol_bder_rel, n_failed)
        call check_vec('hcov', x, hcov_t, hcov_fd, tol_hcov_abs, tol_hcov_rel, n_failed)
        call check_vec('hctr', x, hctr_t, hctr_fd, tol_hctr_abs, tol_hctr_rel, n_failed)
        call check_vec('hcurl', x, hcurl_t, hcurl_fd, tol_hcurl_abs, tol_hcurl_rel, &
                       n_failed)
    end subroutine check_point


    subroutine oracle_field_fd(f, x, phi_period, B0, h0, bder, dh)
        type(splined_field_t), intent(in) :: f
        real(dp), intent(in) :: x(3), phi_period
        real(dp), intent(out) :: B0, h0(3)
        real(dp), intent(out) :: bder(3), dh(3, 3)

        real(dp) :: Acov0(3), Bp1, Bm1, Bp2, Bm2
        real(dp) :: hp1(3), hm1(3), hp2(3), hm2(3)
        integer :: k
        real(dp) :: step, dBdx

        call f%evaluate(x, Acov0, h0, B0)

        do k = 1, 3
            call eval_shift(f, x, k, +1, phi_period, Bp1, hp1)
            call eval_shift(f, x, k, -1, phi_period, Bm1, hm1)
            call eval_shift(f, x, k, +2, phi_period, Bp2, hp2)
            call eval_shift(f, x, k, -2, phi_period, Bm2, hm2)

            if (k == 1) then
                step = h_rho
            else
                step = h_ang
            end if

            dBdx = (-Bp2 + 8.0_dp*Bp1 - 8.0_dp*Bm1 + Bm2)/(12.0_dp*step)
            bder(k) = dBdx/max(B0, 1.0e-30_dp)
            dh(k, :) = (-hp2 + 8.0_dp*hp1 - 8.0_dp*hm1 + hm2)/(12.0_dp*step)
        end do
    end subroutine oracle_field_fd


    subroutine eval_shift(f, x, k, m, phi_period, Bmod, hcov)
        type(splined_field_t), intent(in) :: f
        real(dp), intent(in) :: x(3)
        integer, intent(in) :: k, m
        real(dp), intent(in) :: phi_period
        real(dp), intent(out) :: Bmod, hcov(3)

        real(dp) :: x_shift(3)
        real(dp) :: Acov(3)

        x_shift = x
        select case (k)
        case (1)
            x_shift(1) = x(1) + real(m, dp)*h_rho
        case (2)
            x_shift(2) = wrap_angle(x(2) + real(m, dp)*h_ang, twopi)
        case (3)
            x_shift(3) = wrap_angle(x(3) + real(m, dp)*h_ang, phi_period)
        case default
            error stop 'eval_shift: invalid coordinate index'
        end select

        call f%evaluate(x_shift, Acov, hcov, Bmod)
    end subroutine eval_shift


    subroutine metric_from_cart_fd(cs, x, phi_period, ginv, sqrtg)
        class(*), intent(in) :: cs
        real(dp), intent(in) :: x(3), phi_period
        real(dp), intent(out) :: ginv(3, 3), sqrtg

        real(dp) :: u0(3), u_p1(3), u_m1(3), u_p2(3), u_m2(3)
        real(dp) :: xc_p1(3), xc_m1(3), xc_p2(3), xc_m2(3)
        real(dp) :: e1(3), e2(3), e3(3)
        real(dp) :: g(3, 3), det

        u0 = [x(1)**2, x(2), x(3)]

        u_p1 = [ (x(1) + h_rho)**2, x(2), x(3) ]
        u_m1 = [ (x(1) - h_rho)**2, x(2), x(3) ]
        u_p2 = [ (x(1) + 2.0_dp*h_rho)**2, x(2), x(3) ]
        u_m2 = [ (x(1) - 2.0_dp*h_rho)**2, x(2), x(3) ]
        call eval_cart(cs, u_p1, xc_p1)
        call eval_cart(cs, u_m1, xc_m1)
        call eval_cart(cs, u_p2, xc_p2)
        call eval_cart(cs, u_m2, xc_m2)
        e1 = (-xc_p2 + 8.0_dp*xc_p1 - 8.0_dp*xc_m1 + xc_m2)/(12.0_dp*h_rho)

        u_p1 = [ x(1)**2, wrap_angle(x(2) + h_ang, twopi), x(3) ]
        u_m1 = [ x(1)**2, wrap_angle(x(2) - h_ang, twopi), x(3) ]
        u_p2 = [ x(1)**2, wrap_angle(x(2) + 2.0_dp*h_ang, twopi), x(3) ]
        u_m2 = [ x(1)**2, wrap_angle(x(2) - 2.0_dp*h_ang, twopi), x(3) ]
        call eval_cart(cs, u_p1, xc_p1)
        call eval_cart(cs, u_m1, xc_m1)
        call eval_cart(cs, u_p2, xc_p2)
        call eval_cart(cs, u_m2, xc_m2)
        e2 = (-xc_p2 + 8.0_dp*xc_p1 - 8.0_dp*xc_m1 + xc_m2)/(12.0_dp*h_ang)

        u_p1 = [ x(1)**2, x(2), wrap_angle(x(3) + h_ang, phi_period) ]
        u_m1 = [ x(1)**2, x(2), wrap_angle(x(3) - h_ang, phi_period) ]
        u_p2 = [ x(1)**2, x(2), wrap_angle(x(3) + 2.0_dp*h_ang, phi_period) ]
        u_m2 = [ x(1)**2, x(2), wrap_angle(x(3) - 2.0_dp*h_ang, phi_period) ]
        call eval_cart(cs, u_p1, xc_p1)
        call eval_cart(cs, u_m1, xc_m1)
        call eval_cart(cs, u_p2, xc_p2)
        call eval_cart(cs, u_m2, xc_m2)
        e3 = (-xc_p2 + 8.0_dp*xc_p1 - 8.0_dp*xc_m1 + xc_m2)/(12.0_dp*h_ang)

        sqrtg = dot_product(e1, cross_product(e2, e3))

        g(1, 1) = dot_product(e1, e1)
        g(1, 2) = dot_product(e1, e2)
        g(1, 3) = dot_product(e1, e3)
        g(2, 1) = g(1, 2)
        g(2, 2) = dot_product(e2, e2)
        g(2, 3) = dot_product(e2, e3)
        g(3, 1) = g(1, 3)
        g(3, 2) = g(2, 3)
        g(3, 3) = dot_product(e3, e3)

        det = g(1,1)*(g(2,2)*g(3,3) - g(2,3)*g(3,2)) &
            - g(1,2)*(g(2,1)*g(3,3) - g(2,3)*g(3,1)) &
            + g(1,3)*(g(2,1)*g(3,2) - g(2,2)*g(3,1))

        ginv(1,1) = (g(2,2)*g(3,3) - g(2,3)*g(3,2))/det
        ginv(1,2) = (g(1,3)*g(3,2) - g(1,2)*g(3,3))/det
        ginv(1,3) = (g(1,2)*g(2,3) - g(1,3)*g(2,2))/det
        ginv(2,1) = (g(2,3)*g(3,1) - g(2,1)*g(3,3))/det
        ginv(2,2) = (g(1,1)*g(3,3) - g(1,3)*g(3,1))/det
        ginv(2,3) = (g(1,3)*g(2,1) - g(1,1)*g(2,3))/det
        ginv(3,1) = (g(2,1)*g(3,2) - g(2,2)*g(3,1))/det
        ginv(3,2) = (g(1,2)*g(3,1) - g(1,1)*g(3,2))/det
        ginv(3,3) = (g(1,1)*g(2,2) - g(1,2)*g(2,1))/det
    end subroutine metric_from_cart_fd


    subroutine eval_cart(cs, u, xcart)
        class(*), intent(in) :: cs
        real(dp), intent(in) :: u(3)
        real(dp), intent(out) :: xcart(3)

        select type (cs)
        class is (coordinate_system_t)
            call cs%evaluate_cart(u, xcart)
        class default
            error stop 'eval_cart: invalid coordinate_system_t'
        end select
    end subroutine eval_cart


    subroutine compute_hcurl(sqrtg, dh, hcurl)
        real(dp), intent(in) :: sqrtg
        real(dp), intent(in) :: dh(3, 3)
        real(dp), intent(out) :: hcurl(3)

        if (abs(sqrtg) <= 0.0_dp) error stop 'compute_hcurl: sqrtg must be nonzero'

        hcurl(1) = (dh(2, 3) - dh(3, 2))/sqrtg
        hcurl(2) = (dh(3, 1) - dh(1, 3))/sqrtg
        hcurl(3) = (dh(1, 2) - dh(2, 1))/sqrtg
    end subroutine compute_hcurl


    subroutine check_scalar(name, x, a, b, tol_abs, tol_rel, n_failed)
        character(*), intent(in) :: name
        real(dp), intent(in) :: x(3), a, b, tol_abs, tol_rel
        integer, intent(inout) :: n_failed
        real(dp) :: abs_err

        abs_err = abs(a - b)
        if (abs_err > tol_abs + tol_rel*abs(b)) then
            print *, 'FAIL: ', trim(name), ' at x=(rho,th,ph)=', x
            print *, '  got=', a, ' expected=', b, ' abs=', abs_err
            n_failed = n_failed + 1
        end if
    end subroutine check_scalar


    subroutine check_vec(name, x, a, b, tol_abs, tol_rel, n_failed)
        character(*), intent(in) :: name
        real(dp), intent(in) :: x(3), a(3), b(3), tol_abs, tol_rel
        integer, intent(inout) :: n_failed
        integer :: k
        real(dp) :: abs_err

        do k = 1, 3
            abs_err = abs(a(k) - b(k))
            if (abs_err > tol_abs + tol_rel*abs(b(k))) then
                print *, 'FAIL: ', trim(name), '(', k, ') at x=(rho,th,ph)=', x
                print *, '  got=', a(k), ' expected=', b(k), ' abs=', abs_err
                n_failed = n_failed + 1
            end if
        end do
    end subroutine check_vec


    pure function cross_product(a, b) result(c)
        real(dp), intent(in) :: a(3), b(3)
        real(dp) :: c(3)

        c(1) = a(2)*b(3) - a(3)*b(2)
        c(2) = a(3)*b(1) - a(1)*b(3)
        c(3) = a(1)*b(2) - a(2)*b(1)
    end function cross_product


    pure real(dp) function wrap_angle(angle, period) result(wrapped)
        real(dp), intent(in) :: angle, period

        wrapped = modulo(angle, period)
        if (wrapped < 0.0_dp) wrapped = wrapped + period
    end function wrap_angle

end program test_magfie_refcoords_splines
