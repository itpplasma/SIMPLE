program test_bder_from_splines
    !> Verify bder = d(log Bmod)/dx computed from spline derivatives.
    !> Reference is a high-order finite-difference oracle of Bmod from evaluate().

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use simple, only: init_vmec
    use util, only: twopi
    use new_vmec_stuff_mod, only: nper
    use field_coils, only: coils_field_t, create_coils_field
    use field_splined, only: splined_field_t, create_splined_field
    use reference_coordinates, only: init_reference_coordinates, ref_coords

    implicit none

    type(coils_field_t) :: raw_coils
    type(splined_field_t) :: splined

    real(dp) :: dummy
    real(dp) :: phi_period
    integer :: n_failed

    integer, parameter :: n_points = 240
    real(dp), parameter :: h_rho = 1.0e-4_dp
    real(dp), parameter :: h_ang = 3.0e-5_dp

    real(dp) :: max_ratio(3)
    integer :: i

    n_failed = 0
    max_ratio = 0.0_dp

    ! Derivative consistency needs only a coarse VMEC coordinate fixture.
    ! Keep quintic splines; multiplier 1 bounds the coordinate arrays to 0.43 GiB.
    call init_vmec('wout_ncsx.nc', 5, 5, 1, dummy)
    call init_reference_coordinates('wout_ncsx.nc')

    call create_coils_field('coils.simple', raw_coils)
    call create_splined_field(raw_coils, ref_coords, splined, n_r=32, n_th=33, n_phi=32)

    phi_period = twopi/real(nper, dp)

    do i = 1, n_points
        call check_point(i, phi_period, splined, max_ratio, n_failed)
    end do

    if (n_failed > 0) then
        call print_summary(max_ratio)
        print *, '================================'
        print *, 'FAILED: ', n_failed, ' bder mismatch(es)'
        print *, '================================'
        error stop 1
    end if

    call print_summary(max_ratio)
    print *, '================================'
    print *, 'PASSED: bder from splines matches FD oracle'
    print *, '================================'

contains

    subroutine check_point(i, phi_period, f, max_ratio, n_failed)
        integer, intent(in) :: i
        real(dp), intent(in) :: phi_period
        type(splined_field_t), intent(in) :: f
        real(dp), intent(inout) :: max_ratio(3)
        integer, intent(inout) :: n_failed

        real(dp) :: x(3)
        real(dp) :: Acov(3), hcov(3), Bmod
        real(dp) :: dAcov(3, 3), dhcov(3, 3), dBmod(3)
        real(dp) :: bder_spline(3), bder_fd(3), step, error_bound
        integer :: k

        x(1) = 0.25_dp + 0.55_dp*real(mod(13*i, 997), dp)/997.0_dp
        x(2) = 0.1_dp + (twopi - 0.2_dp)*real(mod(29*i, 991), dp)/991.0_dp
        x(3) = 0.1_dp + (phi_period - 0.2_dp)*real(mod(31*i, 983), dp)/983.0_dp

        call f%evaluate_with_der(x, Acov, hcov, Bmod, dAcov, dhcov, dBmod)
        bder_spline = dBmod/max(Bmod, 1.0e-30_dp)

        do k = 1, 3
            step = fd_step(f, x, k)
            bder_fd(k) = dlogB_dx_fd5(f, x, k, step, phi_period)/ &
                         max(Bmod, 1.0e-30_dp)
            error_bound = fd_error_bound(f, x, k, step)/max(Bmod, 1.0e-30_dp) + &
                          gamma(2)*(abs(bder_spline(k)) + abs(bder_fd(k)))
            call update_max_and_check(k, x, bder_spline(k), bder_fd(k), &
                                      error_bound, max_ratio, n_failed)
        end do
    end subroutine check_point


    real(dp) function dlogB_dx_fd5(f, x, k, step, phi_period) result(dBdx)
        type(splined_field_t), intent(in) :: f
        real(dp), intent(in) :: x(3)
        integer, intent(in) :: k
        real(dp), intent(in) :: step, phi_period

        real(dp) :: Bp1, Bm1, Bp2, Bm2

        call eval_B_only(f, x, k, +1, step, phi_period, Bp1)
        call eval_B_only(f, x, k, -1, step, phi_period, Bm1)
        call eval_B_only(f, x, k, +2, step, phi_period, Bp2)
        call eval_B_only(f, x, k, -2, step, phi_period, Bm2)

        dBdx = (-Bp2 + 8.0_dp*Bp1 - 8.0_dp*Bm1 + Bm2)/(12.0_dp*step)
    end function dlogB_dx_fd5


    subroutine eval_B_only(f, x, k, m, step, phi_period, Bmod)
        type(splined_field_t), intent(in) :: f
        real(dp), intent(in) :: x(3)
        integer, intent(in) :: k, m
        real(dp), intent(in) :: step, phi_period
        real(dp), intent(out) :: Bmod

        real(dp) :: x_shift(3)
        real(dp) :: Acov(3), hcov(3)

        x_shift = x
        select case (k)
        case (1)
            x_shift(1) = x(1) + real(m, dp)*step
        case (2)
            x_shift(2) = wrap_angle(x(2) + real(m, dp)*step, twopi)
        case (3)
            x_shift(3) = wrap_angle(x(3) + real(m, dp)*step, phi_period)
        case default
            error stop 'eval_B_only: invalid coordinate index'
        end select

        call f%evaluate(x_shift, Acov, hcov, Bmod)
    end subroutine eval_B_only

    real(dp) function fd_step(f, x, k) result(step)
        type(splined_field_t), intent(in) :: f
        real(dp), intent(in) :: x(3)
        integer, intent(in) :: k
        real(dp) :: local, width

        ! Keep +/-2h in one polynomial piece. Stored coefficients satisfy C^4
        ! knot continuity only up to roundoff; no knot-jump estimate is needed.
        width = f%splines%h_step(k)
        local = modulo(x(k) - f%splines%x_min(k), width)
        step = min(merge(h_rho, h_ang, k == 1), &
                   0.25_dp*min(local, width - local))
        if (step <= 0.0_dp) error stop 'FD point lies on a spline knot'
    end function fd_step

    real(dp) function fd_error_bound(f, x, k, step) result(bound)
        type(splined_field_t), intent(in) :: f
        real(dp), intent(in) :: x(3), step
        integer, intent(in) :: k

        integer :: lo(3), hi(3), cell(3), orders(3), j
        real(dp) :: extent(3), delta(3), radius(3)
        real(dp) :: value_bound, gradient(3), hessian_row(3), fifth_bound

        ! This tests two evaluations of the SAME tensor-product quintic, not
        ! interpolation accuracy against the coil field. Thus there is no O(H^5)
        ! spline derivative error here (H = [1/31, 2*pi/32, period/31]).
        ! Within one cell the five-point stencil satisfies exactly
        ! D_h B - B' = -h^4 B^(5)/30: Taylor terms through degree four
        ! cancel, degree five gives -4*h^4*c_5, and B^(5) = 120*c_5.
        ! There are no higher terms. See https://dlmf.nist.gov/3.4#i.
        !
        ! For roundoff, u = epsilon/2 and gamma(n) = n*u/(1-n*u).
        ! A coefficient path in 3D quintic Horner has at most 3*5*2 = 30
        ! operations; 40 includes stencil arithmetic, 32 derivative Horner.
        ! The stencil's absolute weight sum is 18/(12*h) = 3/(2*h).
        ! Coordinate formation/normalization takes at most 8 operations; its
        ! absolute error delta propagates through gradient/Hessian bounds.
        ! Use absolute-coefficient polynomial envelopes throughout the cell,
        ! so cancellation in B or B' cannot hide error.
        ! Both bder values share the same B divisor; its two division errors
        ! are accounted for by the caller. No empirical relative floor is used.
        extent = 0.0_dp
        extent(k) = 2.0_dp*step
        delta = gamma(8)*(abs(x) + abs(f%splines%x_min) + &
                          f%splines%h_step + extent)
        lo = floor((x - extent - delta - f%splines%x_min)/f%splines%h_step) + 1
        hi = floor((x + extent + delta - f%splines%x_min)/f%splines%h_step) + 1
        if (any(lo < 1) .or. any(hi >= f%splines%num_points) .or. &
            any(lo /= hi)) error stop 'FD bound requires one interior cell'
        if (any(f%splines%order /= 5)) error stop 'FD bound requires quintics'
        radius = f%splines%h_step + delta
        cell = lo
        orders = 0
        value_bound = envelope(f, cell, radius, orders)
        do j = 1, 3
            orders = 0
            orders(j) = 1
            gradient(j) = envelope(f, cell, radius, orders)
            orders(k) = orders(k) + 1
            hessian_row(j) = envelope(f, cell, radius, orders)
        end do
        orders = 0
        orders(k) = 5
        fifth_bound = envelope(f, cell, radius, orders)
        bound = step**4*fifth_bound/30.0_dp + &
                (1.5_dp/step)*(gamma(40)*value_bound + dot_product(delta, gradient)) + &
                gamma(32)*gradient(k) + dot_product(delta, hessian_row)
    end function fd_error_bound

    real(dp) function envelope(f, cell, radius, orders) result(bound)
        type(splined_field_t), intent(in) :: f
        integer, intent(in) :: cell(3), orders(3)
        real(dp), intent(in) :: radius(3)
        integer, parameter :: factorial(0:5) = [1, 1, 2, 6, 24, 120]
        real(dp) :: weights(0:5, 3)
        integer :: j, p, p1, p2, p3

        ! Sum |c_p| (p)_orders radius**(p-orders), an envelope for the
        ! requested partial derivative throughout a cell. Coefficients are in
        ! physical local coordinates, not normalized cell coordinates.
        weights = 0.0_dp
        do j = 1, 3
            do p = orders(j), 5
                weights(p, j) = real(factorial(p)/factorial(p - orders(j)), dp)* &
                                radius(j)**(p - orders(j))
            end do
        end do
        bound = 0.0_dp
        do p3 = orders(3), 5
            do p2 = orders(2), 5
                do p1 = orders(1), 5
                    bound = bound + &
                            abs(f%splines%coeff(7, p1, p2, p3, cell(1), &
                                                cell(2), cell(3)))* &
                            weights(p1, 1)*weights(p2, 2)*weights(p3, 3)
                end do
            end do
        end do
        ! Outward allowance for evaluating this positive sum (216 terms,
        ! small integer powers and products): fewer than 1024 operations.
        bound = bound/(1.0_dp - gamma(1024))
    end function envelope

    pure real(dp) function gamma(n) result(bound)
        integer, intent(in) :: n
        real(dp), parameter :: u = epsilon(1.0_dp)/2.0_dp

        bound = real(n, dp)*u/(1.0_dp - real(n, dp)*u)
    end function gamma

    subroutine update_max_and_check(k, x, a, b, bound, max_ratio, n_failed)
        integer, intent(in) :: k
        real(dp), intent(in) :: x(3), a, b, bound
        real(dp), intent(inout) :: max_ratio(3)
        integer, intent(inout) :: n_failed

        real(dp) :: error

        error = abs(a - b)
        max_ratio(k) = max(max_ratio(k), error/bound)

        if (.not. (error <= bound)) then
            print *, 'FAIL: bder(', k, ') at x=(rho,th,ph)=', x
            print *, '  got=', a, ' expected=', b, ' error=', error, ' bound=', bound
            n_failed = n_failed + 1
        end if
    end subroutine update_max_and_check


    subroutine print_summary(max_ratio)
        real(dp), intent(in) :: max_ratio(3)

        print *, 'Max fractions of the FD error bound (must be <= 1):'
        print '(a,3(1x,es12.4))', '  bder:', max_ratio
    end subroutine print_summary


    pure real(dp) function wrap_angle(angle, period) result(wrapped)
        real(dp), intent(in) :: angle, period

        wrapped = modulo(angle, period)
        if (wrapped < 0.0_dp) wrapped = wrapped + period
    end function wrap_angle

end program test_bder_from_splines
