program test_boozer_regular_axis_field
    !> The regular near-axis Boozer field (libneo boozer_regular) is built only
    !> on request and, away from the axis, agrees with the polar Boozer splines
    !> that field_can evaluates in (s, theta_B, phi_B).
    !>
    !> Usage: test_boozer_regular_axis_field.x vmec|chartmap
    !>   vmec      internal VMEC -> Boozer conversion of wout.nc
    !>   chartmap  Boozer chartmap exported from the same conversion of wout.nc
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use simple, only: init_vmec
    use magfie_sub, only: BOOZER
    use field_can_base, only: field_can_t, n_field_evaluations
    use field_can_mod, only: init_field_can, set_regular_axis_field, &
                             eval_field => evaluate
    use field_can_boozer, only: evaluate_boozer_regular
    use field_boozer_chartmap, only: boozer_chartmap_field_t, &
                                     create_boozer_chartmap_field
    use boozer_chartmap, only: export_boozer_chartmap
    use boozer_regular, only: regular_field_t, boozer_regular_field
    use new_vmec_stuff_mod, only: nper

    implicit none

    real(dp), parameter :: twopi = 6.283185307179586476925286766559_dp
    !> Spline representations differ by interpolation error: O(1e-6) in the
    !> values and O(1e-3) in the radial derivatives on the low-resolution QA
    !> equilibrium (libneo#437).
    real(dp), parameter :: tol_value = 2.0e-5_dp
    real(dp), parameter :: tol_derivative = 2.0e-2_dp

    character(len=16) :: source
    type(boozer_chartmap_field_t), allocatable :: chartmap
    real(dp) :: fper
    integer :: nfail

    nfail = 0
    call get_command_argument(1, source)

    select case (trim(source))
    case ('vmec')
        call init_vmec('wout.nc', 5, 5, 5, fper)
        call set_regular_axis_field(.false.)
        call init_field_can(BOOZER)
        call check(.not. boozer_regular_field%ready, &
                   'regular field is not built by default')
        call set_regular_axis_field(.true.)
        call init_field_can(BOOZER)
    case ('chartmap')
        ! The analytic fixture test_boozer_chartmap.nc has B_theta constant,
        ! nonzero at the axis, which no regular field represents. Export the
        ! physical QA field instead.
        call init_vmec('wout.nc', 5, 5, 5, fper)
        call init_field_can(BOOZER)
        call export_boozer_chartmap('wout_regular_axis.chartmap.nc')
        call create_boozer_chartmap_field('wout_regular_axis.chartmap.nc', chartmap)
        call set_regular_axis_field(.false.)
        call init_field_can(BOOZER, chartmap)
        call check(.not. boozer_regular_field%ready, &
                   'regular field is not built by default')
        call set_regular_axis_field(.true.)
        call init_field_can(BOOZER, chartmap)
    case default
        error stop 'usage: test_boozer_regular_axis_field.x vmec|chartmap'
    end select
    call set_regular_axis_field(.false.)

    call check(boozer_regular_field%ready, 'regular field built on request')
    call compare_with_polar(twopi/real(nper, dp))
    call check_counter()

    if (nfail > 0) then
        print '(a,i0,a)', 'FAILED: ', nfail, ' check(s)'
        error stop 1
    end if
    print '(a)', 'PASSED: regular Boozer field wiring ('//trim(source)//')'

contains

    subroutine check(condition, label)
        logical, intent(in) :: condition
        character(*), intent(in) :: label

        if (condition) then
            print '(a)', 'PASS: '//label
        else
            print '(a)', 'FAIL: '//label
            nfail = nfail + 1
        end if
    end subroutine check

    subroutine compare_with_polar(phi_period)
        real(dp), intent(in) :: phi_period

        integer, parameter :: n_s = 4, n_th = 5, n_ph = 3
        real(dp), parameter :: s_values(n_s) = [0.05_dp, 0.2_dp, 0.5_dp, 0.8_dp]
        type(field_can_t) :: polar
        type(regular_field_t) :: reg
        real(dp) :: s, th, ph, x, y, err_value, err_derivative
        real(dp) :: reg_value(5), polar_value(5), reg_ds(5), polar_ds(5)
        real(dp) :: reg_dth(5), polar_dth(5), reg_dph(5), polar_dph(5)
        real(dp) :: dscale(5), vscale(5)
        integer :: i_s, i_th, i_ph

        err_value = 0.0_dp
        err_derivative = 0.0_dp
        do i_s = 1, n_s
            do i_th = 1, n_th
                do i_ph = 1, n_ph
                    s = s_values(i_s)
                    th = twopi*(real(i_th, dp) - 0.7_dp)/real(n_th, dp)
                    ph = phi_period*(real(i_ph, dp) - 0.4_dp)/real(n_ph, dp)
                    x = sqrt(s)*cos(th)
                    y = sqrt(s)*sin(th)
                    call eval_field(polar, s, th, ph, 0)
                    call evaluate_boozer_regular(x, y, ph, 1, reg)

                    polar_value = [polar%Ath, polar%Aph, polar%hth, polar%hph, &
                                   polar%Bmod]
                    reg_value = [reg%Ath, reg%Aph, reg%hth, reg%hph, reg%Bmod]
                    polar_ds = [polar%dAth(1), polar%dAph(1), polar%dhth(1), &
                                polar%dhph(1), polar%dBmod(1)]
                    polar_dth = [polar%dAth(2), polar%dAph(2), polar%dhth(2), &
                                 polar%dhph(2), polar%dBmod(2)]
                    polar_dph = [polar%dAth(3), polar%dAph(3), polar%dhth(3), &
                                 polar%dhph(3), polar%dBmod(3)]
                    ! d/ds = (X d/dX + Y d/dY)/(2 s), d/dtheta = X d/dY - Y d/dX
                    reg_ds = polar_from_xy(reg, x, y, s, 1)
                    reg_dth = polar_from_xy(reg, x, y, s, 2)
                    reg_dph = [reg%dAth(3), reg%dAph(3), reg%dhth(3), &
                               reg%dhph(3), reg%dBmod(3)]

                    vscale = max(abs(polar_value), 1.0e-3_dp*maxval(abs(polar_value)))
                    err_value = max(err_value, &
                                    maxval(abs(reg_value - polar_value)/vscale))
                    dscale = max(abs(polar_ds), abs(polar_dth), abs(polar_dph), &
                                 1.0e-3_dp*abs(polar_value))
                    ! h_theta = I/B is a small, noisy flux function on the
                    ! nearly current-free QA equilibrium; measure both h
                    ! components on the common scale of the covariant unit
                    ! vector.
                    dscale(3:4) = maxval(dscale(3:4))
                    err_derivative = max(err_derivative, &
                                         maxval(abs(reg_ds - polar_ds)/dscale), &
                                         maxval(abs(reg_dth - polar_dth)/dscale), &
                                         maxval(abs(reg_dph - polar_dph)/dscale))
                end do
            end do
        end do
        print '(a,es10.3,a,es10.3)', 'max relative value deviation ', err_value, &
            ', derivative deviation ', err_derivative
        call check(err_value < tol_value, 'values agree with polar splines')
        call check(err_derivative < tol_derivative, &
                   'derivatives agree with polar splines')
    end subroutine compare_with_polar

    function polar_from_xy(reg, x, y, s, component) result(d)
        type(regular_field_t), intent(in) :: reg
        real(dp), intent(in) :: x, y, s
        integer, intent(in) :: component
        real(dp) :: d(5)
        real(dp) :: dx(5), dy(5)

        dx = [reg%dAth(1), reg%dAph(1), reg%dhth(1), reg%dhph(1), reg%dBmod(1)]
        dy = [reg%dAth(2), reg%dAph(2), reg%dhth(2), reg%dhph(2), reg%dBmod(2)]
        if (component == 1) then
            d = (x*dx + y*dy)/(2.0_dp*s)
        else
            d = x*dy - y*dx
        end if
    end function polar_from_xy

    subroutine check_counter()
        type(regular_field_t) :: reg
        integer(8) :: before

        before = n_field_evaluations
        call evaluate_boozer_regular(0.0_dp, 0.0_dp, 0.0_dp, 2, reg)
        call check(n_field_evaluations == before + 1, &
                   'regular evaluation counted as one field evaluation')
        call check(reg%Bmod > 0.0_dp .and. abs(reg%hth) <= tiny(1.0_dp), &
                   'axis value finite with h_theta = 0')
    end subroutine check_counter

end program test_boozer_regular_axis_field
