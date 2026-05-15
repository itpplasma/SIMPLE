program test_bminmax_lifecycle
    use bminmax_mod, only: prop
    use find_bminmax_sub, only: get_bminmax, init_bminmax_arrays
    use magfie_sub, only: dp, magfie
    use new_vmec_stuff_mod, only: nper
    use params, only: class_plot, ntcut, num_surf
    use simple_main, only: needs_bminmax_cache

    implicit none

    integer :: errors

    errors = 0
    call test_cache_gate(errors)
    call test_active_backend_cache(errors)

    if (errors /= 0) then
        print *, "test_bminmax_lifecycle failed with", errors, "errors"
        error stop
    end if

    print *, "test_bminmax_lifecycle passed"

contains

    subroutine test_cache_gate(errors)
        integer, intent(inout) :: errors

        ntcut = 0
        class_plot = .false.

        num_surf = 1
        call assert_false(needs_bminmax_cache(), &
            "single-surface normal tracing must use init_starting_surf bmin/bmax", &
            errors)

        num_surf = 0
        call assert_true(needs_bminmax_cache(), &
            "normal volume tracing must initialize the bminmax cache", errors)

        num_surf = 2
        call assert_true(needs_bminmax_cache(), &
            "normal multi-surface tracing must initialize the bminmax cache", errors)

        class_plot = .true.
        num_surf = 0
        call assert_false(needs_bminmax_cache(), &
            "classifier num_surf=0 semantics stay out of the bminmax PR", errors)

        num_surf = 2
        call assert_true(needs_bminmax_cache(), &
            "classifier multi-surface tracing must initialize the bminmax cache", &
            errors)

        class_plot = .false.
        ntcut = 1
        num_surf = 0
        call assert_false(needs_bminmax_cache(), &
            "ntcut classifier num_surf=0 semantics stay out of the bminmax PR", errors)

        ntcut = 0
        num_surf = 1
    end subroutine test_cache_gate

    subroutine test_active_backend_cache(errors)
        integer, intent(inout) :: errors
        real(dp) :: bmin, bmax

        nper = 1
        prop = .true.
        magfie => test_magfie_backend

        call init_bminmax_arrays
        call get_bminmax(0.4_dp, bmin, bmax)

        call assert_close(bmin, 2.25_dp, 1.0e-7_dp, &
            "bminmax cache did not use the active test magfie backend for bmin", errors)
        call assert_close(bmax, 2.55_dp, 1.0e-7_dp, &
            "bminmax cache did not use the active test magfie backend for bmax", errors)

        prop = .true.
    end subroutine test_active_backend_cache

    subroutine test_magfie_backend(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: bmod, sqrtg
        real(dp), intent(out) :: bder(3), hcovar(3), hctrvr(3), hcurl(3)

        real(dp), parameter :: amp_theta = 0.1_dp
        real(dp), parameter :: amp_phi = 0.05_dp
        real(dp), parameter :: theta0 = 0.2_dp
        real(dp), parameter :: phi0 = 0.4_dp
        real(dp) :: phase_theta, phase_phi, kphi

        kphi = real(max(1, nper), dp)
        phase_theta = x(2) - theta0
        phase_phi = kphi*x(3) - phi0

        bmod = 2.0_dp + x(1) + amp_theta*cos(phase_theta) + amp_phi*cos(phase_phi)
        sqrtg = 1.0_dp
        bder(1) = 1.0_dp/bmod
        bder(2) = -amp_theta*sin(phase_theta)/bmod
        bder(3) = -amp_phi*kphi*sin(phase_phi)/bmod
        hcovar = 0.0_dp
        hctrvr = [0.0_dp, 0.0_dp, 1.0_dp]
        hcurl = 0.0_dp
    end subroutine test_magfie_backend

    subroutine assert_true(condition, message, errors)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message
        integer, intent(inout) :: errors

        if (.not. condition) then
            print *, "ERROR:", trim(message)
            errors = errors + 1
        end if
    end subroutine assert_true

    subroutine assert_false(condition, message, errors)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message
        integer, intent(inout) :: errors

        call assert_true(.not. condition, message, errors)
    end subroutine assert_false

    subroutine assert_close(actual, expected, tolerance, message, errors)
        real(dp), intent(in) :: actual, expected, tolerance
        character(len=*), intent(in) :: message
        integer, intent(inout) :: errors

        if (abs(actual - expected) > tolerance) then
            print *, "ERROR:", trim(message)
            print *, "  actual=", actual, " expected=", expected
            errors = errors + 1
        end if
    end subroutine assert_close

end program test_bminmax_lifecycle
