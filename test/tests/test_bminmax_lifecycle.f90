program test_bminmax_lifecycle
    use bminmax_mod, only: prop
    use find_bminmax_sub, only: get_bminmax, init_bminmax_arrays
    use magfie_sub, only: dp, magfie
    use new_vmec_stuff_mod, only: nper
    use params, only: class_plot, ntcut, num_surf
    use simple_main, only: needs_bminmax_cache
    use test_utils, only: check, check_close

    implicit none

    integer :: errors

    errors = 0
    call test_cache_gate_cases(errors)
    call test_active_backend_cache(errors)

    if (errors /= 0) then
        print *, "test_bminmax_lifecycle failed with", errors, "errors"
        error stop
    end if

    print *, "test_bminmax_lifecycle passed"

contains

    subroutine test_cache_gate_cases(errors)
        integer, intent(inout) :: errors
        character(len=256) :: cases_path
        character(len=1) :: class_plot_token, expected_token
        integer :: unit_id, ios, case_num
        logical :: expected, actual

        call get_command_argument(1, cases_path)
        if (len_trim(cases_path) == 0) cases_path = "bminmax_cache_cases.tsv"

        open(newunit=unit_id, file=trim(cases_path), status="old", action="read")
        case_num = 0

        do
            read(unit_id, *, iostat=ios) num_surf, ntcut, class_plot_token, &
                expected_token
            if (ios < 0) exit
            if (ios > 0) then
                print *, "ERROR: failed to read bminmax cache case", case_num + 1
                errors = errors + 1
                exit
            end if

            case_num = case_num + 1
            class_plot = class_plot_token == "T"
            expected = expected_token == "T"
            actual = needs_bminmax_cache()
            call check(actual .eqv. expected, "bminmax cache predicate case failed", &
                errors)
        end do

        close(unit_id)

        call check(case_num > 0, "bminmax cache predicate cases were empty", errors)

        ntcut = 0
        class_plot = .false.
        num_surf = 1
    end subroutine test_cache_gate_cases

    subroutine test_active_backend_cache(errors)
        integer, intent(inout) :: errors
        real(dp) :: bmin, bmax

        nper = 1
        prop = .true.
        magfie => test_magfie_backend

        call init_bminmax_arrays
        call get_bminmax(0.4_dp, bmin, bmax)

        call check_close(real(bmin, kind(1.0d0)), 2.25d0, 1.0d-7, &
            "bminmax cache did not use the active test magfie backend for bmin", errors)
        call check_close(real(bmax, kind(1.0d0)), 2.55d0, 1.0d-7, &
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

end program test_bminmax_lifecycle
