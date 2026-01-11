program test_chartmap_meiss_debug
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use libneo_coordinates, only: coordinate_system_t, make_chartmap_coordinate_system, &
                                  chartmap_coordinate_system_t, RHO_TOR, RHO_POL, &
                                  PSI_TOR_NORM, PSI_POL_NORM, UNKNOWN
    use field_base, only: magnetic_field_t
    use field_splined, only: splined_field_t, create_splined_field
    use field_can_meiss, only: choose_default_scaling
    use coordinate_scaling, only: coordinate_scaling_t, sqrt_s_scaling_t, identity_scaling_t
    implicit none

    character(len=512) :: chartmap_file, coils_file, vmec_file
    class(coordinate_system_t), allocatable :: coords
    integer :: n_tests_passed, n_tests_failed

    n_tests_passed = 0
    n_tests_failed = 0

    call get_command_argument(1, chartmap_file)
    call get_command_argument(2, coils_file)
    call get_command_argument(3, vmec_file)

    if (len_trim(chartmap_file) == 0) then
        print *, "Usage: test_chartmap_meiss_debug <chartmap.nc> [coils_file] [vmec_file]"
        stop 1
    end if

    print *, "Testing chartmap: ", trim(chartmap_file)
    print *, ""

    call test_rho_convention()
    call test_coordinate_evaluation()
    call test_covariant_basis()
    call test_spline_grid_range()
    call test_meiss_scaling_selection()

    print *, ""
    print *, "============================================"
    print '(a,i0,a,i0)', " RESULTS: ", n_tests_passed, " passed, ", n_tests_failed, " failed"
    print *, "============================================"

    if (n_tests_failed > 0) stop 1

contains

    subroutine test_rho_convention()
        print *, "=== Test 1: rho_convention parsing ==="

        call make_chartmap_coordinate_system(coords, trim(chartmap_file))

        select type (cs => coords)
        type is (chartmap_coordinate_system_t)
            print '(a,i0)', "  rho_convention = ", cs%rho_convention

            ! UNKNOWN is acceptable - grid range fix handles it
            if (cs%rho_convention == RHO_TOR .or. cs%rho_convention == UNKNOWN) then
                print *, "  PASS: rho_convention is valid (RHO_TOR or UNKNOWN)"
                n_tests_passed = n_tests_passed + 1
            else
                print *, "  FAIL: unexpected rho_convention"
                n_tests_failed = n_tests_failed + 1
            end if
        class default
            print *, "  FAIL: coords is not chartmap_coordinate_system_t"
            n_tests_failed = n_tests_failed + 1
        end select
        print *, ""
    end subroutine test_rho_convention

    subroutine test_coordinate_evaluation()
        real(dp) :: u(3), x_cart(3), x_cyl(3)
        real(dp) :: R

        print *, "=== Test 2: Coordinate evaluation ==="

        u = [0.5d0, 0d0, 0d0]
        call coords%evaluate_cart(u, x_cart)
        call coords%evaluate_cyl(u, x_cyl)
        R = x_cyl(1)

        print '(a,3f12.4)', "  Input (rho,theta,zeta) = ", u
        print '(a,f12.4)', "  R [m] = ", R/100d0

        if (R > 1000d0 .and. R < 1500d0) then
            print *, "  PASS: R is in expected range (10-15 m)"
            n_tests_passed = n_tests_passed + 1
        else
            print *, "  FAIL: R is out of expected range"
            n_tests_failed = n_tests_failed + 1
        end if
        print *, ""
    end subroutine test_coordinate_evaluation

    subroutine test_covariant_basis()
        real(dp) :: u(3), e_cov(3, 3)
        real(dp) :: norm_rho, norm_theta, norm_zeta

        print *, "=== Test 3: Covariant basis vectors ==="

        u = [0.5d0, 1.5707963267948966d0, 0.1d0]
        call coords%covariant_basis(u, e_cov)

        norm_rho = sqrt(sum(e_cov(:, 1)**2))
        norm_theta = sqrt(sum(e_cov(:, 2)**2))
        norm_zeta = sqrt(sum(e_cov(:, 3)**2))

        print '(a,3f12.4)', "  Input (rho,theta,zeta) = ", u
        print '(a,f12.6)', "  |e_rho|   = ", norm_rho
        print '(a,f12.6)', "  |e_theta| = ", norm_theta
        print '(a,f12.6)', "  |e_zeta|  = ", norm_zeta

        if (norm_rho > 1d0 .and. norm_theta > 1d0 .and. norm_zeta > 100d0) then
            print *, "  PASS: basis vectors are non-degenerate"
            n_tests_passed = n_tests_passed + 1
        else
            print *, "  FAIL: basis vectors may be degenerate"
            n_tests_failed = n_tests_failed + 1
        end if
        print *, ""
    end subroutine test_covariant_basis

    subroutine test_spline_grid_range()
        real(dp) :: rho_min, rho_max

        print *, "=== Test 4: Spline grid range ==="

        select type (cs => coords)
        type is (chartmap_coordinate_system_t)
            if (cs%has_spl_rz) then
                rho_min = cs%spl_rz%x_min(1)
                rho_max = cs%spl_rz%x_min(1) + cs%spl_rz%h_step(1) * &
                          real(cs%spl_rz%num_points(1) - 1, dp)
            else
                rho_min = cs%spl_cart%x_min(1)
                rho_max = cs%spl_cart%x_min(1) + cs%spl_cart%h_step(1) * &
                          real(cs%spl_cart%num_points(1) - 1, dp)
            end if

            print '(a,es12.4)', "  rho_min = ", rho_min
            print '(a,f12.6)', "  rho_max = ", rho_max
            print '(a,es12.4)', "  Meiss needs rho >= ", 1d-6

            ! rho_min should be small enough for Meiss integration
            if (rho_min <= 0.01d0 .and. rho_max >= 0.99d0) then
                print *, "  PASS: grid range covers Meiss integration domain"
                n_tests_passed = n_tests_passed + 1
            else
                print *, "  FAIL: grid range too narrow for Meiss"
                n_tests_failed = n_tests_failed + 1
            end if
        end select
        print *, ""
    end subroutine test_spline_grid_range

    subroutine test_meiss_scaling_selection()
        class(coordinate_scaling_t), allocatable :: scaling
        type(splined_field_t) :: dummy_field
        logical :: expect_sqrt_s

        print *, "=== Test 5: Meiss scaling selection ==="

        ! Create a minimal splined_field_t with our coords
        allocate(dummy_field%coords, source=coords)

        call choose_default_scaling(dummy_field, scaling)

        expect_sqrt_s = .false.
        select type (cs => coords)
        type is (chartmap_coordinate_system_t)
            ! For known rho conventions, chartmap rho is treated as a normalized
            ! flux label s in [0,1], so Meiss integration expects sqrt(s) scaling.
            if (cs%rho_convention /= UNKNOWN) expect_sqrt_s = .true.
        class default
            expect_sqrt_s = .false.
        end select

        select type (s => scaling)
        type is (identity_scaling_t)
            print *, "  Scaling type: identity_scaling"
            if (expect_sqrt_s) then
                print *, "  FAIL: expected sqrt_s scaling for this chartmap rho convention"
                n_tests_failed = n_tests_failed + 1
            else
                print *, "  PASS: identity scaling selected"
                n_tests_passed = n_tests_passed + 1
            end if
        type is (sqrt_s_scaling_t)
            print *, "  Scaling type: sqrt_s_scaling"
            if (expect_sqrt_s) then
                print *, "  PASS: sqrt_s scaling selected"
                n_tests_passed = n_tests_passed + 1
            else
                print *, "  FAIL: unexpected sqrt_s scaling"
                n_tests_failed = n_tests_failed + 1
            end if
        class default
            print *, "  FAIL: unknown scaling type"
            n_tests_failed = n_tests_failed + 1
        end select
        print *, ""
    end subroutine test_meiss_scaling_selection

end program test_chartmap_meiss_debug
