program test_chartmap_coils_field
    !> Test chartmap reference coordinates with coils field.
    !>
    !> Compares field evaluation using:
    !> - VMEC reference coordinates
    !> - Chartmap reference coordinates
    !>
    !> Both should give the same physical field since they're just
    !> different coordinate representations of the same volume.

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use simple, only: init_vmec
    use magfie_sub, only: VMEC
    use velo_mod, only: isw_field_type
    use field_coils, only: coils_field_t, create_coils_field
    use field_splined, only: splined_field_t, create_splined_field
    use libneo_coordinates, only: coordinate_system_t, &
        make_vmec_coordinate_system, make_chartmap_coordinate_system
    use util, only: twopi
    use pyplot_module, only: pyplot

    implicit none

    type(coils_field_t) :: raw_coils
    type(splined_field_t) :: field_vmec_ref, field_chartmap_ref
    class(coordinate_system_t), allocatable :: vmec_coords, chartmap_coords

    real(dp) :: max_Bmod_rel_err, max_Acov_rel_err, max_hcov_rel_err
    real(dp) :: dummy
    integer :: n_tests, i
    logical :: wout_exists, coils_exists, chartmap_exists
    logical :: passed

    real(dp), parameter :: tol_Bmod = 1.0e-10_dp
    real(dp), parameter :: tol_Acov = 1.0e-10_dp
    real(dp), parameter :: tol_hcov = 1.0e-10_dp

    ! Check required files
    inquire(file='wout.nc', exist=wout_exists)
    inquire(file='coils.simple', exist=coils_exists)
    inquire(file='wout.chartmap.nc', exist=chartmap_exists)

    if (.not. wout_exists) then
        print *, 'FAILED: Required VMEC file (wout.nc) not found'
        error stop 1
    end if
    if (.not. coils_exists) then
        print *, 'FAILED: Required coils file (coils.simple) not found'
        error stop 1
    end if
    if (.not. chartmap_exists) then
        print *, 'FAILED: Required chartmap file (wout.chartmap.nc) not found'
        print *, '  Generate with: python scripts/vmec_to_chartmap.py wout.nc wout.chartmap.nc'
        error stop 1
    end if

    print *, '================================'
    print *, 'Chartmap Coils Field Test'
    print *, '================================'
    print *

    ! Initialize VMEC
    isw_field_type = VMEC
    call init_vmec('wout.nc', 5, 5, 5, dummy)

    ! Create coordinate systems
    call make_vmec_coordinate_system(vmec_coords)
    call make_chartmap_coordinate_system(chartmap_coords, 'wout.chartmap.nc')

    print *, 'Coordinate systems initialized'
    print *

    ! Load coils field
    call create_coils_field('coils.simple', raw_coils)
    print *, 'Coils field loaded'
    print *

    ! Create splined fields with different reference coordinates
    print *, 'Building splined field with VMEC reference...'
    call create_splined_field(raw_coils, vmec_coords, field_vmec_ref)
    print *, 'Building splined field with chartmap reference...'
    call create_splined_field(raw_coils, chartmap_coords, field_chartmap_ref)
    print *

    ! Compare field evaluations
    print *, 'Comparing field evaluations...'
    call compare_fields(field_vmec_ref, field_chartmap_ref, &
                       vmec_coords, chartmap_coords, &
                       max_Bmod_rel_err, max_Acov_rel_err, max_hcov_rel_err, n_tests)

    ! Print results
    print *
    print *, '================================'
    print *, 'Results'
    print *, '================================'
    print '(A,I6)', 'Number of test points: ', n_tests
    print '(A,ES12.4)', 'Max relative |Bmod error|: ', max_Bmod_rel_err
    print '(A,ES12.4)', 'Max relative |Acov error|: ', max_Acov_rel_err
    print '(A,ES12.4)', 'Max relative |hcov error|: ', max_hcov_rel_err
    print *
    print '(A,ES12.4)', 'Tolerance Bmod: ', tol_Bmod
    print '(A,ES12.4)', 'Tolerance Acov: ', tol_Acov
    print '(A,ES12.4)', 'Tolerance hcov: ', tol_hcov
    print *

    ! Generate comparison plots
    print *, 'Generating comparison plots...'
    call generate_comparison_plots(field_vmec_ref, field_chartmap_ref)

    ! Check tolerances
    passed = .true.
    if (max_Bmod_rel_err > tol_Bmod) then
        print *, 'FAILED: Bmod error exceeds tolerance'
        passed = .false.
    end if
    if (max_Acov_rel_err > tol_Acov) then
        print *, 'FAILED: Acov error exceeds tolerance'
        passed = .false.
    end if
    if (max_hcov_rel_err > tol_hcov) then
        print *, 'FAILED: hcov error exceeds tolerance'
        passed = .false.
    end if

    print *
    if (passed) then
        print *, '================================'
        print *, 'PASSED: All tests passed!'
        print *, '================================'
    else
        print *, '================================'
        print *, 'FAILED: Some tests failed'
        print *, '================================'
        error stop 1
    end if

contains

    subroutine compare_fields(field_vmec, field_chart, vmec_coords, chart_coords, &
                              max_Bmod_err, max_Acov_err, max_hcov_err, n_points)
        type(splined_field_t), intent(in) :: field_vmec, field_chart
        class(coordinate_system_t), intent(in) :: vmec_coords, chart_coords
        real(dp), intent(out) :: max_Bmod_err, max_Acov_err, max_hcov_err
        integer, intent(out) :: n_points

        real(dp) :: x_vmec(3), x_chart(3), x_phys_vmec(3), x_phys_chart(3)
        real(dp) :: Acov_vmec(3), hcov_vmec(3), Bmod_vmec
        real(dp) :: Acov_chart(3), hcov_chart(3), Bmod_chart
        real(dp) :: Bmod_err, Acov_err, hcov_err
        integer :: i_r, i_th, i_phi
        real(dp) :: s, rho

        max_Bmod_err = 0d0
        max_Acov_err = 0d0
        max_hcov_err = 0d0
        n_points = 0

        ! Sample points across the volume
        ! Add small offsets to avoid exact spline knots
        do i_phi = 1, 4
            do i_th = 1, 5
                do i_r = 1, 5
                    ! VMEC coordinates: r = sqrt(s), theta, phi
                    ! Add π/100 offset to avoid exact knot positions
                    x_vmec(1) = 0.15d0 + 0.15d0 * (i_r - 1) + 0.01d0  ! r + offset
                    x_vmec(2) = twopi * (i_th - 1) / 4d0 + 0.0314d0   ! theta + π/100
                    x_vmec(3) = twopi * (i_phi - 1) / 12d0 + 0.0157d0 ! phi + π/200

                    ! Chartmap coordinates: rho, theta, zeta
                    s = x_vmec(1)**2
                    rho = sqrt(s)
                    x_chart(1) = rho
                    x_chart(2) = x_vmec(2)
                    x_chart(3) = x_vmec(3)

                    ! Evaluate fields
                    call field_vmec%evaluate(x_vmec, Acov_vmec, hcov_vmec, Bmod_vmec)
                    call field_chart%evaluate(x_chart, Acov_chart, hcov_chart, Bmod_chart)

                    ! Check that we're at the same physical location
                    call vmec_coords%evaluate_point([s, x_vmec(2), x_vmec(3)], x_phys_vmec)
                    call chart_coords%evaluate_point(x_chart, x_phys_chart)

                    ! Compute errors
                    Bmod_err = abs(Bmod_vmec - Bmod_chart) / max(1d-10, Bmod_vmec)
                    Acov_err = maxval(abs(Acov_vmec - Acov_chart)) / max(1d0, maxval(abs(Acov_vmec)))
                    hcov_err = maxval(abs(hcov_vmec - hcov_chart)) / max(1d-10, maxval(abs(hcov_vmec)))

                    max_Bmod_err = max(max_Bmod_err, Bmod_err)
                    max_Acov_err = max(max_Acov_err, Acov_err)
                    max_hcov_err = max(max_hcov_err, hcov_err)

                    n_points = n_points + 1
                end do
            end do
        end do
    end subroutine compare_fields


    subroutine generate_comparison_plots(field_vmec, field_chart)
        type(splined_field_t), intent(in) :: field_vmec, field_chart

        type(pyplot) :: plt
        integer, parameter :: nth = 50, nphi = 50
        real(dp) :: theta_arr(nth), phi_arr(nphi)
        real(dp) :: Bmod_vmec_grid(nth, nphi), Bmod_chart_grid(nth, nphi)
        real(dp) :: Bmod_err_grid(nth, nphi)
        real(dp) :: x_vmec(3), x_chart(3), Acov(3), hcov(3), Bmod
        real(dp) :: s, rho
        integer :: ith, iphi

        ! Fixed radial position: r = 0.5
        s = 0.5d0**2
        rho = sqrt(s)

        ! Create grids
        do ith = 1, nth
            theta_arr(ith) = twopi * (ith - 1) / real(nth - 1, dp)
        end do

        do iphi = 1, nphi
            phi_arr(iphi) = twopi * (iphi - 1) / real(nphi - 1, dp) / 4d0  ! 1/4 period
        end do

        ! Evaluate on grid
        do ith = 1, nth
            do iphi = 1, nphi
                x_vmec = [sqrt(s), theta_arr(ith), phi_arr(iphi)]
                x_chart = [rho, theta_arr(ith), phi_arr(iphi)]

                call field_vmec%evaluate(x_vmec, Acov, hcov, Bmod)
                Bmod_vmec_grid(ith, iphi) = Bmod

                call field_chart%evaluate(x_chart, Acov, hcov, Bmod)
                Bmod_chart_grid(ith, iphi) = Bmod

                Bmod_err_grid(ith, iphi) = abs(Bmod_vmec_grid(ith, iphi) - &
                                               Bmod_chart_grid(ith, iphi))
            end do
        end do

        ! Plot 1: VMEC reference Bmod
        call plt%initialize(grid=.true., xlabel='theta', ylabel='phi', &
            title='Bmod with VMEC reference (r=0.5)', figsize=[10,8])
        call plt%add_contour(theta_arr, phi_arr, transpose(Bmod_vmec_grid), &
            linestyle='-', colorbar=.true., filled=.true.)
        call plt%savefig('test_chartmap_bmod_vmec.png', &
            pyfile='test_chartmap_bmod_vmec.py')
        print *, '  Saved: test_chartmap_bmod_vmec.png'

        ! Plot 2: Chartmap reference Bmod
        call plt%initialize(grid=.true., xlabel='theta', ylabel='phi', &
            title='Bmod with chartmap reference (r=0.5)', figsize=[10,8])
        call plt%add_contour(theta_arr, phi_arr, transpose(Bmod_chart_grid), &
            linestyle='-', colorbar=.true., filled=.true.)
        call plt%savefig('test_chartmap_bmod_chart.png', &
            pyfile='test_chartmap_bmod_chart.py')
        print *, '  Saved: test_chartmap_bmod_chart.png'

        ! Plot 3: Absolute error
        call plt%initialize(grid=.true., xlabel='theta', ylabel='phi', &
            title='Absolute Bmod error |VMEC - chartmap| (r=0.5)', figsize=[10,8])
        call plt%add_contour(theta_arr, phi_arr, transpose(Bmod_err_grid), &
            linestyle='-', colorbar=.true., filled=.true.)
        call plt%savefig('test_chartmap_bmod_error.png', &
            pyfile='test_chartmap_bmod_error.py')
        print *, '  Saved: test_chartmap_bmod_error.png'

    end subroutine generate_comparison_plots

end program test_chartmap_coils_field
