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

    implicit none

    type(coils_field_t) :: raw_coils
    type(splined_field_t) :: field_vmec_ref, field_chartmap_ref
    class(coordinate_system_t), allocatable :: vmec_coords, chartmap_coords

    real(dp) :: max_Bmod_rel_err, max_Acov_rel_err, max_hcov_rel_err
    real(dp) :: dummy
    integer :: n_tests, i
    logical :: wout_exists, coils_exists, chartmap_exists
    logical :: passed

    real(dp), parameter :: tol_Bmod = 1.0e-6_dp
    real(dp), parameter :: tol_Acov = 1.0e-5_dp
    real(dp), parameter :: tol_hcov = 1.0e-5_dp

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
        do i_phi = 1, 4
            do i_th = 1, 5
                do i_r = 1, 5
                    ! VMEC coordinates: r = sqrt(s), theta, phi
                    x_vmec(1) = 0.15d0 + 0.15d0 * (i_r - 1)  ! r
                    x_vmec(2) = twopi * (i_th - 1) / 4d0     ! theta
                    x_vmec(3) = twopi * (i_phi - 1) / 12d0   ! phi (1/3 period)

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

end program test_chartmap_coils_field
