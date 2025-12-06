program test_coils_field_compare
    !> Compare splined_field_t vs direct Biot-Savart at specific test points.
    !> This test verifies the coordinate transformation is correct.

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use simple, only: init_vmec
    use magfie_sub, only: VMEC
    use velo_mod, only: isw_field_type
    use field_coils, only: coils_field_t, create_coils_field
    use field_splined, only: splined_field_t, create_splined_field
    use reference_coordinates, only: init_reference_coordinates, ref_coords
    use simple_coordinates, only: transform_vmec_to_cart
    use neo_biotsavart, only: compute_vector_potential, compute_magnetic_field

    implicit none

    type(coils_field_t) :: raw_coils
    type(splined_field_t) :: splined_coils
    real(dp) :: x_test(3), Acov_spl(3), hcov_spl(3), Bmod_spl
    real(dp) :: Acov_direct(3), hcov_direct(3), Bmod_direct
    real(dp) :: max_Acov_rel_diff, max_hcov_rel_diff, max_Bmod_rel_diff
    real(dp) :: x_cart(3), dxcart_dxvmec(3,3), A_cart(3), B_cart(3)
    real(dp) :: s, ds_dr, dummy
    integer :: i_r, i_th, i_phi, n_tests
    logical :: passed

    logical :: wout_exists, coils_exists

    isw_field_type = VMEC

    inquire(file='wout.nc', exist=wout_exists)
    inquire(file='coils.simple', exist=coils_exists)
    if (.not. wout_exists .or. .not. coils_exists) then
        print *, 'SKIP: Missing wout.nc or coils.simple'
        stop 0
    end if

    call init_vmec('wout.nc', 5, 5, 5, dummy)
    call init_reference_coordinates('wout.nc')

    call create_coils_field('coils.simple', raw_coils)
    call create_splined_field(raw_coils, ref_coords, splined_coils)

    max_Acov_rel_diff = 0d0
    max_hcov_rel_diff = 0d0
    max_Bmod_rel_diff = 0d0
    n_tests = 0

    ! Test at various points in the (r, theta, phi) space
    do i_phi = 1, 5
        do i_th = 1, 5
            do i_r = 1, 5
                x_test(1) = 0.1d0 + 0.15d0 * (i_r - 1)
                x_test(2) = 6.28318530718d0 * (i_th - 1) / 4d0
                x_test(3) = 2.0943951024d0 * (i_phi - 1) / 4d0

                call splined_coils%evaluate(x_test, Acov_spl, hcov_spl, Bmod_spl)

                s = x_test(1)**2
                ds_dr = 2d0 * x_test(1)
                call transform_vmec_to_cart([s, x_test(2), x_test(3)], &
                                            x_cart, dxcart_dxvmec)
                A_cart = compute_vector_potential(raw_coils%coils, x_cart)
                B_cart = compute_magnetic_field(raw_coils%coils, x_cart)
                Bmod_direct = sqrt(B_cart(1)**2 + B_cart(2)**2 + B_cart(3)**2)

                Acov_direct = matmul(A_cart, dxcart_dxvmec)
                Acov_direct(1) = Acov_direct(1) * ds_dr
                hcov_direct = matmul(B_cart, dxcart_dxvmec) / Bmod_direct
                hcov_direct(1) = hcov_direct(1) * ds_dr

                ! Debug: print values for first test point
                if (n_tests == 0) then
                    print '(A,3ES15.6)', 'x_test = ', x_test
                    print '(A,3ES15.6)', 'Acov_spl = ', Acov_spl
                    print '(A,3ES15.6)', 'Acov_direct = ', Acov_direct
                    print '(A,3ES15.6)', 'hcov_spl = ', hcov_spl
                    print '(A,3ES15.6)', 'hcov_direct = ', hcov_direct
                    print '(A,ES15.6)', 'Bmod_spl = ', Bmod_spl
                    print '(A,ES15.6)', 'Bmod_direct = ', Bmod_direct
                end if

                max_Acov_rel_diff = max(max_Acov_rel_diff, &
                    maxval(abs(Acov_spl - Acov_direct)) / max(1d0, maxval(abs(Acov_direct))))
                max_hcov_rel_diff = max(max_hcov_rel_diff, &
                    maxval(abs(hcov_spl - hcov_direct)) / max(1d-10, maxval(abs(hcov_direct))))
                max_Bmod_rel_diff = max(max_Bmod_rel_diff, &
                    abs(Bmod_spl - Bmod_direct) / Bmod_direct)

                n_tests = n_tests + 1
            end do
        end do
    end do

    print '(A,I6)', 'Number of test points: ', n_tests
    print '(A,ES12.4)', 'Max relative |Acov_spl - Acov_direct|: ', max_Acov_rel_diff
    print '(A,ES12.4)', 'Max relative |hcov_spl - hcov_direct|: ', max_hcov_rel_diff
    print '(A,ES12.4)', 'Max relative |Bmod_spl - Bmod_direct|: ', max_Bmod_rel_diff

    ! Spline interpolation should be within 1% for reasonably sampled grid
    passed = (max_Acov_rel_diff < 0.01d0) .and. &
             (max_hcov_rel_diff < 0.01d0) .and. &
             (max_Bmod_rel_diff < 0.01d0)

    if (passed) then
        print *, 'PASSED: Splined field matches direct Biot-Savart within 1%'
    else
        print *, 'FAILED: Differences exceed 1% tolerance'
        error stop 1
    end if

end program test_coils_field_compare
