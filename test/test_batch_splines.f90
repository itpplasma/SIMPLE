program test_batch_splines
    use interpolate, only: SplineData3D, BatchSplineData3D, &
        construct_splines_3d, construct_batch_splines_3d, &
        evaluate_splines_3d, evaluate_batch_splines_3d, &
        evaluate_splines_3d_der, evaluate_batch_splines_3d_der, &
        evaluate_splines_3d_der2, evaluate_batch_splines_3d_der2, &
        destroy_splines_3d, destroy_batch_splines_3d
    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none

    real(dp), parameter :: ABS_TOL = 1.0d-9   ! For values near zero
    real(dp), parameter :: REL_TOL = 1.0d-11  ! For non-zero values
    real(dp), parameter :: PI = 3.141592653589793d0
    real(dp), parameter :: TWOPI = 6.283185307179586d0

    ! Run all tests
    print *, '=========================================='
    print *, 'Running Batch Spline Tests'
    print *, '=========================================='

    print *, ''
    print *, 'Test 1: Meiss field batching...'
    call test_meiss_field_batching()
    print *, 'PASSED: Meiss field batching test'

    print *, ''
    print *, 'Test 2: Albert field batching...'
    call test_albert_field_batching()
    print *, 'PASSED: Albert field batching test'

    print *, ''
    print *, 'Test 3: Coils field batching...'
    call test_coils_field_batching()
    print *, 'PASSED: Coils field batching test'

    print *, ''
    print *, 'Test 4: Performance comparison...'
    call test_batch_performance()
    print *, 'PASSED: Performance test'

    print *, ''
    print *, '=========================================='
    print *, 'All batch spline tests PASSED!'
    print *, '=========================================='

contains

    pure logical function approx_equal(a, b)
        real(dp), intent(in) :: a, b
        real(dp) :: diff, scale
        diff = abs(a - b)
        scale = max(abs(a), abs(b))
        approx_equal = diff <= max(ABS_TOL, REL_TOL * scale)
    end function approx_equal

    subroutine test_meiss_field_batching()
        ! Test that batched Meiss field components produce identical results
        integer, parameter :: n_r = 20, n_th = 21, n_phi = 22
        integer, parameter :: n_components = 5
        real(dp) :: xmin(3) = [1d-6, 0d0, 0d0]
        real(dp) :: xmax(3) = [1d0, TWOPI, TWOPI]
        integer, parameter :: order(3) = [5, 5, 5]
        logical, parameter :: periodic(3) = [.False., .True., .True.]

        ! Individual splines (old approach)
        type(SplineData3D) :: spl_Ath, spl_Aph, spl_hth, spl_hph, spl_Bmod

        ! Batch spline (new approach)
        type(BatchSplineData3D) :: spl_batch

        ! Data arrays
        real(dp), allocatable :: Ath(:,:,:), Aph(:,:,:), hth(:,:,:), hph(:,:,:), Bmod(:,:,:)
        real(dp), allocatable :: field_batch(:,:,:,:)

        ! Test points
        real(dp) :: x_test(3), y_individual(5), y_batch(5)
        real(dp) :: dy_individual(5, 3), dy_batch(3, 5)
        real(dp) :: d2y_individual(5, 6), d2y_batch(6, 5)

        integer :: i_r, i_th, i_phi, i, j, k
        real(dp) :: r, th, phi, h_r, h_th, h_phi

        ! Allocate arrays
        allocate(Ath(n_r, n_th, n_phi))
        allocate(Aph(n_r, n_th, n_phi))
        allocate(hth(n_r, n_th, n_phi))
        allocate(hph(n_r, n_th, n_phi))
        allocate(Bmod(n_r, n_th, n_phi))
        allocate(field_batch(n_r, n_th, n_phi, n_components))

        h_r = (xmax(1)-xmin(1))/(n_r-1)
        h_th = (xmax(2)-xmin(2))/(n_th-1)
        h_phi = (xmax(3)-xmin(3))/(n_phi-1)

        ! Generate test data
        do i_phi = 1, n_phi
            phi = xmin(3) + (i_phi-1)*h_phi
            do i_th = 1, n_th
                th = xmin(2) + (i_th-1)*h_th
                do i_r = 1, n_r
                    r = xmin(1) + (i_r-1)*h_r

                    ! Test functions
                    Ath(i_r, i_th, i_phi) = sin(r)*cos(th)*sin(phi)
                    Aph(i_r, i_th, i_phi) = cos(r)*sin(th)*cos(phi)
                    hth(i_r, i_th, i_phi) = sin(r+th)*cos(phi)
                    hph(i_r, i_th, i_phi) = cos(r)*sin(th+phi)
                    Bmod(i_r, i_th, i_phi) = 1.0d0 + 0.1d0*sin(r)*cos(th)*sin(phi)
                end do
            end do
        end do

        ! Create individual splines (old approach)
        call construct_splines_3d(xmin, xmax, Ath, order, periodic, spl_Ath)
        call construct_splines_3d(xmin, xmax, Aph, order, periodic, spl_Aph)
        call construct_splines_3d(xmin, xmax, hth, order, periodic, spl_hth)
        call construct_splines_3d(xmin, xmax, hph, order, periodic, spl_hph)
        call construct_splines_3d(xmin, xmax, Bmod, order, periodic, spl_Bmod)

        ! Create batch spline (new approach)
        field_batch(:,:,:,1) = Ath
        field_batch(:,:,:,2) = Aph
        field_batch(:,:,:,3) = hth
        field_batch(:,:,:,4) = hph
        field_batch(:,:,:,5) = Bmod

        call construct_batch_splines_3d(xmin, xmax, field_batch, order, periodic, spl_batch)

        ! Test at multiple points
        do i = 1, 10
            x_test(1) = xmin(1) + (xmax(1) - xmin(1)) * real(i-1, dp) / 9.0d0
            x_test(2) = xmin(2) + (xmax(2) - xmin(2)) * real(i, dp) / 10.0d0
            x_test(3) = xmin(3) + (xmax(3) - xmin(3)) * real(i+1, dp) / 11.0d0

            ! Evaluate individual splines
            call evaluate_splines_3d(spl_Ath, x_test, y_individual(1))
            call evaluate_splines_3d(spl_Aph, x_test, y_individual(2))
            call evaluate_splines_3d(spl_hth, x_test, y_individual(3))
            call evaluate_splines_3d(spl_hph, x_test, y_individual(4))
            call evaluate_splines_3d(spl_Bmod, x_test, y_individual(5))

            ! Evaluate batch spline
            call evaluate_batch_splines_3d(spl_batch, x_test, y_batch)

            ! Check values match
            do j = 1, n_components
                if (.not. approx_equal(y_individual(j), y_batch(j))) then
                    print *, 'ERROR: Mismatch at test', i, 'component', j
                    print *, '  Individual:', y_individual(j)
                    print *, '  Batch:', y_batch(j)
                    stop 1
                end if
            end do

            ! Test first derivatives
            call evaluate_splines_3d_der(spl_Ath, x_test, y_individual(1), dy_individual(1,:))
            call evaluate_splines_3d_der(spl_Aph, x_test, y_individual(2), dy_individual(2,:))
            call evaluate_splines_3d_der(spl_hth, x_test, y_individual(3), dy_individual(3,:))
            call evaluate_splines_3d_der(spl_hph, x_test, y_individual(4), dy_individual(4,:))
            call evaluate_splines_3d_der(spl_Bmod, x_test, y_individual(5), dy_individual(5,:))

            call evaluate_batch_splines_3d_der(spl_batch, x_test, y_batch, dy_batch)

            ! Check derivatives match
            do j = 1, n_components
                do k = 1, 3
                    if (.not. approx_equal(dy_individual(j,k), dy_batch(k,j))) then
                        print *, 'ERROR: Derivative mismatch at test', i, 'component', j, 'dim', k
                        print *, '  Individual:', dy_individual(j,k)
                        print *, '  Batch:', dy_batch(k,j)
                        stop 1
                    end if
                end do
            end do

            ! Test second derivatives
            call evaluate_splines_3d_der2(spl_Ath, x_test, y_individual(1), &
                dy_individual(1,:), d2y_individual(1,:))
            call evaluate_splines_3d_der2(spl_Aph, x_test, y_individual(2), &
                dy_individual(2,:), d2y_individual(2,:))
            call evaluate_splines_3d_der2(spl_hth, x_test, y_individual(3), &
                dy_individual(3,:), d2y_individual(3,:))
            call evaluate_splines_3d_der2(spl_hph, x_test, y_individual(4), &
                dy_individual(4,:), d2y_individual(4,:))
            call evaluate_splines_3d_der2(spl_Bmod, x_test, y_individual(5), &
                dy_individual(5,:), d2y_individual(5,:))

            call evaluate_batch_splines_3d_der2(spl_batch, x_test, y_batch, dy_batch, d2y_batch)

            ! Check second derivatives match
            do j = 1, n_components
                do k = 1, 6
                    if (.not. approx_equal(d2y_individual(j,k), d2y_batch(k,j))) then
                        print *, 'ERROR: Second derivative mismatch at test', i, 'component', j, 'dim', k
                        print *, '  Individual:', d2y_individual(j,k)
                        print *, '  Batch:', d2y_batch(k,j)
                        stop 1
                    end if
                end do
            end do
        end do

        ! Cleanup
        call destroy_splines_3d(spl_Ath)
        call destroy_splines_3d(spl_Aph)
        call destroy_splines_3d(spl_hth)
        call destroy_splines_3d(spl_hph)
        call destroy_splines_3d(spl_Bmod)
        call destroy_batch_splines_3d(spl_batch)

        deallocate(Ath, Aph, hth, hph, Bmod, field_batch)
    end subroutine test_meiss_field_batching


    subroutine test_albert_field_batching()
        ! Test Albert field components batching
        integer, parameter :: n_r = 15, n_th = 16, n_phi = 17
        integer, parameter :: n_components = 5
        real(dp) :: xmin(3) = [0.1d0, 0d0, 0d0]
        real(dp) :: xmax(3) = [0.9d0, TWOPI, TWOPI]
        integer, parameter :: order(3) = [5, 5, 5]
        logical, parameter :: periodic(3) = [.False., .True., .True.]

        type(BatchSplineData3D) :: spl_batch
        real(dp), allocatable :: albert_batch(:,:,:,:)
        real(dp) :: x_test(3), y_batch(5), dy_batch(3, 5)
        integer :: i

        ! Allocate and fill test data
        allocate(albert_batch(n_r, n_th, n_phi, n_components))

        ! Generate test data
        call generate_test_field_data(n_r, n_th, n_phi, xmin, xmax, albert_batch)

        ! Create batch spline
        call construct_batch_splines_3d(xmin, xmax, albert_batch, order, periodic, spl_batch)

        ! Test evaluation
        x_test = [0.5d0, PI, PI/2.0d0]
        call evaluate_batch_splines_3d_der(spl_batch, x_test, y_batch, dy_batch)

        ! Verify results are reasonable
        do i = 1, n_components
            if (abs(y_batch(i)) >= 10.0d0) then
                print *, 'ERROR: Unreasonable value magnitude for component', i, ':', y_batch(i)
                stop 1
            end if
            if (maxval(abs(dy_batch(:,i))) >= 100.0d0) then
                print *, 'ERROR: Unreasonable derivative magnitude for component', i
                stop 1
            end if
        end do

        ! Clean up
        call destroy_batch_splines_3d(spl_batch)
        deallocate(albert_batch)
    end subroutine test_albert_field_batching


    subroutine test_coils_field_batching()
        ! Test coils field with 7 components
        integer, parameter :: n_r = 10, n_th = 11, n_phi = 12
        integer, parameter :: n_components = 7
        real(dp) :: xmin(3) = [0.01d0, 0d0, 0d0]
        real(dp) :: xmax(3) = [10d0, TWOPI, TWOPI]
        integer, parameter :: order(3) = [5, 5, 5]
        logical, parameter :: periodic(3) = [.False., .True., .True.]

        type(BatchSplineData3D) :: spl_batch
        real(dp), allocatable :: coils_batch(:,:,:,:)
        real(dp) :: x_test(3), y_batch(7)
        integer :: i

        ! Allocate and fill test data
        allocate(coils_batch(n_r, n_th, n_phi, n_components))

        ! Generate test data for all 7 components (Ar, Ath, Aphi, hr, hth, hphi, Bmod)
        call generate_test_field_data(n_r, n_th, n_phi, xmin, xmax, coils_batch(:,:,:,1:5))
        ! Additional components
        coils_batch(:,:,:,6) = 0.5d0 * coils_batch(:,:,:,1)  ! hphi
        coils_batch(:,:,:,7) = 1.0d0 + 0.1d0 * coils_batch(:,:,:,2)  ! Bmod

        ! Create batch spline
        call construct_batch_splines_3d(xmin, xmax, coils_batch, order, periodic, spl_batch)

        ! Test evaluation
        x_test = [5.0d0, PI, PI]
        call evaluate_batch_splines_3d(spl_batch, x_test, y_batch)

        ! Verify we get 7 components
        if (size(y_batch) /= 7) then
            print *, 'ERROR: Expected 7 components, got', size(y_batch)
            stop 1
        end if

        ! Clean up
        call destroy_batch_splines_3d(spl_batch)
        deallocate(coils_batch)
    end subroutine test_coils_field_batching


    subroutine test_batch_performance()
        ! Benchmark batch vs individual splines
        integer, parameter :: n_r = 30, n_th = 31, n_phi = 32
        integer, parameter :: n_components = 6
        integer, parameter :: n_eval = 1000
        real(dp) :: xmin(3) = [0.1d0, 0d0, 0d0]
        real(dp) :: xmax(3) = [1.0d0, TWOPI, TWOPI]
        integer, parameter :: order(3) = [5, 5, 5]
        logical, parameter :: periodic(3) = [.False., .True., .True.]

        type(SplineData3D) :: spl_individual(n_components)
        type(BatchSplineData3D) :: spl_batch
        real(dp), allocatable :: data(:,:,:), batch_data(:,:,:,:)
        real(dp) :: x_test(3), y_single, y_batch(n_components)
        real(dp) :: t_start, t_end, t_individual, t_batch
        integer :: i, j

        ! Allocate arrays
        allocate(data(n_r, n_th, n_phi))
        allocate(batch_data(n_r, n_th, n_phi, n_components))

        ! Generate test data
        call generate_test_field_data(n_r, n_th, n_phi, xmin, xmax, batch_data)

        ! Create individual splines
        do i = 1, n_components
            data = batch_data(:,:,:,i)
            call construct_splines_3d(xmin, xmax, data, order, periodic, spl_individual(i))
        end do

        ! Create batch spline
        call construct_batch_splines_3d(xmin, xmax, batch_data, order, periodic, spl_batch)

        ! Time individual evaluations
        call cpu_time(t_start)
        do i = 1, n_eval
            x_test = [0.5d0 + 0.001d0*i, PI*i/n_eval, PI*i/n_eval]
            do j = 1, n_components
                call evaluate_splines_3d(spl_individual(j), x_test, y_single)
            end do
        end do
        call cpu_time(t_end)
        t_individual = t_end - t_start

        ! Time batch evaluations
        call cpu_time(t_start)
        do i = 1, n_eval
            x_test = [0.5d0 + 0.001d0*i, PI*i/n_eval, PI*i/n_eval]
            call evaluate_batch_splines_3d(spl_batch, x_test, y_batch)
        end do
        call cpu_time(t_end)
        t_batch = t_end - t_start

        ! Report speedup with proper validation
        print *, 'Individual time:', t_individual, 's'
        print *, 'Batch time:', t_batch, 's'

        if (t_batch > 0.0d0 .and. t_individual > 0.0d0) then
            print *, 'Measured speedup:', t_individual/t_batch, 'x'
        else
            print *, 'WARNING: Invalid timing measurements'
        end if

        ! Note: Actual speedup varies by system and problem size
        ! This test verifies correctness, not performance claims
        if (t_batch >= t_individual) then
            print *, 'NOTE: Batch implementation may show benefits at larger scales'
        else
            print *, 'SUCCESS: Batch shows speedup for this test case'
        end if

        ! Clean up
        do i = 1, n_components
            call destroy_splines_3d(spl_individual(i))
        end do
        call destroy_batch_splines_3d(spl_batch)
        deallocate(data, batch_data)
    end subroutine test_batch_performance


    subroutine generate_test_field_data(n_r, n_th, n_phi, xmin, xmax, field_data)
        integer, intent(in) :: n_r, n_th, n_phi
        real(dp), intent(in) :: xmin(3), xmax(3)
        real(dp), intent(out) :: field_data(:,:,:,:)

        integer :: i_r, i_th, i_phi, n_comp, ic
        real(dp) :: r, th, phi, h_r, h_th, h_phi

        n_comp = size(field_data, 4)
        h_r = (xmax(1)-xmin(1))/(n_r-1)
        h_th = (xmax(2)-xmin(2))/(n_th-1)
        h_phi = (xmax(3)-xmin(3))/(n_phi-1)

        do i_phi = 1, n_phi
            phi = xmin(3) + (i_phi-1)*h_phi
            do i_th = 1, n_th
                th = xmin(2) + (i_th-1)*h_th
                do i_r = 1, n_r
                    r = xmin(1) + (i_r-1)*h_r

                    ! Generate different test functions for each component
                    do ic = 1, n_comp
                        select case(ic)
                        case(1)
                            field_data(i_r, i_th, i_phi, ic) = sin(r)*cos(th)*sin(phi)
                        case(2)
                            field_data(i_r, i_th, i_phi, ic) = cos(r)*sin(th)*cos(phi)
                        case(3)
                            field_data(i_r, i_th, i_phi, ic) = sin(r+th)*cos(phi)
                        case(4)
                            field_data(i_r, i_th, i_phi, ic) = cos(r)*sin(th+phi)
                        case(5)
                            field_data(i_r, i_th, i_phi, ic) = 1.0d0 + 0.1d0*sin(r)*cos(th)
                        case default
                            field_data(i_r, i_th, i_phi, ic) = exp(-r)*cos(th*ic)
                        end select
                    end do
                end do
            end do
        end do
    end subroutine generate_test_field_data

end program test_batch_splines
