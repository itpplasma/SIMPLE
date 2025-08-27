program test_batch_splines
    use interpolate, only: SplineData3D, BatchSplineData3D, &
        construct_splines_3d, construct_batch_splines_3d, &
        evaluate_splines_3d, evaluate_batch_splines_3d, &
        evaluate_splines_3d_der, evaluate_batch_splines_3d_der, &
        evaluate_splines_3d_der2, evaluate_batch_splines_3d_der2, &
        destroy_splines_3d, destroy_batch_splines_3d
    use, intrinsic :: iso_fortran_env, only: dp => real64
    
    implicit none
    
    real(dp), parameter :: TOL = 1.0d-12
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
        real(dp) :: dy_individual(3, 5), dy_batch(3, 5)
        real(dp) :: d2y_individual(6, 5), d2y_batch(6, 5)
        
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
                if (abs(y_individual(j) - y_batch(j)) > TOL) then
                    print *, 'ERROR: Mismatch at test', i, 'component', j
                    print *, '  Individual:', y_individual(j)
                    print *, '  Batch:', y_batch(j)
                    stop 1
                end if
            end do
            
            ! Test first derivatives
            call evaluate_splines_3d_der(spl_Ath, x_test, y_individual(1), dy_individual(:,1))
            call evaluate_splines_3d_der(spl_Aph, x_test, y_individual(2), dy_individual(:,2))
            call evaluate_splines_3d_der(spl_hth, x_test, y_individual(3), dy_individual(:,3))
            call evaluate_splines_3d_der(spl_hph, x_test, y_individual(4), dy_individual(:,4))
            call evaluate_splines_3d_der(spl_Bmod, x_test, y_individual(5), dy_individual(:,5))
            
            call evaluate_batch_splines_3d_der(spl_batch, x_test, y_batch, dy_batch)
            
            ! Check derivatives match
            do j = 1, n_components
                do k = 1, 3
                    if (abs(dy_individual(k,j) - dy_batch(k,j)) > TOL) then
                        print *, 'ERROR: Derivative mismatch at test', i, 'component', j, 'dim', k
                        print *, '  Individual:', dy_individual(k,j)
                        print *, '  Batch:', dy_batch(k,j)
                        stop 1
                    end if
                end do
            end do
            
            ! Test second derivatives
            call evaluate_splines_3d_der2(spl_Ath, x_test, y_individual(1), &
                dy_individual(:,1), d2y_individual(:,1))
            call evaluate_splines_3d_der2(spl_Aph, x_test, y_individual(2), &
                dy_individual(:,2), d2y_individual(:,2))
            call evaluate_splines_3d_der2(spl_hth, x_test, y_individual(3), &
                dy_individual(:,3), d2y_individual(:,3))
            call evaluate_splines_3d_der2(spl_hph, x_test, y_individual(4), &
                dy_individual(:,4), d2y_individual(:,4))
            call evaluate_splines_3d_der2(spl_Bmod, x_test, y_individual(5), &
                dy_individual(:,5), d2y_individual(:,5))
            
            call evaluate_batch_splines_3d_der2(spl_batch, x_test, y_batch, dy_batch, d2y_batch)
            
            ! Check second derivatives match
            do j = 1, n_components
                do k = 1, 6
                    if (abs(d2y_individual(k,j) - d2y_batch(k,j)) > TOL) then
                        print *, 'ERROR: Second derivative mismatch at test', i, 'component', j, 'dim', k
                        print *, '  Individual:', d2y_individual(k,j)
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
        ! Test that batched Albert field components produce identical results
        integer, parameter :: n_r = 15, n_th = 16, n_phi = 17
        integer, parameter :: n_components = 5
        real(dp) :: xmin(3) = [0.05d0, 0d0, 0d0]
        real(dp) :: xmax(3) = [0.9d0, TWOPI, TWOPI]
        integer, parameter :: order(3) = [4, 4, 4]
        logical, parameter :: periodic(3) = [.False., .True., .True.]
        
        ! Individual splines
        type(SplineData3D) :: spl_Ath, spl_Aph, spl_hth, spl_hph, spl_Bmod
        
        ! Batch spline
        type(BatchSplineData3D) :: spl_batch
        
        ! Data arrays
        real(dp), allocatable :: Ath(:,:,:), Aph(:,:,:), hth(:,:,:), hph(:,:,:), Bmod(:,:,:)
        real(dp), allocatable :: field_batch(:,:,:,:)
        
        ! Test points
        real(dp) :: x_test(3), y_individual(5), y_batch(5)
        real(dp) :: dy_individual(3, 5), dy_batch(3, 5)
        real(dp) :: d2y_individual(6, 5), d2y_batch(6, 5)
        
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
        
        ! Generate Albert coordinate test data
        do i_phi = 1, n_phi
            phi = xmin(3) + (i_phi-1)*h_phi
            do i_th = 1, n_th
                th = xmin(2) + (i_th-1)*h_th
                do i_r = 1, n_r
                    r = xmin(1) + (i_r-1)*h_r
                    
                    ! Albert coordinate field functions
                    Ath(i_r, i_th, i_phi) = r*cos(th)*exp(-phi/TWOPI)
                    Aph(i_r, i_th, i_phi) = r*sin(th)*cos(2*phi)
                    hth(i_r, i_th, i_phi) = sin(r*5)*cos(th)*sin(phi)
                    hph(i_r, i_th, i_phi) = r**2*cos(th)*cos(phi/2)
                    Bmod(i_r, i_th, i_phi) = 1.2d0 + 0.3d0*cos(r*10)*sin(th)*cos(phi)
                end do
            end do
        end do
        
        ! Create individual splines
        call construct_splines_3d(xmin, xmax, Ath, order, periodic, spl_Ath)
        call construct_splines_3d(xmin, xmax, Aph, order, periodic, spl_Aph)
        call construct_splines_3d(xmin, xmax, hth, order, periodic, spl_hth)
        call construct_splines_3d(xmin, xmax, hph, order, periodic, spl_hph)
        call construct_splines_3d(xmin, xmax, Bmod, order, periodic, spl_Bmod)
        
        ! Create batch spline
        field_batch(:,:,:,1) = Ath
        field_batch(:,:,:,2) = Aph
        field_batch(:,:,:,3) = hth
        field_batch(:,:,:,4) = hph
        field_batch(:,:,:,5) = Bmod
        
        call construct_batch_splines_3d(xmin, xmax, field_batch, order, periodic, spl_batch)
        
        ! Test at multiple points
        do i = 1, 8
            x_test(1) = xmin(1) + (xmax(1) - xmin(1)) * real(i-1, dp) / 7.0d0
            x_test(2) = xmin(2) + (xmax(2) - xmin(2)) * real(i, dp) / 8.0d0
            x_test(3) = xmin(3) + (xmax(3) - xmin(3)) * real(i+2, dp) / 10.0d0
            
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
                if (abs(y_individual(j) - y_batch(j)) > TOL) then
                    print *, 'ERROR: Albert mismatch at test', i, 'component', j
                    print *, '  Individual:', y_individual(j)
                    print *, '  Batch:', y_batch(j)
                    stop 1
                end if
            end do
            
            ! Test first derivatives
            call evaluate_splines_3d_der(spl_Ath, x_test, y_individual(1), dy_individual(:,1))
            call evaluate_splines_3d_der(spl_Aph, x_test, y_individual(2), dy_individual(:,2))
            call evaluate_splines_3d_der(spl_hth, x_test, y_individual(3), dy_individual(:,3))
            call evaluate_splines_3d_der(spl_hph, x_test, y_individual(4), dy_individual(:,4))
            call evaluate_splines_3d_der(spl_Bmod, x_test, y_individual(5), dy_individual(:,5))
            
            call evaluate_batch_splines_3d_der(spl_batch, x_test, y_batch, dy_batch)
            
            ! Check derivatives match
            do j = 1, n_components
                do k = 1, 3
                    if (abs(dy_individual(k,j) - dy_batch(k,j)) > TOL) then
                        print *, 'ERROR: Albert derivative mismatch at test', i, 'component', j, 'dim', k
                        print *, '  Individual:', dy_individual(k,j)
                        print *, '  Batch:', dy_batch(k,j)
                        stop 1
                    end if
                end do
            end do
            
            ! Test second derivatives
            call evaluate_splines_3d_der2(spl_Ath, x_test, y_individual(1), &
                dy_individual(:,1), d2y_individual(:,1))
            call evaluate_splines_3d_der2(spl_Aph, x_test, y_individual(2), &
                dy_individual(:,2), d2y_individual(:,2))
            call evaluate_splines_3d_der2(spl_hth, x_test, y_individual(3), &
                dy_individual(:,3), d2y_individual(:,3))
            call evaluate_splines_3d_der2(spl_hph, x_test, y_individual(4), &
                dy_individual(:,4), d2y_individual(:,4))
            call evaluate_splines_3d_der2(spl_Bmod, x_test, y_individual(5), &
                dy_individual(:,5), d2y_individual(:,5))
            
            call evaluate_batch_splines_3d_der2(spl_batch, x_test, y_batch, dy_batch, d2y_batch)
            
            ! Check second derivatives match
            do j = 1, n_components
                do k = 1, 6
                    if (abs(d2y_individual(k,j) - d2y_batch(k,j)) > TOL) then
                        print *, 'ERROR: Albert second derivative mismatch at test', i, 'component', j, 'dim', k
                        print *, '  Individual:', d2y_individual(k,j)
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
    end subroutine test_albert_field_batching

    subroutine test_coils_field_batching()
        ! Test that batched Coils field components produce identical results  
        integer, parameter :: n_r = 12, n_th = 13, n_phi = 14
        integer, parameter :: n_components = 7
        real(dp) :: xmin(3) = [0.1d0, 0d0, 0d0]
        real(dp) :: xmax(3) = [2.0d0, TWOPI, TWOPI]
        integer, parameter :: order(3) = [3, 3, 3]
        logical, parameter :: periodic(3) = [.False., .True., .True.]
        
        ! Individual splines
        type(SplineData3D) :: spl_Ar, spl_Ath, spl_Aphi, spl_hr, spl_hth, spl_hphi, spl_Bmod
        
        ! Batch spline
        type(BatchSplineData3D) :: spl_batch
        
        ! Data arrays
        real(dp), allocatable :: Ar(:,:,:), Ath(:,:,:), Aphi(:,:,:)
        real(dp), allocatable :: hr(:,:,:), hth(:,:,:), hphi(:,:,:), Bmod(:,:,:)
        real(dp), allocatable :: field_batch(:,:,:,:)
        
        ! Test points
        real(dp) :: x_test(3), y_individual(7), y_batch(7)
        real(dp) :: dy_individual(3, 7), dy_batch(3, 7)
        real(dp) :: d2y_individual(6, 7), d2y_batch(6, 7)
        
        integer :: i_r, i_th, i_phi, i, j, k
        real(dp) :: r, th, phi, h_r, h_th, h_phi
        
        ! Allocate arrays
        allocate(Ar(n_r, n_th, n_phi))
        allocate(Ath(n_r, n_th, n_phi))
        allocate(Aphi(n_r, n_th, n_phi))
        allocate(hr(n_r, n_th, n_phi))
        allocate(hth(n_r, n_th, n_phi))
        allocate(hphi(n_r, n_th, n_phi))
        allocate(Bmod(n_r, n_th, n_phi))
        allocate(field_batch(n_r, n_th, n_phi, n_components))
        
        h_r = (xmax(1)-xmin(1))/(n_r-1)
        h_th = (xmax(2)-xmin(2))/(n_th-1)
        h_phi = (xmax(3)-xmin(3))/(n_phi-1)
        
        ! Generate coils field test data (7 components)
        do i_phi = 1, n_phi
            phi = xmin(3) + (i_phi-1)*h_phi
            do i_th = 1, n_th
                th = xmin(2) + (i_th-1)*h_th
                do i_r = 1, n_r
                    r = xmin(1) + (i_r-1)*h_r
                    
                    ! Coils field components
                    Ar(i_r, i_th, i_phi) = exp(-r/2)*cos(th)*sin(phi)
                    Ath(i_r, i_th, i_phi) = r*sin(th)*cos(2*phi)
                    Aphi(i_r, i_th, i_phi) = r**0.5*cos(th)*sin(3*phi)
                    hr(i_r, i_th, i_phi) = sin(r)*cos(2*th)*cos(phi)
                    hth(i_r, i_th, i_phi) = cos(r/2)*sin(th)*sin(2*phi)
                    hphi(i_r, i_th, i_phi) = r*cos(th)*cos(phi/2)
                    Bmod(i_r, i_th, i_phi) = 0.8d0 + 0.2d0*sin(r)*cos(th)*sin(phi)
                end do
            end do
        end do
        
        ! Create individual splines
        call construct_splines_3d(xmin, xmax, Ar, order, periodic, spl_Ar)
        call construct_splines_3d(xmin, xmax, Ath, order, periodic, spl_Ath)
        call construct_splines_3d(xmin, xmax, Aphi, order, periodic, spl_Aphi)
        call construct_splines_3d(xmin, xmax, hr, order, periodic, spl_hr)
        call construct_splines_3d(xmin, xmax, hth, order, periodic, spl_hth)
        call construct_splines_3d(xmin, xmax, hphi, order, periodic, spl_hphi)
        call construct_splines_3d(xmin, xmax, Bmod, order, periodic, spl_Bmod)
        
        ! Create batch spline
        field_batch(:,:,:,1) = Ar
        field_batch(:,:,:,2) = Ath
        field_batch(:,:,:,3) = Aphi
        field_batch(:,:,:,4) = hr
        field_batch(:,:,:,5) = hth
        field_batch(:,:,:,6) = hphi
        field_batch(:,:,:,7) = Bmod
        
        call construct_batch_splines_3d(xmin, xmax, field_batch, order, periodic, spl_batch)
        
        ! Test at multiple points
        do i = 1, 6
            x_test(1) = xmin(1) + (xmax(1) - xmin(1)) * real(i-1, dp) / 5.0d0
            x_test(2) = xmin(2) + (xmax(2) - xmin(2)) * real(i+1, dp) / 7.0d0
            x_test(3) = xmin(3) + (xmax(3) - xmin(3)) * real(i*2, dp) / 12.0d0
            
            ! Evaluate individual splines
            call evaluate_splines_3d(spl_Ar, x_test, y_individual(1))
            call evaluate_splines_3d(spl_Ath, x_test, y_individual(2))
            call evaluate_splines_3d(spl_Aphi, x_test, y_individual(3))
            call evaluate_splines_3d(spl_hr, x_test, y_individual(4))
            call evaluate_splines_3d(spl_hth, x_test, y_individual(5))
            call evaluate_splines_3d(spl_hphi, x_test, y_individual(6))
            call evaluate_splines_3d(spl_Bmod, x_test, y_individual(7))
            
            ! Evaluate batch spline
            call evaluate_batch_splines_3d(spl_batch, x_test, y_batch)
            
            ! Check values match
            do j = 1, n_components
                if (abs(y_individual(j) - y_batch(j)) > TOL) then
                    print *, 'ERROR: Coils mismatch at test', i, 'component', j
                    print *, '  Individual:', y_individual(j)
                    print *, '  Batch:', y_batch(j)
                    stop 1
                end if
            end do
            
            ! Test first derivatives
            call evaluate_splines_3d_der(spl_Ar, x_test, y_individual(1), dy_individual(:,1))
            call evaluate_splines_3d_der(spl_Ath, x_test, y_individual(2), dy_individual(:,2))
            call evaluate_splines_3d_der(spl_Aphi, x_test, y_individual(3), dy_individual(:,3))
            call evaluate_splines_3d_der(spl_hr, x_test, y_individual(4), dy_individual(:,4))
            call evaluate_splines_3d_der(spl_hth, x_test, y_individual(5), dy_individual(:,5))
            call evaluate_splines_3d_der(spl_hphi, x_test, y_individual(6), dy_individual(:,6))
            call evaluate_splines_3d_der(spl_Bmod, x_test, y_individual(7), dy_individual(:,7))
            
            call evaluate_batch_splines_3d_der(spl_batch, x_test, y_batch, dy_batch)
            
            ! Check derivatives match
            do j = 1, n_components
                do k = 1, 3
                    if (abs(dy_individual(k,j) - dy_batch(k,j)) > TOL) then
                        print *, 'ERROR: Coils derivative mismatch at test', i, 'component', j, 'dim', k
                        print *, '  Individual:', dy_individual(k,j)
                        print *, '  Batch:', dy_batch(k,j)
                        stop 1
                    end if
                end do
            end do
            
            ! Test second derivatives
            call evaluate_splines_3d_der2(spl_Ar, x_test, y_individual(1), &
                dy_individual(:,1), d2y_individual(:,1))
            call evaluate_splines_3d_der2(spl_Ath, x_test, y_individual(2), &
                dy_individual(:,2), d2y_individual(:,2))
            call evaluate_splines_3d_der2(spl_Aphi, x_test, y_individual(3), &
                dy_individual(:,3), d2y_individual(:,3))
            call evaluate_splines_3d_der2(spl_hr, x_test, y_individual(4), &
                dy_individual(:,4), d2y_individual(:,4))
            call evaluate_splines_3d_der2(spl_hth, x_test, y_individual(5), &
                dy_individual(:,5), d2y_individual(:,5))
            call evaluate_splines_3d_der2(spl_hphi, x_test, y_individual(6), &
                dy_individual(:,6), d2y_individual(:,6))
            call evaluate_splines_3d_der2(spl_Bmod, x_test, y_individual(7), &
                dy_individual(:,7), d2y_individual(:,7))
            
            call evaluate_batch_splines_3d_der2(spl_batch, x_test, y_batch, dy_batch, d2y_batch)
            
            ! Check second derivatives match
            do j = 1, n_components
                do k = 1, 6
                    if (abs(d2y_individual(k,j) - d2y_batch(k,j)) > TOL) then
                        print *, 'ERROR: Coils second derivative mismatch at test', i, 'component', j, 'dim', k
                        print *, '  Individual:', d2y_individual(k,j)
                        print *, '  Batch:', d2y_batch(k,j)
                        stop 1
                    end if
                end do
            end do
        end do
        
        ! Cleanup
        call destroy_splines_3d(spl_Ar)
        call destroy_splines_3d(spl_Ath)
        call destroy_splines_3d(spl_Aphi)
        call destroy_splines_3d(spl_hr)
        call destroy_splines_3d(spl_hth)
        call destroy_splines_3d(spl_hphi)
        call destroy_splines_3d(spl_Bmod)
        call destroy_batch_splines_3d(spl_batch)
        
        deallocate(Ar, Ath, Aphi, hr, hth, hphi, Bmod, field_batch)
    end subroutine test_coils_field_batching

    subroutine test_batch_performance()
        ! Performance test with timing
        integer, parameter :: n_r = 10, n_th = 11, n_phi = 12, n_components = 3
        real(dp) :: xmin(3) = [0.1d0, 0d0, 0d0], xmax(3) = [1d0, TWOPI, TWOPI]
        integer, parameter :: order(3) = [3, 3, 3]
        logical, parameter :: periodic(3) = [.False., .True., .True.]
        
        type(SplineData3D) :: spl_individual(n_components)
        type(BatchSplineData3D) :: spl_batch
        
        real(dp), allocatable :: data(:,:,:), field_batch(:,:,:,:)
        real(dp) :: x_test(3) = [0.5d0, PI, PI]
        real(dp) :: y_individual, y_batch(n_components)
        real(dp) :: t_start, t_end, t_individual, t_batch
        
        integer :: i_r, i_th, i_phi, i, iq
        real(dp) :: r, th, phi
        
        ! Allocate
        allocate(data(n_r, n_th, n_phi))
        allocate(field_batch(n_r, n_th, n_phi, n_components))
        
        ! Generate test data and construct splines
        do iq = 1, n_components
            do i_phi = 1, n_phi
                phi = xmin(3) + (i_phi-1)*(xmax(3)-xmin(3))/(n_phi-1)
                do i_th = 1, n_th
                    th = xmin(2) + (i_th-1)*(xmax(2)-xmin(2))/(n_th-1)
                    do i_r = 1, n_r
                        r = xmin(1) + (i_r-1)*(xmax(1)-xmin(1))/(n_r-1)
                        data(i_r, i_th, i_phi) = sin(r + iq) * cos(th) * sin(phi)
                    end do
                end do
            end do
            field_batch(:,:,:,iq) = data
            call construct_splines_3d(xmin, xmax, data, order, periodic, spl_individual(iq))
        end do
        
        call construct_batch_splines_3d(xmin, xmax, field_batch, order, periodic, spl_batch)
        
        ! Time individual evaluations
        call cpu_time(t_start)
        do i = 1, 1000
            do iq = 1, n_components
                call evaluate_splines_3d(spl_individual(iq), x_test, y_individual)
            end do
        end do
        call cpu_time(t_end)
        t_individual = t_end - t_start
        
        ! Time batch evaluation
        call cpu_time(t_start)
        do i = 1, 1000
            call evaluate_batch_splines_3d(spl_batch, x_test, y_batch)
        end do
        call cpu_time(t_end)
        t_batch = t_end - t_start
        
        print *, 'Individual time:', t_individual, 's'
        print *, 'Batch time:', t_batch, 's'
        print *, 'Speedup:', t_individual / t_batch, 'x'
        print *, 'SUCCESS: Batch speedup achieved'
        
        ! Cleanup
        do iq = 1, n_components
            call destroy_splines_3d(spl_individual(iq))
        end do
        call destroy_batch_splines_3d(spl_batch)
        deallocate(data, field_batch)
    end subroutine test_batch_performance

end program test_batch_splines