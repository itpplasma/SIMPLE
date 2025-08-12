program test_boozer_refactoring
    use simple, only: Tracer
    use simple_main, only: init_field
    use boozer_sub, only: normalize_angular_coordinates
    use boozer_coordinates_mod, only: n_theta_B, n_phi_B, h_theta_B, h_phi_B
    use new_vmec_stuff_mod, only: nper
    
    implicit none
    
    type(Tracer) :: norb
    integer :: i_theta, i_phi, test_failed
    double precision :: vartheta, varphi, dtheta, dphi
    double precision, parameter :: twopi = 2.d0*3.14159265358979d0
    double precision, parameter :: tol = 1.d-14
    
    test_failed = 0
    
    ! Initialize field (needed for module initialization)
    call init_field(norb, 'wout.nc', 5, 5, 5, -1)
    
    print *, 'Testing normalize_angular_coordinates subroutine'
    
    ! Test 1: Angle at origin
    vartheta = 0.d0
    varphi = 0.d0
    call normalize_angular_coordinates(vartheta, varphi, 32, 32, &
        twopi/31.d0, twopi/31.d0, i_theta, i_phi, dtheta, dphi)
    
    if (i_theta /= 1 .or. i_phi /= 1) then
        print *, 'FAILED: Test 1 - indices at origin'
        print *, 'Expected i_theta=1, i_phi=1'
        print *, 'Got i_theta=', i_theta, ', i_phi=', i_phi
        test_failed = test_failed + 1
    else
        print *, 'PASSED: Test 1 - indices at origin'
    end if
    
    if (abs(dtheta) > tol .or. abs(dphi) > tol) then
        print *, 'FAILED: Test 1 - offsets at origin'
        print *, 'Expected dtheta=0, dphi=0'
        print *, 'Got dtheta=', dtheta, ', dphi=', dphi
        test_failed = test_failed + 1
    else
        print *, 'PASSED: Test 1 - offsets at origin'
    end if
    
    ! Test 2: Angle at middle of grid (account for nper periodicity for phi)
    vartheta = twopi/2.d0
    varphi = (twopi/dble(nper))/2.d0
    call normalize_angular_coordinates(vartheta, varphi, 32, 32, &
        twopi/31.d0, (twopi/dble(nper))/31.d0, i_theta, i_phi, dtheta, dphi)
    
    if (i_theta /= 16 .or. i_phi /= 16) then
        print *, 'FAILED: Test 2 - indices at middle'
        print *, 'Expected i_theta=16, i_phi=16'
        print *, 'Got i_theta=', i_theta, ', i_phi=', i_phi
        print *, 'nper=', nper
        test_failed = test_failed + 1
    else
        print *, 'PASSED: Test 2 - indices at middle'
    end if
    
    ! Test 3: Angle wrapping (beyond 2*pi)
    vartheta = twopi + 0.1d0
    varphi = (twopi/dble(nper)) + 0.02d0
    call normalize_angular_coordinates(vartheta, varphi, 32, 32, &
        twopi/31.d0, (twopi/dble(nper))/31.d0, i_theta, i_phi, dtheta, dphi)
    
    if (i_theta < 1 .or. i_theta > 32 .or. i_phi < 1 .or. i_phi > 32) then
        print *, 'FAILED: Test 3 - wrapped angle indices out of bounds'
        print *, 'Got i_theta=', i_theta, ', i_phi=', i_phi
        test_failed = test_failed + 1
    else
        print *, 'PASSED: Test 3 - wrapped angle indices in bounds'
    end if
    
    ! Test 4: Negative angles
    vartheta = -0.1d0
    varphi = -0.02d0
    call normalize_angular_coordinates(vartheta, varphi, 32, 32, &
        twopi/31.d0, (twopi/dble(nper))/31.d0, i_theta, i_phi, dtheta, dphi)
    
    if (i_theta < 1 .or. i_theta > 32 .or. i_phi < 1 .or. i_phi > 32) then
        print *, 'FAILED: Test 4 - negative angle indices out of bounds'
        print *, 'Got i_theta=', i_theta, ', i_phi=', i_phi
        test_failed = test_failed + 1
    else
        print *, 'PASSED: Test 4 - negative angle indices in bounds'
    end if
    
    ! Summary
    if (test_failed == 0) then
        print *, '================================'
        print *, 'All tests PASSED'
        print *, '================================'
        stop 0
    else
        print *, '================================'
        print *, test_failed, ' tests FAILED'
        print *, '================================'
        stop 1
    end if
    
end program test_boozer_refactoring