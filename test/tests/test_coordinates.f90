program test_coordinates

use, intrinsic :: iso_fortran_env, only: dp => real64
use simple_coordinates, only: transform_i, get_transform, transform_cyl_to_cart

implicit none

integer :: errors

errors = 0

! Test coordinate transformations and error handling
call test_cyl_to_cart_jacobian(errors)
call test_get_transform_error_handling(errors)

if (errors == 0) then
  print *, "All coordinates tests passed!"
else
  print *, "ERROR: ", errors, " test(s) failed!"
  stop 1
end if

contains

  subroutine test_cyl_to_cart_jacobian(errors)
    integer, intent(inout) :: errors
    real(dp) :: xcyl(3), xcart(3), jacobian(3,3)
    real(dp), parameter :: tol = 1.0e-12_dp
    real(dp), parameter :: pi = 3.14159265358979323846_dp
    real(dp) :: R, phi, Z, cos_phi, sin_phi
    integer :: i
    
    ! Array of test angles: 0°, 30°, 60°, 90°, 120°, 150°, 180°, 270°, -45°, 2π
    real(dp), parameter :: test_angles(10) = [0.0_dp, pi/6.0_dp, pi/3.0_dp, pi/2.0_dp, &
                                             2.0_dp*pi/3.0_dp, 5.0_dp*pi/6.0_dp, pi, &
                                             3.0_dp*pi/2.0_dp, -pi/4.0_dp, 2.0_dp*pi]
    real(dp), parameter :: test_radii(3) = [0.5_dp, 1.0_dp, 3.0_dp]
    real(dp), parameter :: test_heights(3) = [-2.0_dp, 0.0_dp, 1.0_dp]
    
    print *, "Testing cylindrical to cartesian transformation with Jacobian..."
    
    ! Given: A point in cylindrical coordinates (R, phi, Z)
    ! When: We transform to Cartesian with Jacobian
    ! Then: Both coordinates and Jacobian should be mathematically correct
    
    ! Test case 1: Standard point (45 degrees)
    xcyl = [2.0_dp, pi/4.0_dp, 1.5_dp]  ! R=2, phi=45°, Z=1.5
    
    call transform_cyl_to_cart(xcyl, xcart, jacobian)
    
    ! Expected Cartesian coordinates
    ! x = R*cos(phi) = 2*cos(π/4) = 2*√2/2 = √2
    ! y = R*sin(phi) = 2*sin(π/4) = 2*√2/2 = √2
    ! z = Z = 1.5
    if (abs(xcart(1) - sqrt(2.0_dp)) > tol) then
      print *, "ERROR: x-coordinate incorrect for 45° test"
      print *, "Expected:", sqrt(2.0_dp), "Got:", xcart(1)
      errors = errors + 1
    end if
    
    if (abs(xcart(2) - sqrt(2.0_dp)) > tol) then
      print *, "ERROR: y-coordinate incorrect for 45° test"
      print *, "Expected:", sqrt(2.0_dp), "Got:", xcart(2)
      errors = errors + 1
    end if
    
    if (abs(xcart(3) - 1.5_dp) > tol) then
      print *, "ERROR: z-coordinate incorrect for 45° test"
      print *, "Expected: 1.5, Got:", xcart(3)
      errors = errors + 1
    end if
    
    ! Test Jacobian components for 45° case
    ! dx/dR = cos(phi)
    if (abs(jacobian(1,1) - cos(pi/4.0_dp)) > tol) then
      print *, "ERROR: dx/dR incorrect for 45° test"
      print *, "Expected:", cos(pi/4.0_dp), "Got:", jacobian(1,1)
      errors = errors + 1
    end if
    
    ! dx/dphi = -R*sin(phi)
    if (abs(jacobian(1,2) - (-2.0_dp*sin(pi/4.0_dp))) > tol) then
      print *, "ERROR: dx/dphi incorrect for 45° test"
      print *, "Expected:", -2.0_dp*sin(pi/4.0_dp), "Got:", jacobian(1,2)
      errors = errors + 1
    end if
    
    ! dy/dR = sin(phi)
    if (abs(jacobian(2,1) - sin(pi/4.0_dp)) > tol) then
      print *, "ERROR: dy/dR incorrect for 45° test"
      print *, "Expected:", sin(pi/4.0_dp), "Got:", jacobian(2,1)
      errors = errors + 1
    end if
    
    ! dy/dphi = R*cos(phi)
    if (abs(jacobian(2,2) - (2.0_dp*cos(pi/4.0_dp))) > tol) then
      print *, "ERROR: dy/dphi incorrect for 45° test"
      print *, "Expected:", 2.0_dp*cos(pi/4.0_dp), "Got:", jacobian(2,2)
      errors = errors + 1
    end if
    
    ! dz/dZ = 1
    if (abs(jacobian(3,3) - 1.0_dp) > tol) then
      print *, "ERROR: dz/dZ should be 1 for 45° test"
      print *, "Got:", jacobian(3,3)
      errors = errors + 1
    end if
    
    ! Test case 2: Multiple angles including boundary cases
    ! Array of test angles: 0°, 30°, 60°, 90°, 120°, 150°, 180°, 270°, -45°, 2π
    
    ! Test multiple combinations of parameters
    do i = 1, size(test_angles)
      R = 1.5_dp
      phi = test_angles(i) 
      Z = 0.5_dp
      xcyl = [R, phi, Z]
      
      call transform_cyl_to_cart(xcyl, xcart, jacobian)
      
      cos_phi = cos(phi)
      sin_phi = sin(phi)
      
      ! Verify coordinates
      if (abs(xcart(1) - R*cos_phi) > tol) then
        print *, "ERROR: x-coordinate incorrect for angle", phi
        print *, "Expected:", R*cos_phi, "Got:", xcart(1)
        errors = errors + 1
      end if
      
      if (abs(xcart(2) - R*sin_phi) > tol) then
        print *, "ERROR: y-coordinate incorrect for angle", phi
        print *, "Expected:", R*sin_phi, "Got:", xcart(2)
        errors = errors + 1
      end if
      
      if (abs(xcart(3) - Z) > tol) then
        print *, "ERROR: z-coordinate incorrect for angle", phi
        print *, "Expected:", Z, "Got:", xcart(3)
        errors = errors + 1
      end if
      
      ! Verify Jacobian matrix elements
      if (abs(jacobian(1,1) - cos_phi) > tol) then
        print *, "ERROR: dx/dR incorrect for angle", phi
        print *, "Expected:", cos_phi, "Got:", jacobian(1,1)
        errors = errors + 1
      end if
      
      if (abs(jacobian(1,2) - (-R*sin_phi)) > tol) then
        print *, "ERROR: dx/dphi incorrect for angle", phi
        print *, "Expected:", -R*sin_phi, "Got:", jacobian(1,2)
        errors = errors + 1
      end if
      
      if (abs(jacobian(2,1) - sin_phi) > tol) then
        print *, "ERROR: dy/dR incorrect for angle", phi
        print *, "Expected:", sin_phi, "Got:", jacobian(2,1)
        errors = errors + 1
      end if
      
      if (abs(jacobian(2,2) - R*cos_phi) > tol) then
        print *, "ERROR: dy/dphi incorrect for angle", phi
        print *, "Expected:", R*cos_phi, "Got:", jacobian(2,2)
        errors = errors + 1
      end if
    end do
    
    ! Test case 3: Different radii and heights
    do i = 1, size(test_radii)
      R = test_radii(i)
      phi = pi/3.0_dp  ! 60 degrees
      Z = test_heights(i)
      xcyl = [R, phi, Z]
      
      call transform_cyl_to_cart(xcyl, xcart, jacobian)
      
      cos_phi = cos(phi)
      sin_phi = sin(phi)
      
      ! Verify transformation
      if (abs(xcart(1) - R*cos_phi) > tol) then
        print *, "ERROR: x-coordinate incorrect for R=", R
        errors = errors + 1
      end if
      
      if (abs(xcart(2) - R*sin_phi) > tol) then
        print *, "ERROR: y-coordinate incorrect for R=", R
        errors = errors + 1
      end if
      
      if (abs(xcart(3) - Z) > tol) then
        print *, "ERROR: z-coordinate incorrect for Z=", Z
        errors = errors + 1
      end if
      
      ! Verify zero elements remain zero
      if (abs(jacobian(1,3)) > tol .or. abs(jacobian(2,3)) > tol .or. &
          abs(jacobian(3,1)) > tol .or. abs(jacobian(3,2)) > tol) then
        print *, "ERROR: Jacobian zero elements are non-zero for R=", R
        errors = errors + 1
      end if
      
      if (abs(jacobian(3,3) - 1.0_dp) > tol) then
        print *, "ERROR: dz/dZ should be 1 for all cases, R=", R
        errors = errors + 1
      end if
    end do
    
    ! Test case 4: Edge case R=0 (mathematically singular but should not crash)
    xcyl = [0.0_dp, pi/2.0_dp, 0.0_dp]
    call transform_cyl_to_cart(xcyl, xcart, jacobian)
    
    if (abs(xcart(1)) > tol .or. abs(xcart(2)) > tol .or. abs(xcart(3)) > tol) then
      print *, "ERROR: Origin should transform to origin"
      print *, "Got:", xcart
      errors = errors + 1
    end if
    
    ! For R=0, Jacobian should have specific structure
    if (abs(jacobian(1,1) - cos(pi/2.0_dp)) > tol) then  ! cos(π/2) = 0
      print *, "ERROR: dx/dR should be cos(phi) even at R=0"
      errors = errors + 1
    end if
    
    if (abs(jacobian(2,1) - sin(pi/2.0_dp)) > tol) then  ! sin(π/2) = 1
      print *, "ERROR: dy/dR should be sin(phi) even at R=0"
      errors = errors + 1
    end if
    
    ! Test case 5: Negative angles
    xcyl = [1.0_dp, -pi/4.0_dp, 0.0_dp]
    call transform_cyl_to_cart(xcyl, xcart, jacobian)
    
    ! For negative angle, results should be consistent with trigonometry
    if (abs(xcart(1) - cos(-pi/4.0_dp)) > tol) then
      print *, "ERROR: x-coordinate incorrect for negative angle"
      errors = errors + 1
    end if
    
    if (abs(xcart(2) - sin(-pi/4.0_dp)) > tol) then
      print *, "ERROR: y-coordinate incorrect for negative angle"
      errors = errors + 1
    end if
    
    if (errors == 0) then
      print *, "  Cylindrical to Cartesian Jacobian test PASSED"
    end if
    
  end subroutine test_cyl_to_cart_jacobian
  
  subroutine test_get_transform_error_handling(errors)
    integer, intent(inout) :: errors
    procedure(transform_i), pointer :: transform
    
    print *, "Testing get_transform error handling (covers lines 102-131)..."
    
    ! Given: Invalid coordinate system names
    ! When: We call get_transform with invalid names
    ! Then: It should handle errors appropriately
    
    ! Test 1: Valid 'cyl' to 'cart' (should work)
    transform => get_transform('cyl', 'cart')
    if (.not. associated(transform)) then
      print *, "ERROR: Valid 'cyl' to 'cart' transform should be associated"
      errors = errors + 1
    end if
    
    ! Test 2: Valid 'vmec' to 'cart' (should work)
    transform => get_transform('vmec', 'cart')
    if (.not. associated(transform)) then
      print *, "ERROR: Valid 'vmec' to 'cart' transform should be associated"
      errors = errors + 1
    end if
    
    ! Test 3: Valid 'vmec' to 'cyl' (should work)
    transform => get_transform('vmec', 'cyl')
    if (.not. associated(transform)) then
      print *, "ERROR: Valid 'vmec' to 'cyl' transform should be associated"
      errors = errors + 1
    end if
    
    ! Note: The error cases (invalid transform names) cause error stops via handle_transform_error,
    ! so we can't test them directly without causing the program to terminate.
    ! This is a design limitation - error handling uses 'error stop' which is not recoverable.
    ! Alternative approaches would require either:
    ! 1. Changing error handling to use exceptions or return codes (breaking change)
    ! 2. Using external processes or special test frameworks (beyond scope)
    ! 3. Mock/stub testing (requires significant refactoring)
    ! Therefore, we test all valid paths and document that error paths cannot be unit tested.
    
    ! Test 4: Only test cyl->cart transform which doesn't need VMEC data
    transform => get_transform('cyl', 'cart')
    if (associated(transform)) then
      block
        real(dp) :: xcyl(3), xcart(3)
        xcyl = [1.0_dp, 0.5_dp, 1.0_dp]  ! R=1, phi=0.5, Z=1
        
        ! This should not crash and should work without VMEC data
        call transform(xcyl, xcart)
        
        ! Basic sanity check - coordinates should be finite
        if (.not. (xcart(1) > -1.0e10_dp .and. xcart(1) < 1.0e10_dp)) then
          print *, "ERROR: Transform produced non-finite x-coordinate"
          errors = errors + 1
        end if
        
        if (.not. (xcart(2) > -1.0e10_dp .and. xcart(2) < 1.0e10_dp)) then
          print *, "ERROR: Transform produced non-finite y-coordinate"
          errors = errors + 1
        end if
        
        if (.not. (xcart(3) > -1.0e10_dp .and. xcart(3) < 1.0e10_dp)) then
          print *, "ERROR: Transform produced non-finite z-coordinate"
          errors = errors + 1
        end if
      end block
    end if
    
    if (errors == 0) then
      print *, "  get_transform error handling test PASSED"
    end if
    
  end subroutine test_get_transform_error_handling

end program test_coordinates
