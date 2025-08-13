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
    
    print *, "Testing cylindrical to cartesian transformation with Jacobian..."
    
    ! Given: A point in cylindrical coordinates (R, phi, Z)
    ! When: We transform to Cartesian with Jacobian
    ! Then: Both coordinates and Jacobian should be mathematically correct
    
    ! Test case 1: Standard point
    xcyl = [2.0_dp, pi/4.0_dp, 1.5_dp]  ! R=2, phi=45°, Z=1.5
    
    call transform_cyl_to_cart(xcyl, xcart, jacobian)
    
    ! Expected Cartesian coordinates
    ! x = R*cos(phi) = 2*cos(π/4) = 2*√2/2 = √2
    ! y = R*sin(phi) = 2*sin(π/4) = 2*√2/2 = √2
    ! z = Z = 1.5
    if (abs(xcart(1) - sqrt(2.0_dp)) > tol) then
      print *, "ERROR: x-coordinate incorrect"
      print *, "Expected:", sqrt(2.0_dp), "Got:", xcart(1)
      errors = errors + 1
    end if
    
    if (abs(xcart(2) - sqrt(2.0_dp)) > tol) then
      print *, "ERROR: y-coordinate incorrect"
      print *, "Expected:", sqrt(2.0_dp), "Got:", xcart(2)
      errors = errors + 1
    end if
    
    if (abs(xcart(3) - 1.5_dp) > tol) then
      print *, "ERROR: z-coordinate incorrect"
      print *, "Expected: 1.5, Got:", xcart(3)
      errors = errors + 1
    end if
    
    ! Test Jacobian components
    ! dx/dR = cos(phi)
    if (abs(jacobian(1,1) - cos(pi/4.0_dp)) > tol) then
      print *, "ERROR: dx/dR incorrect"
      print *, "Expected:", cos(pi/4.0_dp), "Got:", jacobian(1,1)
      errors = errors + 1
    end if
    
    ! dx/dphi = -R*sin(phi)
    if (abs(jacobian(1,2) - (-2.0_dp*sin(pi/4.0_dp))) > tol) then
      print *, "ERROR: dx/dphi incorrect"
      print *, "Expected:", -2.0_dp*sin(pi/4.0_dp), "Got:", jacobian(1,2)
      errors = errors + 1
    end if
    
    ! dy/dR = sin(phi)
    if (abs(jacobian(2,1) - sin(pi/4.0_dp)) > tol) then
      print *, "ERROR: dy/dR incorrect"
      print *, "Expected:", sin(pi/4.0_dp), "Got:", jacobian(2,1)
      errors = errors + 1
    end if
    
    ! dy/dphi = R*cos(phi)
    if (abs(jacobian(2,2) - (2.0_dp*cos(pi/4.0_dp))) > tol) then
      print *, "ERROR: dy/dphi incorrect"
      print *, "Expected:", 2.0_dp*cos(pi/4.0_dp), "Got:", jacobian(2,2)
      errors = errors + 1
    end if
    
    ! dz/dZ = 1
    if (abs(jacobian(3,3) - 1.0_dp) > tol) then
      print *, "ERROR: dz/dZ should be 1"
      print *, "Got:", jacobian(3,3)
      errors = errors + 1
    end if
    
    ! Test edge case: R=0 (should not crash, though mathematically singular)
    xcyl = [0.0_dp, pi/2.0_dp, 0.0_dp]
    call transform_cyl_to_cart(xcyl, xcart, jacobian)
    
    if (abs(xcart(1)) > tol .or. abs(xcart(2)) > tol .or. abs(xcart(3)) > tol) then
      print *, "ERROR: Origin should transform to origin"
      print *, "Got:", xcart
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
    
    ! Note: The error cases (invalid transform names) cause error stops,
    ! so we can't test them directly without causing the program to terminate.
    ! However, we've tested the valid paths which exercises the select case logic.
    
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
