program test_coordinates
  use simple_coordinates
  implicit none
  
  integer :: errors
  
  errors = 0
  
  ! Test cylindrical to cartesian transformation
  call test_cyl_to_cart_transform(errors)
  
  ! Test coordinate transformation function pointer system
  call test_transform_function_pointers(errors)
  
  if (errors == 0) then
    print *, "All coordinates module tests passed!"
  else
    print *, "ERROR: ", errors, " test(s) failed!"
    stop 1
  end if
  
contains

  subroutine test_cyl_to_cart_transform(errors)
    integer, intent(inout) :: errors
    real(dp) :: xfrom(3), xto(3), dxto_dxfrom(3,3)
    real(dp), parameter :: tolerance = 1.0d-14
    real(dp), parameter :: pi = 3.14159265358979d0
    
    print *, "Testing cylindrical to cartesian transformation..."
    
    ! Given: Cylindrical coordinates (r, phi, z)
    ! When: We transform to cartesian coordinates
    ! Then: The result should follow x = r*cos(phi), y = r*sin(phi), z = z
    
    ! Test case 1: Simple point on x-axis
    xfrom = [2.0_dp, 0.0_dp, 1.0_dp]  ! (r=2, phi=0, z=1)
    call transform_cyl_to_cart(xfrom, xto, dxto_dxfrom)
    
    ! Expected: (x=2, y=0, z=1)
    if (abs(xto(1) - 2.0_dp) > tolerance) then
      print *, "ERROR: Cyl to cart x-coordinate incorrect for phi=0"
      print *, "Expected: 2.0, Got:", xto(1)
      errors = errors + 1
    end if
    
    if (abs(xto(2)) > tolerance) then
      print *, "ERROR: Cyl to cart y-coordinate should be 0 for phi=0"
      print *, "Got:", xto(2)
      errors = errors + 1
    end if
    
    if (abs(xto(3) - 1.0_dp) > tolerance) then
      print *, "ERROR: Cyl to cart z-coordinate should be unchanged"
      print *, "Expected: 1.0, Got:", xto(3)
      errors = errors + 1
    end if
    
    ! Test case 2: Point on y-axis
    xfrom = [3.0_dp, pi/2.0_dp, 2.0_dp]  ! (r=3, phi=pi/2, z=2)
    call transform_cyl_to_cart(xfrom, xto)
    
    ! Expected: (x=0, y=3, z=2)
    if (abs(xto(1)) > tolerance) then
      print *, "ERROR: Cyl to cart x-coordinate should be 0 for phi=pi/2"
      print *, "Got:", xto(1)
      errors = errors + 1
    end if
    
    if (abs(xto(2) - 3.0_dp) > tolerance) then
      print *, "ERROR: Cyl to cart y-coordinate incorrect for phi=pi/2"
      print *, "Expected: 3.0, Got:", xto(2)
      errors = errors + 1
    end if
    
    ! Test case 3: Check Jacobian matrix
    xfrom = [1.0_dp, pi/4.0_dp, 0.0_dp]  ! (r=1, phi=45°, z=0)
    call transform_cyl_to_cart(xfrom, xto, dxto_dxfrom)
    
    ! Check Jacobian elements
    ! dx/dr = cos(phi)
    if (abs(dxto_dxfrom(1,1) - cos(pi/4.0_dp)) > tolerance) then
      print *, "ERROR: Jacobian dx/dr incorrect"
      print *, "Expected:", cos(pi/4.0_dp), "Got:", dxto_dxfrom(1,1)
      errors = errors + 1
    end if
    
    ! dx/dphi = -r*sin(phi)
    if (abs(dxto_dxfrom(1,2) - (-1.0_dp*sin(pi/4.0_dp))) > tolerance) then
      print *, "ERROR: Jacobian dx/dphi incorrect"
      print *, "Expected:", -sin(pi/4.0_dp), "Got:", dxto_dxfrom(1,2)
      errors = errors + 1
    end if
    
    ! dx/dz = 0
    if (abs(dxto_dxfrom(1,3)) > tolerance) then
      print *, "ERROR: Jacobian dx/dz should be 0"
      print *, "Got:", dxto_dxfrom(1,3)
      errors = errors + 1
    end if
    
    ! dz/dz = 1
    if (abs(dxto_dxfrom(3,3) - 1.0_dp) > tolerance) then
      print *, "ERROR: Jacobian dz/dz should be 1"
      print *, "Got:", dxto_dxfrom(3,3)
      errors = errors + 1
    end if
    
    ! Test invariant properties
    call test_cyl_to_cart_invariants(errors)
    
    if (errors == 0) then
      print *, "  Cylindrical to cartesian transformation test PASSED"
    end if
    
  end subroutine test_cyl_to_cart_transform
  
  subroutine test_cyl_to_cart_invariants(errors)
    integer, intent(inout) :: errors
    real(dp) :: xfrom(3), xto(3)
    real(dp), parameter :: tolerance = 1.0d-13
    real(dp) :: r_original, r_computed
    
    ! Given: Cylindrical coordinates
    ! When: We transform to cartesian and compute the radius
    ! Then: The radius should be preserved (r = sqrt(x^2 + y^2))
    
    xfrom = [2.5_dp, 1.2_dp, -0.8_dp]  ! Arbitrary point
    r_original = xfrom(1)
    
    call transform_cyl_to_cart(xfrom, xto)
    r_computed = sqrt(xto(1)**2 + xto(2)**2)
    
    if (abs(r_computed - r_original) > tolerance) then
      print *, "ERROR: Radius not preserved in cyl->cart transformation"
      print *, "Original r:", r_original, "Computed r:", r_computed
      errors = errors + 1
    end if
    
    ! Check that z-coordinate is preserved
    if (abs(xto(3) - xfrom(3)) > tolerance) then
      print *, "ERROR: Z-coordinate not preserved in cyl->cart transformation"
      errors = errors + 1
    end if
    
  end subroutine test_cyl_to_cart_invariants
  
  subroutine test_transform_function_pointers(errors)
    integer, intent(inout) :: errors
    procedure(transform_i), pointer :: transform_ptr
    
    print *, "Testing coordinate transformation function pointers..."
    
    ! Given: The get_transform function
    ! When: We request valid transformations
    ! Then: It should return appropriate function pointers
    
    ! Test valid transformation: cyl to cart
    transform_ptr => get_transform('cyl', 'cart')
    if (.not. associated(transform_ptr)) then
      print *, "ERROR: get_transform should return associated pointer for cyl->cart"
      errors = errors + 1
    else
      ! Test that the function pointer works
      call test_function_pointer_execution(transform_ptr, errors)
    end if
    
    if (errors == 0) then
      print *, "  Transform function pointers test PASSED"
    end if
    
  end subroutine test_transform_function_pointers
  
  subroutine test_function_pointer_execution(transform_ptr, errors)
    procedure(transform_i), pointer, intent(in) :: transform_ptr
    integer, intent(inout) :: errors
    real(dp) :: xfrom(3), xto(3)
    real(dp), parameter :: tolerance = 1.0d-14
    
    ! Given: A function pointer to a coordinate transformation
    ! When: We execute it with test coordinates
    ! Then: It should produce the expected result
    
    xfrom = [1.0_dp, 0.0_dp, 0.0_dp]  ! (r=1, phi=0, z=0)
    call transform_ptr(xfrom, xto)
    
    ! For cyl->cart with (1,0,0), we expect (1,0,0)
    if (abs(xto(1) - 1.0_dp) > tolerance .or. &
        abs(xto(2)) > tolerance .or. &
        abs(xto(3)) > tolerance) then
      print *, "ERROR: Function pointer execution failed"
      print *, "Expected: [1,0,0], Got:", xto
      errors = errors + 1
    end if
    
  end subroutine test_function_pointer_execution

end program test_coordinates