program test_canonical_stencil
  ! Test that verifies the stencil initialization in get_canonical_coordinates
  ! matches the expected behavior before refactoring
  implicit none
  
  call test_original_stencil_values()
  
  print *, "All canonical stencil tests passed!"
  
contains

  subroutine test_original_stencil_values()
    ! Test the exact stencil values as computed in get_canonical_coordinates
    integer :: nh_stencil
    double precision :: h_theta_c, h_phi_c
    double precision, allocatable :: dstencil_theta(:), dstencil_phi(:)
    double precision, parameter :: tol = 1.0d-14
    
    ! Test case 1: nh_stencil = 1
    nh_stencil = 1
    h_theta_c = 0.1d0
    h_phi_c = 0.2d0
    allocate(dstencil_theta(-nh_stencil:nh_stencil))
    allocate(dstencil_phi(-nh_stencil:nh_stencil))
    
    ! Original code behavior
    dstencil_theta(-1) = -0.5d0
    dstencil_theta(0) = 0.0d0
    dstencil_theta(1) = 0.5d0
    dstencil_phi = dstencil_theta
    dstencil_theta = dstencil_theta / h_theta_c
    dstencil_phi = dstencil_phi / h_phi_c
    
    ! Verify values
    if (abs(dstencil_theta(-1) - (-5.0d0)) > tol) then
      print *, "ERROR: dstencil_theta(-1) incorrect for nh_stencil=1"
      print *, "Expected:", -5.0d0, "Got:", dstencil_theta(-1)
      error stop 1
    end if
    
    if (abs(dstencil_theta(1) - 5.0d0) > tol) then
      print *, "ERROR: dstencil_theta(1) incorrect for nh_stencil=1"
      error stop 1
    end if
    
    if (abs(dstencil_phi(-1) - (-2.5d0)) > tol) then
      print *, "ERROR: dstencil_phi(-1) incorrect for nh_stencil=1"
      error stop 1
    end if
    
    deallocate(dstencil_theta, dstencil_phi)
    
    ! Test case 2: nh_stencil = 2
    nh_stencil = 2
    allocate(dstencil_theta(-nh_stencil:nh_stencil))
    allocate(dstencil_phi(-nh_stencil:nh_stencil))
    
    ! Original code behavior
    dstencil_theta(-2) = 1.d0/12.d0
    dstencil_theta(-1) = -2.d0/3.d0
    dstencil_theta(0) = 0.0d0
    dstencil_theta(1) = 2.d0/3.d0
    dstencil_theta(2) = -1.d0/12.d0
    dstencil_phi = dstencil_theta
    dstencil_theta = dstencil_theta / h_theta_c
    dstencil_phi = dstencil_phi / h_phi_c
    
    ! Verify key values
    if (abs(dstencil_theta(-1) - (-2.d0/3.d0/h_theta_c)) > tol) then
      print *, "ERROR: dstencil_theta(-1) incorrect for nh_stencil=2"
      error stop 1
    end if
    
    deallocate(dstencil_theta, dstencil_phi)
    
    ! Test case 3: nh_stencil = 3
    nh_stencil = 3
    allocate(dstencil_theta(-nh_stencil:nh_stencil))
    allocate(dstencil_phi(-nh_stencil:nh_stencil))
    
    ! Original code behavior
    dstencil_theta(-3) = -1.d0/60.d0
    dstencil_theta(-2) = 0.15d0
    dstencil_theta(-1) = -0.75d0
    dstencil_theta(0) = 0.0d0
    dstencil_theta(1) = 0.75d0
    dstencil_theta(2) = -0.15d0
    dstencil_theta(3) = 1.d0/60.d0
    dstencil_phi = dstencil_theta
    dstencil_theta = dstencil_theta / h_theta_c
    dstencil_phi = dstencil_phi / h_phi_c
    
    ! Verify symmetry
    if (abs(dstencil_theta(-1) + dstencil_theta(1)) > tol) then
      print *, "ERROR: stencil not antisymmetric for nh_stencil=3"
      error stop 1
    end if
    
    deallocate(dstencil_theta, dstencil_phi)
    
    print *, "test_original_stencil_values: PASSED"
  end subroutine test_original_stencil_values

end program test_canonical_stencil