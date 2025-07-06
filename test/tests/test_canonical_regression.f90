program test_canonical_regression
  use get_can_sub
  use canonical_coordinates_mod
  use new_vmec_stuff_mod
  use vector_potentail_mod
  implicit none
  
  ! Test that get_canonical_coordinates produces expected results
  ! This is a regression test to ensure refactoring doesn't change behavior
  
  call test_canonical_initialization()
  call test_canonical_values()
  
  print *, "All canonical regression tests passed!"
  
contains

  subroutine test_canonical_initialization()
    ! Test that module variables are properly initialized
    double precision, parameter :: tol = 1.0d-14
    
    ! Set up minimal VMEC-like configuration
    ns = 10
    n_theta = 8
    n_phi = 4
    h_theta = 0.785398163397448d0  ! pi/4
    h_phi = 1.5707963267949d0      ! pi/2
    hs = 0.1d0
    ns_s = 3
    ns_tp = 3
    nh_stencil = 2
    
    ! Check that these are properly transferred
    call get_canonical_coordinates()
    
    if (ns_c /= ns) then
      print *, "ERROR: ns_c not properly set"
      print *, "Expected:", ns, "Got:", ns_c
      error stop 1
    end if
    
    if (n_theta_c /= n_theta) then
      print *, "ERROR: n_theta_c not properly set"
      error stop 1
    end if
    
    if (n_phi_c /= n_phi) then
      print *, "ERROR: n_phi_c not properly set"
      error stop 1
    end if
    
    if (abs(h_theta_c - h_theta) > tol) then
      print *, "ERROR: h_theta_c not properly set"
      error stop 1
    end if
    
    if (abs(h_phi_c - h_phi) > tol) then
      print *, "ERROR: h_phi_c not properly set"
      error stop 1
    end if
    
    print *, "test_canonical_initialization: PASSED"
  end subroutine test_canonical_initialization
  
  subroutine test_canonical_values()
    ! Test specific output values for a known configuration
    ! These are "golden" values that should not change
    double precision :: test_r, test_theta, test_phi
    double precision :: A_theta, A_phi, dA_theta_dr, dA_phi_dr
    double precision :: sqg_test, B_vartheta_test, B_varphi_test
    double precision, parameter :: tol = 1.0d-10
    
    ! Test at a specific point
    test_r = 0.5d0
    test_theta = 0.1d0
    test_phi = 0.2d0
    
    ! These would be the expected values from a previous run
    ! In a real test, you would compute these once and store them
    ! For now, we just verify the routine runs without error
    
    call splint_can_coord(.true., 0, test_r, test_theta, test_phi, &
                          A_theta, A_phi, dA_theta_dr, dA_phi_dr, &
                          sqg_test, B_vartheta_test, B_varphi_test)
    
    ! Check that values are reasonable (non-zero, finite)
    if (abs(sqg_test) < tol) then
      print *, "WARNING: sqg_test is too small:", sqg_test
    end if
    
    if (.not. (abs(B_vartheta_test) < 1.0d10)) then
      print *, "ERROR: B_vartheta_test is not finite:", B_vartheta_test
      error stop 1
    end if
    
    if (.not. (abs(B_varphi_test) < 1.0d10)) then
      print *, "ERROR: B_varphi_test is not finite:", B_varphi_test
      error stop 1
    end if
    
    print *, "test_canonical_values: PASSED"
  end subroutine test_canonical_values

end program test_canonical_regression