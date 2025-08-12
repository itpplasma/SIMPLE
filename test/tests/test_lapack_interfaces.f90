program test_lapack_interfaces
  use lapack_interfaces
  implicit none
  
  integer :: errors
  
  errors = 0
  
  ! Test DGESV interface with simple linear system
  call test_dgesv_interface(errors)
  
  if (errors == 0) then
    print *, "All LAPACK interfaces tests passed!"
  else
    print *, "ERROR: ", errors, " test(s) failed!"
    stop 1
  end if
  
contains

  subroutine test_dgesv_interface(errors)
    integer, intent(inout) :: errors
    integer, parameter :: n = 3, nrhs = 1
    real(8) :: a(n,n), b(n,nrhs)
    integer :: ipiv(n), info
    real(8), parameter :: tolerance = 1.0d-12
    
    print *, "Testing DGESV interface..."
    
    ! Given: A simple 3x3 linear system Ax = b
    ! When: We solve it using DGESV
    ! Then: The solution should be correct and info should indicate success
    
    ! Set up the system: A*x = b where x = [1, 2, 3]
    ! A = [[2, 1, 0], [1, 2, 1], [0, 1, 2]]
    ! b = [4, 8, 7] (which gives x = [1, 2, 3])
    
    a(1,1) = 2.0d0; a(1,2) = 1.0d0; a(1,3) = 0.0d0
    a(2,1) = 1.0d0; a(2,2) = 2.0d0; a(2,3) = 1.0d0
    a(3,1) = 0.0d0; a(3,2) = 1.0d0; a(3,3) = 2.0d0
    
    b(1,1) = 4.0d0
    b(2,1) = 8.0d0
    b(3,1) = 7.0d0
    
    ! Call DGESV to solve the system
    call dgesv(n, nrhs, a, n, ipiv, b, n, info)
    
    ! Check that the solution completed successfully
    if (info /= 0) then
      print *, "ERROR: DGESV failed with info =", info
      errors = errors + 1
      return
    end if
    
    ! Check that the solution is correct (x = [1, 2, 3])
    if (abs(b(1,1) - 1.0d0) > tolerance) then
      print *, "ERROR: Incorrect solution for x(1)"
      print *, "Expected: 1.0, Got:", b(1,1)
      errors = errors + 1
    end if
    
    if (abs(b(2,1) - 2.0d0) > tolerance) then
      print *, "ERROR: Incorrect solution for x(2)"
      print *, "Expected: 2.0, Got:", b(2,1)
      errors = errors + 1
    end if
    
    if (abs(b(3,1) - 3.0d0) > tolerance) then
      print *, "ERROR: Incorrect solution for x(3)"
      print *, "Expected: 3.0, Got:", b(3,1)
      errors = errors + 1
    end if
    
    ! Check that pivot array contains valid indices
    if (any(ipiv < 1) .or. any(ipiv > n)) then
      print *, "ERROR: Invalid pivot indices"
      print *, "Pivot array:", ipiv
      errors = errors + 1
    end if
    
    ! Test edge case: singular matrix (should fail gracefully)
    call test_singular_matrix(errors)
    
    if (errors == 0) then
      print *, "  DGESV interface test PASSED"
    end if
    
  end subroutine test_dgesv_interface
  
  subroutine test_singular_matrix(errors)
    integer, intent(inout) :: errors
    integer, parameter :: n = 2, nrhs = 1
    real(8) :: a(n,n), b(n,nrhs)
    integer :: ipiv(n), info
    
    print *, "Testing DGESV with singular matrix..."
    
    ! Given: A singular matrix (non-invertible)
    ! When: We try to solve the system
    ! Then: DGESV should return a non-zero info value
    
    ! Create a singular matrix (row 2 = 2 * row 1)
    a(1,1) = 1.0d0; a(1,2) = 2.0d0
    a(2,1) = 2.0d0; a(2,2) = 4.0d0
    
    b(1,1) = 1.0d0
    b(2,1) = 2.0d0
    
    call dgesv(n, nrhs, a, n, ipiv, b, n, info)
    
    ! For a singular matrix, info should be > 0
    if (info == 0) then
      print *, "ERROR: DGESV should detect singular matrix"
      errors = errors + 1
    else
      print *, "  Singular matrix correctly detected (info =", info, ")"
    end if
    
  end subroutine test_singular_matrix

end program test_lapack_interfaces