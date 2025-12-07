module mock_field_module
  use field_base
  implicit none
  
  ! Mock implementation of magnetic_field_t for testing
  type, extends(magnetic_field_t) :: MockField
  contains
    procedure :: evaluate => mock_evaluate
  end type MockField
  
contains

  subroutine mock_evaluate(self, x, Acov, hcov, Bmod, sqgBctr)
    class(MockField), intent(in) :: self
    real(dp), intent(in) :: x(3)
    real(dp), intent(out) :: Acov(3)
    real(dp), intent(out) :: hcov(3)
    real(dp), intent(out) :: Bmod
    real(dp), intent(out), optional :: sqgBctr(3)
    
    ! Simple mock implementation that returns predictable values
    Acov = 0.0_dp
    hcov = 1.0_dp
    Bmod = 1.0_dp
    
    if (present(sqgBctr)) then
        sqgBctr = x  ! Just return the input coordinates
    end if
    
  end subroutine mock_evaluate

end module mock_field_module

program test_field_base
  use field_base
  use mock_field_module
  implicit none
  
  integer :: errors
  
  errors = 0
  
  ! Test abstract interface definition
  call test_abstract_interface(errors)
  
  if (errors == 0) then
    print *, "All field_base module tests passed!"
  else
    print *, "ERROR: ", errors, " test(s) failed!"
    stop 1
  end if
  
contains

  subroutine test_abstract_interface(errors)
    integer, intent(inout) :: errors
    
    print *, "Testing abstract magnetic_field_t interface..."
    
    ! Given: The field_base module defines an abstract magnetic_field_t type
    ! When: We create a concrete implementation
    ! Then: It should satisfy the interface requirements
    
    ! Test that we can create a mock field and verify interface
    call test_mock_field_implementation(errors)
    
    if (errors == 0) then
      print *, "  Abstract interface test PASSED"
    end if
    
  end subroutine test_abstract_interface
  
  subroutine test_mock_field_implementation(errors)
    integer, intent(inout) :: errors
    type(MockField) :: mock_field
    real(dp) :: x(3), Acov(3), hcov(3), Bmod, sqgBctr(3)
    real(dp), parameter :: tolerance = 1.0d-14
    
    ! Given: A mock implementation of magnetic_field_t
    ! When: We call the evaluate method
    ! Then: It should return predictable values
    
    x = [0.5_dp, 0.0_dp, 0.0_dp]
    
    call mock_field%evaluate(x, Acov, hcov, Bmod, sqgBctr)
    
    ! Check that the mock field produces expected values
    if (abs(Bmod - 1.0_dp) > tolerance) then
      print *, "ERROR: Mock field Bmod incorrect"
      print *, "Expected: 1.0, Got:", Bmod
      errors = errors + 1
    end if
    
    if (any(abs(Acov - 0.0_dp) > tolerance)) then
      print *, "ERROR: Mock field Acov should be zero"
      errors = errors + 1
    end if
    
    if (any(abs(hcov - 1.0_dp) > tolerance)) then
      print *, "ERROR: Mock field hcov should be unity"
      errors = errors + 1
    end if
    
    ! Test call without optional argument
    call mock_field%evaluate(x, Acov, hcov, Bmod)
    
    if (abs(Bmod - 1.0_dp) > tolerance) then
      print *, "ERROR: Mock field call without optional argument failed"
      errors = errors + 1
    end if
    
  end subroutine test_mock_field_implementation

end program test_field_base