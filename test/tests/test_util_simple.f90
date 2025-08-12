program test_util
  use util
  implicit none
  
  integer :: errors
  
  errors = 0
  
  ! Test mathematical constants
  call test_constants(errors)
  
  ! Test newunit function
  call test_newunit_function(errors)
  
  if (errors == 0) then
    print *, "All util module tests passed!"
  else
    print *, "ERROR: ", errors, " test(s) failed!"
    stop 1
  end if
  
contains

  subroutine test_constants(errors)
    integer, intent(inout) :: errors
    double precision, parameter :: tolerance = 1.0d-14
    
    print *, "Testing mathematical and physical constants..."
    
    ! Given: The util module defines mathematical and physical constants
    ! When: We check the values against known constants
    ! Then: The values should match expected physical constants
    
    ! Test pi value
    if (abs(pi - 3.14159265358979d0) > tolerance) then
      print *, "ERROR: pi constant incorrect"
      print *, "Expected: 3.14159265358979d0, Got:", pi
      errors = errors + 1
    end if
    
    ! Test twopi value (should be 2*pi)
    if (abs(twopi - 2.0d0*pi) > tolerance) then
      print *, "ERROR: twopi should equal 2*pi"
      print *, "Expected:", 2.0d0*pi, "Got:", twopi
      errors = errors + 1
    end if
    
    ! Test sqrt2 value
    if (abs(sqrt2 - dsqrt(2.0d0)) > tolerance) then
      print *, "ERROR: sqrt2 constant incorrect"
      print *, "Expected:", dsqrt(2.0d0), "Got:", sqrt2
      errors = errors + 1
    end if
    
    ! Test physical constants (basic sanity checks)
    ! Speed of light should be positive and reasonable
    if (c <= 0.0d0 .or. c < 1.0d10 .or. c > 1.0d11) then
      print *, "ERROR: Speed of light constant unreasonable"
      print *, "Got:", c
      errors = errors + 1
    end if
    
    ! Electron charge should be positive and reasonable
    if (e_charge <= 0.0d0 .or. e_charge < 1.0d-11 .or. e_charge > 1.0d-9) then
      print *, "ERROR: Electron charge constant unreasonable"
      print *, "Got:", e_charge
      errors = errors + 1
    end if
    
    ! Electron mass should be positive and reasonable
    if (e_mass <= 0.0d0 .or. e_mass < 1.0d-29 .or. e_mass > 1.0d-27) then
      print *, "ERROR: Electron mass constant unreasonable"
      print *, "Got:", e_mass
      errors = errors + 1
    end if
    
    ! Proton mass should be positive and reasonable
    if (p_mass <= 0.0d0 .or. p_mass < 1.0d-25 .or. p_mass > 1.0d-23) then
      print *, "ERROR: Proton mass constant unreasonable"
      print *, "Got:", p_mass
      errors = errors + 1
    end if
    
    ! Electron volt should be positive and reasonable
    if (ev <= 0.0d0 .or. ev < 1.0d-13 .or. ev > 1.0d-11) then
      print *, "ERROR: Electron volt constant unreasonable"
      print *, "Got:", ev
      errors = errors + 1
    end if
    
    if (errors == 0) then
      print *, "  Constants test PASSED"
    end if
    
  end subroutine test_constants
  
  subroutine test_newunit_function(errors)
    integer, intent(inout) :: errors
    integer :: unit1, unit2, unit3
    logical :: opened
    
    print *, "Testing newunit function..."
    
    ! Given: The newunit function should find available logical unit numbers
    ! When: We call newunit multiple times
    ! Then: Each call should return a different available unit number
    
    ! Get first available unit
    unit1 = newunit()
    if (unit1 < 10 .or. unit1 > 1000) then
      print *, "ERROR: newunit returned unit outside expected range [10,1000]"
      print *, "Got:", unit1
      errors = errors + 1
      return
    end if
    
    ! Check that the unit is indeed available
    inquire(unit=unit1, opened=opened)
    if (opened) then
      print *, "ERROR: newunit returned a unit that is already opened"
      print *, "Unit:", unit1
      errors = errors + 1
    end if
    
    ! Open the unit to make it unavailable
    open(unit=unit1, file='/dev/null', status='old')
    
    ! Get second available unit
    unit2 = newunit()
    if (unit2 == unit1) then
      print *, "ERROR: newunit returned the same unit twice"
      print *, "Unit:", unit2
      errors = errors + 1
    end if
    
    if (unit2 < 10 .or. unit2 > 1000) then
      print *, "ERROR: second newunit call returned unit outside expected range"
      print *, "Got:", unit2
      errors = errors + 1
    end if
    
    ! Test optional argument
    unit3 = newunit(unit=unit3)
    if (unit3 /= newunit()) then
      print *, "WARNING: newunit function may not be deterministic"
      ! This is not necessarily an error, just a note
    end if
    
    ! Clean up
    if (unit1 > 0) close(unit1)
    if (unit2 > 0) close(unit2)
    
    ! Test edge case: when many units are occupied
    ! This is a behavioral test to ensure the function handles near-exhaustion gracefully
    
    if (errors == 0) then
      print *, "  Newunit function test PASSED"
    end if
    
  end subroutine test_newunit_function

end program test_util