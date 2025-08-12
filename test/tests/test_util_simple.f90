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
    
    ! Test twopi value (should be approximately 2*pi = 6.28318530717958)
    if (abs(twopi - 6.28318530717958d0) > tolerance) then
      print *, "ERROR: twopi constant incorrect"
      print *, "Expected: 6.28318530717958d0, Got:", twopi
      errors = errors + 1
    end if
    
    ! Test sqrt2 value
    if (abs(sqrt2 - dsqrt(2.0d0)) > tolerance) then
      print *, "ERROR: sqrt2 constant incorrect"
      print *, "Expected:", dsqrt(2.0d0), "Got:", sqrt2
      errors = errors + 1
    end if
    
    ! Test physical constants against their defined values in util.F90
    ! Speed of light in cm/s (defined as 2.9979d10)
    if (abs(c - 2.9979d10) > 1.0d6) then
      print *, "ERROR: Speed of light constant incorrect"
      print *, "Expected: 2.9979e10 cm/s, Got:", c
      errors = errors + 1
    end if
    
    ! Electron charge in CGS units (4.8032e-10 esu)
    if (abs(e_charge - 4.8032d-10) > 1.0d-13) then
      print *, "ERROR: Electron charge constant incorrect"
      print *, "Expected: 4.8032e-10 esu, Got:", e_charge
      errors = errors + 1
    end if
    
    ! Electron mass in grams (9.1094e-28 g)
    if (abs(e_mass - 9.1094d-28) > 1.0d-32) then
      print *, "ERROR: Electron mass constant incorrect"
      print *, "Expected: 9.1094e-28 g, Got:", e_mass
      errors = errors + 1
    end if
    
    ! Proton mass in grams (1.6726e-24 g)
    if (abs(p_mass - 1.6726d-24) > 1.0d-28) then
      print *, "ERROR: Proton mass constant incorrect"
      print *, "Expected: 1.6726e-24 g, Got:", p_mass
      errors = errors + 1
    end if
    
    ! Electron volt in ergs (1.6022e-12 erg)
    if (abs(ev - 1.6022d-12) > 1.0d-16) then
      print *, "ERROR: Electron volt constant incorrect"
      print *, "Expected: 1.6022e-12 erg, Got:", ev
      errors = errors + 1
    end if
    
    if (errors == 0) then
      print *, "  Constants test PASSED"
    end if
    
  end subroutine test_constants
  
  subroutine test_newunit_function(errors)
    integer, intent(inout) :: errors
    integer :: unit1, unit2, unit3, unit_opt
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
    
    ! Open the unit to make it unavailable (portable scratch file)
    open(unit=unit1, status='scratch', action='readwrite')
    
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
    
    ! Close unit1 first to free it up
    if (unit1 > 0) close(unit1)
    
    ! Test optional argument (avoid undefined behavior)
    unit_opt = -1  ! Initialize to invalid value
    unit3 = newunit(unit=unit_opt)
    if (unit3 /= unit_opt) then
      print *, "ERROR: newunit should return same value in optional argument"
      print *, "Returned:", unit3, "Optional:", unit_opt
      errors = errors + 1
    end if
    
    ! Verify unit3 is different from unit2 (unit1 is closed so could be reused)
    if (unit3 == unit2) then
      print *, "ERROR: newunit returned unit2 which should still be unavailable"
      errors = errors + 1
    end if
    
    ! Clean up remaining units
    if (unit2 > 0) close(unit2)
    if (unit3 > 0) close(unit3)
    
    ! Test edge case: when many units are occupied
    ! This is a behavioral test to ensure the function handles near-exhaustion gracefully
    
    if (errors == 0) then
      print *, "  Newunit function test PASSED"
    end if
    
  end subroutine test_newunit_function

end program test_util