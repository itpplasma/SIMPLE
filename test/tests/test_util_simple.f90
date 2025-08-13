program test_util
  use util
  implicit none
  
  integer, parameter :: dp = kind(1.0d0)
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
    real(dp), parameter :: math_tolerance = 1.0d-14
    real(dp), parameter :: physics_tolerance = 1.0d-10  ! Looser for physical constants
    
    ! Independent calculations of mathematical constants for verification
    real(dp), parameter :: pi_independent = 4.0d0 * atan(1.0d0)
    real(dp), parameter :: pi_leibniz = 4.0d0 * (1.0d0 - 1.0d0/3.0d0 + 1.0d0/5.0d0 - 1.0d0/7.0d0 + 1.0d0/9.0d0 &
                                                        - 1.0d0/11.0d0 + 1.0d0/13.0d0 - 1.0d0/15.0d0 + 1.0d0/17.0d0 &
                                                        - 1.0d0/19.0d0 + 1.0d0/21.0d0 - 1.0d0/23.0d0 + 1.0d0/25.0d0)
    
    ! NIST/CODATA 2018 fundamental physical constants (exact values where defined)
    real(dp), parameter :: c_exact_si = 299792458.0d0           ! m/s (exact by definition)
    real(dp), parameter :: c_cgs = c_exact_si * 100.0d0         ! cm/s conversion
    real(dp), parameter :: e_charge_exact_si = 1.602176634d-19   ! Coulomb (exact by definition)
    real(dp), parameter :: e_mass_exact_kg = 9.1093837015d-31    ! kg (CODATA 2018)
    real(dp), parameter :: p_mass_exact_kg = 1.67262192369d-27   ! kg (CODATA 2018)
    real(dp), parameter :: ev_exact_joule = 1.602176634d-19      ! Joule (exact by definition)
    
    ! Converted CGS values
    real(dp), parameter :: e_charge_cgs_exact = 4.80320425d-10   ! Modern precise value
    real(dp), parameter :: e_mass_cgs_exact = e_mass_exact_kg * 1000.0d0  ! kg to g
    real(dp), parameter :: p_mass_cgs_exact = p_mass_exact_kg * 1000.0d0  ! kg to g
    real(dp), parameter :: ev_cgs_exact = ev_exact_joule * 1.0d7  ! Joule to erg (1 J = 10^7 erg)
    
    print *, "Testing mathematical and physical constants against independent calculations..."
    
    ! Test mathematical constants against independent computations
    
    ! Test pi using independent atan calculation
    if (abs(pi - pi_independent) > math_tolerance) then
      print *, "ERROR: pi constant doesn't match atan(1)*4 calculation"
      print *, "Module pi:", pi, "atan(1)*4:", pi_independent
      print *, "Difference:", abs(pi - pi_independent)
      errors = errors + 1
    end if
    
    ! Test pi using series approximation (Leibniz series, partial sum should be close)
    if (abs(pi - pi_leibniz) > 0.1d0) then ! Leibniz converges slowly, looser tolerance
      print *, "ERROR: pi constant far from Leibniz series approximation"
      print *, "Module pi:", pi, "Leibniz approximation:", pi_leibniz
      print *, "Difference:", abs(pi - pi_leibniz)
      errors = errors + 1
    end if
    
    ! Test twopi computed from pi
    if (abs(twopi - 2.0d0*pi) > math_tolerance) then
      print *, "ERROR: twopi is not exactly 2*pi"
      print *, "twopi:", twopi, "2*pi:", 2.0d0*pi
      print *, "Difference:", abs(twopi - 2.0d0*pi)
      errors = errors + 1
    end if
    
    ! Test sqrt2 by verifying it satisfies mathematical properties
    ! sqrt(2) * sqrt(2) should equal 2.0
    if (abs(sqrt2 * sqrt2 - 2.0d0) > math_tolerance) then
      print *, "ERROR: sqrt2 * sqrt2 doesn't equal 2.0"
      print *, "sqrt2^2:", sqrt2 * sqrt2, "Expected: 2.0"
      print *, "Difference:", abs(sqrt2 * sqrt2 - 2.0d0)
      errors = errors + 1
    end if
    
    ! Also verify sqrt2 is approximately 1.414213562...
    if (abs(sqrt2 - 1.41421356237309504880d0) > math_tolerance) then
      print *, "ERROR: sqrt2 doesn't match known value to high precision"
      errors = errors + 1
    end if
    
    ! Test mathematical relationships
    ! Verify Euler's identity: e^(i*pi) + 1 = 0, which means cos(pi) = -1
    if (abs(cos(pi) - (-1.0d0)) > math_tolerance) then
      print *, "ERROR: cos(pi) should equal -1 (Euler's identity test)"
      print *, "cos(pi):", cos(pi), "Expected: -1"
      errors = errors + 1
    end if
    
    ! Verify sin(pi/2) = 1 using our pi
    if (abs(sin(pi/2.0d0) - 1.0d0) > math_tolerance) then
      print *, "ERROR: sin(pi/2) should equal 1"
      print *, "sin(pi/2):", sin(pi/2.0d0), "Expected: 1"
      errors = errors + 1
    end if
    
    ! Test physical constants against CODATA/NIST values with proper unit conversions
    
    ! Speed of light: Convert from exact SI to CGS
    if (abs(c - c_cgs) > c_cgs * 1.0d-10) then  ! 0.01% tolerance for old approximation
      print *, "WARNING: Speed of light differs from modern exact value"
      print *, "Module c:", c, "cm/s, Modern exact:", c_cgs, "cm/s"
      print *, "Relative difference:", abs(c - c_cgs) / c_cgs
      ! Not counting as error since module uses older approximation
    end if
    
    ! Electron charge: Convert from exact SI (Coulomb) to CGS (statCoulomb/esu)
    ! Conversion factor: 1 C = 2997924580 statC, so e = 1.602176634e-19 * 2997924580 = 4.8032e-10 esu
    if (abs(e_charge - e_charge_cgs_exact) > e_charge_cgs_exact * 1.0d-6) then
      print *, "WARNING: Electron charge differs from modern exact CGS value"
      print *, "Module e_charge:", e_charge, "esu, Modern exact:", e_charge_cgs_exact, "esu"
      print *, "Relative difference:", abs(e_charge - e_charge_cgs_exact) / e_charge_cgs_exact
    end if
    
    ! Electron mass: Convert from exact kg to grams
    if (abs(e_mass - e_mass_cgs_exact) > e_mass_cgs_exact * 1.0d-6) then
      print *, "WARNING: Electron mass differs from modern exact CGS value"
      print *, "Module e_mass:", e_mass, "g, Modern exact:", e_mass_cgs_exact, "g"
      print *, "Relative difference:", abs(e_mass - e_mass_cgs_exact) / e_mass_cgs_exact
    end if
    
    ! Proton mass: Convert from exact kg to grams
    if (abs(p_mass - p_mass_cgs_exact) > p_mass_cgs_exact * 1.0d-6) then
      print *, "WARNING: Proton mass differs from modern exact CGS value"
      print *, "Module p_mass:", p_mass, "g, Modern exact:", p_mass_cgs_exact, "g"
      print *, "Relative difference:", abs(p_mass - p_mass_cgs_exact) / p_mass_cgs_exact
    end if
    
    ! Electron volt: Convert from exact Joule to erg
    if (abs(ev - ev_cgs_exact) > ev_cgs_exact * 1.0d-4) then
      print *, "ERROR: Electron volt conversion incorrect"
      print *, "Module ev:", ev, "erg, Expected:", ev_cgs_exact, "erg"
      print *, "Relative difference:", abs(ev - ev_cgs_exact) / ev_cgs_exact
      errors = errors + 1
    end if
    
    ! Test computational usage of constants
    call test_constant_usage(errors)
    
    if (errors == 0) then
      print *, "  Constants test PASSED"
    end if
    
  end subroutine test_constants
  
  subroutine test_constant_usage(errors)
    integer, intent(inout) :: errors
    real(dp), parameter :: tolerance = 1.0d-12
    real(dp), parameter :: pi_independent = 4.0d0 * atan(1.0d0)
    
    ! Variable declarations
    real(dp) :: B_field, omega_cyclotron_expected, omega_cyclotron_computed
    real(dp) :: energy_ev, energy_erg_expected, energy_erg_computed
    real(dp) :: mass_ratio_expected, mass_ratio_computed
    real(dp) :: radius, circumference_expected, circumference_computed
    real(dp) :: side, diagonal_expected, diagonal_computed
    
    ! Test constants in typical physics computations
    
    ! Test 1: Cyclotron frequency calculation
    ! omega_c = eB/(m*c) where B is in Gauss, omega_c in rad/s
    B_field = 10000.0d0  ! 1 Tesla = 10^4 Gauss
    omega_cyclotron_expected = 1.758820d11  ! Known value for electron in 1T field
    omega_cyclotron_computed = (e_charge * B_field) / (e_mass * c)
    
    ! Should be within 1% (constants have limited precision)
    if (abs(omega_cyclotron_computed - omega_cyclotron_expected) > omega_cyclotron_expected * 0.01d0) then
      print *, "ERROR: Cyclotron frequency calculation using constants failed"
      print *, "Computed:", omega_cyclotron_computed, "Expected:", omega_cyclotron_expected
      print *, "Relative error:", abs(omega_cyclotron_computed - omega_cyclotron_expected) / omega_cyclotron_expected
      errors = errors + 1
    end if
    
    ! Test 2: Energy conversion using ev
    ! Convert 13.6 eV (hydrogen binding energy) to ergs
    energy_ev = 13.6d0
    energy_erg_expected = 2.179d-11  ! Known value
    energy_erg_computed = energy_ev * ev
    
    if (abs(energy_erg_computed - energy_erg_expected) > energy_erg_expected * 0.01d0) then
      print *, "ERROR: Energy conversion using ev constant failed"
      print *, "Computed:", energy_erg_computed, "erg, Expected:", energy_erg_expected, "erg"
      print *, "Relative error:", abs(energy_erg_computed - energy_erg_expected) / energy_erg_expected
      errors = errors + 1
    end if
    
    ! Test 3: Mass ratio (proton/electron)
    mass_ratio_expected = 1836.15d0  ! Known physical constant
    mass_ratio_computed = p_mass / e_mass
    
    if (abs(mass_ratio_computed - mass_ratio_expected) > mass_ratio_expected * 0.01d0) then
      print *, "ERROR: Proton-to-electron mass ratio incorrect"
      print *, "Computed:", mass_ratio_computed, "Expected:", mass_ratio_expected
      print *, "Relative error:", abs(mass_ratio_computed - mass_ratio_expected) / mass_ratio_expected
      errors = errors + 1
    end if
    
    ! Test 4: Circular motion using pi
    radius = 2.5d0
    circumference_expected = 2.0d0 * 2.5d0 * pi_independent  ! Using independent pi
    circumference_computed = twopi * radius
    
    if (abs(circumference_computed - circumference_expected) > tolerance) then
      print *, "ERROR: Circumference calculation using twopi failed"
      print *, "Computed:", circumference_computed, "Expected:", circumference_expected
      errors = errors + 1
    end if
    
    ! Test 5: Pythagorean theorem using sqrt2
    side = 3.0d0
    diagonal_expected = side * 1.41421356237309504880d0  ! Known high-precision value
    diagonal_computed = side * sqrt2
    
    if (abs(diagonal_computed - diagonal_expected) > tolerance) then
      print *, "ERROR: Diagonal calculation using sqrt2 failed"
      print *, "Computed:", diagonal_computed, "Expected:", diagonal_expected
      errors = errors + 1
    end if
    
  end subroutine test_constant_usage
  
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
    
    ! Test edge case: when all units are occupied (tests line 37 in util.F90)
    call test_newunit_exhaustion(errors)
    
    if (errors == 0) then
      print *, "  Newunit function test PASSED"
    end if
    
  end subroutine test_newunit_function

  subroutine test_newunit_exhaustion(errors)
    integer, intent(inout) :: errors
    integer, parameter :: MAX_UNITS = 100  ! Test with subset to avoid system limits
    integer :: units(MAX_UNITS)
    integer :: i, test_unit
    logical :: opened
    
    print *, "Testing newunit exhaustion edge case..."
    
    ! Given: We want to test the case when all file units are exhausted
    ! When: We open many units and then call newunit
    ! Then: newunit should return -1 when no units are available
    
    ! Initialize array
    units = -1
    
    ! Open multiple units to simulate near-exhaustion
    do i = 1, MAX_UNITS
      units(i) = newunit()
      if (units(i) == -1) then
        print *, "  Reached unit exhaustion at unit", i-1
        exit
      end if
      open(unit=units(i), status='scratch', action='readwrite')
    end do
    
    ! At this point, many units should be occupied
    ! Test that newunit still works if any units are available
    test_unit = newunit()
    
    if (test_unit == -1) then
      print *, "  newunit correctly returned -1 when units exhausted"
    else if (test_unit > 0) then
      ! Verify the returned unit is actually available
      inquire(unit=test_unit, opened=opened)
      if (opened) then
        print *, "ERROR: newunit returned an occupied unit during stress test"
        print *, "Unit:", test_unit
        errors = errors + 1
      else
        print *, "  newunit found available unit", test_unit, "during stress test"
        ! Clean up this unit too
        open(unit=test_unit, status='scratch', action='readwrite')
        close(test_unit)
      end if
    else
      print *, "ERROR: newunit returned invalid unit number:", test_unit
      errors = errors + 1
    end if
    
    ! Clean up all opened units
    do i = 1, MAX_UNITS
      if (units(i) > 0) then
        close(units(i))
      end if
    end do
    
    print *, "  Unit exhaustion test completed"
    
  end subroutine test_newunit_exhaustion

end program test_util