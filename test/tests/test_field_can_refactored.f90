program test_field_can_refactored
! Unit tests for refactored field_can module
! Tests follow behavior-driven design with Given-When-Then structure

use field_can_mod, only: dp, optional_or_default, FieldCan_init, &
                         compute_d2vpar, compute_d2H, compute_d2pth, &
                         name_from_id, id_from_name, TEST, CANFLUX, BOOZER, MEISS, ALBERT
use field_can_base, only: FieldCan
implicit none

integer :: n_tests_passed, n_tests_failed
real(dp), parameter :: tol = 1.0e-12_dp

n_tests_passed = 0
n_tests_failed = 0

call test_optional_or_default()
call test_fieldcan_init()
call test_name_id_conversion()
call test_d2vpar_computation()
call test_d2H_computation()
call test_d2pth_computation()

print *, '============================================'
print *, 'Test Summary:'
print *, 'Passed:', n_tests_passed
print *, 'Failed:', n_tests_failed
print *, '============================================'

if (n_tests_failed > 0) then
  error stop 'Some tests failed!'
end if

contains

  !===========================================================================
  ! Test: Optional parameter helper function
  !===========================================================================
  subroutine test_optional_or_default()
    ! Given: Optional and default values
    ! When: We call optional_or_default with present and absent values
    ! Then: It should return the correct value
    
    real(dp) :: result
    real(dp) :: test_val
    logical :: test_passed
    
    print *, 'Testing optional_or_default helper...'
    
    test_passed = .true.
    
    ! Test with value present
    test_val = 5.0_dp
    result = optional_or_default(test_val, 2.0_dp)
    if (abs(result - 5.0_dp) > tol) then
      test_passed = .false.
      print *, '  Failed: Expected 5.0, got', result
    end if
    
    ! Test with value absent (simulate with null call)
    result = optional_or_default_test_absent(2.0_dp)
    if (abs(result - 2.0_dp) > tol) then
      test_passed = .false.
      print *, '  Failed: Expected 2.0 (default), got', result
    end if
    
    if (test_passed) then
      print *, '  PASSED: optional_or_default works correctly'
      n_tests_passed = n_tests_passed + 1
    else
      print *, '  FAILED: optional_or_default validation'
      n_tests_failed = n_tests_failed + 1
    end if
    
  end subroutine test_optional_or_default
  
  ! Helper function to test absent optional
  function optional_or_default_test_absent(default_val) result(val)
    real(dp), intent(in) :: default_val
    real(dp) :: val
    real(dp) :: dummy_optional
    ! Call without providing the optional argument
    val = call_optional_helper(default_val)
  end function optional_or_default_test_absent
  
  ! Additional helper to test absent optional
  function call_optional_helper(default_val, opt_val) result(val)
    real(dp), intent(in) :: default_val
    real(dp), intent(in), optional :: opt_val
    real(dp) :: val
    val = optional_or_default(opt_val, default_val)
  end function call_optional_helper

  !===========================================================================
  ! Test: FieldCan initialization
  !===========================================================================
  subroutine test_fieldcan_init()
    ! Given: A FieldCan structure
    ! When: We initialize it with different optional parameters
    ! Then: The values should be set correctly
    
    type(FieldCan) :: f
    logical :: test_passed
    
    print *, 'Testing FieldCan_init...'
    
    test_passed = .true.
    
    ! Test with all defaults
    call FieldCan_init(f)
    if (abs(f%mu) > tol .or. abs(f%ro0) > tol .or. abs(f%vpar) > tol) then
      test_passed = .false.
      print *, '  Failed: Default initialization not all zeros'
    end if
    
    ! Test with specific values
    call FieldCan_init(f, mu=1.5_dp, ro0=2.5_dp, vpar=3.5_dp)
    if (abs(f%mu - 1.5_dp) > tol .or. abs(f%ro0 - 2.5_dp) > tol .or. abs(f%vpar - 3.5_dp) > tol) then
      test_passed = .false.
      print *, '  Failed: Specific value initialization incorrect'
    end if
    
    ! Test with partial values
    call FieldCan_init(f, mu=4.0_dp)
    if (abs(f%mu - 4.0_dp) > tol .or. abs(f%ro0) > tol .or. abs(f%vpar) > tol) then
      test_passed = .false.
      print *, '  Failed: Partial initialization incorrect'
    end if
    
    if (test_passed) then
      print *, '  PASSED: FieldCan_init works correctly'
      n_tests_passed = n_tests_passed + 1
    else
      print *, '  FAILED: FieldCan_init validation'
      n_tests_failed = n_tests_failed + 1
    end if
    
  end subroutine test_fieldcan_init

  !===========================================================================
  ! Test: Name to ID and ID to name conversion
  !===========================================================================
  subroutine test_name_id_conversion()
    ! Given: Field names and IDs
    ! When: We convert between them
    ! Then: The conversions should be consistent
    
    integer :: id
    character(128) :: name
    logical :: test_passed
    
    print *, 'Testing name/ID conversion functions...'
    
    test_passed = .true.
    
    ! Test ID to name to ID round trip
    id = BOOZER
    name = name_from_id(id)
    if (id_from_name(name) /= BOOZER) then
      test_passed = .false.
      print *, '  Failed: BOOZER round trip conversion'
    end if
    
    ! Test name to ID to name round trip
    name = "flux"
    id = id_from_name(name)
    if (trim(name_from_id(id)) /= "flux") then
      test_passed = .false.
      print *, '  Failed: flux round trip conversion'
    end if
    
    ! Test all field types
    if (id_from_name("test") /= TEST) test_passed = .false.
    if (id_from_name("flux") /= CANFLUX) test_passed = .false.
    if (id_from_name("boozer") /= BOOZER) test_passed = .false.
    if (id_from_name("meiss") /= MEISS) test_passed = .false.
    if (id_from_name("albert") /= ALBERT) test_passed = .false.
    
    if (test_passed) then
      print *, '  PASSED: Name/ID conversions are consistent'
      n_tests_passed = n_tests_passed + 1
    else
      print *, '  FAILED: Name/ID conversion validation'
      n_tests_failed = n_tests_failed + 1
    end if
    
  end subroutine test_name_id_conversion

  !===========================================================================
  ! Test: d2vpar computation
  !===========================================================================
  subroutine test_d2vpar_computation()
    ! Given: Input arrays for d2vpar calculation
    ! When: We compute d2vpar
    ! Then: The result should follow the expected pattern
    
    real(dp) :: d2vpar(6)
    real(dp) :: d2Aph(6), d2hph(6), dhph(3), dvpar(3)
    real(dp) :: vpar, ro0, hph
    logical :: test_passed
    
    print *, 'Testing compute_d2vpar...'
    
    test_passed = .true.
    
    ! Set up test data
    d2Aph = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp]
    d2hph = [0.1_dp, 0.2_dp, 0.3_dp, 0.4_dp, 0.5_dp, 0.6_dp]
    dhph = [0.5_dp, 0.6_dp, 0.7_dp]
    dvpar = [0.2_dp, 0.3_dp, 0.4_dp]
    vpar = 2.0_dp
    ro0 = 1.0_dp
    hph = 2.0_dp
    
    call compute_d2vpar(d2vpar, d2Aph, d2hph, dhph, dvpar, vpar, ro0, hph)
    
    ! Check that computation completes and returns finite values
    if (any(abs(d2vpar) > 1000.0_dp)) then
      test_passed = .false.
      print *, '  Failed: d2vpar values out of reasonable range'
    end if
    
    ! Test with zero inputs
    d2Aph = 0.0_dp
    d2hph = 0.0_dp
    dhph = 0.0_dp
    dvpar = 0.0_dp
    vpar = 0.0_dp
    ro0 = 1.0_dp
    hph = 1.0_dp
    
    call compute_d2vpar(d2vpar, d2Aph, d2hph, dhph, dvpar, vpar, ro0, hph)
    
    if (any(abs(d2vpar) > tol)) then
      test_passed = .false.
      print *, '  Failed: Zero inputs should give zero output'
    end if
    
    if (test_passed) then
      print *, '  PASSED: compute_d2vpar works correctly'
      n_tests_passed = n_tests_passed + 1
    else
      print *, '  FAILED: compute_d2vpar validation'
      n_tests_failed = n_tests_failed + 1
    end if
    
  end subroutine test_d2vpar_computation

  !===========================================================================
  ! Test: d2H computation
  !===========================================================================
  subroutine test_d2H_computation()
    ! Given: Input parameters for d2H calculation
    ! When: We compute d2H
    ! Then: The result should be physically reasonable
    
    real(dp) :: d2H(6)
    real(dp) :: d2vpar(6), d2Bmod(6), dvpar(3)
    real(dp) :: vpar, mu
    logical :: test_passed
    integer :: i
    
    print *, 'Testing compute_d2H...'
    
    test_passed = .true.
    
    ! Set up test data
    d2vpar = [0.1_dp, 0.2_dp, 0.3_dp, 0.4_dp, 0.5_dp, 0.6_dp]
    d2Bmod = [0.01_dp, 0.02_dp, 0.03_dp, 0.04_dp, 0.05_dp, 0.06_dp]
    dvpar = [1.0_dp, 2.0_dp, 3.0_dp]
    vpar = 2.0_dp
    mu = 0.5_dp
    
    call compute_d2H(d2H, vpar, d2vpar, mu, d2Bmod, dvpar)
    
    ! Check basic computation
    ! d2H(1) should include dvpar(1)**2 term
    if (d2H(1) < dvpar(1)**2) then
      test_passed = .false.
      print *, '  Failed: d2H(1) should include dvpar(1)^2 contribution'
    end if
    
    ! Check symmetry-like properties
    ! d2H values should be finite and reasonable
    do i = 1, 6
      if (.not. (d2H(i) > -1000.0_dp .and. d2H(i) < 1000.0_dp)) then
        test_passed = .false.
        print *, '  Failed: d2H(', i, ') out of reasonable range:', d2H(i)
      end if
    end do
    
    if (test_passed) then
      print *, '  PASSED: compute_d2H works correctly'
      n_tests_passed = n_tests_passed + 1
    else
      print *, '  FAILED: compute_d2H validation'
      n_tests_failed = n_tests_failed + 1
    end if
    
  end subroutine test_d2H_computation

  !===========================================================================
  ! Test: d2pth computation
  !===========================================================================
  subroutine test_d2pth_computation()
    ! Given: Input parameters for d2pth calculation
    ! When: We compute d2pth
    ! Then: The result should follow expected patterns
    
    real(dp) :: d2pth(6)
    real(dp) :: d2vpar(6), d2hth(6), d2Ath(6)
    real(dp) :: dvpar(3), dhth(3)
    real(dp) :: hth, vpar, ro0
    logical :: test_passed
    
    print *, 'Testing compute_d2pth...'
    
    test_passed = .true.
    
    ! Set up test data
    d2vpar = [0.1_dp, 0.2_dp, 0.3_dp, 0.4_dp, 0.5_dp, 0.6_dp]
    d2hth = [0.01_dp, 0.02_dp, 0.03_dp, 0.04_dp, 0.05_dp, 0.06_dp]
    d2Ath = [0.001_dp, 0.002_dp, 0.003_dp, 0.004_dp, 0.005_dp, 0.006_dp]
    dvpar = [0.5_dp, 0.6_dp, 0.7_dp]
    dhth = [0.05_dp, 0.06_dp, 0.07_dp]
    hth = 1.5_dp
    vpar = 2.0_dp
    ro0 = 1.0_dp
    
    call compute_d2pth(d2pth, d2vpar, hth, vpar, d2hth, d2Ath, ro0, dvpar, dhth)
    
    ! Check that d2pth(1) includes the 2*dvpar(1)*dhth(1) term
    ! The base value should be at least d2vpar(1)*hth
    if (d2pth(1) < d2vpar(1)*hth) then
      test_passed = .false.
      print *, '  Failed: d2pth(1) missing base contribution'
    end if
    
    ! Test with all zeros except one term to verify contribution
    d2vpar = 0.0_dp
    d2hth = 0.0_dp
    d2Ath = 0.0_dp
    dvpar = [1.0_dp, 0.0_dp, 0.0_dp]
    dhth = [1.0_dp, 0.0_dp, 0.0_dp]
    
    call compute_d2pth(d2pth, d2vpar, hth, vpar, d2hth, d2Ath, ro0, dvpar, dhth)
    
    ! d2pth(1) should be 2*dvpar(1)*dhth(1) = 2.0
    if (abs(d2pth(1) - 2.0_dp) > tol) then
      test_passed = .false.
      print *, '  Failed: d2pth(1) cross term incorrect, expected 2.0, got', d2pth(1)
    end if
    
    if (test_passed) then
      print *, '  PASSED: compute_d2pth works correctly'
      n_tests_passed = n_tests_passed + 1
    else
      print *, '  FAILED: compute_d2pth validation'
      n_tests_failed = n_tests_failed + 1
    end if
    
  end subroutine test_d2pth_computation

end program test_field_can_refactored