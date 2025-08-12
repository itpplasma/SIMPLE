program test_find_bminmax_refactored
  use find_bminmax_sub
  use bminmax_mod
  implicit none

  integer :: num_passed, num_failed
  
  num_passed = 0
  num_failed = 0
  
  print *, "========================================="
  print *, "Testing Refactored Find B Min/Max Module"
  print *, "========================================="
  
  call test_scan_grid_for_extrema(num_passed, num_failed)
  call test_compute_b_second_derivatives(num_passed, num_failed)
  call test_initialize_bminmax_arrays(num_passed, num_failed)
  call test_interpolate_bminmax(num_passed, num_failed)
  call test_get_bminmax_caching(num_passed, num_failed)
  
  print *, "========================================="
  print *, "Test Summary:"
  print *, "  Passed: ", num_passed
  print *, "  Failed: ", num_failed
  print *, "========================================="
  
  if (num_failed > 0) then
    error stop "Some tests failed"
  end if

contains

  subroutine test_scan_grid_for_extrema(num_passed, num_failed)
    integer, intent(inout) :: num_passed, num_failed
    real(kind(1.0d0)) :: s, bmin, bmax, tmin, tmax, pmin, pmax
    logical :: test_passed
    
    print *, ""
    print *, "Test: scan_grid_for_extrema"
    print *, "  Given: Flux surface coordinate and grid parameters"
    print *, "  When: Scanning grid for B field extrema"
    print *, "  Then: Should find min and max B values with locations"
    
    test_passed = .true.
    
    ! Test at s = 0.5
    s = 0.5d0
    call scan_grid_for_extrema(s, 10, 10, bmin, bmax, tmin, tmax, pmin, pmax)
    
    ! Check that bmin <= bmax
    if (bmin > bmax) then
      test_passed = .false.
      print *, "    ERROR: bmin > bmax: ", bmin, " > ", bmax
    end if
    
    ! Check that angles are within [0, 2π]
    if (tmin < 0.0d0 .or. tmin > 8.0d0*atan(1.0d0) .or. &
        tmax < 0.0d0 .or. tmax > 8.0d0*atan(1.0d0)) then
      test_passed = .false.
      print *, "    ERROR: Theta angles out of range"
    end if
    
    if (pmin < 0.0d0 .or. pmin > 8.0d0*atan(1.0d0) .or. &
        pmax < 0.0d0 .or. pmax > 8.0d0*atan(1.0d0)) then
      test_passed = .false.
      print *, "    ERROR: Phi angles out of range"
    end if
    
    if (test_passed) then
      print *, "    ✓ Grid scan finds valid extrema"
      print *, "      bmin=", bmin, " at (θ=", tmin, ", φ=", pmin, ")"
      print *, "      bmax=", bmax, " at (θ=", tmax, ", φ=", pmax, ")"
      num_passed = num_passed + 1
    else
      print *, "    ✗ Grid scan failed"
      num_failed = num_failed + 1
    end if
  end subroutine test_scan_grid_for_extrema
  
  subroutine test_compute_b_second_derivatives(num_passed, num_failed)
    integer, intent(inout) :: num_passed, num_failed
    real(kind(1.0d0)) :: s, theta, phi, hdt, hdp
    real(kind(1.0d0)) :: btt, btp, bpt, bpp
    logical :: test_passed
    
    print *, ""
    print *, "Test: compute_b_second_derivatives"
    print *, "  Given: Position and finite difference steps"
    print *, "  When: Computing second derivatives of B field"
    print *, "  Then: Should return Hessian matrix elements"
    
    test_passed = .true.
    
    s = 0.5d0
    theta = 0.5d0
    phi = 0.5d0
    hdt = 1.0d-3
    hdp = 1.0d-3
    
    call compute_b_second_derivatives(s, theta, phi, hdt, hdp, &
                                      btt, btp, bpt, bpp)
    
    ! Check that mixed derivatives are symmetric (within numerical tolerance)
    if (abs(btp - bpt) > 1.0d-8 * max(abs(btp), abs(bpt), 1.0d0)) then
      test_passed = .false.
      print *, "    ERROR: Mixed derivatives not symmetric: btp=", btp, " bpt=", bpt
    else
      print *, "    ✓ Mixed derivatives are symmetric"
    end if
    
    ! Check that derivatives are finite
    if (.not. (abs(btt) < huge(btt) .and. abs(bpp) < huge(bpp))) then
      test_passed = .false.
      print *, "    ERROR: Non-finite derivatives"
    else
      print *, "    ✓ All derivatives are finite"
    end if
    
    if (test_passed) then
      print *, "      Hessian: btt=", btt, " btp=", btp, " bpp=", bpp
      num_passed = num_passed + 1
    else
      num_failed = num_failed + 1
    end if
  end subroutine test_compute_b_second_derivatives
  
  subroutine test_initialize_bminmax_arrays(num_passed, num_failed)
    integer, intent(inout) :: num_passed, num_failed
    logical :: test_passed
    integer :: k
    real(kind(1.0d0)) :: expected_hsbmnx
    
    print *, ""
    print *, "Test: initialize_bminmax_arrays"
    print *, "  Given: Uninitialized B min/max arrays"
    print *, "  When: Initializing arrays"
    print *, "  Then: Should compute and store B min/max for all flux surfaces"
    
    test_passed = .true.
    
    ! Reset prop flag to force initialization
    prop = .true.
    
    call initialize_bminmax_arrays
    
    ! Check that prop flag is now false
    if (prop) then
      test_passed = .false.
      print *, "    ERROR: prop flag not set to false"
    else
      print *, "    ✓ Initialization flag correctly set"
    end if
    
    ! Check that hsbmnx is correctly set
    expected_hsbmnx = 1.0d0 / dble(nsbmnx)
    if (abs(hsbmnx - expected_hsbmnx) > 1.0d-14) then
      test_passed = .false.
      print *, "    ERROR: Incorrect hsbmnx: ", hsbmnx, " expected ", expected_hsbmnx
    else
      print *, "    ✓ Grid spacing correctly calculated"
    end if
    
    ! Check that arrays are filled with valid values
    do k = 0, nsbmnx
      if (bmin_arr(k) > bmax_arr(k)) then
        test_passed = .false.
        print *, "    ERROR: bmin > bmax at k=", k
        exit
      end if
    end do
    
    if (test_passed) then
      print *, "    ✓ All array values are valid (bmin <= bmax)"
      print *, "      Array size: ", nsbmnx + 1, " points"
      num_passed = num_passed + 1
    else
      num_failed = num_failed + 1
    end if
  end subroutine test_initialize_bminmax_arrays
  
  subroutine test_interpolate_bminmax(num_passed, num_failed)
    integer, intent(inout) :: num_passed, num_failed
    real(kind(1.0d0)) :: s, bmin, bmax
    real(kind(1.0d0)) :: bmin_test1, bmax_test1, bmin_test2, bmax_test2
    logical :: test_passed
    
    print *, ""
    print *, "Test: interpolate_bminmax"
    print *, "  Given: Initialized arrays and flux surface coordinate"
    print *, "  When: Interpolating B min/max values"
    print *, "  Then: Should return smoothly interpolated values"
    
    test_passed = .true.
    
    ! Ensure arrays are initialized
    call initialize_bminmax_arrays
    
    ! Test at s = 0.25
    s = 0.25d0
    call interpolate_bminmax(s, bmin_test1, bmax_test1)
    
    ! Test at s = 0.26 (should be close to previous)
    s = 0.26d0
    call interpolate_bminmax(s, bmin_test2, bmax_test2)
    
    ! Check continuity - nearby points should have similar values
    if (abs(bmin_test2 - bmin_test1) > 0.1d0 * abs(bmin_test1)) then
      test_passed = .false.
      print *, "    ERROR: Large discontinuity in bmin interpolation"
    else
      print *, "    ✓ Bmin interpolation is continuous"
    end if
    
    if (abs(bmax_test2 - bmax_test1) > 0.1d0 * abs(bmax_test1)) then
      test_passed = .false.
      print *, "    ERROR: Large discontinuity in bmax interpolation"
    else
      print *, "    ✓ Bmax interpolation is continuous"
    end if
    
    ! Test boundary case s = 0
    s = 0.0d0
    call interpolate_bminmax(s, bmin, bmax)
    
    if (bmin > bmax) then
      test_passed = .false.
      print *, "    ERROR: Invalid interpolation at s=0"
    else
      print *, "    ✓ Valid interpolation at boundary s=0"
    end if
    
    ! Test boundary case s = 1
    s = 1.0d0
    call interpolate_bminmax(s, bmin, bmax)
    
    if (bmin > bmax) then
      test_passed = .false.
      print *, "    ERROR: Invalid interpolation at s=1"
    else
      print *, "    ✓ Valid interpolation at boundary s=1"
    end if
    
    if (test_passed) then
      num_passed = num_passed + 1
    else
      num_failed = num_failed + 1
    end if
  end subroutine test_interpolate_bminmax
  
  subroutine test_get_bminmax_caching(num_passed, num_failed)
    integer, intent(inout) :: num_passed, num_failed
    real(kind(1.0d0)) :: s, bmin1, bmax1, bmin2, bmax2
    logical :: test_passed
    
    print *, ""
    print *, "Test: get_bminmax with caching"
    print *, "  Given: Multiple calls to get_bminmax"
    print *, "  When: Getting B min/max for same flux surface"
    print *, "  Then: Should use cached values for efficiency"
    
    test_passed = .true.
    
    ! Reset prop to test initialization
    prop = .true.
    
    ! First call - should initialize arrays
    s = 0.5d0
    call get_bminmax(s, bmin1, bmax1)
    
    ! Check that prop is now false (arrays initialized)
    if (prop) then
      test_passed = .false.
      print *, "    ERROR: Arrays not initialized on first call"
    else
      print *, "    ✓ Arrays initialized on first call"
    end if
    
    ! Second call - should use cached values
    call get_bminmax(s, bmin2, bmax2)
    
    ! Results should be identical
    if (abs(bmin1 - bmin2) > 1.0d-14 .or. abs(bmax1 - bmax2) > 1.0d-14) then
      test_passed = .false.
      print *, "    ERROR: Inconsistent results between calls"
    else
      print *, "    ✓ Consistent results from cached arrays"
    end if
    
    ! Check that values are valid
    if (bmin1 > bmax1) then
      test_passed = .false.
      print *, "    ERROR: Invalid B field extrema: bmin > bmax"
    else
      print *, "    ✓ Valid B field extrema returned"
      print *, "      s=", s, " bmin=", bmin1, " bmax=", bmax1
    end if
    
    if (test_passed) then
      num_passed = num_passed + 1
    else
      num_failed = num_failed + 1
    end if
  end subroutine test_get_bminmax_caching
  
end program test_find_bminmax_refactored