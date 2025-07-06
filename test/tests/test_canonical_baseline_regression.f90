program test_canonical_baseline_regression
  ! Regression test that compares current output against known baseline
  ! This test runs simple.x and verifies debug output matches expected values
  
  implicit none
  
  ! Expected baseline output key lines  
  character(len=80), parameter :: expected_lines(6) = [ &
    " DEBUG: get_canonical_coordinates starting                                    ", &
    " DEBUG: ns_c=          50 n_theta_c=          36 n_phi_c=          81        ", &
    " DEBUG: nh_stencil=           3                                               ", &
    " DEBUG: G_c(1,1,1)=   1.0000000000000000E-008                               ", &
    " DEBUG: spline_can_coord starting, fullset= T                               ", &
    " DEBUG: spline_can_coord completed                                           " &
  ]
  
  integer, parameter :: max_lines = 100
  character(len=200) :: actual_lines(max_lines)
  integer :: num_lines, i, iostat
  logical :: all_match
  
  ! Run simple.x and capture debug output
  call system('./build/simple.x | grep "DEBUG:" > current_debug_output.txt 2>&1')
  
  ! Read the current output
  open(unit=10, file='current_debug_output.txt', status='old', action='read')
  num_lines = 0
  do i = 1, max_lines
    read(10, '(A)', iostat=iostat) actual_lines(i)
    if (iostat /= 0) exit
    num_lines = i
  end do
  close(10)
  
  ! Compare with expected baseline
  all_match = .true.
  
  if (num_lines < size(expected_lines)) then
    print *, "ERROR: Too few output lines. Expected:", size(expected_lines), "Got:", num_lines
    all_match = .false.
  end if
  
  ! Compare the first key lines that should match exactly
  do i = 1, min(size(expected_lines), num_lines)
    if (trim(actual_lines(i)) /= trim(expected_lines(i))) then
      print *, "ERROR: Line", i, "does not match"
      print *, "Expected: '", trim(expected_lines(i)), "'"
      print *, "Actual:   '", trim(actual_lines(i)), "'"
      all_match = .false.
    end if
  end do
  
  if (all_match) then
    print *, "test_canonical_baseline_regression: PASSED"
    print *, "All", size(expected_lines), "baseline lines match current output"
  else
    print *, "test_canonical_baseline_regression: FAILED"
    error stop 1
  end if
  
  ! Clean up
  call system('rm -f current_debug_output.txt')
  
end program test_canonical_baseline_regression