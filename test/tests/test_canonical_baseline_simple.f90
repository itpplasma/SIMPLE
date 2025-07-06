program test_canonical_baseline_simple
  ! Simple baseline regression test that verifies key debug output lines
  
  implicit none
  
  integer :: iostat, line_count
  character(len=200) :: line
  logical :: found_start, found_ns, found_stencil, found_spline_start, found_spline_end
  
  found_start = .false.
  found_ns = .false.
  found_stencil = .false.
  found_spline_start = .false.
  found_spline_end = .false.
  line_count = 0
  
  ! Copy required input files and run simple.x 
  call system('cp ../../../simple.in .')
  call system('cp ../../../wout.nc .')
  call system('../../simple.x > temp_output.txt 2>&1')
  
  ! Read and check for key lines
  open(unit=10, file='temp_output.txt', status='old', action='read')
  
  do
    read(10, '(A)', iostat=iostat) line
    if (iostat /= 0) exit
    line_count = line_count + 1
    
    if (index(line, 'DEBUG: get_canonical_coordinates starting') > 0) then
      found_start = .true.
    end if
    
    if (index(line, 'DEBUG: ns_c=          50') > 0 .and. &
        index(line, 'n_theta_c=          36') > 0 .and. &
        index(line, 'n_phi_c=          81') > 0) then
      found_ns = .true.
    end if
    
    if (index(line, 'DEBUG: nh_stencil=           3') > 0) then
      found_stencil = .true.
    end if
    
    if (index(line, 'DEBUG: spline_can_coord starting') > 0) then
      found_spline_start = .true.
    end if
    
    if (index(line, 'DEBUG: spline_can_coord completed') > 0) then
      found_spline_end = .true.
    end if
  end do
  
  close(10)
  
  ! Check all required lines were found
  if (.not. found_start) then
    print *, "ERROR: Did not find 'get_canonical_coordinates starting'"
    error stop 1
  end if
  
  if (.not. found_ns) then
    print *, "ERROR: Did not find correct ns_c, n_theta_c, n_phi_c line"
    error stop 1
  end if
  
  if (.not. found_stencil) then
    print *, "ERROR: Did not find 'nh_stencil=3'"
    error stop 1
  end if
  
  if (.not. found_spline_start) then
    print *, "ERROR: Did not find 'spline_can_coord starting'"
    error stop 1
  end if
  
  if (.not. found_spline_end) then
    print *, "ERROR: Did not find 'spline_can_coord completed'"
    error stop 1
  end if
  
  print *, "test_canonical_baseline_simple: PASSED"
  print *, "Found all expected debug output lines"
  print *, "Total output lines:", line_count
  
  ! Clean up
  call system('rm -f temp_output.txt')
  
end program test_canonical_baseline_simple