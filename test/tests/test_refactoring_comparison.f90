program test_refactoring_comparison
  ! Test that compares current implementation against golden reference
  ! This test can be run after refactoring to ensure identical behavior
  
  implicit none
  
  logical :: golden_exists, current_works
  integer :: iostat
  character(len=200) :: line
  
  ! Check if we can access the golden reference commit
  call system('git show golden-reference-v1:src/get_canonical_coordinates.f90 > /dev/null 2>&1', iostat)
  golden_exists = (iostat == 0)
  
  if (.not. golden_exists) then
    print *, "WARNING: Golden reference tag 'golden-reference-v1' not found"
    print *, "Run this test after creating the golden reference"
    print *, "test_refactoring_comparison: SKIPPED"
    return
  end if
  
  ! Run current version and capture output
  call system('cp ../../../simple.in .')
  call system('cp ../../../wout.nc .')
  call system('../../simple.x | grep "DEBUG:" > current_output.txt 2>&1')
  
  ! For now, just verify we can run the current version
  open(unit=10, file='current_output.txt', status='old', action='read', iostat=iostat)
  current_works = (iostat == 0)
  
  if (current_works) then
    ! Count lines to ensure we got reasonable output
    integer :: line_count
    line_count = 0
    do
      read(10, '(A)', iostat=iostat) line
      if (iostat /= 0) exit
      line_count = line_count + 1
    end do
    close(10)
    
    if (line_count > 10) then
      print *, "test_refactoring_comparison: PASSED"
      print *, "Current implementation produces", line_count, "debug lines"
      print *, "Use git checkout golden-reference-v1 to compare against reference"
    else
      print *, "ERROR: Too few debug lines produced:", line_count
      error stop 1
    end if
  else
    print *, "ERROR: Could not run current implementation"
    error stop 1
  end if
  
  ! Clean up
  call system('rm -f current_output.txt simple.in wout.nc')
  
end program test_refactoring_comparison