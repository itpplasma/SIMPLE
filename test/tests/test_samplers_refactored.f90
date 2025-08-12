program test_samplers_refactored
  use samplers
  use util, only: pi, twopi
  implicit none

  integer :: num_passed, num_failed
  
  num_passed = 0
  num_failed = 0
  
  print *, "========================================="
  print *, "Testing Refactored Samplers Module"
  print *, "========================================="
  
  call test_random_in_range(num_passed, num_failed)
  call test_random_pitch(num_passed, num_failed)
  call test_random_angles(num_passed, num_failed)
  call test_initialize_particle_velocity(num_passed, num_failed)
  call test_write_read_particles(num_passed, num_failed)
  call test_compute_b_extrema(num_passed, num_failed)
  call test_calculate_grid_params(num_passed, num_failed)
  call test_reallocate_grid_arrays(num_passed, num_failed)
  
  print *, "========================================="
  print *, "Test Summary:"
  print *, "  Passed: ", num_passed
  print *, "  Failed: ", num_failed
  print *, "========================================="
  
  if (num_failed > 0) then
    error stop "Some tests failed"
  end if

contains

  subroutine test_random_in_range(num_passed, num_failed)
    integer, intent(inout) :: num_passed, num_failed
    real(dp) :: value
    real(dp) :: min_val, max_val
    integer :: i
    logical :: test_passed
    
    print *, ""
    print *, "Test: random_in_range"
    print *, "  Given: Min and max values for range"
    print *, "  When: Generating random values"
    print *, "  Then: Values should be within specified range"
    
    min_val = 0.2_dp
    max_val = 0.8_dp
    test_passed = .true.
    
    ! Generate multiple samples to test range
    do i = 1, 100
      value = random_in_range(min_val, max_val)
      if (value < min_val .or. value > max_val) then
        test_passed = .false.
        print *, "    ERROR: Value out of range: ", value
        exit
      end if
    end do
    
    if (test_passed) then
      print *, "    ✓ Random values within range"
      num_passed = num_passed + 1
    else
      print *, "    ✗ Random values out of range"
      num_failed = num_failed + 1
    end if
  end subroutine test_random_in_range
  
  subroutine test_random_pitch(num_passed, num_failed)
    integer, intent(inout) :: num_passed, num_failed
    real(dp) :: pitch
    integer :: i
    logical :: test_passed
    
    print *, ""
    print *, "Test: random_pitch"
    print *, "  Given: No input parameters"
    print *, "  When: Generating random pitch angles"
    print *, "  Then: Values should be in range [-1, 1]"
    
    test_passed = .true.
    
    ! Generate multiple samples to test range
    do i = 1, 100
      pitch = random_pitch()
      if (pitch < -1.0_dp .or. pitch > 1.0_dp) then
        test_passed = .false.
        print *, "    ERROR: Pitch out of range: ", pitch
        exit
      end if
    end do
    
    if (test_passed) then
      print *, "    ✓ Pitch values within [-1, 1]"
      num_passed = num_passed + 1
    else
      print *, "    ✗ Pitch values out of range"
      num_failed = num_failed + 1
    end if
  end subroutine test_random_pitch
  
  subroutine test_random_angles(num_passed, num_failed)
    integer, intent(inout) :: num_passed, num_failed
    real(dp) :: theta, phi
    integer :: i
    logical :: test_passed
    
    print *, ""
    print *, "Test: random_angles"
    print *, "  Given: No input parameters"
    print *, "  When: Generating random angles"
    print *, "  Then: Both angles should be in range [0, 2π]"
    
    test_passed = .true.
    
    ! Generate multiple samples to test range
    do i = 1, 100
      call random_angles(theta, phi)
      if (theta < 0.0_dp .or. theta > twopi .or. &
          phi < 0.0_dp .or. phi > twopi) then
        test_passed = .false.
        print *, "    ERROR: Angles out of range: theta=", theta, " phi=", phi
        exit
      end if
    end do
    
    if (test_passed) then
      print *, "    ✓ Angles within [0, 2π]"
      num_passed = num_passed + 1
    else
      print *, "    ✗ Angles out of range"
      num_failed = num_failed + 1
    end if
  end subroutine test_random_angles
  
  subroutine test_initialize_particle_velocity(num_passed, num_failed)
    integer, intent(inout) :: num_passed, num_failed
    real(dp), dimension(5) :: z
    logical :: test_passed
    
    print *, ""
    print *, "Test: initialize_particle_velocity"
    print *, "  Given: Particle state vector"
    print *, "  When: Initializing velocity components"
    print *, "  Then: Should set normalized velocity and pitch correctly"
    
    test_passed = .true.
    
    ! Test with default values
    z = 0.0_dp
    call initialize_particle_velocity(z)
    
    if (abs(z(4) - 1.0_dp) > 1e-10) then
      test_passed = .false.
      print *, "    ERROR: Default normalized velocity not 1.0: ", z(4)
    else if (z(5) < -1.0_dp .or. z(5) > 1.0_dp) then
      test_passed = .false.
      print *, "    ERROR: Random pitch out of range: ", z(5)
    else
      print *, "    ✓ Default initialization correct"
    end if
    
    ! Test with specified values
    z = 0.0_dp
    call initialize_particle_velocity(z, normalized_v=2.0_dp, pitch=0.5_dp)
    
    if (abs(z(4) - 2.0_dp) > 1e-10) then
      test_passed = .false.
      print *, "    ERROR: Specified normalized velocity incorrect: ", z(4)
    else if (abs(z(5) - 0.5_dp) > 1e-10) then
      test_passed = .false.
      print *, "    ERROR: Specified pitch incorrect: ", z(5)
    else
      print *, "    ✓ Specified initialization correct"
    end if
    
    if (test_passed) then
      num_passed = num_passed + 1
    else
      num_failed = num_failed + 1
    end if
  end subroutine test_initialize_particle_velocity
  
  subroutine test_write_read_particles(num_passed, num_failed)
    integer, intent(inout) :: num_passed, num_failed
    real(dp), dimension(5, 3) :: particles_out, particles_in
    integer :: i, j
    logical :: test_passed
    
    print *, ""
    print *, "Test: write_particles_to_file and read_particles_from_file"
    print *, "  Given: Particle array to write"
    print *, "  When: Writing and reading back from file"
    print *, "  Then: Should recover original data"
    
    test_passed = .true.
    
    ! Initialize test data
    do i = 1, 5
      do j = 1, 3
        particles_out(i, j) = real(i, dp) * 10.0_dp + real(j, dp)
      end do
    end do
    
    ! Write to file
    call write_particles_to_file('test_particles.dat', particles_out)
    
    ! Read back from file
    call read_particles_from_file('test_particles.dat', particles_in)
    
    ! Compare
    do i = 1, 5
      do j = 1, 3
        if (abs(particles_in(i, j) - particles_out(i, j)) > 1e-10) then
          test_passed = .false.
          print *, "    ERROR: Mismatch at (", i, ",", j, "): ", &
                   particles_in(i, j), " != ", particles_out(i, j)
        end if
      end do
    end do
    
    if (test_passed) then
      print *, "    ✓ File I/O preserves particle data"
      num_passed = num_passed + 1
    else
      print *, "    ✗ File I/O corrupts particle data"
      num_failed = num_failed + 1
    end if
    
    ! Clean up
    call execute_command_line('rm -f test_particles.dat')
  end subroutine test_write_read_particles
  
  subroutine test_compute_b_extrema(num_passed, num_failed)
    integer, intent(inout) :: num_passed, num_failed
    real(dp), dimension(5) :: bstart_test
    real(dp) :: bmin_test, bmax_test
    logical :: test_passed
    
    print *, ""
    print *, "Test: compute_b_extrema"
    print *, "  Given: Array of B field values"
    print *, "  When: Computing min and max"
    print *, "  Then: Should correctly identify extrema"
    
    test_passed = .true.
    
    bstart_test = [1.0_dp, 3.0_dp, 0.5_dp, 2.0_dp, 1.5_dp]
    
    call compute_b_extrema(bstart_test, bmin_test, bmax_test)
    
    if (abs(bmin_test - 0.5_dp) > 1e-10) then
      test_passed = .false.
      print *, "    ERROR: Incorrect bmin: ", bmin_test, " expected 0.5"
    else if (abs(bmax_test - 3.0_dp) > 1e-10) then
      test_passed = .false.
      print *, "    ERROR: Incorrect bmax: ", bmax_test, " expected 3.0"
    else
      print *, "    ✓ B field extrema computed correctly"
    end if
    
    if (test_passed) then
      num_passed = num_passed + 1
    else
      num_failed = num_failed + 1
    end if
  end subroutine test_compute_b_extrema
  
  subroutine test_calculate_grid_params(num_passed, num_failed)
    integer, intent(inout) :: num_passed, num_failed
    real(dp) :: grid_density, xsize
    integer :: ngrid, ntestpart_calc
    logical :: test_passed
    
    print *, ""
    print *, "Test: calculate_grid_params"
    print *, "  Given: Grid density parameter"
    print *, "  When: Calculating grid parameters"
    print *, "  Then: Should compute correct grid size and particle count"
    
    test_passed = .true.
    
    grid_density = 0.25_dp
    
    call calculate_grid_params(grid_density, xsize, ngrid, ntestpart_calc)
    
    ! Check expected values
    if (abs(xsize - 2.0_dp * pi * grid_density) > 1e-10) then
      test_passed = .false.
      print *, "    ERROR: Incorrect xsize: ", xsize
    else if (ngrid /= 3) then  ! int(1.0/0.25 - 1.0) = 3
      test_passed = .false.
      print *, "    ERROR: Incorrect ngrid: ", ngrid, " expected 3"
    else if (ntestpart_calc /= 9) then  ! 3 * 3 = 9
      test_passed = .false.
      print *, "    ERROR: Incorrect ntestpart: ", ntestpart_calc, " expected 9"
    else
      print *, "    ✓ Grid parameters calculated correctly"
    end if
    
    if (test_passed) then
      num_passed = num_passed + 1
    else
      num_failed = num_failed + 1
    end if
  end subroutine test_calculate_grid_params
  
  subroutine test_reallocate_grid_arrays(num_passed, num_failed)
    use params, only: zstart_dim1, zend, times_lost, trap_par, perp_inv, iclass
    integer, intent(inout) :: num_passed, num_failed
    real(dp), dimension(:,:), allocatable :: zstart_test
    integer :: new_size
    logical :: test_passed
    
    print *, ""
    print *, "Test: reallocate_grid_arrays"
    print *, "  Given: New particle count"
    print *, "  When: Reallocating arrays"
    print *, "  Then: Should allocate arrays with correct dimensions"
    
    test_passed = .true.
    new_size = 16
    
    ! Ensure arrays are deallocated first
    if (allocated(zstart_test)) deallocate(zstart_test)
    if (allocated(zend)) deallocate(zend)
    if (allocated(times_lost)) deallocate(times_lost)
    if (allocated(trap_par)) deallocate(trap_par)
    if (allocated(perp_inv)) deallocate(perp_inv)
    if (allocated(iclass)) deallocate(iclass)
    
    ! Set zstart_dim1 for the test
    zstart_dim1 = 5
    
    call reallocate_grid_arrays(zstart_test, new_size)
    
    ! Check allocations
    if (.not. allocated(zstart_test)) then
      test_passed = .false.
      print *, "    ERROR: zstart not allocated"
    else if (size(zstart_test, 1) /= 5 .or. size(zstart_test, 2) /= new_size) then
      test_passed = .false.
      print *, "    ERROR: zstart dimensions incorrect: ", &
               size(zstart_test, 1), "x", size(zstart_test, 2)
    end if
    
    if (.not. allocated(zend)) then
      test_passed = .false.
      print *, "    ERROR: zend not allocated"
    else if (size(zend, 2) /= new_size) then
      test_passed = .false.
      print *, "    ERROR: zend size incorrect: ", size(zend, 2)
    end if
    
    if (.not. allocated(times_lost)) then
      test_passed = .false.
      print *, "    ERROR: times_lost not allocated"
    else if (size(times_lost) /= new_size) then
      test_passed = .false.
      print *, "    ERROR: times_lost size incorrect: ", size(times_lost)
    end if
    
    if (test_passed) then
      print *, "    ✓ Arrays reallocated correctly"
      num_passed = num_passed + 1
    else
      print *, "    ✗ Array reallocation failed"
      num_failed = num_failed + 1
    end if
    
    ! Clean up
    if (allocated(zstart_test)) deallocate(zstart_test)
    if (allocated(zend)) deallocate(zend)
    if (allocated(times_lost)) deallocate(times_lost)
    if (allocated(trap_par)) deallocate(trap_par)
    if (allocated(perp_inv)) deallocate(perp_inv)
    if (allocated(iclass)) deallocate(iclass)
  end subroutine test_reallocate_grid_arrays
  
end program test_samplers_refactored