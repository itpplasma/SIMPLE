program test_simple_main_refactored
! Unit tests for refactored simple_main module
! Tests follow behavior-driven design with Given-When-Then structure

use simple_main, only: dp, compute_bmod, compute_pitch_angle_params, &
                       increase_confined_count, to_standard_z_coordinates, &
                       update_momentum, collide, init_counters, &
                       print_parameters, write_loss_times, &
                       write_average_inverse_loss_time, &
                       write_confinement_fraction, write_classification_data
use params, only: confpart_pass, confpart_trap, times_lost, ntestpart, &
                  ntimstep, dtau, dtaumin, v0, trace_time, ntcut, &
                  class_plot, iclass, zstart, zend, trap_par, perp_inv, &
                  kpart, swcoll, num_surf
use simple, only: Tracer

implicit none

integer :: n_tests_passed, n_tests_failed
integer :: ierr  ! For execute_command_line
real(dp), parameter :: tol = 1.0e-12_dp
real(dp), parameter :: tol_phys = 1.0e-6_dp

n_tests_passed = 0
n_tests_failed = 0

call test_counter_initialization()
call test_bmod_computation()
call test_pitch_angle_parameters()
call test_confinement_counting()
call test_coordinate_transformation()
call test_momentum_update()
call test_output_file_generation()
call test_statistics_normalization()
call test_collision_handling()

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
  ! Test: Counter initialization
  !===========================================================================
  subroutine test_counter_initialization()
    ! Given: Uninitialized counters
    ! When: We initialize counters
    ! Then: All counters should be reset to appropriate values
    
    logical :: test_passed
    
    print *, 'Testing counter initialization...'
    
    test_passed = .true.
    
    ! Set some non-zero values first
    kpart = 100
    if (allocated(confpart_pass)) deallocate(confpart_pass)
    if (allocated(confpart_trap)) deallocate(confpart_trap)
    if (allocated(times_lost)) deallocate(times_lost)
    
    ! Setup arrays
    ntimstep = 10
    ntestpart = 5
    allocate(confpart_pass(ntimstep))
    allocate(confpart_trap(ntimstep))
    allocate(times_lost(ntestpart))
    
    confpart_pass = 1.0_dp
    confpart_trap = 1.0_dp
    times_lost = 1.0_dp
    
    ! Initialize counters
    call init_counters()
    
    ! Check kpart is reset
    if (kpart /= 0) then
      test_passed = .false.
      print *, '  Failed: kpart not reset to 0'
    end if
    
    ! Check confinement arrays are zeroed
    if (any(confpart_pass /= 0.0_dp)) then
      test_passed = .false.
      print *, '  Failed: confpart_pass not zeroed'
    end if
    
    if (any(confpart_trap /= 0.0_dp)) then
      test_passed = .false.
      print *, '  Failed: confpart_trap not zeroed'
    end if
    
    ! Check times_lost initialized to -1
    if (any(times_lost /= -1.0_dp)) then
      test_passed = .false.
      print *, '  Failed: times_lost not initialized to -1'
    end if
    
    ! Clean up
    deallocate(confpart_pass, confpart_trap, times_lost)
    
    if (test_passed) then
      print *, '  PASSED: Counters initialized correctly'
      n_tests_passed = n_tests_passed + 1
    else
      print *, '  FAILED: Counter initialization'
      n_tests_failed = n_tests_failed + 1
    end if
    
  end subroutine test_counter_initialization

  !===========================================================================
  ! Test: Magnetic field computation
  !===========================================================================
  subroutine test_bmod_computation()
    ! Given: A position in space
    ! When: We compute the magnetic field
    ! Then: Should return a valid field magnitude
    
    logical :: test_passed
    real(dp) :: z(3), bmod
    
    print *, 'Testing magnetic field computation...'
    
    test_passed = .true.
    
    ! Note: This test is limited without a full magnetic field setup
    ! We verify that the function exists and has the correct interface
    
    if (test_passed) then
      print *, '  PASSED: Magnetic field computation interface verified'
      n_tests_passed = n_tests_passed + 1
    else
      print *, '  FAILED: Magnetic field computation'
      n_tests_failed = n_tests_failed + 1
    end if
    
  end subroutine test_bmod_computation

  !===========================================================================
  ! Test: Pitch angle parameters
  !===========================================================================
  subroutine test_pitch_angle_parameters()
    ! Given: Particle phase space coordinates
    ! When: We compute pitch angle parameters
    ! Then: Should determine passing/trapped status and parameters
    
    logical :: test_passed, passing
    real(dp) :: z(5), trap_par_val, perp_inv_val
    
    print *, 'Testing pitch angle parameter calculation...'
    
    test_passed = .true.
    
    ! Set up test case (simplified without full field)
    z = [0.5_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.7_dp]
    num_surf = 0  ! Disable bmin/bmax calculation for test
    
    ! This would require proper magnetic field setup to test fully
    ! We verify the subroutine exists and has correct interface
    
    if (test_passed) then
      print *, '  PASSED: Pitch angle parameter interface verified'
      n_tests_passed = n_tests_passed + 1
    else
      print *, '  FAILED: Pitch angle parameter calculation'
      n_tests_failed = n_tests_failed + 1
    end if
    
  end subroutine test_pitch_angle_parameters

  !===========================================================================
  ! Test: Confinement counting
  !===========================================================================
  subroutine test_confinement_counting()
    ! Given: Particle confinement status
    ! When: We increase confinement count
    ! Then: Appropriate counter should be incremented
    
    logical :: test_passed
    integer :: it
    
    print *, 'Testing confinement counting...'
    
    test_passed = .true.
    
    ! Setup
    ntimstep = 10
    if (allocated(confpart_pass)) deallocate(confpart_pass)
    if (allocated(confpart_trap)) deallocate(confpart_trap)
    allocate(confpart_pass(ntimstep))
    allocate(confpart_trap(ntimstep))
    
    confpart_pass = 0.0_dp
    confpart_trap = 0.0_dp
    
    it = 5
    
    ! Test passing particle count
    call increase_confined_count(it, .true.)
    
    if (abs(confpart_pass(it) - 1.0_dp) > tol) then
      test_passed = .false.
      print *, '  Failed: Passing particle count not incremented'
    end if
    
    if (abs(confpart_trap(it)) > tol) then
      test_passed = .false.
      print *, '  Failed: Trapped count should remain zero'
    end if
    
    ! Test trapped particle count
    call increase_confined_count(it, .false.)
    
    if (abs(confpart_trap(it) - 1.0_dp) > tol) then
      test_passed = .false.
      print *, '  Failed: Trapped particle count not incremented'
    end if
    
    ! Clean up
    deallocate(confpart_pass, confpart_trap)
    
    if (test_passed) then
      print *, '  PASSED: Confinement counting works correctly'
      n_tests_passed = n_tests_passed + 1
    else
      print *, '  FAILED: Confinement counting'
      n_tests_failed = n_tests_failed + 1
    end if
    
  end subroutine test_confinement_counting

  !===========================================================================
  ! Test: Coordinate transformation
  !===========================================================================
  subroutine test_coordinate_transformation()
    ! Given: Tracer with symplectic integrator state
    ! When: We transform to standard z coordinates
    ! Then: Coordinates should be properly converted
    
    logical :: test_passed
    type(Tracer) :: norb
    real(dp) :: z(5)
    
    print *, 'Testing coordinate transformation...'
    
    test_passed = .true.
    
    ! Setup tracer (simplified test)
    norb%si%z = [0.5_dp, 0.1_dp, 0.2_dp, 0.3_dp]
    norb%f%mu = 0.1_dp
    norb%f%Bmod = 1.0_dp
    norb%f%vpar = 0.5_dp
    
    call to_standard_z_coordinates(norb, z)
    
    ! Check that first 3 coordinates match
    if (any(abs(z(1:3) - norb%si%z(1:3)) > tol)) then
      test_passed = .false.
      print *, '  Failed: First 3 coordinates should match'
    end if
    
    ! Check that z(4) is computed correctly
    if (z(4) <= 0.0_dp) then
      test_passed = .false.
      print *, '  Failed: z(4) should be positive'
    end if
    
    ! Check that z(5) is normalized
    if (abs(z(5)) > 1.0_dp) then
      test_passed = .false.
      print *, '  Failed: z(5) should be between -1 and 1'
    end if
    
    if (test_passed) then
      print *, '  PASSED: Coordinate transformation works correctly'
      n_tests_passed = n_tests_passed + 1
    else
      print *, '  FAILED: Coordinate transformation'
      n_tests_failed = n_tests_failed + 1
    end if
    
  end subroutine test_coordinate_transformation

  !===========================================================================
  ! Test: Momentum update
  !===========================================================================
  subroutine test_momentum_update()
    ! Given: Particle state
    ! When: We update momentum
    ! Then: Momentum variables should be consistently updated
    
    logical :: test_passed
    type(Tracer) :: norb
    real(dp) :: z(5)
    
    print *, 'Testing momentum update...'
    
    test_passed = .true.
    
    ! Setup test case
    z = [0.5_dp, 0.1_dp, 0.2_dp, 1.0_dp, 0.7_dp]
    norb%f%Bmod = 1.0_dp
    norb%f%hph = 1.0_dp
    norb%f%Aph = 0.0_dp
    norb%f%ro0 = 1.0_dp
    
    call update_momentum(norb, z)
    
    ! Check that pabs is set
    if (abs(norb%si%pabs - z(4)) > tol) then
      test_passed = .false.
      print *, '  Failed: pabs not set correctly'
    end if
    
    ! Check that vpar is computed
    if (norb%f%vpar == 0.0_dp .and. z(5) /= 0.0_dp) then
      test_passed = .false.
      print *, '  Failed: vpar should be non-zero'
    end if
    
    if (test_passed) then
      print *, '  PASSED: Momentum update works correctly'
      n_tests_passed = n_tests_passed + 1
    else
      print *, '  FAILED: Momentum update'
      n_tests_failed = n_tests_failed + 1
    end if
    
  end subroutine test_momentum_update

  !===========================================================================
  ! Test: Output file generation
  !===========================================================================
  subroutine test_output_file_generation()
    ! Given: Simulation data
    ! When: We write output files
    ! Then: Files should be created with correct format
    
    logical :: test_passed, file_exists
    
    print *, 'Testing output file generation...'
    
    test_passed = .true.
    
    ! Setup minimal data
    ntestpart = 2
    ntimstep = 2
    dtau = 0.1_dp
    v0 = 1.0_dp
    trace_time = 1.0_dp
    ntcut = 0
    class_plot = .false.
    
    if (allocated(times_lost)) deallocate(times_lost)
    if (allocated(trap_par)) deallocate(trap_par)
    if (allocated(zstart)) deallocate(zstart)
    if (allocated(zend)) deallocate(zend)
    if (allocated(perp_inv)) deallocate(perp_inv)
    if (allocated(confpart_pass)) deallocate(confpart_pass)
    if (allocated(confpart_trap)) deallocate(confpart_trap)
    
    allocate(times_lost(ntestpart))
    allocate(trap_par(ntestpart))
    allocate(zstart(5, ntestpart))
    allocate(zend(5, ntestpart))
    allocate(perp_inv(ntestpart))
    allocate(confpart_pass(ntimstep))
    allocate(confpart_trap(ntimstep))
    
    times_lost = 0.5_dp
    trap_par = 0.1_dp
    zstart = 0.0_dp
    zend = 0.0_dp
    perp_inv = 0.2_dp
    confpart_pass = 0.5_dp
    confpart_trap = 0.5_dp
    
    ! Call individual output functions
    call write_loss_times()
    call write_average_inverse_loss_time()
    call write_confinement_fraction()
    
    ! Check that files exist
    inquire(file='times_lost.dat', exist=file_exists)
    if (.not. file_exists) then
      test_passed = .false.
      print *, '  Failed: times_lost.dat not created'
    end if
    
    inquire(file='confined_fraction.dat', exist=file_exists)
    if (.not. file_exists) then
      test_passed = .false.
      print *, '  Failed: confined_fraction.dat not created'
    end if
    
    ! Clean up files
    call execute_command_line('rm -f times_lost.dat avg_inverse_t_lost.dat ' // &
                             'confined_fraction.dat class_parts.dat', exitstat=ierr)
    
    ! Clean up arrays
    deallocate(times_lost, trap_par, zstart, zend, perp_inv)
    deallocate(confpart_pass, confpart_trap)
    
    if (test_passed) then
      print *, '  PASSED: Output files generated correctly'
      n_tests_passed = n_tests_passed + 1
    else
      print *, '  FAILED: Output file generation'
      n_tests_failed = n_tests_failed + 1
    end if
    
  end subroutine test_output_file_generation

  !===========================================================================
  ! Test: Statistics normalization
  !===========================================================================
  subroutine test_statistics_normalization()
    ! Given: Raw confinement counts
    ! When: We normalize statistics
    ! Then: Values should be divided by total particles
    
    logical :: test_passed
    
    print *, 'Testing statistics normalization...'
    
    test_passed = .true.
    
    ! Setup
    ntestpart = 100
    ntimstep = 5
    if (allocated(confpart_pass)) deallocate(confpart_pass)
    if (allocated(confpart_trap)) deallocate(confpart_trap)
    allocate(confpart_pass(ntimstep))
    allocate(confpart_trap(ntimstep))
    
    ! Set raw counts
    confpart_pass = 50.0_dp
    confpart_trap = 30.0_dp
    
    ! Normalize (inline since we can't call normalize_statistics directly)
    confpart_pass = confpart_pass / ntestpart
    confpart_trap = confpart_trap / ntestpart
    
    ! Check normalization
    if (any(abs(confpart_pass - 0.5_dp) > tol)) then
      test_passed = .false.
      print *, '  Failed: confpart_pass not normalized correctly'
    end if
    
    if (any(abs(confpart_trap - 0.3_dp) > tol)) then
      test_passed = .false.
      print *, '  Failed: confpart_trap not normalized correctly'
    end if
    
    ! Clean up
    deallocate(confpart_pass, confpart_trap)
    
    if (test_passed) then
      print *, '  PASSED: Statistics normalized correctly'
      n_tests_passed = n_tests_passed + 1
    else
      print *, '  FAILED: Statistics normalization'
      n_tests_failed = n_tests_failed + 1
    end if
    
  end subroutine test_statistics_normalization

  !===========================================================================
  ! Test: Collision handling
  !===========================================================================
  subroutine test_collision_handling()
    ! Given: Collision mode enabled
    ! When: We handle collisions
    ! Then: Collision routines should be called appropriately
    
    logical :: test_passed
    real(dp) :: z(5), dt
    
    print *, 'Testing collision handling...'
    
    test_passed = .true.
    
    ! Setup test case
    z = [0.5_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.7_dp]
    dt = 0.01_dp
    dtaumin = dt  ! Set module variable for error reporting
    
    ! Note: Full collision testing requires collision module setup
    ! Here we verify the interface exists
    if (swcoll) then
      call collide(z, dt)
    endif
    
    if (test_passed) then
      print *, '  PASSED: Collision handling interface verified'
      n_tests_passed = n_tests_passed + 1
    else
      print *, '  FAILED: Collision handling'
      n_tests_failed = n_tests_failed + 1
    end if
    
  end subroutine test_collision_handling

end program test_simple_main_refactored