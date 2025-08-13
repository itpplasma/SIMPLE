program test_timing
  use timing
  implicit none
  
  integer :: errors
  
  errors = 0
  
  ! Test timer initialization
  call test_timer_initialization(errors)
  
  ! Test time measurement functions
  call test_time_measurement(errors)
  
  ! Test elapsed time with real operations
  call test_elapsed_time_real_ops(errors)
  
  ! Test phase time printing (covers lines 34-49)
  call test_phase_time_printing(errors)
  
  if (errors == 0) then
    print *, "All timing module tests passed!"
  else
    print *, "ERROR: ", errors, " test(s) failed!"
    stop 1
  end if
  
contains

  subroutine test_timer_initialization(errors)
    integer, intent(inout) :: errors
    
    print *, "Testing timer initialization..."
    
    ! Given: The timing module provides timer initialization
    ! When: We call init_timer
    ! Then: The timer should be properly initialized with reasonable values
    
    call init_timer()
    
    ! Check that clock_rate is positive (system should have a working clock)
    if (clock_rate <= 0) then
      print *, "ERROR: Clock rate should be positive after initialization"
      print *, "Got:", clock_rate
      errors = errors + 1
    end if
    
    ! Check that clock_max is positive
    if (clock_max <= 0) then
      print *, "ERROR: Clock max should be positive after initialization"
      print *, "Got:", clock_max
      errors = errors + 1
    end if
    
    ! Check that program_start_time is non-negative
    if (program_start_time < 0) then
      print *, "ERROR: Program start time should be non-negative"
      print *, "Got:", program_start_time
      errors = errors + 1
    end if
    
    ! Check that phase_start_time is initialized to program_start_time
    if (phase_start_time /= program_start_time) then
      print *, "ERROR: Phase start time should equal program start time initially"
      print *, "Program start:", program_start_time
      print *, "Phase start:", phase_start_time
      errors = errors + 1
    end if
    
    if (errors == 0) then
      print *, "  Timer initialization test PASSED"
    end if
    
  end subroutine test_timer_initialization
  
  subroutine test_time_measurement(errors)
    integer, intent(inout) :: errors
    real(dp) :: time1, time2, elapsed
    real(dp), parameter :: max_elapsed = 1.0_dp   ! Maximum expected time (1 second)
    real(dp) :: dummy_work
    integer :: i, j
    
    print *, "Testing time measurement functions..."
    
    ! Given: The timing module provides wall-clock time measurement
    ! When: We measure time during known operations
    ! Then: The elapsed time should be within expected bounds
    
    call init_timer()
    
    ! Test 1: Basic time measurement with substantial workload
    time1 = get_wtime()
    
    ! Do substantial work that cannot be optimized away
    dummy_work = 0.0_dp
    do i = 1, 10000
      do j = 1, 100
        ! Use transcendental functions that are expensive and can't be optimized away
        dummy_work = dummy_work + sin(real(i, dp) * 0.001_dp) + cos(real(j, dp) * 0.001_dp)
      end do
    end do
    ! Force use of result to prevent optimization
    if (dummy_work > 1.0e20_dp .or. dummy_work < -1.0e20_dp) then
      print *, "Note: dummy_work=", dummy_work
    end if
    
    time2 = get_wtime()
    elapsed = time2 - time1
    
    ! Check that time advances properly
    ! With the substantial workload (10000x100 operations), we must detect non-zero elapsed time
    if (elapsed <= 0.0_dp) then
      print *, "ERROR: Timer failed to measure elapsed time for substantial workload"
      print *, "This indicates a serious timing system problem"
      print *, "Elapsed:", elapsed, "Time1:", time1, "Time2:", time2
      errors = errors + 1
    else if (elapsed < 1.0d-6) then
      print *, "WARNING: Very small elapsed time detected:", elapsed
      print *, "This may indicate insufficient workload or high-resolution timer"
    end if
    
    ! Check that elapsed time is within reasonable bounds
    if (elapsed > max_elapsed) then
      print *, "ERROR: Elapsed time too large for simple loop"
      print *, "Elapsed:", elapsed, "seconds (expected <", max_elapsed, ")"
      errors = errors + 1
    end if
    
    ! Test 2: Resolution test - rapid successive calls
    time1 = get_wtime()
    time2 = get_wtime()
    
    ! Successive calls should either be equal or show small increment
    if (time2 < time1) then
      print *, "ERROR: Time went backwards!"
      print *, "Time1:", time1, "Time2:", time2
      errors = errors + 1
    else if (time2 - time1 > 0.01_dp) then
      print *, "WARNING: Large gap between successive get_wtime calls:", time2 - time1
    end if
    
    ! Test 3: Additional stress test for timing precision
    if (elapsed > 0.0_dp .and. elapsed < 1.0d-4) then
      print *, "Running additional precision test with even more expensive operations..."
      time1 = get_wtime()
      
      ! Even more expensive operation to test timer precision limits
      dummy_work = 0.0_dp
      do i = 1, 50000
        do j = 1, 200
          dummy_work = dummy_work + sqrt(abs(real(i*j, dp)) + 1.0_dp) 
          dummy_work = dummy_work + log(max(1.0_dp, real(i+j, dp)))
          dummy_work = dummy_work + exp(min(1.0_dp, real(i, dp)*1.0d-6))
        end do
      end do
      ! Force use of result
      if (dummy_work > 1.0e30_dp .or. dummy_work < -1.0e30_dp) then
        print *, "Note: Extended precision test result:", dummy_work
      end if
      
      time2 = get_wtime()
      elapsed = time2 - time1
      
      if (elapsed <= 0.0_dp) then
        print *, "ERROR: Extended precision test still shows zero elapsed time"
        print *, "This indicates a fundamental timing system problem"
        print *, "Extended elapsed:", elapsed, "Time1:", time1, "Time2:", time2
        errors = errors + 1
      else
        print *, "Success: Extended precision test measured elapsed time:", elapsed
        
        ! Verify extended operation took longer than basic operation
        if (elapsed < 1.0d-5) then
          print *, "WARNING: Extended operation elapsed time unexpectedly small:", elapsed
        end if
      end if
    end if
    
    if (errors == 0) then
      print *, "  Time measurement test PASSED"
    end if
    
  end subroutine test_time_measurement
  
  subroutine test_elapsed_time_real_ops(errors)
    integer, intent(inout) :: errors
    real(dp) :: time1, time2, time3, elapsed, actual_elapsed
    
    print *, "Testing elapsed time with real operations..."
    
    ! Given: The timing module provides elapsed time calculation during real operations
    ! When: We measure elapsed time during actual computation
    ! Then: The elapsed time should match real wall-clock time progression
    
    call init_timer()
    
    ! Test 1: Elapsed time during real computation
    time1 = get_wtime()
    ! Do some actual work that takes measurable time
    call cpu_intensive_work()
    time2 = get_wtime()
    
    elapsed = get_elapsed_time(time1, time2)
    actual_elapsed = time2 - time1
    
    ! Check that get_elapsed_time returns consistent result with direct calculation
    if (abs(elapsed - actual_elapsed) > 1.0d-12) then
      print *, "ERROR: get_elapsed_time not consistent with direct calculation"
      print *, "Direct:", actual_elapsed, "Function:", elapsed
      errors = errors + 1
    end if
    
    ! Test 2: Monotonic time progression during sequential operations
    time1 = get_wtime()
    call cpu_intensive_work()
    time2 = get_wtime()
    call cpu_intensive_work()
    time3 = get_wtime()
    
    ! Time should be strictly increasing
    if (time1 >= time2 .or. time2 >= time3) then
      print *, "ERROR: Time not monotonically increasing"
      print *, "Times:", time1, time2, time3
      errors = errors + 1
    end if
    
    ! Test 3: Check that elapsed time is positive and reasonable
    elapsed = get_elapsed_time(time1, time3)
    if (elapsed <= 0.0_dp) then
      print *, "ERROR: Total elapsed time should be positive"
      print *, "Elapsed:", elapsed
      errors = errors + 1
    else if (elapsed > 1.0_dp) then
      print *, "ERROR: Elapsed time unreasonably large"
      print *, "Elapsed:", elapsed
      errors = errors + 1
    end if
    
    if (errors == 0) then
      print *, "  Elapsed time real ops test PASSED"
    end if
    
  end subroutine test_elapsed_time_real_ops
  
  ! Helper subroutine to introduce a measurable but bounded delay
  subroutine cpu_intensive_work()
    integer :: i, j
    real(dp) :: dummy
    dummy = 0.0_dp
    ! Perform floating point operations that can't be easily optimized away
    ! Increased by 10x to make timing more visible (was 1000 x 100)
    do i = 1, 10000
      do j = 1, 100
        dummy = dummy + sin(real(i*j, dp)) * cos(real(i+j, dp))
      end do
    end do
    ! Use the result to prevent optimization
    if (dummy > 1.0e10_dp) print *, "Note: dummy=", dummy
  end subroutine cpu_intensive_work

  subroutine test_phase_time_printing(errors)
    integer, intent(inout) :: errors
    integer(kind=8) :: start_phase, end_phase, saved_start, saved_phase
    
    print *, "Testing phase time printing (covers lines 34-49)..."
    
    ! Given: The timing module provides phase timing output
    ! When: We call print_phase_time with different phase names
    ! Then: It should print timing information and update phase_start_time
    
    call init_timer()
    
    ! Save original timing state
    saved_start = program_start_time
    saved_phase = phase_start_time
    
    ! Test 1: Basic phase timing with work between phases
    print *, "  Testing basic phase timing output..."
    
    call cpu_intensive_work()  ! Do some work before first phase
    call print_phase_time("Phase A")
    
    ! Check that phase_start_time was updated (should be different from program start)
    if (phase_start_time == saved_start) then
      print *, "ERROR: phase_start_time should be updated after print_phase_time"
      errors = errors + 1
    end if
    
    start_phase = phase_start_time
    call cpu_intensive_work()  ! Do more work
    call print_phase_time("Phase B")
    
    ! Check that phase_start_time was updated again
    if (phase_start_time == start_phase) then
      print *, "ERROR: phase_start_time should be updated after second print_phase_time"
      errors = errors + 1
    end if
    
    ! Test 2: Phase timing with longer phase name
    call print_phase_time("Long Phase Name With Spaces")
    
    ! Test 3: Phase timing with special characters (but not problematic ones)
    call print_phase_time("Phase_3-Final")
    
    ! Test 4: Minimal phase timing (empty work)
    call print_phase_time("Quick Phase")
    
    ! Test 5: Check phase progression is monotonic
    start_phase = phase_start_time
    call cpu_intensive_work()
    call print_phase_time("Final Phase")
    end_phase = phase_start_time
    
    if (end_phase <= start_phase) then
      print *, "ERROR: Phase times should be monotonic"
      print *, "Start phase:", start_phase, "End phase:", end_phase
      errors = errors + 1
    end if
    
    if (errors == 0) then
      print *, "  Phase time printing test PASSED"
    end if
    
  end subroutine test_phase_time_printing

end program test_timing