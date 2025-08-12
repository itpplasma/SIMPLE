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
    real(dp), parameter :: max_elapsed = 0.1_dp   ! Maximum expected time (100ms)
    integer :: loop_count
    
    print *, "Testing time measurement functions..."
    
    ! Given: The timing module provides wall-clock time measurement
    ! When: We measure time during known operations
    ! Then: The elapsed time should be within expected bounds
    
    call init_timer()
    
    ! Test 1: Basic time measurement with controlled workload
    time1 = get_wtime()
    
    ! Do a calibrated amount of work
    loop_count = 0
    do while (loop_count < 100000)
      loop_count = loop_count + 1
    end do
    
    time2 = get_wtime()
    elapsed = time2 - time1
    
    ! Check that time advances
    if (elapsed <= 0.0_dp) then
      print *, "ERROR: Time should advance during computation"
      print *, "Elapsed:", elapsed
      errors = errors + 1
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
    else if (time2 - time1 > 0.001_dp) then
      print *, "WARNING: Large gap between successive get_wtime calls:", time2 - time1
    end if
    
    if (errors == 0) then
      print *, "  Time measurement test PASSED"
    end if
    
  end subroutine test_time_measurement
  
  subroutine test_elapsed_time_real_ops(errors)
    integer, intent(inout) :: errors
    real(dp) :: time1, time2, time3, elapsed, actual_elapsed
    integer :: i
    
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
    do i = 1, 1000
      do j = 1, 100
        dummy = dummy + sin(real(i*j, dp)) * cos(real(i+j, dp))
      end do
    end do
    ! Use the result to prevent optimization
    if (dummy > 1.0e10_dp) print *, "Note: dummy=", dummy
  end subroutine cpu_intensive_work

end program test_timing