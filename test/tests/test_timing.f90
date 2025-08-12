program test_timing
  use timing
  implicit none
  
  integer :: errors
  
  errors = 0
  
  ! Test timer initialization
  call test_timer_initialization(errors)
  
  ! Test time measurement functions
  call test_time_measurement(errors)
  
  ! Test elapsed time calculation
  call test_elapsed_time_calculation(errors)
  
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
    real(dp) :: time1, time2
    real(dp), parameter :: tolerance = 1.0d-6  ! 1 microsecond tolerance
    real(dp), parameter :: sleep_time = 0.01_dp  ! Target sleep time in seconds
    
    print *, "Testing time measurement functions..."
    
    ! Given: The timing module provides wall-clock time measurement
    ! When: We measure time before and after a delay
    ! Then: The elapsed time should be reasonable
    
    call init_timer()
    
    ! Test get_wtime function
    time1 = get_wtime()
    
    ! Introduce a small delay by doing some work
    call cpu_intensive_work()
    
    time2 = get_wtime()
    
    ! Check that time advances
    if (time2 <= time1) then
      print *, "ERROR: Time should advance between measurements"
      print *, "Time1:", time1, "Time2:", time2
      errors = errors + 1
    end if
    
    ! Check that elapsed time is reasonable (should be small but positive)
    if (time2 - time1 > 1.0_dp) then
      print *, "ERROR: Elapsed time too large for simple operation"
      print *, "Elapsed:", time2 - time1
      errors = errors + 1
    end if
    
    ! Test that times are positive (wall-clock time should be positive)
    if (time1 < 0.0_dp .or. time2 < 0.0_dp) then
      print *, "ERROR: Wall-clock times should be positive"
      print *, "Time1:", time1, "Time2:", time2
      errors = errors + 1
    end if
    
    if (errors == 0) then
      print *, "  Time measurement test PASSED"
    end if
    
  end subroutine test_time_measurement
  
  subroutine test_elapsed_time_calculation(errors)
    integer, intent(inout) :: errors
    real(dp) :: start_time, end_time, elapsed
    real(dp), parameter :: tolerance = 1.0d-12
    
    print *, "Testing elapsed time calculation..."
    
    ! Given: The timing module provides elapsed time calculation
    ! When: We calculate elapsed time between two time points
    ! Then: The result should equal the simple difference
    
    start_time = 100.0_dp
    end_time = 150.5_dp
    
    elapsed = get_elapsed_time(start_time, end_time)
    
    ! Check that elapsed time equals simple difference
    if (abs(elapsed - (end_time - start_time)) > tolerance) then
      print *, "ERROR: Elapsed time calculation incorrect"
      print *, "Expected:", end_time - start_time
      print *, "Got:", elapsed
      errors = errors + 1
    end if
    
    ! Test with negative elapsed time (end before start)
    elapsed = get_elapsed_time(end_time, start_time)
    if (abs(elapsed - (start_time - end_time)) > tolerance) then
      print *, "ERROR: Negative elapsed time calculation incorrect"
      print *, "Expected:", start_time - end_time
      print *, "Got:", elapsed
      errors = errors + 1
    end if
    
    ! Test with zero elapsed time
    elapsed = get_elapsed_time(start_time, start_time)
    if (abs(elapsed) > tolerance) then
      print *, "ERROR: Zero elapsed time should be exactly zero"
      print *, "Got:", elapsed
      errors = errors + 1
    end if
    
    if (errors == 0) then
      print *, "  Elapsed time calculation test PASSED"
    end if
    
  end subroutine test_elapsed_time_calculation
  
  ! Helper subroutine to introduce a predictable delay
  subroutine cpu_intensive_work()
    integer :: i, j, dummy
    dummy = 0
    do i = 1, 1000
      do j = 1, 100
        dummy = dummy + i * j
      end do
    end do
    ! Prevent compiler optimization
    if (dummy < 0) print *, "Unexpected result"
  end subroutine cpu_intensive_work

end program test_timing