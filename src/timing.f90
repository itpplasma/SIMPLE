module timing
  implicit none
  
  ! Define real(dp) kind parameter
  integer, parameter :: dp = kind(1.0d0)
  
  ! Timer variables (now storing clock counts)
  integer(kind=8) :: program_start_time
  integer(kind=8) :: phase_start_time
  integer(kind=8) :: clock_rate
  integer(kind=8) :: clock_max
  
contains

  subroutine init_timer()
    real(dp) :: resolution_ms, max_seconds
    
    call system_clock(program_start_time, clock_rate, clock_max)
    phase_start_time = program_start_time
    
    ! Calculate resolution in milliseconds
    resolution_ms = 1000.0_dp / real(clock_rate, dp)
    
    ! Calculate maximum time in seconds before wraparound
    max_seconds = real(clock_max, dp) / real(clock_rate, dp)
    
    ! Print timing system information
    write(*,'(A)') 'INFO: Timing system initialized'
    write(*,'(A,I12)') '  Clock rate (ticks per second): ', clock_rate
    write(*,'(A,F12.6,A)') '  Resolution per tick: ', resolution_ms, ' ms'
    write(*,'(A,ES12.3,A)') '  Maximum time before wraparound: ', max_seconds, ' seconds'
  end subroutine init_timer
  
  subroutine print_phase_time(phase_name)
    character(*), intent(in) :: phase_name
    integer(kind=8) :: current_time
    real(dp) :: elapsed_ms, total_ms
    
    call system_clock(current_time, clock_rate, clock_max)
    
    ! Convert clock counts to milliseconds
    elapsed_ms = real(current_time - phase_start_time, dp) / real(clock_rate, dp) * 1000.0_dp
    total_ms = real(current_time - program_start_time, dp) / real(clock_rate, dp) * 1000.0_dp
    
    write(*,'(A,A,F12.3,A,F12.3,A)') 'INFO: ', phase_name, elapsed_ms, &
      ' ms (Total: ', total_ms, ' ms)'
    
    phase_start_time = current_time
  end subroutine print_phase_time
  
  ! Get current wall-clock time in seconds (replacement for omp_get_wtime)
  function get_wtime() result(wtime)
    real(dp) :: wtime
    integer(kind=8) :: current_time
    integer(kind=8) :: rate, max
    
    call system_clock(current_time, rate, max)
    wtime = real(current_time, dp) / real(rate, dp)
  end function get_wtime
  
  ! Get elapsed time in seconds between two time points
  function get_elapsed_time(start_time, end_time) result(elapsed)
    real(dp), intent(in) :: start_time, end_time
    real(dp) :: elapsed
    
    elapsed = end_time - start_time
  end function get_elapsed_time

end module timing