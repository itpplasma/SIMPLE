module timing
  implicit none
  
  ! Define real(dp) kind parameter
  integer, parameter :: dp = kind(1.0d0)
  
  ! Timer variables
  real(dp) :: program_start_time
  real(dp) :: phase_start_time
  
contains

  subroutine init_timer()
    call cpu_time(program_start_time)
    phase_start_time = program_start_time
  end subroutine init_timer
  
  subroutine print_phase_time(phase_name)
    character(*), intent(in) :: phase_name
    real(dp) :: current_time, elapsed_ms, total_ms
    
    call cpu_time(current_time)
    elapsed_ms = (current_time - phase_start_time) * 1000.0_dp
    total_ms = (current_time - program_start_time) * 1000.0_dp
    
    write(*,'(A,A,F12.3,A,F12.3,A)') 'INFO: ', phase_name, elapsed_ms, &
      ' ms (Total: ', total_ms, ' ms)'
    
    phase_start_time = current_time
  end subroutine print_phase_time

end module timing