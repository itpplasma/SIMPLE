program neo_orb_main

  use params, only : read_config, netcdffile, ns_s, ns_tp, multharm, &
    integmode, params_init
  use simple_main, only : init_field, run, write_output
  use simple, only : Tracer
  use timing, only : init_timer, print_phase_time

  implicit none

  character(256) :: config_file
  type(Tracer) :: norb

  ! Initialize timing
  call init_timer()

  ! read configuration file name from command line arguments
  if (command_argument_count() == 0) then
    config_file = 'simple.in'
  else
    call get_command_argument(1, config_file)
  end if
  call print_phase_time('Command line parsing completed')

  ! Must be called in this order. TODO: Fix
  call read_config(config_file)
  call print_phase_time('Configuration reading completed')
  
  call init_field(norb, netcdffile, ns_s, ns_tp, multharm, integmode)
  call print_phase_time('Field initialization completed')
  
  call params_init
  call print_phase_time('Parameter initialization completed')
  
  call run(norb)
  call print_phase_time('Main simulation run completed')
  
  call write_output
  call print_phase_time('Output writing completed')

end program neo_orb_main
