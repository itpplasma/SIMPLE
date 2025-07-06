program neo_orb_main

  use params, only : read_config, netcdffile, ns_s, ns_tp, multharm, &
    integmode, params_init
  use simple_main, only : init_field, run, write_output
  use simple, only : Tracer

  implicit none

  character(256) :: config_file
  type(Tracer) :: norb

  ! read configuration file name from command line arguments
  if (command_argument_count() == 0) then
    config_file = 'simple.in'
  else
    call get_command_argument(1, config_file)
  end if

  ! Must be called in this order. TODO: Fix
  call read_config(config_file)
  call init_field(norb, netcdffile, ns_s, ns_tp, multharm, integmode)
  call params_init
  call run(norb)
  call write_output

end program neo_orb_main
