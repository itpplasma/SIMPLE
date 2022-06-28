program neo_orb_main

  use params, only : read_config, netcdffile, ns_s, ns_tp, multharm, &
    integmode, params_init
  use simple_main, only : run, write_output
  use simple, only : Tracer, init_field

  implicit none

  type(Tracer) :: norb

  call read_config
  call init_field(norb, netcdffile, ns_s, ns_tp, multharm, integmode)
  call params_init
  call run(norb)
  call write_output

end program neo_orb_main
