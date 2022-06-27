program neo_orb_main
  use params, only : read_config
  use simple_main, only : run, write_output

  call read_config
  call run
  call write_output

end program neo_orb_main
