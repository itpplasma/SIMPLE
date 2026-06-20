program simple_diag_cp_gc
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use params, only: read_config, netcdffile, ns_s, ns_tp, multharm, integmode, &
    params_init, dtaumin, ntimstep, ntestpart, zstart, startmode, grid_density, &
    special_ants_file, reuse_batch, num_surf, sbeg, v0, ntau_macro
  use simple, only: tracer_t, init_cp, orbit_timestep_cp_canonical
  use simple_main, only: init_field
  use magfie_sub, only: init_magfie, VMEC
  use samplers, only: init_starting_surf, sample, START_FILE
  use orbit_cpp_canonical, only: cpp_canon_boozer_guiding_center, cpp_canon_to_gc

  implicit none

  character(256) :: config_file, arg
  type(tracer_t) :: norb
  real(dp) :: z(5), xgc(3), t, r, th, ph, vpar
  integer(8) :: kt
  integer :: particle_number, ierr, unit, it, ktau

  config_file = 'simple.in'
  particle_number = 1
  select case (command_argument_count())
  case (0)
  case (1)
    call get_command_argument(1, arg)
    read(arg, *) particle_number
  case (2)
    call get_command_argument(1, config_file)
    call get_command_argument(2, arg)
    read(arg, *) particle_number
  case default
    print *, 'Usage: simple_diag_cp_gc [config_file] [particle_number]'
    stop 1
  end select

  call read_config(config_file)
  block
    use boozer_coordinates_mod, only: use_B_r, use_del_tp_B
    use_B_r = .true.
    use_del_tp_B = .true.
  end block
  call init_field(norb, netcdffile, ns_s, ns_tp, multharm, integmode)
  block
    use boozer_coordinates_mod, only: use_B_r, use_del_tp_B
    use boozer_sub, only: get_boozer_coordinates
    use_B_r = .true.
    use_del_tp_B = .true.
    call get_boozer_coordinates
  end block
  call params_init
  call init_magfie(VMEC)
  call init_starting_surf
  call sample_start_points

  if (particle_number < 1 .or. particle_number > ntestpart) then
    print *, 'particle out of range', particle_number, ntestpart
    error stop
  end if

  z = zstart(:, particle_number)
  call init_cp(norb%cp, norb%f, z, dtaumin)

  open(newunit=unit, file='cp_gc_trace.dat', status='replace')
  write(unit, '(A)') '# t_s s_full theta_full phi_full s_gc theta_gc phi_gc p_abs v_par'
  call write_state(unit, 0.0_dp, norb%cp)

  ierr = 0
  kt = 0
  do it = 2, ntimstep
    do ktau = 1, ntau_macro(it)
      call orbit_timestep_cp_canonical(norb%cp, norb%f, z, ierr)
      if (ierr /= 0) exit
      kt = kt + 1
    end do
    if (ierr /= 0) exit
    t = kt*dtaumin/v0
    call write_state(unit, t, norb%cp)
  end do
  close(unit)

  print '(A,A)', 'trace written: ', 'cp_gc_trace.dat'
  print '(A,I0)', 'rows written: ', it - 1

contains

  subroutine sample_start_points
    if (startmode == 1) then
      if ((0.0_dp < grid_density) .and. (1.0_dp > grid_density)) then
        call sample(zstart, grid_density)
      else
        call sample(zstart)
      end if
    else if (startmode == 2) then
      call sample(zstart, START_FILE)
    else if (startmode == 3) then
      call sample(special_ants_file)
    else if (startmode == 4) then
      call sample(zstart, reuse_batch)
    else if (startmode == 5) then
      if (num_surf == 1) then
        call sample(zstart, 0.0_dp, sbeg(1))
      else
        call sample(zstart, sbeg(1), sbeg(num_surf))
      end if
    else
      print *, 'Invalid startmode: ', startmode
      error stop
    end if
  end subroutine sample_start_points

  subroutine write_state(unit, t, st)
    use orbit_cpp_canonical, only: cpp_canon_state_t
    integer, intent(in) :: unit
    real(dp), intent(in) :: t
    type(cpp_canon_state_t), intent(in) :: st

    call cpp_canon_to_gc(st, r, th, ph, vpar)
    call cpp_canon_boozer_guiding_center(st, xgc)
    write(unit, '(9ES24.16)') t, r, th, ph, xgc(1), xgc(2), xgc(3), st%pabs, vpar
  end subroutine write_state
end program simple_diag_cp_gc
