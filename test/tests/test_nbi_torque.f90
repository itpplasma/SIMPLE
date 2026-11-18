program test_nbi_torque
    !> Beam torque in an axisymmetric tokamak. Along every collisionless orbit
    !> step P = m v_par h_phi + (q/c) A_phi is conserved, so the torque left for
    !> a non-axisymmetric field must vanish, and the angular momentum given to
    !> the bulk, lost and still carried by fast ions adds up to the injected one.
    use, intrinsic :: iso_fortran_env, only: dp => real64
  use simple_api, only: simple_reset_config, simple_init, simple_sample, simple_trace, &
                          simple_write_output
    use params, only: field_input, coord_input, startmode, isw_field_type, integmode, &
                     deterministic, n_d, n_e, particle_energy_eV, nbi_file, ntestpart, &
                  swcoll, thermal_stop_factor, trace_time, ntimstep, am1, am2, Z1, Z2, &
                      contr_pp, relerr
    use simple_profiles, only: Te_scale, Ti1_scale, Ti2_scale, ni1_scale, ni2_scale
    use nbi_torque, only: nbi_tally_total
    use timing, only: init_timer
    implicit none
    real(dp) :: total(5)
    real(dp), allocatable :: z(:, :)
    character(1024) :: filename, mode
    integer :: unit

    call get_command_argument(1, filename)
    call get_command_argument(2, mode)
    call init_timer
    call simple_reset_config
    field_input = trim(filename)
    coord_input = trim(filename)
    isw_field_type = 1
    integmode = 0
    startmode = 8
    deterministic = .true.
    n_d = 2
    n_e = 1
    particle_energy_eV = 1.0e5_dp
    ntestpart = 32
    swcoll = .true.
    contr_pp = -1.0e10_dp
    thermal_stop_factor = 1.5_dp
    trace_time = 1.0e-3_dp
    if (trim(mode) == 'collisionless') swcoll = .false.
    ntimstep = 100
    am1 = 2.0_dp
    am2 = 3.0_dp
    Z1 = 1.0_dp
    Z2 = 1.0_dp
    Te_scale = 1.0e3_dp
    Ti1_scale = 1.0e3_dp
    Ti2_scale = 1.0e3_dp
    ni1_scale = 0.3e20_dp
    ni2_scale = 0.0_dp
    nbi_file = 'beamlets.dat'
    open (newunit=unit, file='beamlets.dat', status='replace')
    ! Tangential to R = 6 m in the midplane, co-current direction undetermined.
    write (unit, '(a)') 'T 2 1.0e5 1.0e6 1.0 0.0 0.0 6.0 -9.0 0.0 0.0 1.0 0.0'
    close (unit)
    call simple_init
    allocate (z(5, ntestpart))
    call simple_sample(ntestpart, z)
    call simple_trace
    call simple_write_output
    total = nbi_tally_total()
    print '(a,5es14.5)', 'bulk field lost fast injected [N m] ', total
    if (abs(sum(total(1:4))/total(5) - 1.0_dp) > 1.0e-10_dp) &
        error stop 'angular momentum budget does not close'
    if (abs(total(2)) > 1.0e-3_dp*abs(total(1))) &
        error stop 'field torque does not vanish in an axisymmetric field'
    print '(a)', 'PASS'
end program test_nbi_torque
