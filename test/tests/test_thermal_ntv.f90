program test_thermal_ntv
    !> Thermal delta-f NTV tally in an axisymmetric tokamak. P = (q/c) psi* is
    !> conserved along every collisionless step, so the torque vanishes with and
    !> without collisions, and without collisions the weight is
    !> w = -(A_1 + A_2 E/T) (s_end - s_start) at constant energy.
    use, intrinsic :: iso_fortran_env, only: dp => real64
  use simple_api, only: simple_reset_config, simple_init, simple_sample, simple_trace, &
                          simple_write_output
    use params, only: field_input, coord_input, startmode, isw_field_type, integmode, &
                      deterministic, n_d, n_e, particle_energy_eV, ntestpart, swcoll, &
                     trace_time, ntimstep, am1, am2, Z1, Z2, contr_pp, num_surf, sbeg, &
                      ntv_tally, ntv_drive, ntv_temperature_eV, ntv_bin, zend, &
                      particle_energy_eV_effective
    use simple_profiles, only: Te_scale, Ti1_scale, Ti2_scale, ni1_scale, ni2_scale
    use thermal_ntv, only: ntv_tally_weight, ntv_tally_marker
    use timing, only: init_timer
    implicit none
    real(dp), parameter :: a1 = 1.0_dp, a2 = 0.5_dp, s0 = 0.4_dp
    real(dp), allocatable :: z(:, :)
    real(dp) :: c_i, d_i, s_end, x2, err, torque, dwell, wmax
    character(1024) :: filename, mode
    integer :: i

    call get_command_argument(1, filename)
    call get_command_argument(2, mode)
    call init_timer
    call simple_reset_config
    field_input = trim(filename)
    coord_input = trim(filename)
    isw_field_type = 2
    integmode = 1
    startmode = 1
    num_surf = 1
    sbeg(1) = s0
    deterministic = .true.
    n_d = 2
    n_e = 1
    particle_energy_eV = 1.0e3_dp
    ntestpart = 32
    contr_pp = -1.0e10_dp
    trace_time = 1.0e-4_dp
    ntimstep = 100
    swcoll = trim(mode) == 'collisional'
    am1 = 2.0_dp
    am2 = 3.0_dp
    Z1 = 1.0_dp
    Z2 = 1.0_dp
    Te_scale = 1.0e3_dp
    Ti1_scale = 1.0e3_dp
    Ti2_scale = 1.0e3_dp
    ni1_scale = 1.0e20_dp
    ni2_scale = 0.0_dp
    ntv_tally = .true.
    ntv_drive = [a1, a2]
    ntv_temperature_eV = 1.0e3_dp
    ntv_bin = [0.0_dp, 1.0_dp]
    call simple_init
    allocate (z(5, ntestpart))
    call simple_sample(ntestpart, z)
    call simple_trace
    call simple_write_output

    torque = 0
    dwell = 0
    err = 0
    wmax = 0
    do i = 1, ntestpart
        call ntv_tally_marker(i, c_i, d_i)
        torque = max(torque, abs(c_i))
        dwell = dwell + d_i
        if (.not. swcoll) then
            s_end = zend(1, i)
            x2 = particle_energy_eV_effective/ntv_temperature_eV*zend(4, i)**2
            err = max(err, abs(ntv_tally_weight(i) + (a1 + a2*x2)*(s_end - z(1, i))))
            wmax = max(wmax, abs(ntv_tally_weight(i)))
        end if
    end do
    print '(a,3es14.5)', 'max |c_i| [kg m^2/s], dwell [s], weight error ', &
        torque, dwell, err
    if (dwell <= 0) error stop 'no marker time in the tally bin'
    ! A thermal NTV of 1e-21 N m per particle over the 1e-4 s run gives
    ! |c_i| ~ 1e-25 kg m^2/s.
   if (torque > 1.0e-30_dp) error stop 'torque does not vanish in an axisymmetric field'
    ! The symplectic energy error (|dp/p| < 1e-3) enters through E/T.
    if (.not. swcoll .and. err > 1.0e-2_dp*wmax) &
        error stop 'weight differs from -(A_1 + A_2 E/T) ds'
    print '(a)', 'PASS'
end program test_thermal_ntv
