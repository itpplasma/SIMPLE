program test_thermal_ntv
    !> Thermal delta-f NTV tally in an axisymmetric tokamak. P = (q/c) psi* is
    !> conserved along every collisionless step, so the torque vanishes with and
    !> without collisions, and every weight stays zero: the canonical Maxwellian
    !> solves the unperturbed equation. A drive by ds would move the weights
    !> with the bounce displacement of the orbits. The bin selects markers by
    !> birth, so a bin narrower than the orbits still holds their whole trace.
    use, intrinsic :: iso_fortran_env, only: dp => real64
  use simple_api, only: simple_reset_config, simple_init, simple_sample, simple_trace, &
                          simple_write_output
    use params, only: field_input, coord_input, startmode, isw_field_type, integmode, &
                      deterministic, n_d, n_e, particle_energy_eV, ntestpart, swcoll, &
                     trace_time, ntimstep, am1, am2, Z1, Z2, contr_pp, num_surf, sbeg, &
                      ntv_tally, ntv_drive, ntv_temperature_eV, ntv_bin
    use simple_profiles, only: Te_scale, Ti1_scale, Ti2_scale, ni1_scale, ni2_scale
    use thermal_ntv, only: ntv_tally_weight, ntv_tally_marker
    use timing, only: init_timer
    implicit none
    real(dp), parameter :: a1 = 1.0_dp, a2 = 0.5_dp, s0 = 0.4_dp
    real(dp), allocatable :: z(:, :)
    real(dp) :: c_i, d_i, torque, dwell, wmax
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
    ntv_bin = [s0 - 1.0e-6_dp, s0 + 1.0e-6_dp]
    call simple_init
    allocate (z(5, ntestpart))
    call simple_sample(ntestpart, z)
    call simple_trace
    call simple_write_output

    torque = 0
    dwell = 0
    wmax = 0
    do i = 1, ntestpart
        call ntv_tally_marker(i, c_i, d_i)
        torque = max(torque, abs(c_i))
        dwell = dwell + d_i
        wmax = max(wmax, abs(ntv_tally_weight(i)))
    end do
    print '(a,3es14.5)', 'max |c_i| [kg m^2/s], dwell [s], max |w| ', &
        torque, dwell, wmax
    ! Every marker starts on s0 and stays confined for the whole trace.
    if (abs(dwell/(ntestpart*trace_time) - 1) > 1.0e-6_dp) &
        error stop 'tally misses marker time outside the birth bin'
    ! A thermal NTV of 1e-21 N m per particle over the 1e-4 s run gives
    ! |c_i| ~ 1e-25 kg m^2/s.
   if (torque > 1.0e-30_dp) error stop 'torque does not vanish in an axisymmetric field'
    ! A ds drive gives |w| ~ (A_1 + A_2 E/T) times the orbit width in s, 1e-3.
    if (wmax > 1.0e-10_dp) error stop 'weights move in an axisymmetric field'
    print '(a)', 'PASS'
end program test_thermal_ntv
