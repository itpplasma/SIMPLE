program test_nbi_torque
    !> Beam torque in an axisymmetric tokamak. Along every collisionless orbit
    !> step P = m v_par h_phi + (q/c) A_phi is conserved, so the torque left for
    !> a non-axisymmetric field must vanish, and the angular momentum given to
    !> the bulk, lost and still carried by fast ions adds up to the injected one.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use simple_api, only: simple_reset_config, simple_init, simple_sample, &
                          simple_set_start, simple_trace, simple_write_output
    use params, only: field_input, coord_input, startmode, isw_field_type, integmode, &
                     deterministic, n_d, n_e, particle_energy_eV, nbi_file, ntestpart, &
                  swcoll, thermal_stop_factor, trace_time, ntimstep, am1, am2, Z1, Z2, &
                      contr_pp, nbi_path_length, nbi_torque_bins, npoiper2, &
                      ntau_macro, orbit_exit_code, &
                      ORBIT_EXIT_LCFS, ORBIT_EXIT_COMPLETED, v0
    use simple_profiles, only: Te_scale, Ti1_scale, Ti2_scale, ni1_scale, ni2_scale
    use nbi_torque, only: nbi_tally_total
    use timing, only: init_timer
    implicit none
    real(dp) :: total(6)
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
    if (trim(mode) == 'fractional_loss') then
        swcoll = .false.
        trace_time = 5.0e-7_dp
        ntimstep = 2
        npoiper2 = 16
        nbi_torque_bins = 1
    end if
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
    ! From 9 m before tangency through the far side of the plasma.
    nbi_path_length = 20.0_dp
    open (newunit=unit, file='beamlets.dat', status='replace')
    ! Tangential to R = 6 m in the midplane, co-current direction undetermined.
    write (unit, '(a)') 'T 2 1.0e5 1.0e6 1.0 0.0 0.0 6.0 -9.0 0.0 0.0 1.0 0.0'
    close (unit)
    call simple_init
    allocate (z(5, ntestpart))
    call simple_sample(ntestpart, z)
    if (trim(mode) == 'fractional_loss') call start_near_lcfs(z)
    call simple_trace
    call simple_write_output
    total = nbi_tally_total()
    print '(a,5es14.5)', 'bulk field lost fast injected [N m] ', total(1:5)
    print '(a,es14.5)', 'beam-target neutrons [1/s] ', total(6)
    if (abs(sum(total(1:4))/total(5) - 1.0_dp) > 1.0e-10_dp) &
        error stop 'angular momentum budget does not close'
    if (abs(total(2)) > 1.0e-3_dp*abs(total(1))) &
        error stop 'field torque does not vanish in an axisymmetric field'
    if (trim(mode) == 'fractional_loss') call check_fractional_loss(total)
    print '(a)', 'PASS'

contains

    subroutine start_near_lcfs(start)
        real(dp), intent(out) :: start(:, :)
        real(dp), parameter :: twopi = 2.0_dp*acos(-1.0_dp)
        integer :: i

        ! Poloidal coverage gives both outward and inward drifts.
        do i = 1, size(start, 2)
            start(:, i) = [0.99999_dp, (i - 0.5_dp)*twopi/size(start, 2), &
                           0.2_dp, 1.0_dp, 0.7_dp]
        end do
        call simple_set_start(size(start, 2), start)
    end subroutine start_near_lcfs

    subroutine check_fractional_loss(tally)
        use field_can_mod, only: ref_to_integ
        use fusion_source, only: dd_beam_maxwellian_rate
        use nbi_birth, only: nbi_marker_rate
        use nbi_torque, only: nbi_gc_momentum
        use simple_profiles, only: get_plasma_params
        use util, only: p_mass
        real(dp), intent(in) :: tally(6)
        real(dp) :: expected_lost, expected_neutron, l, a, s, z_end(5), z_integ(5)
        real(dp) :: loss_time, trapped, birth_s, perp, te, ti1, ti2, n1, n2, rate
        integer :: row, i, n_lost, unit

        if (ntau_macro(2) /= 1) error stop 'fixture must have one microstep'
        expected_lost = 0.0_dp
        expected_neutron = 0.0_dp
        n_lost = 0
        open (newunit=unit, file='times_lost.dat', status='old')
        do i = 1, ntestpart
            read (unit, *) row, loss_time, trapped, birth_s, perp, z_end
            if (row /= i) error stop 'unexpected endpoint marker'
            call ref_to_integ(z_end(1:3), z_integ(1:3))
            z_integ(4:5) = z_end(4:5)
            call nbi_gc_momentum(z_integ, l, a, s)
            select case (orbit_exit_code(i))
            case (ORBIT_EXIT_LCFS)
                if (.not. (0.0_dp < loss_time .and. loss_time < trace_time)) &
                    error stop 'loss must occur within the first microstep'
                expected_lost = expected_lost + l
                n_lost = n_lost + 1
            case (ORBIT_EXIT_COMPLETED)
                if (abs(loss_time - trace_time) > 1.0e-12_dp*trace_time) &
                    error stop 'still-fast marker did not complete the microstep'
            case default
                error stop 'unexpected exit in fractional-loss fixture'
            end select
            call get_plasma_params(s, te, ti1, ti2, n1, n2)
            rate = dd_beam_maxwellian_rate(1.0e-2_dp*v0*z_end(4), ti1, 2.0e-3_dp*p_mass)
            expected_neutron = expected_neutron + loss_time*n1*rate
        end do
        close (unit)
        if (n_lost == 0 .or. n_lost == ntestpart) &
            error stop 'fixture must have lost and still-fast markers'
        expected_lost = nbi_marker_rate*expected_lost
        expected_neutron = nbi_marker_rate*expected_neutron
        print '(a,2es22.12)', 'lost torque, tally and endpoints [N m] ', &
            tally(3), expected_lost
        print '(a,2es22.12)', 'neutron rate, tally and endpoints [1/s] ', &
            tally(6), expected_neutron
        if (abs(tally(3) - expected_lost) > 1.0e-12_dp*abs(expected_lost)) &
            error stop 'lost rate does not use saved endpoints'
        if (abs(tally(6) - expected_neutron) > 1.0e-12_dp*expected_neutron) &
            error stop 'neutron rate does not use accepted durations'
    end subroutine check_fractional_loss
end program test_nbi_torque
