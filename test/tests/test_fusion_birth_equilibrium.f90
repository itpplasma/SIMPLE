program test_fusion_birth_equilibrium
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use simple_api, only: simple_reset_config, simple_init, simple_sample, &
        simple_trace, simple_write_output
    use params, only: field_input, coord_input, startmode, num_surf, &
        isw_field_type, integmode, nper, npoiper, ntimstep, deterministic, &
        fusion_grid, facE_al, sbeg, trace_time, orbit_exit_code
    use magfie_sub, only: magfie
    use field_can_mod, only: ref_to_integ
    use fusion_source, only: alpha_energy_joule, dt_reactivity
    use timing, only: init_timer
    implicit none
    integer, parameter :: n = 100000, nb = 20, nr = 200, na = 48
    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), allocatable :: z(:, :)
    real(dp) :: expected(nb), x(3), jac, bmod, bder(3), hcov(3), hctr(3), hcurl(3)
    real(dp) :: volume, power, reported, marker_power, chi2, radius, slo, shi
    integer :: ir, it, iz, b, i, counts(nb), unit, ios, np
    character(1024) :: filename, mode, line, key

    call get_command_argument(1, filename)
    call get_command_argument(2, mode)
    call init_timer
    call simple_reset_config
    field_input = trim(filename)
    coord_input = trim(filename)
    isw_field_type = 1
    integmode = 0
    if (trim(mode) == 'chartmap') isw_field_type = 2
    startmode = 7
    num_surf = 0
    deterministic = .true.
    nper = 1
    npoiper = 4
    ntimstep = 2
    fusion_grid = [128, 32, 32]
    slo = 0.0_dp
    shi = 1.0_dp
    if (trim(mode) == 'annulus') then
        num_surf = 2
        slo = 0.2_dp
        shi = 0.8_dp
        sbeg(1:2) = [slo, shi]
    end if
    np = n
    if (trim(mode) == 'trace') then
        np = 8
        trace_time = 1.0e-8_dp
    end if
    if (trim(mode) == 'bad_energy') facE_al = 2.0_dp
    call simple_init
    allocate (z(5, np))
    call simple_sample(np, z)
    if (trim(mode) == 'bad_energy') stop
    if (trim(mode) == 'trace') then
        call simple_trace
        if (any(orbit_exit_code >= 100)) error stop 'fusion trace numerical failure'
        call simple_write_output
        print *, 'PASS eight fusion markers traced and results written'
        stop
    end if
    counts = 0
    do i = 1, n
        if (trim(mode) == 'chartmap') then
            call ref_to_integ(z(1:3, i), x)
            radius = x(1)
        else
            radius = z(1, i)
        end if
        if (radius < slo .or. radius >= shi) error stop 'equilibrium radial bounds'
        b = min(nb, int((radius - slo)/(shi - slo)*nb) + 1)
        counts(b) = counts(b) + 1
    end do

    ! Independent full-torus volume quadrature, with different nodes/resolution.
    ! This uses the physical equilibrium, not the source CDF or its cell masses.
    expected = 0.0_dp
    do ir = 1, nr
        x(1) = slo + (shi - slo)*(real(ir, dp) - 0.5_dp)/nr
        b = (ir - 1)/(nr/nb) + 1
        do it = 1, na
            x(2) = 2*pi*(real(it, dp) - 0.5_dp)/na
            do iz = 1, na
                x(3) = 2*pi*(real(iz, dp) - 0.5_dp)/na
                call magfie(x, bmod, jac, bder, hcov, hctr, hcurl)
                expected(b) = expected(b) + abs(jac)
            end do
        end do
    end do
    volume = sum(expected)*(shi - slo)*(2*pi)**2/(nr*na*na)*1.0e-6_dp
    power = volume*(5.0e19_dp)**2*dt_reactivity(1.0e4_dp)*alpha_energy_joule
    expected = expected/sum(expected)
    chi2 = sum((counts - n*expected)**2/(n*expected))
    print *, trim(mode), ' flat dV/ds: chi2(19) = ', chi2, ' limit = 65'
    if (chi2 > 65) error stop 'equilibrium radial chi-squared'
    reported = -1.0_dp
    marker_power = -1.0_dp
    open (newunit=unit, file='alpha_source.dat', status='old')
    do
        read (unit, '(a)', iostat=ios) line
        if (ios /= 0) exit
        if (line(1:1) == '#') cycle
        read (line, *) key
        if (trim(key) == 'alpha_power_W') read (line, *) key, reported
        if (trim(key) == 'marker_power_W') read (line, *) key, marker_power
    end do
    close (unit)
    print *, 'power: reported, independent [W], relative error:', &
        reported, power, abs(reported/power - 1)
    if (abs(reported/power - 1) > 0.002_dp) error stop 'equilibrium alpha power'
    if (abs(n*marker_power/reported - 1) > 1e-12_dp) error stop 'marker weights'
    if (maxval(abs(z(4, :) - 1)) > 1e-14_dp) error stop 'birth energy'
    open (newunit=unit, file='equilibrium_radial.csv', status='replace')
    write (unit, '(a)') 's_left,s_right,count,reference_probability'
    do b = 1, nb
        write (unit, '(2(f8.4,","),i8,",",es24.16)') &
            slo + (shi - slo)*real(b - 1, dp)/nb, &
            slo + (shi - slo)*real(b, dp)/nb, counts(b), expected(b)
    end do
    close (unit)
    print *, 'PASS equilibrium distribution, total alpha power, marker normalization'
end program test_fusion_birth_equilibrium
