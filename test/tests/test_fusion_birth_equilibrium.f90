program test_fusion_birth_equilibrium
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use simple_api, only: simple_reset_config, simple_init, simple_sample, &
        simple_trace, simple_write_output
    use params, only: field_input, coord_input, startmode, num_surf, &
        isw_field_type, integmode, nper, npoiper, ntimstep, deterministic, &
        fusion_grid, facE_al, sbeg, trace_time, orbit_exit_code
    use magfie_sub, only: magfie
    use field_can_mod, only: ref_to_integ
    use fusion_source, only: alpha_energy_joule, dt_reactivity, &
        build_birth_cdf, draw_birth
    use simple_profiles, only: ni1_coef
    use netcdf, only: nf90_open, nf90_nowrite, nf90_inq_varid, nf90_get_var, &
        nf90_inq_dimid, nf90_inquire_dimension, nf90_close, nf90_noerr, nf90_strerror
    use timing, only: init_timer
    !$ use omp_lib, only: omp_set_num_threads, omp_set_dynamic, omp_get_num_threads
    implicit none
    integer, parameter :: nb = 20, nr = 200, na = 48
    integer :: n = 100000
    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), allocatable :: z(:, :), cdf(:)
    real(dp) :: birth(5)
    real(dp) :: expected(nb), x(3), jac, bmod, bder(3), hcov(3), hctr(3), hcurl(3)
    real(dp) :: volume, power, reported, marker_power, chi2, radius, slo, shi
    real(dp) :: bin_volume(nb), no_jac(nb), source_factor, cell_volume
    real(dp) :: wout_volume(nb), wout_source(nb), volume_p, sensitivity
    logical :: chartmap
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
    chartmap = index(trim(mode), 'chartmap') == 1
    if (chartmap) isw_field_type = 2
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
    if (trim(mode) == 'vmec_peaked') then
        ! n_D = 5e19 (1-s), n_T = 5e19, T_D = T_T = 10 keV.
        ! Enough markers to resolve this wout's weak radial Jacobian variation.
        ni1_coef(2) = -1.0_dp
        n = 3000000
        fusion_grid(1) = 400
    end if
    ! Exercise API output on 100000 markers; stream additional draws into the
    ! sensitive histogram without allocating/writing millions of start rows.
    np = min(n, 100000)
    if (trim(mode) == 'trace') then
        np = 8
        trace_time = 1.0e-8_dp
    end if
    if (trim(mode) == 'bad_energy') facE_al = 2.0_dp
    call simple_init
    if (index(trim(mode), '_parallel') > 0) then
        call check_parallel
        stop
    end if
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
    if (n > np) call build_birth_cdf(slo, shi, fusion_grid, cdf, power)
    counts = 0
    do i = 1, n
        if (i > np) then
            call draw_birth(slo, shi, fusion_grid, cdf, birth)
            radius = birth(1)
        else if (chartmap) then
            call ref_to_integ(z(1:3, i), x)
            radius = x(1)
        else
            radius = z(1, i)
        end if
        if (radius < slo .or. radius >= shi) error stop 'equilibrium radial bounds'
        b = min(nb, int((radius - slo)/(shi - slo)*nb) + 1)
        counts(b) = counts(b) + 1
    end do

    ! Separate midpoint quadrature checks the sampler's discretization. For
    ! chartmaps it uses the SAME field-derived Boozer Jacobian, not geometry.
    expected = 0.0_dp
    bin_volume = 0.0_dp
    no_jac = 0.0_dp
    cell_volume = (shi - slo)*(2*pi)**2/(nr*na*na)*1.0e-6_dp
    do ir = 1, nr
        x(1) = slo + (shi - slo)*(real(ir, dp) - 0.5_dp)/nr
        source_factor = radial_source(x(1))
        b = (ir - 1)/(nr/nb) + 1
        no_jac(b) = no_jac(b) + source_factor
        do it = 1, na
            x(2) = 2*pi*(real(it, dp) - 0.5_dp)/na
            do iz = 1, na
                x(3) = 2*pi*(real(iz, dp) - 0.5_dp)/na
                call magfie(x, bmod, jac, bder, hcov, hctr, hcurl)
                bin_volume(b) = bin_volume(b) + abs(jac)*cell_volume
                expected(b) = expected(b) + source_factor*abs(jac)*cell_volume
            end do
        end do
    end do
    volume = sum(bin_volume)
    if (.not. chartmap) then
        call read_wout_reference(wout_volume, wout_source, volume_p)
        print *, 'VMEC volume: magfie, wout interval, volume_p [m^3]:', &
            volume, sum(wout_volume), volume_p
        print *, 'max relative bin dV/ds error:', &
            maxval(abs(bin_volume/wout_volume - 1))
        if (maxval(abs(bin_volume/wout_volume - 1)) > 5e-4_dp) &
            error stop 'VMEC per-bin volume disagrees with wout vp'
        if (slo == 0 .and. shi == 1) then
            if (abs(volume/volume_p - 1) > 5e-4_dp) &
                error stop 'VMEC total volume disagrees with wout volume_p'
        end if
        ! The VMEC distribution and power oracle use wout vp, not magfie.
        bin_volume = wout_volume
        expected = wout_source
    end if
    power = sum(expected)*(5.0e19_dp)**2* &
        dt_reactivity(1.0e4_dp)*alpha_energy_joule
    expected = expected/sum(expected)
    no_jac = no_jac/sum(no_jac)
    if (trim(mode) == 'vmec_peaked') then
        sensitivity = n*sum((no_jac - expected)**2/expected)
        print *, 'omitted-Jacobian expected chi2 excess:', sensitivity
        if (sensitivity <= 65) error stop 'insensitive equilibrium test'
    end if
    chi2 = sum((counts - n*expected)**2/(n*expected))
    print *, trim(mode), ' source dV/ds: chi2(19) = ', chi2, ' limit = 65'
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
    if (abs(np*marker_power/reported - 1) > 1e-12_dp) error stop 'marker weights'
    if (maxval(abs(z(4, :) - 1)) > 1e-14_dp) error stop 'birth energy'
    open (newunit=unit, file='equilibrium_radial.csv', status='replace')
    write (unit, '(a)') 's_left,s_right,count,reference_probability,'// &
        'without_jacobian_probability,reference_dV_ds_m3'
    do b = 1, nb
        write (unit, '(2(f8.4,","),i8,3(",",es24.16))') &
            slo + (shi - slo)*real(b - 1, dp)/nb, &
            slo + (shi - slo)*real(b, dp)/nb, counts(b), expected(b), &
            no_jac(b), bin_volume(b)*nb/(shi - slo)
    end do
    close (unit)
    print *, 'PASS equilibrium distribution, total alpha power, marker normalization'
contains

    function radial_source(s) result(f)
        real(dp), intent(in) :: s
        real(dp) :: f
        f = 1.0_dp
        if (trim(mode) == 'vmec_peaked') f = 1.0_dp - s
    end function radial_source

    subroutine nc_check(status)
        integer, intent(in) :: status
        if (status /= nf90_noerr) then
            print *, nf90_strerror(status)
            error stop 'reading independent wout volume reference'
        end if
    end subroutine nc_check

    subroutine read_wout_reference(volumes, sources, total_volume)
        real(dp), intent(out) :: volumes(nb), sources(nb), total_volume
        real(dp), allocatable :: vp(:)
        real(dp) :: s, h, t, dv, file_volume
        integer :: ncid, varid, dimid, ns, j, l, bin
        integer, parameter :: nq = 20000

        call nc_check(nf90_open(trim(filename), nf90_nowrite, ncid))
        call nc_check(nf90_inq_dimid(ncid, 'radius', dimid))
        call nc_check(nf90_inquire_dimension(ncid, dimid, len=ns))
        allocate (vp(ns))
        call nc_check(nf90_inq_varid(ncid, 'vp', varid))
        call nc_check(nf90_get_var(ncid, varid, vp))
        call nc_check(nf90_inq_varid(ncid, 'volume_p', varid))
        call nc_check(nf90_get_var(ncid, varid, total_volume))
        call nc_check(nf90_close(ncid))
        ! In this wout vp(1)=0 is padding, vp(2:ns) lives at half mesh
        ! s=(j-1.5)/(ns-1). vp is in m^3, with dV/ds=4*pi^2*vp.
        ! Verify the convention against the independently stored volume_p.
        h = 1.0_dp/(ns - 1)
        file_volume = 4*pi**2*h*sum(vp(2:ns))
        if (abs(file_volume/total_volume - 1) > 5e-6_dp) &
            error stop 'wout vp half-mesh/volume_p convention'
        volumes = 0.0_dp
        sources = 0.0_dp
        do j = 1, nq
            s = slo + (shi - slo)*(j - 0.5_dp)/nq
            ! Linear half-mesh interpolation, with endpoint extrapolation.
            t = s/h + 1.5_dp
            l = max(2, min(ns - 1, int(t)))
            t = t - l
            dv = 4*pi**2*((1 - t)*vp(l) + t*vp(l + 1))*(shi - slo)/nq
            bin = (j - 1)/(nq/nb) + 1
            volumes(bin) = volumes(bin) + dv
            sources(bin) = sources(bin) + radial_source(s)*dv
        end do
    end subroutine read_wout_reference

    subroutine check_parallel
        real(dp), allocatable :: serial_cdf(:), parallel_cdf(:)
        real(dp) :: serial_power, parallel_power, elapsed(3, 2)
        real(dp) :: draw(5), cdf_error
        integer :: trial, threads, column, clock_start, clock_end, clock_rate
        integer :: bins(nb, 2), seed_size, sample, bin, actual_threads
        integer, allocatable :: seed(:)
        integer, parameter :: grid(3) = [128, 64, 64]

        !$ call omp_set_dynamic(.false.)
        call random_seed(size=seed_size)
        allocate (seed(seed_size))
        seed = 543
        bins = 0
        do column = 1, 2
            threads = 1
            if (column == 2) threads = 8
            !$ call omp_set_num_threads(threads)
            actual_threads = 1
            !$omp parallel
            !$omp single
            !$ actual_threads = omp_get_num_threads()
            !$omp end single
            !$omp end parallel
            if (actual_threads /= threads) error stop 'OpenMP team size'
            do trial = 1, 3
                call system_clock(clock_start, clock_rate)
                if (column == 1) then
                    call build_birth_cdf(0.0_dp, 1.0_dp, grid, &
                        serial_cdf, serial_power)
                else
                    call build_birth_cdf(0.0_dp, 1.0_dp, grid, &
                        parallel_cdf, parallel_power)
                end if
                call system_clock(clock_end)
                elapsed(trial, column) = real(clock_end - clock_start, dp)/clock_rate
                if (column == 2) then
                    cdf_error = maxval(abs(parallel_cdf - serial_cdf))
                    if (cdf_error > 5e-14_dp) error stop 'parallel CDF changed'
                    if (abs(parallel_power/serial_power - 1) > 5e-14_dp) &
                        error stop 'parallel power changed'
                end if
            end do
            call random_seed(put=seed)
            do sample = 1, 100000
                if (column == 1) then
                    call draw_birth(0.0_dp, 1.0_dp, grid, serial_cdf, draw)
                else
                    call draw_birth(0.0_dp, 1.0_dp, grid, parallel_cdf, draw)
                end if
                bin = min(nb, int(draw(1)*nb) + 1)
                bins(bin, column) = bins(bin, column) + 1
            end do
            print *, trim(mode), ' CDF seconds, threads=', threads, elapsed(:, column)
        end do
        if (any(bins(:, 1) /= bins(:, 2))) error stop 'parallel histogram changed'
        print *, 'PASS 1/8 threads: CDF max error=', cdf_error, &
            ' histogram max count difference=', maxval(abs(bins(:, 1) - bins(:, 2)))
    end subroutine check_parallel
end program test_fusion_birth_equilibrium
