program test_geoflux_rk45_orbit_plot

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use params, only: read_config, params_init, netcdffile, ns_s, ns_tp, multharm, &
        integmode, isw_field_type, dtaumin, relerr, ntimstep, v0, zstart
    use simple, only: tracer_t
    use simple_main, only: init_field
    use magfie_sub, only: init_magfie
    use alpha_lifetime_sub, only: orbit_timestep_axis
    use samplers, only: load_starting_points
    use geoflux_coordinates, only: geoflux_to_cyl
    use fortplot, only: figure, plot, savefig, xlabel, ylabel, title
    use util, only: twopi

    implicit none

    type(tracer_t) :: norb
    character(len=1024) :: out_root, out_dir, config_file, start_file, geqdsk_file
    character(len=1024) :: pdf_rz, pdf_s, cmd
    integer :: status, ierr, i, nmax, n_used, mkdir_stat
    real(dp), allocatable :: s_traj(:), theta_traj(:), phi_traj(:), time_traj(:)
    real(dp), allocatable :: r_traj(:), z_traj(:)
    real(dp) :: z(5), xcyl(3)
    logical :: exists

    out_root = ''
    call get_environment_variable('SIMPLE_VISUAL_DIR', value=out_root, status=status)
    if (status /= 0 .or. len_trim(out_root) == 0) then
        out_root = '/tmp/SIMPLE_visual_artifacts'
    end if

    out_dir = trim(out_root)//'/geoflux_rk45_orbit'
    cmd = 'mkdir -p '//trim(out_dir)
    call execute_command_line(trim(cmd), exitstat=mkdir_stat)
    if (mkdir_stat /= 0) then
        error stop 'test_geoflux_rk45_orbit_plot: failed to create output directory'
    end if

    geqdsk_file = 'EQDSK_I.geqdsk'
    call get_environment_variable('LIBNEO_TEST_GEQDSK', value=geqdsk_file, status=status)
    if (status /= 0 .or. len_trim(geqdsk_file) == 0) then
        geqdsk_file = 'EQDSK_I.geqdsk'
    end if

    config_file = trim(out_dir)//'/simple_geoflux_rk45_plot.in'
    start_file = trim(out_dir)//'/start.dat'

    call write_config(trim(config_file), trim(geqdsk_file))
    call write_start(trim(start_file))

    call read_config(trim(config_file))
    call init_field(norb, netcdffile, ns_s, ns_tp, multharm, integmode)
    call params_init
    call init_magfie(isw_field_type)

    call load_starting_points(zstart, trim(start_file))
    z = zstart(:, 1)

    nmax = ntimstep
    allocate(s_traj(nmax), theta_traj(nmax), phi_traj(nmax), time_traj(nmax))
    allocate(r_traj(nmax), z_traj(nmax))

    ierr = 0
    n_used = 0
    do i = 1, nmax
        s_traj(i) = z(1)
        theta_traj(i) = z(2)
        phi_traj(i) = z(3)
        time_traj(i) = real(i - 1, dp) * dtaumin / max(v0, 1.0d-12)

        call geoflux_to_cyl((/ z(1), z(2), z(3) /), xcyl)
        r_traj(i) = xcyl(1)
        z_traj(i) = xcyl(3)

        n_used = i
        call orbit_timestep_axis(z, dtaumin, dtaumin, relerr, ierr)
        if (ierr /= 0) exit
        if (z(1) < 0.0_dp .or. z(1) > 1.0_dp) exit
    end do

    pdf_rz = trim(out_dir)//'/geoflux_rk45_orbit_RZ.pdf'
    pdf_s = trim(out_dir)//'/geoflux_rk45_orbit_s.pdf'

    call figure()
    call plot(r_traj(1:n_used), z_traj(1:n_used), linestyle='b-')
    call xlabel('R (cm)')
    call ylabel('Z (cm)')
    call title('GEQDSK geoflux RK45 orbit projection (R,Z)')
    call savefig(trim(pdf_rz))

    call figure()
    call plot(time_traj(1:n_used), s_traj(1:n_used), linestyle='k-')
    call xlabel('t (s) [scaled]')
    call ylabel('s')
    call title('GEQDSK geoflux RK45 orbit: s(t)')
    call savefig(trim(pdf_s))

    inquire(file=trim(pdf_rz), exist=exists)
    if (.not. exists) then
        error stop 'test_geoflux_rk45_orbit_plot: missing RZ plot output'
    end if
    inquire(file=trim(pdf_s), exist=exists)
    if (.not. exists) then
        error stop 'test_geoflux_rk45_orbit_plot: missing s plot output'
    end if

    print *, 'VISUAL_ARTIFACT: ', trim(pdf_rz)
    print *, 'VISUAL_ARTIFACT: ', trim(pdf_s)

contains

    subroutine write_config(path, geqdsk_path)
        character(len=*), intent(in) :: path
        character(len=*), intent(in) :: geqdsk_path
        integer :: unit

        open(newunit=unit, file=trim(path), status='replace', action='write')
        write(unit, '(A)') '&config'
        write(unit, '(A)') 'netcdffile = '''//trim(geqdsk_path)//''''
        write(unit, '(A)') 'integmode = 0'
        write(unit, '(A)') 'isw_field_type = 1'
        write(unit, '(A)') 'ntestpart = 1'
        write(unit, '(A)') 'ntimstep = 400'
        write(unit, '(A)') 'npoiper2 = 64'
        write(unit, '(A)') 'multharm = 3'
        write(unit, '(A)') 'trace_time = 1d-10'
        write(unit, '(A)') 'relerr = 1d-10'
        write(unit, '(A)') 'deterministic = .True.'
        write(unit, '(A)') '/'
        close(unit)
    end subroutine write_config

    subroutine write_start(path)
        character(len=*), intent(in) :: path
        integer :: unit
        real(dp) :: s0, th0, ph0, p0, lam0

        s0 = 0.25_dp
        th0 = 0.1_dp * twopi
        ph0 = 0.0_dp
        p0 = 1.0_dp
        lam0 = 0.2_dp

        open(newunit=unit, file=trim(path), status='replace', action='write')
        write(unit, *) s0, th0, ph0, p0, lam0
        close(unit)
    end subroutine write_start

end program test_geoflux_rk45_orbit_plot

