program test_geoflux_rk45_orbit_plot_chartmap_compare

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use params, only: read_config, params_init, netcdffile, ns_s, ns_tp, multharm, &
        integmode, isw_field_type, dtaumin, relerr, ntimstep, v0, zstart
    use simple, only: tracer_t
    use simple_main, only: init_field
    use magfie_sub, only: init_magfie, VMEC, magfie
    use alpha_lifetime_sub, only: orbit_timestep_axis
    use samplers, only: load_starting_points
    use reference_coordinates, only: ref_coords
    use pyplot_module, only: pyplot
    use util, only: twopi
    use geoflux_coordinates, only: geoflux_to_cart
    use netcdf

    implicit none

    type(tracer_t) :: norb
    type(pyplot) :: plt

    character(len=1024) :: out_root, out_dir
    character(len=1024) :: out_geoflux, out_chartmap, out_compare
    character(len=1024) :: cfg_geoflux, cfg_chartmap
    character(len=1024) :: start_pass, start_trap
    character(len=1024) :: geqdsk_file, chartmap_file, chartmap_env
    integer :: status

    integer, parameter :: norbits = 2
    integer, parameter :: ntheta_plot = 361
    integer, parameter :: nsurf_plot = 6
    integer :: nmax, i, iorb
    integer :: n_used_geo(norbits), n_used_cm(norbits)
    real(dp) :: color_geo(3), color_cm(3)

    real(dp), allocatable :: time_traj(:)

    real(dp), allocatable :: s_geo(:, :), th_geo(:, :), ph_geo(:, :)
    real(dp), allocatable :: R_geo(:, :), Z_geo(:, :), B_geo(:, :)

    real(dp), allocatable :: s_cm(:, :), th_cm(:, :), ph_cm(:, :)
    real(dp), allocatable :: R_cm(:, :), Z_cm(:, :), B_cm(:, :)

    real(dp) :: theta_grid(ntheta_plot), surf_s(nsurf_plot)
    real(dp) :: R_surf_geo(ntheta_plot, nsurf_plot), Z_surf_geo(ntheta_plot, nsurf_plot)
    real(dp) :: R_surf_cm(ntheta_plot, nsurf_plot), Z_surf_cm(ntheta_plot, nsurf_plot)

    out_root = ''
    call get_environment_variable('SIMPLE_ARTIFACT_DIR', value=out_root, status=status)
    if (status /= 0 .or. len_trim(out_root) == 0) then
        out_root = '/tmp/SIMPLE_artifacts'
    end if

    out_dir = trim(out_root)//'/plot/test_geoflux_rk45_orbit_plot_chartmap_compare'
    out_geoflux = trim(out_dir)//'/geoflux'
    out_chartmap = trim(out_dir)//'/chartmap'
    out_compare = trim(out_dir)//'/compare'

    call mkdir_p(trim(out_geoflux)//'/orbit')
    call mkdir_p(trim(out_geoflux)//'/flux_surfaces')
    call mkdir_p(trim(out_chartmap)//'/orbit')
    call mkdir_p(trim(out_chartmap)//'/flux_surfaces')
    call mkdir_p(trim(out_compare))

    geqdsk_file = 'EQDSK_I.geqdsk'
    call get_environment_variable('LIBNEO_TEST_GEQDSK', value=geqdsk_file, status=status)
    if (status /= 0 .or. len_trim(geqdsk_file) == 0) geqdsk_file = 'EQDSK_I.geqdsk'

    chartmap_file = trim(out_dir)//'/chartmap_from_geoflux.nc'
    call get_environment_variable('SIMPLE_CHARTMAP_FILE', value=chartmap_env, status=status)
    if (status == 0 .and. len_trim(chartmap_env) > 0) chartmap_file = trim(chartmap_env)

    cfg_geoflux = trim(out_dir)//'/simple_geoflux_ref.in'
    cfg_chartmap = trim(out_dir)//'/simple_chartmap_ref.in'
    start_pass = trim(out_dir)//'/start_passing.dat'
    start_trap = trim(out_dir)//'/start_trapped.dat'

    call write_start(start_pass, 0.25_dp, 0.1_dp*twopi, 0.0_dp, 1.0_dp, 0.7_dp)
    call write_start(start_trap, 0.25_dp, 0.25_dp*twopi, 0.0_dp, 1.0_dp, 0.0_dp)

    call write_config_geoflux(cfg_geoflux, geqdsk_file)
    call write_config_chartmap(cfg_chartmap, geqdsk_file, chartmap_file)

    call load_case(cfg_geoflux)
    if (.not. (status == 0 .and. len_trim(chartmap_env) > 0)) then
        call write_chartmap_from_geoflux(trim(chartmap_file), 65, 129, 65)
    end if

    nmax = ntimstep
    allocate(time_traj(nmax))
    do i = 1, nmax
        time_traj(i) = real(i - 1, dp) * dtaumin / max(v0, 1.0d-12)
    end do

    allocate(s_geo(nmax, norbits), th_geo(nmax, norbits), ph_geo(nmax, norbits))
    allocate(R_geo(nmax, norbits), Z_geo(nmax, norbits), B_geo(nmax, norbits))
    allocate(s_cm(nmax, norbits), th_cm(nmax, norbits), ph_cm(nmax, norbits))
    allocate(R_cm(nmax, norbits), Z_cm(nmax, norbits), B_cm(nmax, norbits))

    call integrate_orbit_from_start(start_pass, 1, n_used_geo, s_geo, th_geo, ph_geo, R_geo, Z_geo, B_geo)
    call integrate_orbit_from_start(start_trap, 2, n_used_geo, s_geo, th_geo, ph_geo, R_geo, Z_geo, B_geo)

    color_geo = [0.0_dp, 0.0_dp, 1.0_dp]
    color_cm = [1.0_dp, 0.0_dp, 0.0_dp]

    call build_theta_surf_grids(theta_grid, surf_s)
    call build_flux_surfaces(theta_grid, surf_s, R_surf_geo, Z_surf_geo)
    call plot_case(out_geoflux, 'geoflux reference', time_traj, n_used_geo, R_geo, Z_geo, s_geo, th_geo, ph_geo, B_geo, &
        R_surf_geo, Z_surf_geo, theta_grid, surf_s)

    call load_case(cfg_chartmap)
    call integrate_orbit_from_start(start_pass, 1, n_used_cm, s_cm, th_cm, ph_cm, R_cm, Z_cm, B_cm)
    call integrate_orbit_from_start(start_trap, 2, n_used_cm, s_cm, th_cm, ph_cm, R_cm, Z_cm, B_cm)
    call build_flux_surfaces(theta_grid, surf_s, R_surf_cm, Z_surf_cm)
    call plot_case(out_chartmap, 'chartmap reference', time_traj, n_used_cm, R_cm, Z_cm, s_cm, th_cm, ph_cm, B_cm, &
        R_surf_cm, Z_surf_cm, theta_grid, surf_s)

    do iorb = 1, norbits
        if (min(n_used_geo(iorb), n_used_cm(iorb)) < 50) then
            error stop 'test_geoflux_rk45_orbit_plot_chartmap_compare: orbit produced too few points'
        end if
    end do

    call plot_compare(out_compare, time_traj, n_used_geo, n_used_cm, R_geo, Z_geo, R_cm, Z_cm, color_geo, color_cm)

    print *, 'ARTIFACT_DIR: ', trim(out_dir)
    print *, 'ARTIFACT: ', trim(out_compare)//'/orbit_RZ_overlay.png'
    print *, 'ARTIFACT: ', trim(out_compare)//'/orbit_deltaRZ_t.png'
    print *, 'ARTIFACT: ', trim(out_compare)//'/flux_surfaces_RZ_phi0_overlay.png'

contains

    subroutine mkdir_p(path)
        character(len=*), intent(in) :: path
        character(len=2048) :: cmd
        integer :: stat

        cmd = 'mkdir -p '//trim(path)
        call execute_command_line(trim(cmd), exitstat=stat)
        if (stat /= 0) error stop 'mkdir_p failed'
    end subroutine mkdir_p

    logical function file_exists(path)
        character(len=*), intent(in) :: path
        inquire(file=trim(path), exist=file_exists)
    end function file_exists

    subroutine write_config_geoflux(path, geqdsk_path)
        character(len=*), intent(in) :: path
        character(len=*), intent(in) :: geqdsk_path
        integer :: unit

        open(newunit=unit, file=trim(path), status='replace', action='write')
        write(unit, '(A)') '&config'
        write(unit, '(A)') 'netcdffile = '''//trim(geqdsk_path)//''''
        write(unit, '(A)') 'field_input = '''//trim(geqdsk_path)//''''
        write(unit, '(A)') 'coord_input = '''//trim(geqdsk_path)//''''
        write(unit, '(A)') 'integmode = 0'
        write(unit, '(A)') 'isw_field_type = 1'
        write(unit, '(A)') 'ntestpart = 1'
        write(unit, '(A)') 'ntimstep = 800'
        write(unit, '(A)') 'npoiper2 = 64'
        write(unit, '(A)') 'multharm = 3'
        write(unit, '(A)') 'trace_time = 1d-6'
        write(unit, '(A)') 'relerr = 1d-11'
        write(unit, '(A)') 'deterministic = .True.'
        write(unit, '(A)') '/'
        close(unit)
    end subroutine write_config_geoflux

    subroutine write_config_chartmap(path, geqdsk_path, chartmap_path)
        character(len=*), intent(in) :: path
        character(len=*), intent(in) :: geqdsk_path
        character(len=*), intent(in) :: chartmap_path
        integer :: unit

        open(newunit=unit, file=trim(path), status='replace', action='write')
        write(unit, '(A)') '&config'
        write(unit, '(A)') 'netcdffile = '''//trim(geqdsk_path)//''''
        write(unit, '(A)') 'field_input = '''//trim(geqdsk_path)//''''
        write(unit, '(A)') 'coord_input = '''//trim(chartmap_path)//''''
        write(unit, '(A)') 'integmode = 0'
        write(unit, '(A)') 'isw_field_type = 1'
        write(unit, '(A)') 'ntestpart = 1'
        write(unit, '(A)') 'ntimstep = 800'
        write(unit, '(A)') 'npoiper2 = 64'
        write(unit, '(A)') 'multharm = 3'
        write(unit, '(A)') 'trace_time = 1d-6'
        write(unit, '(A)') 'relerr = 1d-11'
        write(unit, '(A)') 'deterministic = .True.'
        write(unit, '(A)') '/'
        close(unit)
    end subroutine write_config_chartmap

    subroutine load_case(config_path)
        character(len=*), intent(in) :: config_path
        character(len=256) :: cfg_local

        cfg_local = trim(config_path)
        call read_config(cfg_local)
        call init_field(norb, netcdffile, ns_s, ns_tp, multharm, integmode)
        call params_init
        call init_magfie(VMEC)
    end subroutine load_case

    subroutine integrate_orbit_from_start(start_path, orbit_index, n_used, s_arr, th_arr, ph_arr, R_arr, Z_arr, B_arr)
        character(len=*), intent(in) :: start_path
        integer, intent(in) :: orbit_index
        integer, intent(inout) :: n_used(norbits)
        real(dp), intent(inout) :: s_arr(:, :), th_arr(:, :), ph_arr(:, :)
        real(dp), intent(inout) :: R_arr(:, :), Z_arr(:, :), B_arr(:, :)

        integer :: i_local, ierr_local
        real(dp) :: z_local(5)
        real(dp) :: xyz(3), xcyl(3)
        real(dp) :: bmod_local, sqrtg_local
        real(dp) :: bder(3), hcov(3), hctrvr(3), hcurl(3)

        call load_starting_points(zstart, trim(start_path))
        z_local = zstart(:, 1)

        ierr_local = 0
        n_used(orbit_index) = 0
        do i_local = 1, nmax
            s_arr(i_local, orbit_index) = z_local(1)
            th_arr(i_local, orbit_index) = z_local(2)
            ph_arr(i_local, orbit_index) = z_local(3)

            call ref_coords%evaluate_point((/ z_local(1), z_local(2), z_local(3) /), xyz)
            call cart_to_cyl(xyz, xcyl)
            R_arr(i_local, orbit_index) = xcyl(1)
            Z_arr(i_local, orbit_index) = xcyl(3)

            call magfie((/ z_local(1), z_local(2), z_local(3) /), bmod_local, sqrtg_local, bder, hcov, hctrvr, hcurl)
            B_arr(i_local, orbit_index) = bmod_local

            n_used(orbit_index) = i_local
            call orbit_timestep_axis(z_local, dtaumin, dtaumin, relerr, ierr_local)
            if (ierr_local /= 0) exit
            if (z_local(1) < 0.0_dp .or. z_local(1) > 1.0_dp) exit
        end do
    end subroutine integrate_orbit_from_start

    subroutine cart_to_cyl(xyz, xcyl)
        real(dp), intent(in) :: xyz(3)
        real(dp), intent(out) :: xcyl(3)

        xcyl(1) = sqrt(xyz(1)*xyz(1) + xyz(2)*xyz(2))
        xcyl(2) = atan2(xyz(2), xyz(1))
        xcyl(3) = xyz(3)
    end subroutine cart_to_cyl

    subroutine build_theta_surf_grids(theta, svals)
        real(dp), intent(out) :: theta(ntheta_plot)
        real(dp), intent(out) :: svals(nsurf_plot)
        integer :: it

        do it = 1, ntheta_plot
            theta(it) = (real(it - 1, dp) / real(ntheta_plot - 1, dp)) * twopi
        end do
        svals = [0.10_dp, 0.25_dp, 0.40_dp, 0.60_dp, 0.80_dp, 0.95_dp]
    end subroutine build_theta_surf_grids

    subroutine build_flux_surfaces(theta, svals, R_surf, Z_surf)
        real(dp), intent(in) :: theta(ntheta_plot)
        real(dp), intent(in) :: svals(nsurf_plot)
        real(dp), intent(out) :: R_surf(ntheta_plot, nsurf_plot)
        real(dp), intent(out) :: Z_surf(ntheta_plot, nsurf_plot)

        integer :: isurf, it
        real(dp) :: xyz(3), xcyl(3)

        do isurf = 1, nsurf_plot
            do it = 1, ntheta_plot
                call ref_coords%evaluate_point((/ svals(isurf), theta(it), 0.0_dp /), xyz)
                call cart_to_cyl(xyz, xcyl)
                R_surf(it, isurf) = xcyl(1)
                Z_surf(it, isurf) = xcyl(3)
            end do
        end do
    end subroutine build_flux_surfaces

    subroutine plot_case(out_case, title_prefix, t, n_used, R, Z, s, th, ph, B, R_surf, Z_surf, theta, svals)
        character(len=*), intent(in) :: out_case, title_prefix
        real(dp), intent(in) :: t(:)
        integer, intent(in) :: n_used(norbits)
        real(dp), intent(in) :: R(:, :), Z(:, :), s(:, :), th(:, :), ph(:, :), B(:, :)
        real(dp), intent(in) :: R_surf(ntheta_plot, nsurf_plot), Z_surf(ntheta_plot, nsurf_plot)
        real(dp), intent(in) :: theta(ntheta_plot), svals(nsurf_plot)

        integer :: isurf

        call plt%initialize(grid=.true., xlabel='R (cm)', ylabel='Z (cm)', &
            title=trim(title_prefix)//': orbit projection (R,Z)', legend=.true., figsize=[10, 8])
        call plt%add_plot(R(1:n_used(1), 1), Z(1:n_used(1), 1), label='passing', linestyle='-')
        call plt%add_plot(R(1:n_used(2), 2), Z(1:n_used(2), 2), label='trapped', linestyle='-')
        call plt%savefig(trim(out_case)//'/orbit/orbit_RZ.png', pyfile=trim(out_case)//'/orbit/orbit_RZ.py')

        call plt%initialize(grid=.true., xlabel='t (s) [scaled]', ylabel='s', &
            title=trim(title_prefix)//': s(t)', figsize=[10, 6])
        call plt%add_plot(t(1:n_used(1)), s(1:n_used(1), 1), label='passing', linestyle='-')
        call plt%add_plot(t(1:n_used(2)), s(1:n_used(2), 2), label='trapped', linestyle='-')
        call plt%savefig(trim(out_case)//'/orbit/orbit_s_t.png', pyfile=trim(out_case)//'/orbit/orbit_s_t.py')

        call plt%initialize(grid=.true., xlabel='t (s) [scaled]', ylabel='Bmod (G)', &
            title=trim(title_prefix)//': Bmod(t)', figsize=[10, 6])
        call plt%add_plot(t(1:n_used(1)), B(1:n_used(1), 1), label='passing', linestyle='-')
        call plt%add_plot(t(1:n_used(2)), B(1:n_used(2), 2), label='trapped', linestyle='-')
        call plt%savefig(trim(out_case)//'/orbit/orbit_Bmod_t.png', pyfile=trim(out_case)//'/orbit/orbit_Bmod_t.py')

        call plt%initialize(grid=.true., xlabel='R (cm)', ylabel='Z (cm)', &
            title=trim(title_prefix)//': flux surfaces at phi=0', legend=.true., figsize=[10, 8])
        do isurf = 1, nsurf_plot
            call plt%add_plot(R_surf(:, isurf), Z_surf(:, isurf), label='surface', linestyle='-')
        end do
        call plt%add_plot(R(1:n_used(1), 1), Z(1:n_used(1), 1), label='passing', linestyle='-')
        call plt%add_plot(R(1:n_used(2), 2), Z(1:n_used(2), 2), label='trapped', linestyle='-')
        call plt%savefig(trim(out_case)//'/flux_surfaces/flux_surfaces_RZ_phi0.png', pyfile=trim(out_case)//'/flux_surfaces/flux_surfaces_RZ_phi0.py')
    end subroutine plot_case

    subroutine plot_compare(out_cmp, t, n_geo, n_cm, Rg, Zg, Rc, Zc, cgeo, ccm)
        character(len=*), intent(in) :: out_cmp
        real(dp), intent(in) :: t(:)
        integer, intent(in) :: n_geo(norbits), n_cm(norbits)
        real(dp), intent(in) :: Rg(:, :), Zg(:, :), Rc(:, :), Zc(:, :)
        real(dp), intent(in) :: cgeo(3), ccm(3)

        integer :: n1, n2
        real(dp), allocatable :: dR_pass(:), dZ_pass(:), dR_trap(:), dZ_trap(:)

        call plt%initialize(grid=.true., xlabel='R (cm)', ylabel='Z (cm)', &
            title='GEQDSK RK45: orbit overlay (geoflux vs chartmap)', legend=.true., figsize=[10, 8])
        call plt%add_plot(Rg(1:n_geo(1), 1), Zg(1:n_geo(1), 1), label='geoflux passing', linestyle='-', color=cgeo)
        call plt%add_plot(Rg(1:n_geo(2), 2), Zg(1:n_geo(2), 2), label='geoflux trapped', linestyle='--', color=cgeo)
        call plt%add_plot(Rc(1:n_cm(1), 1), Zc(1:n_cm(1), 1), label='chartmap passing', linestyle='-', color=ccm)
        call plt%add_plot(Rc(1:n_cm(2), 2), Zc(1:n_cm(2), 2), label='chartmap trapped', linestyle='--', color=ccm)
        call plt%savefig(trim(out_cmp)//'/orbit_RZ_overlay.png', pyfile=trim(out_cmp)//'/orbit_RZ_overlay.py')

        n1 = min(n_geo(1), n_cm(1))
        n2 = min(n_geo(2), n_cm(2))
        allocate(dR_pass(n1), dZ_pass(n1), dR_trap(n2), dZ_trap(n2))
        dR_pass = Rc(1:n1, 1) - Rg(1:n1, 1)
        dZ_pass = Zc(1:n1, 1) - Zg(1:n1, 1)
        dR_trap = Rc(1:n2, 2) - Rg(1:n2, 2)
        dZ_trap = Zc(1:n2, 2) - Zg(1:n2, 2)

        call plt%initialize(grid=.true., xlabel='t (s) [scaled]', ylabel='delta (cm)', &
            title='GEQDSK RK45: chartmap - geoflux deltaR, deltaZ', legend=.true., figsize=[10, 6])
        call plt%add_plot(t(1:n1), dR_pass, label='passing deltaR', linestyle='-')
        call plt%add_plot(t(1:n1), dZ_pass, label='passing deltaZ', linestyle='-')
        call plt%add_plot(t(1:n2), dR_trap, label='trapped deltaR', linestyle='--')
        call plt%add_plot(t(1:n2), dZ_trap, label='trapped deltaZ', linestyle='--')
        call plt%savefig(trim(out_cmp)//'/orbit_deltaRZ_t.png', pyfile=trim(out_cmp)//'/orbit_deltaRZ_t.py')

        deallocate(dR_pass, dZ_pass, dR_trap, dZ_trap)

        call plt%initialize(grid=.true., xlabel='R (cm)', ylabel='Z (cm)', &
            title='GEQDSK: flux surfaces overlay (geoflux vs chartmap)', legend=.true., figsize=[10, 8])
        call plt%add_plot(R_surf_geo(:, nsurf_plot), Z_surf_geo(:, nsurf_plot), label='geoflux outer', linestyle='-', color=cgeo)
        call plt%add_plot(R_surf_cm(:, nsurf_plot), Z_surf_cm(:, nsurf_plot), label='chartmap outer', linestyle='-', color=ccm)
        call plt%savefig(trim(out_cmp)//'/flux_surfaces_RZ_phi0_overlay.png', pyfile=trim(out_cmp)//'/flux_surfaces_RZ_phi0_overlay.py')
    end subroutine plot_compare

    subroutine write_start(path, s0, th0, ph0, p0, lam0)
        character(len=*), intent(in) :: path
        real(dp), intent(in) :: s0, th0, ph0, p0, lam0
        integer :: unit

        open(newunit=unit, file=trim(path), status='replace', action='write')
        write(unit, *) s0, th0, ph0, p0, lam0
        close(unit)
    end subroutine write_start

    subroutine write_chartmap_from_geoflux(filename, nrho, ntheta, nzeta)
        character(len=*), intent(in) :: filename
        integer, intent(in) :: nrho, ntheta, nzeta

        integer :: ncid, dim_rho, dim_theta, dim_zeta
        integer :: var_rho, var_theta, var_zeta
        integer :: var_x, var_y, var_z
        integer :: dimids_xyz(3)
        integer :: ierr
        integer :: i_r, i_t, i_z
        real(dp), allocatable :: rho(:), theta(:), zeta(:)
        real(dp), allocatable :: x(:, :, :), y(:, :, :), z(:, :, :)
        real(dp) :: u(3), xyz(3)

        allocate(rho(nrho), theta(ntheta), zeta(nzeta))
        allocate(x(nrho, ntheta, nzeta), y(nrho, ntheta, nzeta), z(nrho, ntheta, nzeta))

        do i_r = 1, nrho
            rho(i_r) = real(i_r - 1, dp)/real(nrho - 1, dp)
        end do
        do i_t = 1, ntheta
            theta(i_t) = (real(i_t - 1, dp)/real(ntheta - 1, dp))*twopi
        end do
        do i_z = 1, nzeta
            zeta(i_z) = (real(i_z - 1, dp)/real(nzeta - 1, dp))*twopi
        end do

        do i_z = 1, nzeta
            do i_t = 1, ntheta
                do i_r = 1, nrho
                    u = [rho(i_r), theta(i_t), zeta(i_z)]
                    call geoflux_to_cart(u, xyz)
                    x(i_r, i_t, i_z) = xyz(1)
                    y(i_r, i_t, i_z) = xyz(2)
                    z(i_r, i_t, i_z) = xyz(3)
                end do
            end do
        end do

        ierr = nf90_create(trim(filename), nf90_clobber, ncid)
        call nc_check(ierr, 'nf90_create')

        ierr = nf90_def_dim(ncid, 'rho', nrho, dim_rho)
        call nc_check(ierr, 'def_dim rho')
        ierr = nf90_def_dim(ncid, 'theta', ntheta, dim_theta)
        call nc_check(ierr, 'def_dim theta')
        ierr = nf90_def_dim(ncid, 'zeta', nzeta, dim_zeta)
        call nc_check(ierr, 'def_dim zeta')

        ierr = nf90_def_var(ncid, 'rho', nf90_double, (/dim_rho/), var_rho)
        call nc_check(ierr, 'def_var rho')
        ierr = nf90_def_var(ncid, 'theta', nf90_double, (/dim_theta/), var_theta)
        call nc_check(ierr, 'def_var theta')
        ierr = nf90_def_var(ncid, 'zeta', nf90_double, (/dim_zeta/), var_zeta)
        call nc_check(ierr, 'def_var zeta')

        dimids_xyz = (/dim_rho, dim_theta, dim_zeta/)
        ierr = nf90_def_var(ncid, 'x', nf90_double, dimids_xyz, var_x)
        call nc_check(ierr, 'def_var x')
        ierr = nf90_def_var(ncid, 'y', nf90_double, dimids_xyz, var_y)
        call nc_check(ierr, 'def_var y')
        ierr = nf90_def_var(ncid, 'z', nf90_double, dimids_xyz, var_z)
        call nc_check(ierr, 'def_var z')

        ierr = nf90_enddef(ncid)
        call nc_check(ierr, 'enddef')

        ierr = nf90_put_var(ncid, var_rho, rho)
        call nc_check(ierr, 'put rho')
        ierr = nf90_put_var(ncid, var_theta, theta)
        call nc_check(ierr, 'put theta')
        ierr = nf90_put_var(ncid, var_zeta, zeta)
        call nc_check(ierr, 'put zeta')
        ierr = nf90_put_var(ncid, var_x, x)
        call nc_check(ierr, 'put x')
        ierr = nf90_put_var(ncid, var_y, y)
        call nc_check(ierr, 'put y')
        ierr = nf90_put_var(ncid, var_z, z)
        call nc_check(ierr, 'put z')

        ierr = nf90_close(ncid)
        call nc_check(ierr, 'close')

        deallocate(rho, theta, zeta, x, y, z)
    end subroutine write_chartmap_from_geoflux

    subroutine nc_check(ierr, what)
        integer, intent(in) :: ierr
        character(len=*), intent(in) :: what

        if (ierr /= nf90_noerr) then
            print *, 'NetCDF error in ', trim(what), ': ', trim(nf90_strerror(ierr))
            error stop 'NetCDF failure'
        end if
    end subroutine nc_check

end program test_geoflux_rk45_orbit_plot_chartmap_compare
