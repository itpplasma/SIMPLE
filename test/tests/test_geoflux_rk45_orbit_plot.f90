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
    use geoflux_field, only: splint_geoflux_field
    use field_sub, only: field_eq, psif
    use pyplot_module, only: pyplot
    use util, only: twopi

    implicit none

    type(tracer_t) :: norb
    type(pyplot) :: plt
    character(len=1024) :: out_root, out_dir, out_orbit, out_flux, out_field
    character(len=1024) :: config_file, start_file, geqdsk_file
    character(len=1024) :: png_orbit_rz, png_s_t, png_theta_t, png_phi_t, png_bmod_t
    character(len=1024) :: png_flux_rz, png_psi_rz, png_bmod_st
    character(len=1024) :: png_bmod_rz, png_br_rz, png_bphi_rz, png_bz_rz
    character(len=1024) :: traj_dat
    character(len=1024) :: cmd
    integer :: status, ierr, i, itheta, isurf, nmax, n_used, mkdir_stat
    integer :: ntheta_plot, nsurf_plot, ns_plot
    real(dp), allocatable :: s_traj(:), theta_traj(:), phi_traj(:), time_traj(:)
    real(dp), allocatable :: r_traj(:), z_traj(:), bmod_traj(:)
    real(dp), allocatable :: theta_grid(:), s_grid(:)
    real(dp), allocatable :: bmod_grid(:, :)
    real(dp), allocatable :: r_surf(:, :), z_surf(:, :)
    real(dp), allocatable :: surf_s(:)
    real(dp), allocatable :: r_grid(:), z_grid(:)
    real(dp), allocatable :: psi_rz(:, :), bmod_rz(:, :)
    real(dp), allocatable :: br_rz(:, :), bphi_rz(:, :), bz_rz(:, :)
    real(dp) :: z(5), xcyl(3)
    real(dp) :: acov_tmp(3), hcov_tmp(3), sqg_tmp(3)
    real(dp) :: smin, smax, rmin, rmax, zmin, zmax
    logical :: exists

    out_root = ''
    call get_environment_variable('SIMPLE_ARTIFACT_DIR', value=out_root, status=status)
    if (status /= 0 .or. len_trim(out_root) == 0) then
        out_root = '/tmp/SIMPLE_artifacts'
    end if

    out_dir = trim(out_root)//'/plot/test_geoflux_rk45_orbit_plot'
    out_orbit = trim(out_dir)//'/orbit'
    out_flux = trim(out_dir)//'/flux_surfaces'
    out_field = trim(out_dir)//'/fields'

    cmd = 'mkdir -p '//trim(out_orbit)//' '//trim(out_flux)//' '//trim(out_field)
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
    allocate(r_traj(nmax), z_traj(nmax), bmod_traj(nmax))

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
        call splint_geoflux_field(z(1), z(2), z(3), acov_tmp, hcov_tmp, bmod_traj(i), sqg_tmp)

        n_used = i
        call orbit_timestep_axis(z, dtaumin, dtaumin, relerr, ierr)
        if (ierr /= 0) exit
        if (z(1) < 0.0_dp .or. z(1) > 1.0_dp) exit
    end do

    if (n_used < 50) then
        error stop 'test_geoflux_rk45_orbit_plot: orbit produced too few points'
    end if

    call compute_ranges(s_traj(1:n_used), r_traj(1:n_used), z_traj(1:n_used), smin, smax, rmin, rmax, zmin, zmax)
    if (.not. (rmax > rmin .and. zmax > zmin .and. smax > smin)) then
        error stop 'test_geoflux_rk45_orbit_plot: orbit appears degenerate'
    end if

    png_orbit_rz = trim(out_orbit)//'/orbit_RZ.png'
    png_s_t = trim(out_orbit)//'/orbit_s_t.png'
    png_theta_t = trim(out_orbit)//'/orbit_theta_t.png'
    png_phi_t = trim(out_orbit)//'/orbit_phi_t.png'
    png_bmod_t = trim(out_orbit)//'/orbit_Bmod_t.png'
    traj_dat = trim(out_orbit)//'/trajectory.dat'

    call plt%initialize(grid=.true., xlabel='R (cm)', ylabel='Z (cm)', &
        title='GEQDSK geoflux RK45 orbit projection (R,Z)', legend=.true., figsize=[10, 8])
    call plt%add_plot(r_traj(1:n_used), z_traj(1:n_used), label='orbit', linestyle='-')
    call plt%savefig(trim(png_orbit_rz), pyfile=trim(out_orbit)//'/orbit_RZ.py')

    call write_trajectory_table(traj_dat, time_traj(1:n_used), s_traj(1:n_used), theta_traj(1:n_used), &
        phi_traj(1:n_used), r_traj(1:n_used), z_traj(1:n_used), bmod_traj(1:n_used))

    call plt%initialize(grid=.true., xlabel='t (s) [scaled]', ylabel='s', &
        title='GEQDSK geoflux RK45 orbit: s(t)', figsize=[10, 6])
    call plt%add_plot(time_traj(1:n_used), s_traj(1:n_used), label='s', linestyle='-')
    call plt%savefig(trim(png_s_t), pyfile=trim(out_orbit)//'/orbit_s_t.py')

    call plt%initialize(grid=.true., xlabel='t (s) [scaled]', ylabel='theta (rad)', &
        title='GEQDSK geoflux RK45 orbit: theta(t)', figsize=[10, 6])
    call plt%add_plot(time_traj(1:n_used), theta_traj(1:n_used), label='theta', linestyle='-')
    call plt%savefig(trim(png_theta_t), pyfile=trim(out_orbit)//'/orbit_theta_t.py')

    call plt%initialize(grid=.true., xlabel='t (s) [scaled]', ylabel='phi (rad)', &
        title='GEQDSK geoflux RK45 orbit: phi(t)', figsize=[10, 6])
    call plt%add_plot(time_traj(1:n_used), phi_traj(1:n_used), label='phi', linestyle='-')
    call plt%savefig(trim(png_phi_t), pyfile=trim(out_orbit)//'/orbit_phi_t.py')

    call plt%initialize(grid=.true., xlabel='t (s) [scaled]', ylabel='Bmod (G)', &
        title='GEQDSK geoflux RK45 orbit: Bmod(t)', figsize=[10, 6])
    call plt%add_plot(time_traj(1:n_used), bmod_traj(1:n_used), label='Bmod', linestyle='-')
    call plt%savefig(trim(png_bmod_t), pyfile=trim(out_orbit)//'/orbit_Bmod_t.py')

    png_flux_rz = trim(out_flux)//'/flux_surfaces_RZ_phi0.png'

    nsurf_plot = 6
    ntheta_plot = 361
    ns_plot = 128
    allocate(surf_s(nsurf_plot))
    surf_s = [0.10_dp, 0.25_dp, 0.40_dp, 0.60_dp, 0.80_dp, 0.95_dp]
    allocate(r_surf(ntheta_plot, nsurf_plot), z_surf(ntheta_plot, nsurf_plot))
    allocate(theta_grid(ntheta_plot))

    do itheta = 1, ntheta_plot
        theta_grid(itheta) = (real(itheta - 1, dp) / real(ntheta_plot - 1, dp)) * twopi
    end do

    do isurf = 1, nsurf_plot
        do itheta = 1, ntheta_plot
            call geoflux_to_cyl((/ surf_s(isurf), theta_grid(itheta), 0.0_dp /), xcyl)
            r_surf(itheta, isurf) = xcyl(1)
            z_surf(itheta, isurf) = xcyl(3)
        end do
    end do

    call plt%initialize(grid=.true., xlabel='R (cm)', ylabel='Z (cm)', &
        title='GEQDSK flux surfaces at phi=0 with orbit overlay', legend=.true., figsize=[10, 8])
    do isurf = 1, nsurf_plot
        call plt%add_plot(r_surf(:, isurf), z_surf(:, isurf), label='surface', linestyle='-')
    end do
    call plt%add_plot(r_traj(1:n_used), z_traj(1:n_used), label='orbit', linestyle='-')
    call plt%savefig(trim(png_flux_rz), pyfile=trim(out_flux)//'/flux_surfaces_overlay.py')

    call build_rz_field_maps(plt, out_field, out_flux, r_surf(:, nsurf_plot), z_surf(:, nsurf_plot), &
        r_traj(1:n_used), z_traj(1:n_used))

    allocate(s_grid(ns_plot))
    allocate(bmod_grid(nsurf_plot, ntheta_plot))
    do isurf = 1, nsurf_plot
        do itheta = 1, ntheta_plot
            call splint_geoflux_field(surf_s(isurf), theta_grid(itheta), 0.0_dp, acov_tmp, hcov_tmp, bmod_grid(isurf, itheta), sqg_tmp)
        end do
    end do

    png_bmod_st = trim(out_field)//'/Bmod_s_theta_phi0.png'
    call plt%initialize(grid=.true., xlabel='theta (rad)', ylabel='surface index', &
        title='GEQDSK Bmod(s,theta) at phi=0 (sampled on flux surfaces)', figsize=[10, 6])
    call plt%add_imshow(bmod_grid)
    call plt%savefig(trim(png_bmod_st), pyfile=trim(out_field)//'/Bmod_s_theta_phi0.py')

    inquire(file=trim(png_orbit_rz), exist=exists)
    if (.not. exists) then
        error stop 'test_geoflux_rk45_orbit_plot: missing RZ plot output'
    end if
    inquire(file=trim(png_flux_rz), exist=exists)
    if (.not. exists) then
        error stop 'test_geoflux_rk45_orbit_plot: missing flux surfaces plot output'
    end if

    print *, 'ARTIFACT_DIR: ', trim(out_dir)
    print *, 'ARTIFACT: ', trim(png_orbit_rz)
    print *, 'ARTIFACT: ', trim(png_s_t)
    print *, 'ARTIFACT: ', trim(png_theta_t)
    print *, 'ARTIFACT: ', trim(png_phi_t)
    print *, 'ARTIFACT: ', trim(png_bmod_t)
    print *, 'ARTIFACT: ', trim(png_flux_rz)
    print *, 'ARTIFACT: ', trim(png_bmod_st)
    print *, 'ARTIFACT: ', trim(traj_dat)

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
        write(unit, '(A)') 'ntimstep = 800'
        write(unit, '(A)') 'npoiper2 = 64'
        write(unit, '(A)') 'multharm = 3'
        write(unit, '(A)') 'trace_time = 1d-6'
        write(unit, '(A)') 'relerr = 1d-11'
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
        lam0 = 0.7_dp

        open(newunit=unit, file=trim(path), status='replace', action='write')
        write(unit, *) s0, th0, ph0, p0, lam0
        close(unit)
    end subroutine write_start

    subroutine compute_ranges(s_arr, r_arr, z_arr, smin, smax, rmin, rmax, zmin, zmax)
        real(dp), intent(in) :: s_arr(:), r_arr(:), z_arr(:)
        real(dp), intent(out) :: smin, smax, rmin, rmax, zmin, zmax

        smin = minval(s_arr)
        smax = maxval(s_arr)
        rmin = minval(r_arr)
        rmax = maxval(r_arr)
        zmin = minval(z_arr)
        zmax = maxval(z_arr)
    end subroutine compute_ranges

    subroutine write_trajectory_table(path, t, s, theta, phi, r, zc, bmod)
        character(len=*), intent(in) :: path
        real(dp), intent(in) :: t(:), s(:), theta(:), phi(:), r(:), zc(:), bmod(:)
        integer :: unit, i

        open(newunit=unit, file=trim(path), status='replace', action='write')
        write(unit, '(A)') '# t_scaled  s  theta  phi  R_cm  Z_cm  Bmod_G'
        do i = 1, size(t)
            write(unit, '(7ES22.14)') t(i), s(i), theta(i), phi(i), r(i), zc(i), bmod(i)
        end do
        close(unit)
    end subroutine write_trajectory_table

    subroutine build_rz_field_maps(plt, out_field_dir, out_flux_dir, r_lcfs, z_lcfs, r_orbit, z_orbit)
        type(pyplot), intent(inout) :: plt
        character(len=*), intent(in) :: out_field_dir, out_flux_dir
        real(dp), intent(in) :: r_lcfs(:), z_lcfs(:)
        real(dp), intent(in) :: r_orbit(:), z_orbit(:)

        integer, parameter :: nr = 220, nz = 220
        real(dp) :: rmin_g, rmax_g, zmin_g, zmax_g, dr, dz
        real(dp) :: phi0
        real(dp), allocatable :: rgrid(:), zgrid(:)
        real(dp), allocatable :: psi_map(:, :), bmod_map(:, :)
        real(dp), allocatable :: br_map(:, :), bphi_map(:, :), bz_map(:, :)
        integer :: i, j
        real(dp) :: br, bphi, bz, d1, d2, d3, d4, d5, d6, d7, d8, d9

        phi0 = 0.0_dp

        rmin_g = min(minval(r_lcfs), minval(r_orbit))
        rmax_g = max(maxval(r_lcfs), maxval(r_orbit))
        zmin_g = min(minval(z_lcfs), minval(z_orbit))
        zmax_g = max(maxval(z_lcfs), maxval(z_orbit))

        rmin_g = rmin_g - 0.1_dp * (rmax_g - rmin_g)
        rmax_g = rmax_g + 0.1_dp * (rmax_g - rmin_g)
        zmin_g = zmin_g - 0.1_dp * (zmax_g - zmin_g)
        zmax_g = zmax_g + 0.1_dp * (zmax_g - zmin_g)

        allocate(rgrid(nr), zgrid(nz))
        allocate(psi_map(nr, nz), bmod_map(nr, nz))
        allocate(br_map(nr, nz), bphi_map(nr, nz), bz_map(nr, nz))

        dr = (rmax_g - rmin_g) / real(nr - 1, dp)
        dz = (zmax_g - zmin_g) / real(nz - 1, dp)

        do i = 1, nr
            rgrid(i) = rmin_g + real(i - 1, dp) * dr
        end do
        do j = 1, nz
            zgrid(j) = zmin_g + real(j - 1, dp) * dz
        end do

        do j = 1, nz
            do i = 1, nr
                call field_eq(rgrid(i), phi0, zgrid(j), br, bphi, bz, d1, d2, d3, d4, d5, d6, d7, d8, d9)
                psi_map(i, j) = psif
                br_map(i, j) = br
                bphi_map(i, j) = bphi
                bz_map(i, j) = bz
                bmod_map(i, j) = sqrt(br * br + bphi * bphi + bz * bz)
            end do
        end do

        call plt%initialize(grid=.true., xlabel='R (cm)', ylabel='Z (cm)', &
            title='GEQDSK psi(R,Z) contours at phi=0 with orbit overlay', legend=.true., figsize=[10, 8])
        call plt%add_contour(rgrid, zgrid, transpose(psi_map), linestyle='-', colorbar=.false.)
        call plt%add_plot(r_lcfs, z_lcfs, label='LCFS approx', linestyle='-')
        call plt%add_plot(r_orbit, z_orbit, label='orbit', linestyle='-')
        call plt%savefig(trim(out_flux_dir)//'/psi_contours_RZ_phi0.png', pyfile=trim(out_flux_dir)//'/psi_contours_RZ_phi0.py')

        call plt%initialize(grid=.true., xlabel='R (cm)', ylabel='Z (cm)', &
            title='GEQDSK |B|(R,Z) at phi=0 with orbit overlay', legend=.true., figsize=[10, 8])
        call plt%add_contour(rgrid, zgrid, transpose(bmod_map), linestyle='-', colorbar=.false.)
        call plt%add_plot(r_lcfs, z_lcfs, label='LCFS approx', linestyle='-')
        call plt%add_plot(r_orbit, z_orbit, label='orbit', linestyle='-')
        call plt%savefig(trim(out_field_dir)//'/Bmod_RZ_phi0.png', pyfile=trim(out_field_dir)//'/Bmod_RZ_phi0.py')

        call plt%initialize(grid=.true., xlabel='R (cm)', ylabel='Z (cm)', &
            title='GEQDSK Br(R,Z) at phi=0', legend=.false., figsize=[10, 8])
        call plt%add_contour(rgrid, zgrid, transpose(br_map), linestyle='-', colorbar=.false.)
        call plt%savefig(trim(out_field_dir)//'/Br_RZ_phi0.png', pyfile=trim(out_field_dir)//'/Br_RZ_phi0.py')

        call plt%initialize(grid=.true., xlabel='R (cm)', ylabel='Z (cm)', &
            title='GEQDSK Bphi(R,Z) at phi=0', legend=.false., figsize=[10, 8])
        call plt%add_contour(rgrid, zgrid, transpose(bphi_map), linestyle='-', colorbar=.false.)
        call plt%savefig(trim(out_field_dir)//'/Bphi_RZ_phi0.png', pyfile=trim(out_field_dir)//'/Bphi_RZ_phi0.py')

        call plt%initialize(grid=.true., xlabel='R (cm)', ylabel='Z (cm)', &
            title='GEQDSK Bz(R,Z) at phi=0', legend=.false., figsize=[10, 8])
        call plt%add_contour(rgrid, zgrid, transpose(bz_map), linestyle='-', colorbar=.false.)
        call plt%savefig(trim(out_field_dir)//'/Bz_RZ_phi0.png', pyfile=trim(out_field_dir)//'/Bz_RZ_phi0.py')

        deallocate(rgrid, zgrid, psi_map, bmod_map, br_map, bphi_map, bz_map)
    end subroutine build_rz_field_maps

end program test_geoflux_rk45_orbit_plot
