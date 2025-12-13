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
    character(len=1024) :: config_file, start_pass_file, start_trap_file, geqdsk_file
    character(len=256) :: config_file_256
    character(len=1024) :: png_orbit_rz, png_s_t, png_theta_t, png_phi_t, png_bmod_t
    character(len=1024) :: png_flux_rz, png_psi_rz, png_bmod_st
    character(len=1024) :: png_bmod_rz, png_br_rz, png_bphi_rz, png_bz_rz
    character(len=1024) :: traj_dat
    character(len=1024) :: cmd
    integer, parameter :: norbits = 2
    integer :: status, ierr, i, itheta, isurf, nmax, mkdir_stat
    integer :: iorb
    integer :: n_used(norbits)
    integer :: ntheta_plot, nsurf_plot, ns_plot
    real(dp), allocatable :: time_traj(:)
    real(dp), allocatable :: s_traj(:, :), theta_traj(:, :), phi_traj(:, :)
    real(dp), allocatable :: r_traj(:, :), z_traj(:, :), bmod_traj(:, :)
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
    real(dp) :: color_pass(3), color_trap(3)
    character(len=32) :: orbit_label(norbits)
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
    start_pass_file = trim(out_dir)//'/start_passing.dat'
    start_trap_file = trim(out_dir)//'/start_trapped.dat'

    call write_config(trim(config_file), trim(geqdsk_file))
    call write_start(trim(start_pass_file), 0.25_dp, 0.1_dp*twopi, 0.0_dp, 1.0_dp, 0.7_dp)
    call write_start(trim(start_trap_file), 0.25_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp)

    config_file_256 = trim(config_file)
    call read_config(config_file_256)
    call init_field(norb, netcdffile, ns_s, ns_tp, multharm, integmode)
    call params_init
    call init_magfie(isw_field_type)

    orbit_label(1) = 'passing'
    orbit_label(2) = 'trapped'
    color_pass = [0.0_dp, 0.0_dp, 1.0_dp]
    color_trap = [1.0_dp, 0.0_dp, 0.0_dp]

    nmax = ntimstep
    allocate(time_traj(nmax))
    allocate(s_traj(nmax, norbits), theta_traj(nmax, norbits), phi_traj(nmax, norbits))
    allocate(r_traj(nmax, norbits), z_traj(nmax, norbits), bmod_traj(nmax, norbits))

    do i = 1, nmax
        time_traj(i) = real(i - 1, dp) * dtaumin / max(v0, 1.0d-12)
    end do

    call integrate_orbit_from_start(trim(start_pass_file), 1)
    call integrate_orbit_from_start(trim(start_trap_file), 2)

    do iorb = 1, norbits
        if (n_used(iorb) < 50) then
            error stop 'test_geoflux_rk45_orbit_plot: orbit produced too few points'
        end if
    end do

    do iorb = 1, norbits
        call compute_ranges(s_traj(1:n_used(iorb), iorb), r_traj(1:n_used(iorb), iorb), z_traj(1:n_used(iorb), iorb), &
            smin, smax, rmin, rmax, zmin, zmax)
        if (.not. (rmax > rmin .and. zmax > zmin .and. smax > smin)) then
            error stop 'test_geoflux_rk45_orbit_plot: orbit appears degenerate'
        end if
    end do

    png_orbit_rz = trim(out_orbit)//'/orbit_RZ.png'
    png_s_t = trim(out_orbit)//'/orbit_s_t.png'
    png_theta_t = trim(out_orbit)//'/orbit_theta_t.png'
    png_phi_t = trim(out_orbit)//'/orbit_phi_t.png'
    png_bmod_t = trim(out_orbit)//'/orbit_Bmod_t.png'
    traj_dat = trim(out_orbit)//'/trajectory.dat'

    call plt%initialize(grid=.true., xlabel='R (cm)', ylabel='Z (cm)', &
        title='GEQDSK geoflux RK45 orbit projection (R,Z)', legend=.true., figsize=[10, 8])
    call plt%add_plot(r_traj(1:n_used(1), 1), z_traj(1:n_used(1), 1), label='passing', linestyle='-', color=color_pass)
    call plt%add_plot(r_traj(1:n_used(2), 2), z_traj(1:n_used(2), 2), label='trapped', linestyle='-', color=color_trap)
    call plt%savefig(trim(png_orbit_rz), pyfile=trim(out_orbit)//'/orbit_RZ.py')

    call write_trajectory_table(trim(out_orbit)//'/trajectory_passing.dat', time_traj(1:n_used(1)), s_traj(1:n_used(1), 1), &
        theta_traj(1:n_used(1), 1), phi_traj(1:n_used(1), 1), r_traj(1:n_used(1), 1), z_traj(1:n_used(1), 1), bmod_traj(1:n_used(1), 1))
    call write_trajectory_table(trim(out_orbit)//'/trajectory_trapped.dat', time_traj(1:n_used(2)), s_traj(1:n_used(2), 2), &
        theta_traj(1:n_used(2), 2), phi_traj(1:n_used(2), 2), r_traj(1:n_used(2), 2), z_traj(1:n_used(2), 2), bmod_traj(1:n_used(2), 2))

    call plt%initialize(grid=.true., xlabel='t (s) [scaled]', ylabel='s', &
        title='GEQDSK geoflux RK45 orbit: s(t)', figsize=[10, 6])
    call plt%add_plot(time_traj(1:n_used(1)), s_traj(1:n_used(1), 1), label='passing', linestyle='-', color=color_pass)
    call plt%add_plot(time_traj(1:n_used(2)), s_traj(1:n_used(2), 2), label='trapped', linestyle='-', color=color_trap)
    call plt%savefig(trim(png_s_t), pyfile=trim(out_orbit)//'/orbit_s_t.py')

    call plt%initialize(grid=.true., xlabel='t (s) [scaled]', ylabel='theta (rad)', &
        title='GEQDSK geoflux RK45 orbit: theta(t)', figsize=[10, 6])
    call plt%add_plot(time_traj(1:n_used(1)), theta_traj(1:n_used(1), 1), label='passing', linestyle='-', color=color_pass)
    call plt%add_plot(time_traj(1:n_used(2)), theta_traj(1:n_used(2), 2), label='trapped', linestyle='-', color=color_trap)
    call plt%savefig(trim(png_theta_t), pyfile=trim(out_orbit)//'/orbit_theta_t.py')

    call plt%initialize(grid=.true., xlabel='t (s) [scaled]', ylabel='phi (rad)', &
        title='GEQDSK geoflux RK45 orbit: phi(t)', figsize=[10, 6])
    call plt%add_plot(time_traj(1:n_used(1)), phi_traj(1:n_used(1), 1), label='passing', linestyle='-', color=color_pass)
    call plt%add_plot(time_traj(1:n_used(2)), phi_traj(1:n_used(2), 2), label='trapped', linestyle='-', color=color_trap)
    call plt%savefig(trim(png_phi_t), pyfile=trim(out_orbit)//'/orbit_phi_t.py')

    call plt%initialize(grid=.true., xlabel='t (s) [scaled]', ylabel='Bmod (G)', &
        title='GEQDSK geoflux RK45 orbit: Bmod(t)', figsize=[10, 6])
    call plt%add_plot(time_traj(1:n_used(1)), bmod_traj(1:n_used(1), 1), label='passing', linestyle='-', color=color_pass)
    call plt%add_plot(time_traj(1:n_used(2)), bmod_traj(1:n_used(2), 2), label='trapped', linestyle='-', color=color_trap)
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
    call plt%add_plot(r_traj(1:n_used(1), 1), z_traj(1:n_used(1), 1), label='passing', linestyle='-', color=color_pass)
    call plt%add_plot(r_traj(1:n_used(2), 2), z_traj(1:n_used(2), 2), label='trapped', linestyle='-', color=color_trap)
    call plt%savefig(trim(png_flux_rz), pyfile=trim(out_flux)//'/flux_surfaces_overlay.py')

    call build_rz_field_maps(plt, out_field, out_flux, r_surf(:, nsurf_plot), z_surf(:, nsurf_plot), &
        r_traj(1:n_used(1), 1), z_traj(1:n_used(1), 1), r_traj(1:n_used(2), 2), z_traj(1:n_used(2), 2), &
        color_pass, color_trap)

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
    print *, 'ARTIFACT: ', trim(out_orbit)//'/trajectory_passing.dat'
    print *, 'ARTIFACT: ', trim(out_orbit)//'/trajectory_trapped.dat'

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

    subroutine write_start(path, s0, th0, ph0, p0, lam0)
        character(len=*), intent(in) :: path
        real(dp), intent(in) :: s0, th0, ph0, p0, lam0
        integer :: unit

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

    subroutine build_rz_field_maps(plt, out_field_dir, out_flux_dir, r_lcfs, z_lcfs, &
        r_orbit_pass, z_orbit_pass, r_orbit_trap, z_orbit_trap, color_pass, color_trap)
        type(pyplot), intent(inout) :: plt
        character(len=*), intent(in) :: out_field_dir, out_flux_dir
        real(dp), intent(in) :: r_lcfs(:), z_lcfs(:)
        real(dp), intent(in) :: r_orbit_pass(:), z_orbit_pass(:)
        real(dp), intent(in) :: r_orbit_trap(:), z_orbit_trap(:)
        real(dp), intent(in) :: color_pass(3), color_trap(3)

        integer, parameter :: nr = 220, nz = 220
        real(dp) :: rmin_g, rmax_g, zmin_g, zmax_g, dr, dz
        real(dp) :: phi0
        real(dp), allocatable :: rgrid(:), zgrid(:)
        real(dp), allocatable :: psi_map(:, :), bmod_map(:, :)
        real(dp), allocatable :: br_map(:, :), bphi_map(:, :), bz_map(:, :)
        integer :: i, j
        real(dp) :: br, bphi, bz, d1, d2, d3, d4, d5, d6, d7, d8, d9

        phi0 = 0.0_dp

        rmin_g = min(minval(r_lcfs), min(minval(r_orbit_pass), minval(r_orbit_trap)))
        rmax_g = max(maxval(r_lcfs), max(maxval(r_orbit_pass), maxval(r_orbit_trap)))
        zmin_g = min(minval(z_lcfs), min(minval(z_orbit_pass), minval(z_orbit_trap)))
        zmax_g = max(maxval(z_lcfs), max(maxval(z_orbit_pass), maxval(z_orbit_trap)))

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
        call plt%add_contour(rgrid, zgrid, psi_map, linestyle='-', colorbar=.false.)
        call plt%add_plot(r_lcfs, z_lcfs, label='LCFS approx', linestyle='-')
        call plt%add_plot(r_orbit_pass, z_orbit_pass, label='passing', linestyle='-', color=color_pass)
        call plt%add_plot(r_orbit_trap, z_orbit_trap, label='trapped', linestyle='-', color=color_trap)
        call plt%savefig(trim(out_flux_dir)//'/psi_contours_RZ_phi0.png', pyfile=trim(out_flux_dir)//'/psi_contours_RZ_phi0.py')

        call plt%initialize(grid=.true., xlabel='R (cm)', ylabel='Z (cm)', &
            title='GEQDSK |B|(R,Z) at phi=0 with orbit overlay', legend=.true., figsize=[10, 8])
        call plt%add_contour(rgrid, zgrid, bmod_map, linestyle='-', colorbar=.false.)
        call plt%add_plot(r_lcfs, z_lcfs, label='LCFS approx', linestyle='-')
        call plt%add_plot(r_orbit_pass, z_orbit_pass, label='passing', linestyle='-', color=color_pass)
        call plt%add_plot(r_orbit_trap, z_orbit_trap, label='trapped', linestyle='-', color=color_trap)
        call plt%savefig(trim(out_field_dir)//'/Bmod_RZ_phi0.png', pyfile=trim(out_field_dir)//'/Bmod_RZ_phi0.py')

        call plt%initialize(grid=.true., xlabel='R (cm)', ylabel='Z (cm)', &
            title='GEQDSK Br(R,Z) at phi=0', legend=.false., figsize=[10, 8])
        call plt%add_contour(rgrid, zgrid, br_map, linestyle='-', colorbar=.false.)
        call plt%savefig(trim(out_field_dir)//'/Br_RZ_phi0.png', pyfile=trim(out_field_dir)//'/Br_RZ_phi0.py')

        call plt%initialize(grid=.true., xlabel='R (cm)', ylabel='Z (cm)', &
            title='GEQDSK Bphi(R,Z) at phi=0', legend=.false., figsize=[10, 8])
        call plt%add_contour(rgrid, zgrid, bphi_map, linestyle='-', colorbar=.false.)
        call plt%savefig(trim(out_field_dir)//'/Bphi_RZ_phi0.png', pyfile=trim(out_field_dir)//'/Bphi_RZ_phi0.py')

        call plt%initialize(grid=.true., xlabel='R (cm)', ylabel='Z (cm)', &
            title='GEQDSK Bz(R,Z) at phi=0', legend=.false., figsize=[10, 8])
        call plt%add_contour(rgrid, zgrid, bz_map, linestyle='-', colorbar=.false.)
        call plt%savefig(trim(out_field_dir)//'/Bz_RZ_phi0.png', pyfile=trim(out_field_dir)//'/Bz_RZ_phi0.py')

        deallocate(rgrid, zgrid, psi_map, bmod_map, br_map, bphi_map, bz_map)
    end subroutine build_rz_field_maps

    subroutine integrate_orbit_from_start(start_path, orbit_index)
        character(len=*), intent(in) :: start_path
        integer, intent(in) :: orbit_index

        integer :: i_local
        real(dp) :: z_local(5)
        integer :: ierr_local

        call load_starting_points(zstart, trim(start_path))
        z_local = zstart(:, 1)

        ierr_local = 0
        n_used(orbit_index) = 0
        do i_local = 1, nmax
            s_traj(i_local, orbit_index) = z_local(1)
            theta_traj(i_local, orbit_index) = z_local(2)
            phi_traj(i_local, orbit_index) = z_local(3)

            call geoflux_to_cyl((/ z_local(1), z_local(2), z_local(3) /), xcyl)
            r_traj(i_local, orbit_index) = xcyl(1)
            z_traj(i_local, orbit_index) = xcyl(3)
            call splint_geoflux_field(z_local(1), z_local(2), z_local(3), acov_tmp, hcov_tmp, bmod_traj(i_local, orbit_index), sqg_tmp)

            n_used(orbit_index) = i_local
            call orbit_timestep_axis(z_local, dtaumin, dtaumin, relerr, ierr_local)
            if (ierr_local /= 0) exit
            if (z_local(1) < 0.0_dp .or. z_local(1) > 1.0_dp) exit
        end do
    end subroutine integrate_orbit_from_start

end program test_geoflux_rk45_orbit_plot
