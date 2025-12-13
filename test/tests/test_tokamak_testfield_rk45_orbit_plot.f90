program test_tokamak_testfield_rk45_orbit_plot

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use params, only: read_config, params_init, netcdffile, ns_s, ns_tp, multharm, &
        integmode, isw_field_type, dtaumin, relerr, ntimstep, v0, zstart
    use simple, only: tracer_t
    use simple_main, only: init_field
    use magfie_sub, only: init_magfie, TEST
    use alpha_lifetime_sub, only: orbit_timestep_axis, velo_can
    use samplers, only: load_starting_points
    use pyplot_module, only: pyplot
    use util, only: twopi
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

    implicit none

    type(tracer_t) :: norb
    type(pyplot) :: plt
    character(len=1024) :: out_root, out_dir, out_orbit, out_flux, out_field
    character(len=1024) :: config_file, start_pass_file, start_trap_file
    character(len=256) :: config_file_256
    character(len=1024) :: cmd
    integer :: status, mkdir_stat

    integer, parameter :: norbits = 2
    integer :: nmax, i, iorb, ierr
    integer :: n_used(norbits)
    character(len=16) :: orbit_label(norbits)
    real(dp) :: color_pass(3), color_trap(3)

    real(dp), allocatable :: time_traj(:)
    real(dp), allocatable :: s_traj(:, :), theta_traj(:, :), phi_traj(:, :)
    real(dp), allocatable :: r_traj(:, :), z_traj(:, :), bmod_traj(:, :)
    real(dp), allocatable :: p_traj(:, :), lam_traj(:, :), mu_traj(:, :)

    out_root = ''
    call get_environment_variable('SIMPLE_ARTIFACT_DIR', value=out_root, status=status)
    if (status /= 0 .or. len_trim(out_root) == 0) then
        out_root = '/tmp/SIMPLE_artifacts'
    end if

    out_dir = trim(out_root)//'/plot/test_tokamak_testfield_rk45_orbit_plot'
    out_orbit = trim(out_dir)//'/orbit'
    out_flux = trim(out_dir)//'/flux_surfaces'
    out_field = trim(out_dir)//'/fields'

    cmd = 'mkdir -p '//trim(out_orbit)//' '//trim(out_flux)//' '//trim(out_field)
    call execute_command_line(trim(cmd), exitstat=mkdir_stat)
    if (mkdir_stat /= 0) then
        error stop 'test_tokamak_testfield_rk45_orbit_plot: failed to create output directory'
    end if

    config_file = trim(out_dir)//'/simple_tokamak_testfield_rk45_plot.in'
    start_pass_file = trim(out_dir)//'/start_passing.dat'
    start_trap_file = trim(out_dir)//'/start_trapped.dat'

    call write_config(trim(config_file))

    config_file_256 = trim(config_file)
    call read_config(config_file_256)
    call init_field(norb, netcdffile, ns_s, ns_tp, multharm, integmode)
    call params_init
    call init_magfie(TEST)

    call write_passing_start(trim(start_pass_file))
    call write_trapped_start_at_bmax(trim(start_trap_file))

    orbit_label(1) = 'passing'
    orbit_label(2) = 'trapped'
    color_pass = [0.0_dp, 0.0_dp, 1.0_dp]
    color_trap = [1.0_dp, 0.0_dp, 0.0_dp]

    nmax = ntimstep
    allocate(time_traj(nmax))
    allocate(s_traj(nmax, norbits), theta_traj(nmax, norbits), phi_traj(nmax, norbits))
    allocate(r_traj(nmax, norbits), z_traj(nmax, norbits), bmod_traj(nmax, norbits))
    allocate(p_traj(nmax, norbits), lam_traj(nmax, norbits), mu_traj(nmax, norbits))

    do i = 1, nmax
        time_traj(i) = real(i - 1, dp) * dtaumin / max(v0, 1.0d-12)
    end do

    call integrate_orbit_from_start(trim(start_pass_file), 1)
    call integrate_orbit_from_start(trim(start_trap_file), 2)

    do iorb = 1, norbits
        if (n_used(iorb) < 50) then
            error stop 'test_tokamak_testfield_rk45_orbit_plot: orbit produced too few points'
        end if
    end do

    call plot_orbit_and_diagnostics()
    call plot_flux_surfaces_and_fields()

    print *, 'ARTIFACT_DIR: ', trim(out_dir)
    call print_artifacts()

contains

    subroutine write_config(path)
        character(len=*), intent(in) :: path
        integer :: unit

        open(newunit=unit, file=trim(path), status='replace', action='write')
        write(unit, '(A)') '&config'
        write(unit, '(A)') 'netcdffile = ''wout.nc'''
        write(unit, '(A)') 'integmode = 0'
        write(unit, '(A)') 'isw_field_type = -1'
        write(unit, '(A)') 'ntestpart = 1'
        write(unit, '(A)') 'ntimstep = 800'
        write(unit, '(A)') 'npoiper2 = 128'
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

    subroutine write_passing_start(path)
        character(len=*), intent(in) :: path
        call write_start(path, 0.25_dp, 0.1_dp*twopi, 0.0_dp, 1.0_dp, 0.7_dp)
    end subroutine write_passing_start

    subroutine write_trapped_start_at_bmax(path)
        character(len=*), intent(in) :: path
        real(dp) :: s0, phi0, theta0

        s0 = 0.25_dp
        phi0 = 0.0_dp
        theta0 = 0.25_dp*twopi

        call write_start(path, s0, theta0, phi0, 1.0_dp, 0.0_dp)
    end subroutine write_trapped_start_at_bmax

    subroutine integrate_orbit_from_start(start_path, orbit_index)
        character(len=*), intent(in) :: start_path
        integer, intent(in) :: orbit_index

        real(dp) :: z_local(5)
        real(dp) :: vz_local(5)
        real(dp) :: xcyl(3)
        integer :: i_local, ierr_local
        real(dp) :: bmod_val

        call load_starting_points(zstart, trim(start_path))
        z_local = zstart(:, 1)

        call velo_can(0.0_dp, z_local, vz_local)
        if (.not. all(ieee_is_finite(vz_local))) then
            print *, 'Non-finite velo_can at start for orbit ', orbit_index
            print *, 'z = ', z_local
            print *, 'vz = ', vz_local
            error stop 'test_tokamak_testfield_rk45_orbit_plot: non-finite initial velocity'
        end if

        ierr_local = 0
        n_used(orbit_index) = 0
        do i_local = 1, nmax
            s_traj(i_local, orbit_index) = z_local(1)
            theta_traj(i_local, orbit_index) = z_local(2)
            phi_traj(i_local, orbit_index) = z_local(3)
            p_traj(i_local, orbit_index) = z_local(4)
            lam_traj(i_local, orbit_index) = z_local(5)

            call testfield_to_cyl((/z_local(1), z_local(2), z_local(3)/), xcyl)
            r_traj(i_local, orbit_index) = xcyl(1)
            z_traj(i_local, orbit_index) = xcyl(3)

            call eval_testfield_bmod(z_local(1), z_local(2), z_local(3), bmod_val)
            bmod_traj(i_local, orbit_index) = bmod_val
            mu_traj(i_local, orbit_index) = 0.5_dp*z_local(4)**2 * max(0.0_dp, 1.0_dp - z_local(5)**2) / max(bmod_val, 1.0d-12)

            n_used(orbit_index) = i_local
            call orbit_timestep_axis(z_local, dtaumin, dtaumin, relerr, ierr_local)
            if (ierr_local /= 0) exit
            if (z_local(1) < 0.0_dp .or. z_local(1) > 1.0_dp) exit
        end do
    end subroutine integrate_orbit_from_start

    subroutine testfield_to_cyl(x_geo, x_cyl)
        real(dp), intent(in) :: x_geo(3)
        real(dp), intent(out) :: x_cyl(3)

        real(dp), parameter :: R0 = 1.0_dp, a = 0.5_dp
        real(dp) :: s, theta, phi, r

        s = max(0.0_dp, min(1.0_dp, x_geo(1)))
        theta = x_geo(2)
        phi = x_geo(3)

        r = a*sqrt(s)
        x_cyl(1) = R0 + r*cos(theta)
        x_cyl(2) = phi
        x_cyl(3) = r*sin(theta)
    end subroutine testfield_to_cyl

    subroutine cyl_to_testfield(x_cyl, x_geo, inside)
        real(dp), intent(in) :: x_cyl(3)
        real(dp), intent(out) :: x_geo(3)
        logical, intent(out) :: inside

        real(dp), parameter :: R0 = 1.0_dp, a = 0.5_dp
        real(dp) :: dR, Z, r, theta

        dR = x_cyl(1) - R0
        Z = x_cyl(3)
        r = sqrt(dR*dR + Z*Z)

        inside = (r <= a)
        if (inside) then
            theta = atan2(Z, dR)
            x_geo(1) = (r/a)**2
            x_geo(2) = theta
            x_geo(3) = x_cyl(2)
        else
            x_geo = 0.0_dp
        end if
    end subroutine cyl_to_testfield

    subroutine eval_testfield_bmod(s, theta, phi, bmod)
        real(dp), intent(in) :: s, theta, phi
        real(dp), intent(out) :: bmod

        real(dp), parameter :: B0 = 1.0_dp, R0 = 1.0_dp, a = 0.5_dp, iota0 = 1.0_dp
        real(dp) :: r, R_cyl, Bphi, Bpol

        r = a*sqrt(max(0.0_dp, min(1.0_dp, s)))
        R_cyl = R0 + r*cos(theta)

        Bphi = B0*(1.0_dp - r/R0*cos(theta))
        Bpol = B0*iota0*(1.0_dp - (r/a)**2)*r/max(R_cyl, 1.0d-12)

        bmod = sqrt(Bphi*Bphi + Bpol*Bpol)
    end subroutine eval_testfield_bmod

    subroutine eval_testfield_cyl(s, theta, phi, Br, Bphi, Bz, Bmod)
        real(dp), intent(in) :: s, theta, phi
        real(dp), intent(out) :: Br, Bphi, Bz, Bmod

        real(dp), parameter :: B0 = 1.0_dp, R0 = 1.0_dp, a = 0.5_dp, iota0 = 1.0_dp
        real(dp) :: r, R_cyl, Bpol

        r = a*sqrt(max(0.0_dp, min(1.0_dp, s)))
        R_cyl = R0 + r*cos(theta)

        Bphi = B0*(1.0_dp - r/R0*cos(theta))
        Bpol = B0*iota0*(1.0_dp - (r/a)**2)*r/max(R_cyl, 1.0d-12)
        Br = -Bpol*sin(theta)
        Bz = Bpol*cos(theta)
        Bmod = sqrt(Br*Br + Bphi*Bphi + Bz*Bz)
    end subroutine eval_testfield_cyl

    subroutine plot_orbit_and_diagnostics()
        character(len=1024) :: png_orbit_rz

        png_orbit_rz = trim(out_orbit)//'/orbit_RZ.png'
        call plt%initialize(grid=.true., xlabel='R', ylabel='Z', &
            title='TEST tokamak RK45 orbit projection (R,Z)', legend=.true., figsize=[10, 8])
        call plt%add_plot(r_traj(1:n_used(1), 1), z_traj(1:n_used(1), 1), label='passing', linestyle='-', color=color_pass)
        call plt%add_plot(r_traj(1:n_used(2), 2), z_traj(1:n_used(2), 2), label='trapped', linestyle='-', color=color_trap)
        call plt%savefig(trim(png_orbit_rz), pyfile=trim(out_orbit)//'/orbit_RZ.py')

        call plt%initialize(grid=.true., xlabel='t (s) [scaled]', ylabel='s', &
            title='TEST tokamak RK45 orbit: s(t)', figsize=[10, 6])
        call plt%add_plot(time_traj(1:n_used(1)), s_traj(1:n_used(1), 1), label='passing', linestyle='-', color=color_pass)
        call plt%add_plot(time_traj(1:n_used(2)), s_traj(1:n_used(2), 2), label='trapped', linestyle='-', color=color_trap)
        call plt%savefig(trim(out_orbit)//'/orbit_s_t.png', pyfile=trim(out_orbit)//'/orbit_s_t.py')

        call plt%initialize(grid=.true., xlabel='t (s) [scaled]', ylabel='theta (rad)', &
            title='TEST tokamak RK45 orbit: theta(t)', figsize=[10, 6])
        call plt%add_plot(time_traj(1:n_used(1)), theta_traj(1:n_used(1), 1), label='passing', linestyle='-', color=color_pass)
        call plt%add_plot(time_traj(1:n_used(2)), theta_traj(1:n_used(2), 2), label='trapped', linestyle='-', color=color_trap)
        call plt%savefig(trim(out_orbit)//'/orbit_theta_t.png', pyfile=trim(out_orbit)//'/orbit_theta_t.py')

        call plt%initialize(grid=.true., xlabel='t (s) [scaled]', ylabel='phi (rad)', &
            title='TEST tokamak RK45 orbit: phi(t)', figsize=[10, 6])
        call plt%add_plot(time_traj(1:n_used(1)), phi_traj(1:n_used(1), 1), label='passing', linestyle='-', color=color_pass)
        call plt%add_plot(time_traj(1:n_used(2)), phi_traj(1:n_used(2), 2), label='trapped', linestyle='-', color=color_trap)
        call plt%savefig(trim(out_orbit)//'/orbit_phi_t.png', pyfile=trim(out_orbit)//'/orbit_phi_t.py')

        call plt%initialize(grid=.true., xlabel='t (s) [scaled]', ylabel='Bmod', &
            title='TEST tokamak RK45 orbit: Bmod(t)', figsize=[10, 6])
        call plt%add_plot(time_traj(1:n_used(1)), bmod_traj(1:n_used(1), 1), label='passing', linestyle='-', color=color_pass)
        call plt%add_plot(time_traj(1:n_used(2)), bmod_traj(1:n_used(2), 2), label='trapped', linestyle='-', color=color_trap)
        call plt%savefig(trim(out_orbit)//'/orbit_Bmod_t.png', pyfile=trim(out_orbit)//'/orbit_Bmod_t.py')

        call plt%initialize(grid=.true., xlabel='t (s) [scaled]', ylabel='p (normalized)', &
            title='TEST tokamak RK45 orbit: p(t)', figsize=[10, 6])
        call plt%add_plot(time_traj(1:n_used(1)), p_traj(1:n_used(1), 1), label='passing', linestyle='-', color=color_pass)
        call plt%add_plot(time_traj(1:n_used(2)), p_traj(1:n_used(2), 2), label='trapped', linestyle='-', color=color_trap)
        call plt%savefig(trim(out_orbit)//'/orbit_p_t.png', pyfile=trim(out_orbit)//'/orbit_p_t.py')

        call plt%initialize(grid=.true., xlabel='t (s) [scaled]', ylabel='lambda = v_par/v', &
            title='TEST tokamak RK45 orbit: lambda(t)', figsize=[10, 6])
        call plt%add_plot(time_traj(1:n_used(1)), lam_traj(1:n_used(1), 1), label='passing', linestyle='-', color=color_pass)
        call plt%add_plot(time_traj(1:n_used(2)), lam_traj(1:n_used(2), 2), label='trapped', linestyle='-', color=color_trap)
        call plt%savefig(trim(out_orbit)//'/orbit_lambda_t.png', pyfile=trim(out_orbit)//'/orbit_lambda_t.py')

        call plt%initialize(grid=.true., xlabel='t (s) [scaled]', ylabel='mu ~ p^2 (1-lambda^2) / (2 B)', &
            title='TEST tokamak RK45 orbit: mu(t) diagnostic', figsize=[10, 6])
        call plt%add_plot(time_traj(1:n_used(1)), mu_traj(1:n_used(1), 1), label='passing', linestyle='-', color=color_pass)
        call plt%add_plot(time_traj(1:n_used(2)), mu_traj(1:n_used(2), 2), label='trapped', linestyle='-', color=color_trap)
        call plt%savefig(trim(out_orbit)//'/orbit_mu_t.png', pyfile=trim(out_orbit)//'/orbit_mu_t.py')
    end subroutine plot_orbit_and_diagnostics

    subroutine plot_flux_surfaces_and_fields()
        integer, parameter :: nsurf_plot = 6, ntheta_plot = 361
        integer, parameter :: nr = 240, nz = 240
        real(dp) :: surf_s(nsurf_plot), theta_grid(ntheta_plot)
        real(dp) :: r_surf(ntheta_plot, nsurf_plot), z_surf(ntheta_plot, nsurf_plot)
        real(dp) :: rmin_g, rmax_g, zmin_g, zmax_g, dr, dz
        real(dp) :: rgrid(nr), zgrid(nz)
        real(dp) :: bmod_map(nr, nz), br_map(nr, nz), bphi_map(nr, nz), bz_map(nr, nz)
        real(dp) :: x_geo(3), x_cyl(3)
        logical :: inside
        integer :: isurf, itheta, iR, iZ
        real(dp) :: Br, Bphi, Bz, Bmod

        surf_s = [0.05_dp, 0.15_dp, 0.25_dp, 0.45_dp, 0.70_dp, 0.95_dp]
        do itheta = 1, ntheta_plot
            theta_grid(itheta) = (real(itheta - 1, dp) / real(ntheta_plot - 1, dp)) * twopi
        end do

        do isurf = 1, nsurf_plot
            do itheta = 1, ntheta_plot
                call testfield_to_cyl((/surf_s(isurf), theta_grid(itheta), 0.0_dp/), x_cyl)
                r_surf(itheta, isurf) = x_cyl(1)
                z_surf(itheta, isurf) = x_cyl(3)
            end do
        end do

        call plt%initialize(grid=.true., xlabel='R', ylabel='Z', &
            title='TEST tokamak flux surfaces at phi=0 with orbit overlay', legend=.true., figsize=[10, 8])
        do isurf = 1, nsurf_plot
            call plt%add_plot(r_surf(:, isurf), z_surf(:, isurf), label='surface', linestyle='-')
        end do
        call plt%add_plot(r_traj(1:n_used(1), 1), z_traj(1:n_used(1), 1), label='passing', linestyle='-', color=color_pass)
        call plt%add_plot(r_traj(1:n_used(2), 2), z_traj(1:n_used(2), 2), label='trapped', linestyle='-', color=color_trap)
        call plt%savefig(trim(out_flux)//'/flux_surfaces_RZ_phi0.png', pyfile=trim(out_flux)//'/flux_surfaces_overlay.py')

        rmin_g = minval(r_surf(:, nsurf_plot))
        rmax_g = maxval(r_surf(:, nsurf_plot))
        zmin_g = minval(z_surf(:, nsurf_plot))
        zmax_g = maxval(z_surf(:, nsurf_plot))

        rmin_g = rmin_g - 0.1_dp*(rmax_g - rmin_g)
        rmax_g = rmax_g + 0.1_dp*(rmax_g - rmin_g)
        zmin_g = zmin_g - 0.1_dp*(zmax_g - zmin_g)
        zmax_g = zmax_g + 0.1_dp*(zmax_g - zmin_g)

        dr = (rmax_g - rmin_g)/real(nr - 1, dp)
        dz = (zmax_g - zmin_g)/real(nz - 1, dp)
        do iR = 1, nr
            rgrid(iR) = rmin_g + real(iR - 1, dp)*dr
        end do
        do iZ = 1, nz
            zgrid(iZ) = zmin_g + real(iZ - 1, dp)*dz
        end do

        do iZ = 1, nz
            do iR = 1, nr
                x_cyl = [rgrid(iR), 0.0_dp, zgrid(iZ)]
                call cyl_to_testfield(x_cyl, x_geo, inside)
                if (inside) then
                    call eval_testfield_cyl(x_geo(1), x_geo(2), 0.0_dp, Br, Bphi, Bz, Bmod)
                    br_map(iR, iZ) = Br
                    bphi_map(iR, iZ) = Bphi
                    bz_map(iR, iZ) = Bz
                    bmod_map(iR, iZ) = Bmod
                else
                    br_map(iR, iZ) = 0.0_dp
                    bphi_map(iR, iZ) = 0.0_dp
                    bz_map(iR, iZ) = 0.0_dp
                    bmod_map(iR, iZ) = 0.0_dp
                end if
            end do
        end do

        call plt%initialize(grid=.true., xlabel='R', ylabel='Z', &
            title='TEST tokamak |B|(R,Z) at phi=0', legend=.false., figsize=[10, 8])
        call plt%add_contour(rgrid, zgrid, bmod_map, linestyle='-', colorbar=.false.)
        call plt%savefig(trim(out_field)//'/Bmod_RZ_phi0.png', pyfile=trim(out_field)//'/Bmod_RZ_phi0.py')

        call plt%initialize(grid=.true., xlabel='R', ylabel='Z', &
            title='TEST tokamak Br(R,Z) at phi=0', legend=.false., figsize=[10, 8])
        call plt%add_contour(rgrid, zgrid, br_map, linestyle='-', colorbar=.false.)
        call plt%savefig(trim(out_field)//'/Br_RZ_phi0.png', pyfile=trim(out_field)//'/Br_RZ_phi0.py')

        call plt%initialize(grid=.true., xlabel='R', ylabel='Z', &
            title='TEST tokamak Bphi(R,Z) at phi=0', legend=.false., figsize=[10, 8])
        call plt%add_contour(rgrid, zgrid, bphi_map, linestyle='-', colorbar=.false.)
        call plt%savefig(trim(out_field)//'/Bphi_RZ_phi0.png', pyfile=trim(out_field)//'/Bphi_RZ_phi0.py')

        call plt%initialize(grid=.true., xlabel='R', ylabel='Z', &
            title='TEST tokamak Bz(R,Z) at phi=0', legend=.false., figsize=[10, 8])
        call plt%add_contour(rgrid, zgrid, bz_map, linestyle='-', colorbar=.false.)
        call plt%savefig(trim(out_field)//'/Bz_RZ_phi0.png', pyfile=trim(out_field)//'/Bz_RZ_phi0.py')
    end subroutine plot_flux_surfaces_and_fields

    subroutine print_artifacts()
        print *, 'ARTIFACT: ', trim(out_orbit)//'/orbit_RZ.png'
        print *, 'ARTIFACT: ', trim(out_orbit)//'/orbit_s_t.png'
        print *, 'ARTIFACT: ', trim(out_orbit)//'/orbit_theta_t.png'
        print *, 'ARTIFACT: ', trim(out_orbit)//'/orbit_phi_t.png'
        print *, 'ARTIFACT: ', trim(out_orbit)//'/orbit_Bmod_t.png'
        print *, 'ARTIFACT: ', trim(out_orbit)//'/orbit_p_t.png'
        print *, 'ARTIFACT: ', trim(out_orbit)//'/orbit_lambda_t.png'
        print *, 'ARTIFACT: ', trim(out_orbit)//'/orbit_mu_t.png'
        print *, 'ARTIFACT: ', trim(out_flux)//'/flux_surfaces_RZ_phi0.png'
        print *, 'ARTIFACT: ', trim(out_field)//'/Bmod_RZ_phi0.png'
        print *, 'ARTIFACT: ', trim(out_field)//'/Br_RZ_phi0.png'
        print *, 'ARTIFACT: ', trim(out_field)//'/Bphi_RZ_phi0.png'
        print *, 'ARTIFACT: ', trim(out_field)//'/Bz_RZ_phi0.png'
    end subroutine print_artifacts

end program test_tokamak_testfield_rk45_orbit_plot
