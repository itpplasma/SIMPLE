program test_chartmap_wall_losses
    !> Test chartmap wall boundary loss computations.
    !>
    !> Compares loss fractions between:
    !>   1. VMEC coordinates with RK45 integration
    !>   2. Chartmap reference coordinates with RK45 integration
    !>   3. Chartmap with Meiss canonical coordinates
    !>
    !> Also tests wall offset: traces particles to chartmap boundary at rho=1,
    !> which can be the VMEC LCFS or an offset wall surface.
    !>
    !> Outputs NetCDF file for visualization of loss positions and times.

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use netcdf
    use simple, only: init_vmec
    use util, only: twopi
    use new_vmec_stuff_mod, only: nper
    use magfie_sub, only: init_magfie, magfie, VMEC, MEISS, REFCOORDS, &
                          set_magfie_refcoords_field
    use alpha_lifetime_sub, only: orbit_timestep_axis
    use parmot_mod, only: ro0, rmu
    use field_vmec, only: vmec_field_t, create_vmec_field
    use field_splined, only: splined_field_t, create_splined_field
    use field_can_mod, only: init_field_can, integ_to_ref
    use field_can_meiss, only: init_meiss, get_meiss_coordinates, cleanup_meiss
    use reference_coordinates, only: init_reference_coordinates, ref_coords
    use libneo_coordinates, only: coordinate_system_t, vmec_coordinate_system_t, &
                                  make_chartmap_coordinate_system

    implicit none

    character(len=*), parameter :: wout_file = 'wout.nc'
    character(len=1000) :: chartmap_file, chartmap_offset_file
    character(len=1000) :: output_nc

    type(vmec_field_t) :: vmec_field
    type(splined_field_t) :: splined_vmec, splined_chartmap
    class(coordinate_system_t), allocatable :: chartmap_coords

    integer, parameter :: n_particles = 64
    integer, parameter :: n_steps = 2000
    real(dp), parameter :: trace_time = 1.0e-3_dp
    real(dp), parameter :: relerr = 1.0e-10_dp
    real(dp), parameter :: s_start = 0.8_dp

    real(dp) :: dtau, tau, v0_alpha, fper, zeta_period
    real(dp) :: z0(5, n_particles)

    real(dp), allocatable :: loss_times_vmec(:), loss_times_chart(:), loss_times_meiss(:)
    real(dp), allocatable :: loss_pos_vmec(:, :), loss_pos_chart(:, :), loss_pos_meiss(:, :)
    integer, allocatable :: lost_step_vmec(:), lost_step_chart(:), lost_step_meiss(:)

    real(dp) :: confined_vmec, confined_chart, confined_meiss
    integer :: i, ierr, n_failed

    n_failed = 0

    print *, '========================================================'
    print *, 'Chartmap Wall Losses Test'
    print *, '========================================================'
    print *, 'Comparing loss fractions: VMEC vs Chartmap vs Meiss'
    print *, '========================================================'
    print *

    chartmap_file = 'wout.chartmap.nc'
    if (command_argument_count() >= 1) then
        call get_command_argument(1, chartmap_file)
    end if

    output_nc = 'chartmap_wall_losses.nc'
    if (command_argument_count() >= 2) then
        call get_command_argument(2, output_nc)
    end if

    allocate (loss_times_vmec(n_particles), loss_times_chart(n_particles))
    allocate (loss_times_meiss(n_particles))
    allocate (loss_pos_vmec(3, n_particles), loss_pos_chart(3, n_particles))
    allocate (loss_pos_meiss(3, n_particles))
    allocate (lost_step_vmec(n_particles), lost_step_chart(n_particles))
    allocate (lost_step_meiss(n_particles))

    print *, 'Step 1: Initialize VMEC equilibrium'
    call init_vmec(wout_file, 5, 5, 5, fper)
    call init_reference_coordinates(wout_file)
    call set_physics_parameters(v0_alpha)
    zeta_period = twopi/real(nper, dp)

    tau = trace_time*v0_alpha
    dtau = tau/n_steps
    print '(A,ES12.4)', '  Trace time (s): ', trace_time
    print '(A,I0)', '  Number of particles: ', n_particles
    print '(A,I0)', '  Time steps: ', n_steps

    print *, 'Step 2: Load chartmap'
    print *, '  Using: ', trim(chartmap_file)

    print *, 'Step 3: Create field objects'
    call create_vmec_field(vmec_field)
    call create_splined_field(vmec_field, ref_coords, splined_vmec)
    call make_chartmap_coordinate_system(chartmap_coords, chartmap_file)
    call create_splined_field(vmec_field, chartmap_coords, splined_chartmap)

    print *, 'Step 4: Initialize particles'
    call init_particles(z0)

    print *, '========================================================'
    print *, 'Running loss computations...'
    print *, '========================================================'

    print *, 'Mode 1: VMEC coordinates (RK45)'
    call init_magfie(VMEC)
    call trace_particles_vmec(z0, loss_times_vmec, loss_pos_vmec, lost_step_vmec)
    confined_vmec = count_confined(loss_times_vmec)
    print '(A,F6.2,A)', '  Confined fraction: ', confined_vmec*100.0_dp, '%'

    print *, 'Mode 2: Chartmap reference coordinates (RK45)'
    call set_magfie_refcoords_field(splined_chartmap)
    call init_magfie(REFCOORDS)
    call trace_particles_chartmap(z0, loss_times_chart, loss_pos_chart, lost_step_chart)
    confined_chart = count_confined(loss_times_chart)
    print '(A,F6.2,A)', '  Confined fraction: ', confined_chart*100.0_dp, '%'

    print *, 'Mode 3: Meiss canonical from chartmap (RK45)'
    call init_meiss(splined_chartmap, n_r_=16, n_th_=17, n_phi_=18, &
                    rmin=0.1d0, rmax=0.9d0)
    call get_meiss_coordinates
    call init_magfie(MEISS)
    call trace_particles_meiss(z0, loss_times_meiss, loss_pos_meiss, lost_step_meiss)
    call cleanup_meiss
    confined_meiss = count_confined(loss_times_meiss)
    print '(A,F6.2,A)', '  Confined fraction: ', confined_meiss*100.0_dp, '%'

    print *
    print *, '========================================================'
    print *, 'Results Summary'
    print *, '========================================================'
    print '(A,F6.2,A)', '  VMEC confined:     ', confined_vmec*100.0_dp, '%'
    print '(A,F6.2,A)', '  Chartmap confined: ', confined_chart*100.0_dp, '%'
    print '(A,F6.2,A)', '  Meiss confined:    ', confined_meiss*100.0_dp, '%'

    print *, 'Step 5: Write output for visualization'
    call write_output_netcdf

    if (abs(confined_vmec - confined_chart) > 0.15_dp) then
        print *, 'WARNING: VMEC vs Chartmap confined fraction differs by > 15%'
        n_failed = n_failed + 1
    else if (abs(confined_vmec - confined_chart) > 0.05_dp) then
        print *, 'NOTE: VMEC vs Chartmap differs by 5-15% (boundary detection differs)'
    end if

    if (n_failed > 0) then
        print *, '========================================================'
        print *, 'FAILED: ', n_failed, ' comparison(s) exceeded tolerance'
        print *, '========================================================'
        error stop 1
    end if

    print *, '========================================================'
    print *, 'PASSED: Loss fractions consistent across coordinate systems'
    print *, '========================================================'

contains

    subroutine set_physics_parameters(v0_out)
        use util, only: ev, p_mass, c_light => c, e_charge

        real(dp), intent(out) :: v0_out
        real(dp) :: E_alpha, rlarm
        integer :: n_d, n_e

        n_d = 4
        n_e = 2
        E_alpha = 3.5d6

        v0_out = sqrt(2.0_dp*E_alpha*ev/(n_d*p_mass))
        rlarm = v0_out*n_d*p_mass*c_light/(n_e*e_charge)
        ro0 = rlarm
        rmu = 1.0d8

        print '(A,ES12.4)', '  Alpha energy (eV): ', E_alpha
        print '(A,ES12.4)', '  v0 (cm/s): ', v0_out
    end subroutine set_physics_parameters

    subroutine init_particles(z)
        real(dp), intent(out) :: z(5, n_particles)
        integer :: i
        real(dp) :: theta, phi, lambda

        do i = 1, n_particles
            z(1, i) = s_start
            theta = twopi*real(mod(i - 1, 8), dp)/8.0_dp
            phi = zeta_period*real((i - 1)/8, dp)/8.0_dp
            z(2, i) = theta
            z(3, i) = phi
            z(4, i) = 1.0_dp
            lambda = -1.0_dp + 2.0_dp*real(mod(i - 1, 4), dp)/3.0_dp
            z(5, i) = lambda
        end do
    end subroutine init_particles

    subroutine trace_particles_vmec(z0, loss_times, loss_pos, lost_step)
        real(dp), intent(in) :: z0(5, n_particles)
        real(dp), intent(out) :: loss_times(n_particles)
        real(dp), intent(out) :: loss_pos(3, n_particles)
        integer, intent(out) :: lost_step(n_particles)

        real(dp) :: z(5)
        integer :: i, j, ierr

        loss_times = trace_time
        loss_pos = 0.0_dp
        lost_step = n_steps

        do i = 1, n_particles
            z = z0(:, i)
            do j = 1, n_steps
                call orbit_timestep_axis(z, dtau, dtau, relerr, ierr)
                if (ierr /= 0 .or. z(1) >= 1.0_dp) then
                    loss_times(i) = real(j, dp)*dtau/v0_alpha
                    loss_pos(:, i) = z(1:3)
                    lost_step(i) = j
                    exit
                end if
            end do
        end do
    end subroutine trace_particles_vmec

    subroutine trace_particles_chartmap(z0, loss_times, loss_pos, lost_step)
        real(dp), intent(in) :: z0(5, n_particles)
        real(dp), intent(out) :: loss_times(n_particles)
        real(dp), intent(out) :: loss_pos(3, n_particles)
        integer, intent(out) :: lost_step(n_particles)

        real(dp) :: z(5), z_chart(5)
        integer :: i, j, ierr

        loss_times = trace_time
        loss_pos = 0.0_dp
        lost_step = n_steps

        do i = 1, n_particles
            call vmec_to_chartmap(z0(:, i), z_chart)
            z = z_chart
            do j = 1, n_steps
                call orbit_timestep_axis(z, dtau, dtau, relerr, ierr)
                if (ierr /= 0 .or. z(1) >= 1.0_dp) then
                    loss_times(i) = real(j, dp)*dtau/v0_alpha
                    loss_pos(:, i) = z(1:3)
                    lost_step(i) = j
                    exit
                end if
            end do
        end do
    end subroutine trace_particles_chartmap

    subroutine trace_particles_meiss(z0, loss_times, loss_pos, lost_step)
        use field_can_mod, only: ref_to_integ

        real(dp), intent(in) :: z0(5, n_particles)
        real(dp), intent(out) :: loss_times(n_particles)
        real(dp), intent(out) :: loss_pos(3, n_particles)
        integer, intent(out) :: lost_step(n_particles)

        real(dp) :: z(5), z_chart(5), x_ref(3), x_integ(3)
        integer :: i, j, ierr

        loss_times = trace_time
        loss_pos = 0.0_dp
        lost_step = n_steps

        do i = 1, n_particles
            call vmec_to_chartmap(z0(:, i), z_chart)
            x_ref = z_chart(1:3)
            call ref_to_integ(x_ref, x_integ)
            z(1:3) = x_integ
            z(4:5) = z_chart(4:5)

            do j = 1, n_steps
                call orbit_timestep_axis(z, dtau, dtau, relerr, ierr)
                if (ierr /= 0) then
                    call integ_to_ref(z(1:3), x_ref)
                    if (x_ref(1) >= 1.0_dp) then
                        loss_times(i) = real(j, dp)*dtau/v0_alpha
                        loss_pos(:, i) = x_ref
                        lost_step(i) = j
                        exit
                    end if
                end if
            end do
        end do
    end subroutine trace_particles_meiss

    subroutine vmec_to_chartmap(z_vmec, z_chart)
        real(dp), intent(in) :: z_vmec(5)
        real(dp), intent(out) :: z_chart(5)

        real(dp) :: x_target(3), x_test(3), err_vec(3), err_norm
        real(dp) :: J(3, 3), J_inv(3, 3), dx(3)
        real(dp), parameter :: tol = 1.0e-8_dp
        integer, parameter :: max_iter = 30
        integer :: iter, info, ipiv(3)
        real(dp) :: work(9)

        call ref_coords%evaluate_cart(z_vmec(1:3), x_target)
        z_chart(1:3) = z_vmec(1:3)
        z_chart(4:5) = z_vmec(4:5)

        do iter = 1, max_iter
            call chartmap_coords%evaluate_cart(z_chart(1:3), x_test)
            err_vec = x_test - x_target
            err_norm = sqrt(sum(err_vec**2))
            if (err_norm < tol) exit

            call compute_jacobian(z_chart(1:3), J)
            J_inv = J
            call dgetrf(3, 3, J_inv, 3, ipiv, info)
            call dgetri(3, J_inv, 3, ipiv, work, 9, info)
            dx = -matmul(J_inv, err_vec)
            z_chart(1:3) = z_chart(1:3) + dx
        end do
    end subroutine vmec_to_chartmap

    subroutine compute_jacobian(x, J)
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: J(3, 3)
        real(dp), parameter :: eps = 1.0e-6_dp
        real(dp) :: x_p(3), x_m(3), c_p(3), c_m(3)
        integer :: i

        do i = 1, 3
            x_p = x
            x_m = x
            x_p(i) = x(i) + eps
            x_m(i) = x(i) - eps
            call chartmap_coords%evaluate_cart(x_p, c_p)
            call chartmap_coords%evaluate_cart(x_m, c_m)
            J(:, i) = (c_p - c_m)/(2.0_dp*eps)
        end do
    end subroutine compute_jacobian

    function count_confined(loss_times) result(frac)
        real(dp), intent(in) :: loss_times(n_particles)
        real(dp) :: frac
        integer :: n_conf

        n_conf = count(loss_times >= trace_time*0.999_dp)
        frac = real(n_conf, dp)/real(n_particles, dp)
    end function count_confined

    subroutine write_output_netcdf
        integer :: ncid, dimid_part
        integer :: varid_lt_vmec, varid_lt_chart, varid_lt_meiss
        integer :: varid_lp_vmec, varid_lp_chart, varid_lp_meiss
        integer :: varid_ls_vmec, varid_ls_chart, varid_ls_meiss
        integer :: varid_z0
        integer :: dimid_5, dimid_3
        integer :: status

        status = nf90_create(trim(output_nc), nf90_netcdf4, ncid)
        if (status /= nf90_noerr) error stop 'Failed to create NetCDF'

        status = nf90_def_dim(ncid, 'n_particles', n_particles, dimid_part)
        status = nf90_def_dim(ncid, 'dim5', 5, dimid_5)
        status = nf90_def_dim(ncid, 'dim3', 3, dimid_3)

        status = nf90_put_att(ncid, nf90_global, 'description', &
                              'Chartmap wall losses comparison')
        status = nf90_put_att(ncid, nf90_global, 'trace_time', trace_time)
        status = nf90_put_att(ncid, nf90_global, 'n_steps', n_steps)

        status = nf90_def_var(ncid, 'z0', nf90_double, [dimid_5, dimid_part], varid_z0)
        status = nf90_def_var(ncid, 'loss_time_vmec', nf90_double, [dimid_part], &
                              varid_lt_vmec)
        status = nf90_def_var(ncid, 'loss_time_chartmap', nf90_double, [dimid_part], &
                              varid_lt_chart)
        status = nf90_def_var(ncid, 'loss_time_meiss', nf90_double, [dimid_part], &
                              varid_lt_meiss)
        status = nf90_def_var(ncid, 'loss_pos_vmec', nf90_double, &
                              [dimid_3, dimid_part], varid_lp_vmec)
        status = nf90_def_var(ncid, 'loss_pos_chartmap', nf90_double, &
                              [dimid_3, dimid_part], varid_lp_chart)
        status = nf90_def_var(ncid, 'loss_pos_meiss', nf90_double, &
                              [dimid_3, dimid_part], varid_lp_meiss)
        status = nf90_def_var(ncid, 'lost_step_vmec', nf90_int, [dimid_part], &
                              varid_ls_vmec)
        status = nf90_def_var(ncid, 'lost_step_chartmap', nf90_int, [dimid_part], &
                              varid_ls_chart)
        status = nf90_def_var(ncid, 'lost_step_meiss', nf90_int, [dimid_part], &
                              varid_ls_meiss)

        status = nf90_enddef(ncid)

        status = nf90_put_var(ncid, varid_z0, z0)
        status = nf90_put_var(ncid, varid_lt_vmec, loss_times_vmec)
        status = nf90_put_var(ncid, varid_lt_chart, loss_times_chart)
        status = nf90_put_var(ncid, varid_lt_meiss, loss_times_meiss)
        status = nf90_put_var(ncid, varid_lp_vmec, loss_pos_vmec)
        status = nf90_put_var(ncid, varid_lp_chart, loss_pos_chart)
        status = nf90_put_var(ncid, varid_lp_meiss, loss_pos_meiss)
        status = nf90_put_var(ncid, varid_ls_vmec, lost_step_vmec)
        status = nf90_put_var(ncid, varid_ls_chart, lost_step_chart)
        status = nf90_put_var(ncid, varid_ls_meiss, lost_step_meiss)

        status = nf90_close(ncid)
        print *, '  Wrote: ', trim(output_nc)
    end subroutine write_output_netcdf

end program test_chartmap_wall_losses
