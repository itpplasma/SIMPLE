program test_orbit_chartmap_comparison
    !> Visual comparison test: orbit integration in 4 different coordinate systems.
    !>
    !> This test validates that orbit integration produces identical results in
    !> real-space (Cartesian) coordinates regardless of the coordinate system used:
    !>
    !>   a) VMEC coordinates (direct magfie_vmec)
    !>   b) Meiss canonical coordinates based on VMEC field
    !>   c) Chartmap reference coordinates (generated from VMEC)
    !>   d) Meiss canonical coordinates based on chartmap-splined field
    !>
    !> All trajectories are converted to Cartesian (X, Y, Z) for comparison.
    !> Integration runs for 1e-4 seconds (normalized time).

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
    character(len=*), parameter :: chartmap_file = 'wout.chartmap.nc'

    type(vmec_field_t) :: vmec_field
    type(splined_field_t) :: splined_vmec, splined_chartmap
    type(vmec_coordinate_system_t) :: vmec_coords
    class(coordinate_system_t), allocatable :: chartmap_coords

    integer, parameter :: n_steps = 1000
    real(dp), parameter :: total_time = 1.0e-4_dp
    real(dp) :: dtau
    real(dp), parameter :: relerr = 1.0e-10_dp

    real(dp), parameter :: tol_cartesian = 1.0e-4_dp  ! 1 micrometer tolerance

    real(dp) :: z0(5), z0_chartmap(5), fper, bmod_0
    real(dp) :: z_vmec(5), z_meiss_vmec(5), z_chartmap(5), z_meiss_chartmap(5)

    real(dp), allocatable :: traj_vmec(:, :), traj_meiss_vmec(:, :)
    real(dp), allocatable :: traj_chartmap(:, :), traj_meiss_chartmap(:, :)
    real(dp), allocatable :: cart_vmec(:, :), cart_meiss_vmec(:, :)
    real(dp), allocatable :: cart_chartmap(:, :), cart_meiss_chartmap(:, :)
    real(dp), allocatable :: time_arr(:)

    integer :: i, ierr, n_failed
    real(dp) :: max_dev_meiss_vmec, max_dev_chartmap, max_dev_meiss_chartmap

    n_failed = 0

    print *, '========================================================'
    print *, 'Orbit Integration Comparison: 4 Coordinate Systems'
    print *, '========================================================'
    print *, 'a) VMEC coordinates (magfie_vmec)'
    print *, 'b) Meiss canonical from VMEC'
    print *, 'c) Chartmap reference coordinates'
    print *, 'd) Meiss canonical from chartmap'
    print *, '========================================================'
    print *

    allocate (traj_vmec(5, n_steps + 1), traj_meiss_vmec(5, n_steps + 1))
    allocate (traj_chartmap(5, n_steps + 1), traj_meiss_chartmap(5, n_steps + 1))
    allocate (cart_vmec(3, n_steps + 1), cart_meiss_vmec(3, n_steps + 1))
    allocate (cart_chartmap(3, n_steps + 1), cart_meiss_chartmap(3, n_steps + 1))
    allocate (time_arr(n_steps + 1))

    dtau = total_time/n_steps

    print *, 'Step 1: Initialize VMEC equilibrium'
    call init_vmec(wout_file, 5, 5, 5, fper)
    call init_reference_coordinates(wout_file)
    call set_physics_parameters

    print *, 'Step 2: Load chartmap (pre-generated from VMEC)'
    print *, '  Using: ', chartmap_file

    print *, 'Step 3: Create field objects'
    call create_vmec_field(vmec_field)
    call create_splined_field(vmec_field, ref_coords, splined_vmec)

    call make_chartmap_coordinate_system(chartmap_coords, chartmap_file)
    call create_splined_field(vmec_field, chartmap_coords, splined_chartmap)

    print *, 'Step 4: Set initial conditions'
    call set_initial_conditions(z0, bmod_0)
    print '(A,3ES14.6)', '  VMEC initial (s, theta, phi): ', z0(1:3)
    call find_chartmap_coords_for_vmec_point(z0(1:3), z0_chartmap(1:3))
    z0_chartmap(4:5) = z0(4:5)
    print '(A,3ES14.6)', '  Chartmap initial (same point): ', z0_chartmap(1:3)
    print '(A,ES14.6)', '  Initial p: ', z0(4)
    print '(A,ES14.6)', '  Initial lambda: ', z0(5)
    print *

    print *, '========================================================'
    print *, 'Integrating orbits...'
    print *, '========================================================'
    print *

    print *, 'Path A: VMEC coordinates (magfie_vmec)'
    call init_magfie(VMEC)
    z_vmec = z0
    traj_vmec(:, 1) = z_vmec
    time_arr(1) = 0.0_dp
    do i = 1, n_steps
        call orbit_timestep_axis(z_vmec, dtau, dtau, relerr, ierr)
        if (ierr /= 0) then
            print *, '  Particle left domain at step ', i
            exit
        end if
        traj_vmec(:, i + 1) = z_vmec
        time_arr(i + 1) = i*dtau
    end do
    print '(A,3ES14.6)', '  Final position: ', z_vmec(1:3)

    print *, 'Path B: Meiss canonical from VMEC'
    call init_meiss(vmec_field)
    call get_meiss_coordinates
    call init_magfie(MEISS)
    z_meiss_vmec = z0
    call convert_vmec_to_meiss(z_meiss_vmec)
    traj_meiss_vmec(:, 1) = z_meiss_vmec
    do i = 1, n_steps
        call orbit_timestep_axis(z_meiss_vmec, dtau, dtau, relerr, ierr)
        if (ierr /= 0) then
            print *, '  Particle left domain at step ', i
            exit
        end if
        traj_meiss_vmec(:, i + 1) = z_meiss_vmec
    end do
    print '(A,3ES14.6)', '  Final position (canonical): ', z_meiss_vmec(1:3)
    call cleanup_meiss

    print *, 'Path C: Chartmap reference coordinates'
    call set_magfie_refcoords_field(splined_chartmap)
    call init_magfie(REFCOORDS)
    z_chartmap = z0_chartmap
    traj_chartmap(:, 1) = z_chartmap
    do i = 1, n_steps
        call orbit_timestep_axis(z_chartmap, dtau, dtau, relerr, ierr)
        if (ierr /= 0) then
            print *, '  Particle left domain at step ', i
            exit
        end if
        traj_chartmap(:, i + 1) = z_chartmap
    end do
    print '(A,3ES14.6)', '  Final position: ', z_chartmap(1:3)

    print *, 'Path D: Meiss canonical from chartmap'
    print *, '  Using coarse grid (16x17x18) for numerical stability'
    call init_meiss(splined_chartmap, n_r_=16, n_th_=17, n_phi_=18, &
                    rmin=0.1d0, rmax=0.9d0)
    call get_meiss_coordinates
    call init_magfie(MEISS)
    z_meiss_chartmap = z0_chartmap
    call convert_vmec_to_meiss(z_meiss_chartmap)
    traj_meiss_chartmap(:, 1) = z_meiss_chartmap
    do i = 1, n_steps
        call orbit_timestep_axis(z_meiss_chartmap, dtau, dtau, relerr, ierr)
        if (ierr /= 0) then
            print *, '  Particle left domain at step ', i
            exit
        end if
        traj_meiss_chartmap(:, i + 1) = z_meiss_chartmap
    end do
    print '(A,3ES14.6)', '  Final position (canonical): ', z_meiss_chartmap(1:3)
    call cleanup_meiss
    print *

    print *, '========================================================'
    print *, 'Converting to Cartesian coordinates...'
    print *, '========================================================'
    call convert_trajectories_to_cartesian

    print *, '========================================================'
    print *, 'Comparing trajectories in Cartesian space...'
    print *, '========================================================'
    call compute_deviations(max_dev_meiss_vmec, max_dev_chartmap, &
                            max_dev_meiss_chartmap)

    print '(A,ES12.4)', '  Max deviation (Meiss-VMEC vs VMEC):     ', max_dev_meiss_vmec
    print '(A,ES12.4)', '  Max deviation (Chartmap vs VMEC):       ', max_dev_chartmap
    print '(A,ES12.4)', '  Max deviation (Meiss-Chartmap vs VMEC): ', &
        max_dev_meiss_chartmap
    print *

    print *, '========================================================'
    print *, 'Writing NetCDF output for visualization...'
    print *, '========================================================'
    call write_comparison_netcdf
    print *, '  Wrote: orbit_chartmap_comparison.nc'
    print *

    if (max_dev_meiss_vmec > tol_cartesian) then
        print *, 'FAIL: Meiss-VMEC deviation exceeds tolerance'
        n_failed = n_failed + 1
    end if
    if (max_dev_chartmap > tol_cartesian) then
        print *, 'FAIL: Chartmap deviation exceeds tolerance'
        n_failed = n_failed + 1
    end if
    if (max_dev_meiss_chartmap > tol_cartesian) then
        print *, 'FAIL: Meiss-Chartmap deviation exceeds tolerance'
        n_failed = n_failed + 1
    end if

    if (n_failed > 0) then
        print *, '========================================================'
        print *, 'FAILED: ', n_failed, ' comparison(s) exceeded tolerance'
        print *, '========================================================'
        error stop 1
    end if

    print *, '========================================================'
    print *, 'PASSED: All trajectories match in Cartesian space'
    print *, '========================================================'

contains

    subroutine set_physics_parameters
        use vector_potentail_mod, only: torflux

        ro0 = 1.0d-5
        rmu = 1.0d8
    end subroutine set_physics_parameters

    subroutine set_initial_conditions(z, bmod)
        real(dp), intent(out) :: z(5), bmod

        real(dp) :: sqrtg, bder(3), hcov(3), hctr(3), hcurl(3)

        call init_magfie(VMEC)

        z(1) = 0.25_dp
        z(2) = 0.5_dp
        z(3) = 0.1_dp
        z(4) = 1.0_dp
        z(5) = 0.5_dp

        call magfie(z(1:3), bmod, sqrtg, bder, hcov, hctr, hcurl)
    end subroutine set_initial_conditions

    subroutine find_chartmap_coords_for_vmec_point(z_vmec, z_chart)
        real(dp), intent(in) :: z_vmec(3)
        real(dp), intent(out) :: z_chart(3)

        real(dp) :: x_target(3), x_test(3), err_vec(3), err_norm
        real(dp) :: J(3, 3), J_inv(3, 3), dx(3)
        real(dp), parameter :: tol = 1.0e-10_dp
        integer, parameter :: max_iter = 50
        integer :: iter, info
        integer :: ipiv(3)
        real(dp) :: work(9)
        logical :: converged

        call ref_coords%evaluate_cart(z_vmec, x_target)
        z_chart = z_vmec
        converged = .false.

        do iter = 1, max_iter
            call chartmap_coords%evaluate_cart(z_chart, x_test)
            err_vec = x_test - x_target
            err_norm = sqrt(sum(err_vec**2))

            if (err_norm < tol) then
                converged = .true.
                exit
            end if

            call compute_jacobian_cart(chartmap_coords, z_chart, J)

            J_inv = J
            call dgetrf(3, 3, J_inv, 3, ipiv, info)
            if (info /= 0) then
                print *, 'ERROR: Jacobian LU factorization failed, info =', info
                error stop 'Newton iteration: singular Jacobian'
            end if
            call dgetri(3, J_inv, 3, ipiv, work, 9, info)
            if (info /= 0) then
                print *, 'ERROR: Jacobian inversion failed, info =', info
                error stop 'Newton iteration: matrix inversion failed'
            end if

            dx = -matmul(J_inv, err_vec)
            z_chart = z_chart + dx
        end do

        if (.not. converged) then
            print *, 'ERROR: Newton iteration did not converge'
            print *, '  Final error norm:', err_norm
            error stop 'Chartmap coordinate search failed to converge'
        end if
    end subroutine find_chartmap_coords_for_vmec_point

    subroutine compute_jacobian_cart(coords, x, J)
        class(coordinate_system_t), intent(in) :: coords
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: J(3, 3)

        real(dp), parameter :: eps = 1.0e-6_dp
        real(dp) :: x_plus(3), x_minus(3), cart_plus(3), cart_minus(3)
        integer :: i

        do i = 1, 3
            x_plus = x
            x_minus = x
            x_plus(i) = x(i) + eps
            x_minus(i) = x(i) - eps
            call coords%evaluate_cart(x_plus, cart_plus)
            call coords%evaluate_cart(x_minus, cart_minus)
            J(:, i) = (cart_plus - cart_minus)/(2.0_dp*eps)
        end do
    end subroutine compute_jacobian_cart

    subroutine convert_vmec_to_meiss(z)
        use field_can_mod, only: ref_to_integ
        real(dp), intent(inout) :: z(5)
        real(dp) :: x_ref(3), x_integ(3)

        x_ref = z(1:3)
        call ref_to_integ(x_ref, x_integ)
        z(1:3) = x_integ
    end subroutine convert_vmec_to_meiss

    subroutine convert_meiss_to_vmec(z_meiss, z_vmec)
        use field_can_mod, only: integ_to_ref
        real(dp), intent(in) :: z_meiss(3)
        real(dp), intent(out) :: z_vmec(3)

        call integ_to_ref(z_meiss, z_vmec)
    end subroutine convert_meiss_to_vmec

    subroutine convert_trajectories_to_cartesian
        integer :: j
        real(dp) :: x_ref(3), x_cart(3)

        do j = 1, n_steps + 1
            x_ref = traj_vmec(1:3, j)
            call ref_coords%evaluate_cart(x_ref, x_cart)
            cart_vmec(:, j) = x_cart

            call convert_meiss_to_vmec(traj_meiss_vmec(1:3, j), x_ref)
            call ref_coords%evaluate_cart(x_ref, x_cart)
            cart_meiss_vmec(:, j) = x_cart

            x_ref = traj_chartmap(1:3, j)
            call chartmap_coords%evaluate_cart(x_ref, x_cart)
            cart_chartmap(:, j) = x_cart

            call convert_meiss_to_vmec(traj_meiss_chartmap(1:3, j), x_ref)
            call chartmap_coords%evaluate_cart(x_ref, x_cart)
            cart_meiss_chartmap(:, j) = x_cart
        end do
    end subroutine convert_trajectories_to_cartesian

    subroutine compute_deviations(dev_meiss, dev_chart, dev_meiss_chart)
        real(dp), intent(out) :: dev_meiss, dev_chart, dev_meiss_chart
        integer :: j
        real(dp) :: d

        dev_meiss = 0.0_dp
        dev_chart = 0.0_dp
        dev_meiss_chart = 0.0_dp

        do j = 1, n_steps + 1
            d = sqrt(sum((cart_meiss_vmec(:, j) - cart_vmec(:, j))**2))
            dev_meiss = max(dev_meiss, d)

            d = sqrt(sum((cart_chartmap(:, j) - cart_vmec(:, j))**2))
            dev_chart = max(dev_chart, d)

            d = sqrt(sum((cart_meiss_chartmap(:, j) - cart_vmec(:, j))**2))
            dev_meiss_chart = max(dev_meiss_chart, d)
        end do
    end subroutine compute_deviations

    subroutine write_comparison_netcdf
        integer :: ncid, dimid_time
        integer :: varid_time, varid_x_vmec, varid_y_vmec, varid_z_vmec
        integer :: varid_x_meiss, varid_y_meiss, varid_z_meiss
        integer :: varid_x_chart, varid_y_chart, varid_z_chart
        integer :: varid_x_meiss_chart, varid_y_meiss_chart, varid_z_meiss_chart
        integer :: status

        status = nf90_create('orbit_chartmap_comparison.nc', nf90_netcdf4, ncid)
        call check_nc(status, 'create')

        status = nf90_def_dim(ncid, 'time', n_steps + 1, dimid_time)
        call check_nc(status, 'def_dim time')

        status = nf90_put_att(ncid, nf90_global, 'description', &
                   'Orbit comparison: VMEC vs Meiss-VMEC vs Chartmap vs Meiss-Chartmap')
        call check_nc(status, 'put_att')

        status = nf90_def_var(ncid, 'time', nf90_double, [dimid_time], varid_time)
        call check_nc(status, 'def_var time')

        status = nf90_def_var(ncid, 'x_vmec', nf90_double, [dimid_time], varid_x_vmec)
        call check_nc(status, 'def_var x_vmec')
        status = nf90_def_var(ncid, 'y_vmec', nf90_double, [dimid_time], varid_y_vmec)
        call check_nc(status, 'def_var y_vmec')
        status = nf90_def_var(ncid, 'z_vmec', nf90_double, [dimid_time], varid_z_vmec)
        call check_nc(status, 'def_var z_vmec')

        status = nf90_def_var(ncid, 'x_meiss_vmec', nf90_double, [dimid_time], &
                              varid_x_meiss)
        call check_nc(status, 'def_var x_meiss')
        status = nf90_def_var(ncid, 'y_meiss_vmec', nf90_double, [dimid_time], &
                              varid_y_meiss)
        call check_nc(status, 'def_var y_meiss')
        status = nf90_def_var(ncid, 'z_meiss_vmec', nf90_double, [dimid_time], &
                              varid_z_meiss)
        call check_nc(status, 'def_var z_meiss')

        status = nf90_def_var(ncid, 'x_chartmap', nf90_double, [dimid_time], &
                              varid_x_chart)
        call check_nc(status, 'def_var x_chart')
        status = nf90_def_var(ncid, 'y_chartmap', nf90_double, [dimid_time], &
                              varid_y_chart)
        call check_nc(status, 'def_var y_chart')
        status = nf90_def_var(ncid, 'z_chartmap', nf90_double, [dimid_time], &
                              varid_z_chart)
        call check_nc(status, 'def_var z_chart')

        status = nf90_def_var(ncid, 'x_meiss_chartmap', nf90_double, &
                              [dimid_time], varid_x_meiss_chart)
        call check_nc(status, 'def_var x_meiss_chart')
        status = nf90_def_var(ncid, 'y_meiss_chartmap', nf90_double, &
                              [dimid_time], varid_y_meiss_chart)
        call check_nc(status, 'def_var y_meiss_chart')
        status = nf90_def_var(ncid, 'z_meiss_chartmap', nf90_double, &
                              [dimid_time], varid_z_meiss_chart)
        call check_nc(status, 'def_var z_meiss_chart')

        status = nf90_enddef(ncid)
        call check_nc(status, 'enddef')

        status = nf90_put_var(ncid, varid_time, time_arr)
        call check_nc(status, 'put_var time')

        status = nf90_put_var(ncid, varid_x_vmec, cart_vmec(1, :))
        call check_nc(status, 'put_var x_vmec')
        status = nf90_put_var(ncid, varid_y_vmec, cart_vmec(2, :))
        call check_nc(status, 'put_var y_vmec')
        status = nf90_put_var(ncid, varid_z_vmec, cart_vmec(3, :))
        call check_nc(status, 'put_var z_vmec')

        status = nf90_put_var(ncid, varid_x_meiss, cart_meiss_vmec(1, :))
        call check_nc(status, 'put_var x_meiss')
        status = nf90_put_var(ncid, varid_y_meiss, cart_meiss_vmec(2, :))
        call check_nc(status, 'put_var y_meiss')
        status = nf90_put_var(ncid, varid_z_meiss, cart_meiss_vmec(3, :))
        call check_nc(status, 'put_var z_meiss')

        status = nf90_put_var(ncid, varid_x_chart, cart_chartmap(1, :))
        call check_nc(status, 'put_var x_chart')
        status = nf90_put_var(ncid, varid_y_chart, cart_chartmap(2, :))
        call check_nc(status, 'put_var y_chart')
        status = nf90_put_var(ncid, varid_z_chart, cart_chartmap(3, :))
        call check_nc(status, 'put_var z_chart')

        status = nf90_put_var(ncid, varid_x_meiss_chart, cart_meiss_chartmap(1, :))
        call check_nc(status, 'put_var x_meiss_chart')
        status = nf90_put_var(ncid, varid_y_meiss_chart, cart_meiss_chartmap(2, :))
        call check_nc(status, 'put_var y_meiss_chart')
        status = nf90_put_var(ncid, varid_z_meiss_chart, cart_meiss_chartmap(3, :))
        call check_nc(status, 'put_var z_meiss_chart')

        status = nf90_close(ncid)
        call check_nc(status, 'close')
    end subroutine write_comparison_netcdf

    subroutine check_nc(status, location)
        integer, intent(in) :: status
        character(len=*), intent(in) :: location

        if (status /= nf90_noerr) then
            print *, 'NetCDF error at ', trim(location), ': ', &
                trim(nf90_strerror(status))
            error stop 'NetCDF operation failed'
        end if
    end subroutine check_nc

end program test_orbit_chartmap_comparison
