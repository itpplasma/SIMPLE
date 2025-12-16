program test_orbit_refcoords_rk45
    !> Integration test: compare RK45 orbit trajectories using magfie_vmec vs
    !> magfie_refcoords.
    !>
    !> This validates that the assembled magfie quantities from the refcoords-style
    !> interface are consistent enough to reproduce trajectories computed with
    !> the legacy magfie_vmec implementation.
    !>
    !> Test approach:
    !>   1. Initialize a single particle at a mid-radius location
    !>   2. Integrate for N timesteps using both magfie implementations
    !>   3. Compare trajectories at final position
    !>   4. Check mu conservation along each trajectory separately
    !>   5. Write trajectories to NetCDF for plotting
    !>
    !> The test uses the VMEC equilibrium as source for both paths - the difference
    !> is whether we use direct VMEC splines (magfie_vmec) or re-splined field
    !> with analytic derivatives (magfie_refcoords).

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use netcdf
    use simple, only: init_vmec
    use util, only: twopi
    use new_vmec_stuff_mod, only: nper
    use magfie_sub, only: init_magfie, magfie, VMEC, REFCOORDS, &
                          set_magfie_refcoords_field
    use alpha_lifetime_sub, only: orbit_timestep_axis
    use parmot_mod, only: ro0, rmu
    use field_vmec, only: vmec_field_t, create_vmec_field
    use field_splined, only: splined_field_t, create_splined_field
    use reference_coordinates, only: init_reference_coordinates, ref_coords

    implicit none

    type(vmec_field_t) :: vmec_field
    type(splined_field_t) :: splined_from_vmec
    real(dp) :: fper

    integer, parameter :: n_steps = 100
    real(dp), parameter :: dtaumin = 1.0e-3_dp
    real(dp), parameter :: relerr = 1.0e-10_dp

    real(dp), parameter :: tol_pos_rel = 1.0e-4_dp
    real(dp), parameter :: tol_mu_conservation = 1.0e-3_dp

    real(dp) :: z0(5), z_vmec(5), z_refcoords(5)
    real(dp) :: mu0_vmec, mu0_refcoords, mu_vmec_final, mu_refcoords_final
    real(dp) :: bmod_0
    integer :: i, ierr, n_failed
    real(dp) :: dev_s, dev_th, dev_phi
    real(dp) :: mu_drift_vmec, mu_drift_refcoords

    real(dp) :: traj_vmec(5, n_steps+1), traj_refcoords(5, n_steps+1)
    real(dp) :: time_arr(n_steps+1), mu_vmec_arr(n_steps+1), mu_refcoords_arr(n_steps+1)

    n_failed = 0

    print *, '================================'
    print *, 'RK45 Orbit Integration Test: magfie_vmec vs magfie_refcoords'
    print *, '================================'
    print *

    call init_vmec('wout.nc', 5, 5, 5, fper)
    call init_reference_coordinates('wout.nc')

    call set_physics_parameters

    call create_vmec_field(vmec_field)
    call create_splined_field(vmec_field, ref_coords, splined_from_vmec)

    call set_initial_conditions(z0, bmod_0)

    print '(A,3ES14.6)', 'Initial position (s, theta, phi): ', z0(1:3)
    print '(A,ES14.6)', 'Initial momentum p: ', z0(4)
    print '(A,ES14.6)', 'Initial pitch lambda: ', z0(5)
    print '(A,ES14.6)', 'Initial Bmod (vmec): ', bmod_0
    print *

    print *, 'Integrating with magfie_vmec...'
    call init_magfie(VMEC)
    mu0_vmec = compute_mu_at_point(z0)
    z_vmec = z0
    traj_vmec(:, 1) = z_vmec
    time_arr(1) = 0.0_dp
    mu_vmec_arr(1) = mu0_vmec

    do i = 1, n_steps
        call orbit_timestep_axis(z_vmec, dtaumin, dtaumin, relerr, ierr)
        if (ierr /= 0) then
            print *, 'magfie_vmec: particle left domain at step ', i
            exit
        end if
        traj_vmec(:, i+1) = z_vmec
        time_arr(i+1) = i * dtaumin
        mu_vmec_arr(i+1) = compute_mu_at_point(z_vmec)
    end do
    mu_vmec_final = compute_mu_at_point(z_vmec)
    mu_drift_vmec = abs(mu_vmec_final - mu0_vmec)/mu0_vmec
    print '(A,3ES14.6)', 'Final position (vmec): ', z_vmec(1:3)
    print '(A,ES14.6)', 'Initial mu (vmec): ', mu0_vmec
    print '(A,ES14.6)', 'Final mu (vmec): ', mu_vmec_final
    print '(A,ES14.6)', 'mu drift (vmec): ', mu_drift_vmec
    print *

    print *, 'Integrating with magfie_refcoords...'
    call set_magfie_refcoords_field(splined_from_vmec)
    call init_magfie(REFCOORDS)
    mu0_refcoords = compute_mu_at_point(z0)
    z_refcoords = z0
    traj_refcoords(:, 1) = z_refcoords
    mu_refcoords_arr(1) = mu0_refcoords

    do i = 1, n_steps
        call orbit_timestep_axis(z_refcoords, dtaumin, dtaumin, relerr, ierr)
        if (ierr /= 0) then
            print *, 'magfie_refcoords: particle left domain at step ', i
            exit
        end if
        traj_refcoords(:, i+1) = z_refcoords
        mu_refcoords_arr(i+1) = compute_mu_at_point(z_refcoords)
    end do
    mu_refcoords_final = compute_mu_at_point(z_refcoords)
    mu_drift_refcoords = abs(mu_refcoords_final - mu0_refcoords)/mu0_refcoords
    print '(A,3ES14.6)', 'Final position (refcoords): ', z_refcoords(1:3)
    print '(A,ES14.6)', 'Initial mu (refcoords): ', mu0_refcoords
    print '(A,ES14.6)', 'Final mu (refcoords): ', mu_refcoords_final
    print '(A,ES14.6)', 'mu drift (refcoords): ', mu_drift_refcoords
    print *

    dev_s = abs(z_vmec(1) - z_refcoords(1))
    dev_th = abs(z_vmec(2) - z_refcoords(2))
    dev_phi = abs(z_vmec(3) - z_refcoords(3))

    print *, 'Trajectory comparison (vmec vs refcoords):'
    print '(A,ES12.4)', '  Final deviation in s:     ', dev_s
    print '(A,ES12.4)', '  Final deviation in theta: ', dev_th
    print '(A,ES12.4)', '  Final deviation in phi:   ', dev_phi
    print '(A,ES12.4)', '  Relative s deviation:     ', dev_s/max(z0(1), 1d-30)
    print '(A,ES12.4)', '  Relative theta deviation: ', dev_th/twopi
    print '(A,ES12.4)', '  Relative phi deviation:   ', dev_phi/twopi
    print *

    print *, 'Invariant (mu) conservation:'
    print '(A,ES12.4)', '  vmec mu drift:      ', mu_drift_vmec
    print '(A,ES12.4)', '  refcoords mu drift: ', mu_drift_refcoords
    print *

    call write_orbits_netcdf(traj_vmec, traj_refcoords, time_arr, &
                             mu_vmec_arr, mu_refcoords_arr, n_steps+1)
    print *, 'Wrote orbit comparison to orbit_refcoords_comparison.nc'
    print *

    if (dev_s/max(z0(1), 1d-30) > tol_pos_rel) then
        print *, 'FAIL: relative s deviation exceeds tolerance'
        n_failed = n_failed + 1
    end if

    if (dev_th/twopi > tol_pos_rel) then
        print *, 'FAIL: relative theta deviation exceeds tolerance'
        n_failed = n_failed + 1
    end if

    if (dev_phi/twopi > tol_pos_rel) then
        print *, 'FAIL: relative phi deviation exceeds tolerance'
        n_failed = n_failed + 1
    end if

    if (mu_drift_vmec > tol_mu_conservation) then
        print *, 'FAIL: vmec mu drift exceeds tolerance'
        n_failed = n_failed + 1
    end if

    if (mu_drift_refcoords > tol_mu_conservation) then
        print *, 'FAIL: refcoords mu drift exceeds tolerance'
        n_failed = n_failed + 1
    end if

    if (n_failed > 0) then
        print *, '================================'
        print *, 'FAILED: ', n_failed, ' check(s) failed'
        print *, '================================'
        error stop 1
    end if

    print *, '================================'
    print *, 'PASSED: RK45 orbit trajectories match between magfie_vmec and magfie_refcoords'
    print *, '================================'

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


    function compute_mu(z, bmod) result(mu)
        real(dp), intent(in) :: z(5), bmod
        real(dp) :: mu

        real(dp) :: p, lambda, coala

        p = z(4)
        lambda = z(5)
        coala = 1.0_dp - lambda**2
        mu = 0.5_dp*p**2*coala/bmod
    end function compute_mu


    function compute_mu_at_point(z) result(mu)
        real(dp), intent(in) :: z(5)
        real(dp) :: mu

        real(dp) :: bmod, sqrtg, bder(3), hcov(3), hctr(3), hcurl(3)

        call magfie(z(1:3), bmod, sqrtg, bder, hcov, hctr, hcurl)
        mu = compute_mu(z, bmod)
    end function compute_mu_at_point


    subroutine check_nc(status, location)
        integer, intent(in) :: status
        character(len=*), intent(in) :: location

        if (status /= nf90_noerr) then
            print *, 'NetCDF error at ', trim(location), ': ', trim(nf90_strerror(status))
            error stop 'NetCDF operation failed'
        end if
    end subroutine check_nc


    subroutine write_orbits_netcdf(traj_vmec, traj_refcoords, time_arr, &
                                   mu_vmec, mu_refcoords, n_points)
        real(dp), intent(in) :: traj_vmec(5, n_points), traj_refcoords(5, n_points)
        real(dp), intent(in) :: time_arr(n_points)
        real(dp), intent(in) :: mu_vmec(n_points), mu_refcoords(n_points)
        integer, intent(in) :: n_points

        integer :: ncid, dimid_time, status
        integer :: varid_time, varid_s_vmec, varid_th_vmec, varid_phi_vmec
        integer :: varid_s_ref, varid_th_ref, varid_phi_ref
        integer :: varid_mu_vmec, varid_mu_ref
        integer :: varid_p_vmec, varid_lambda_vmec
        integer :: varid_p_ref, varid_lambda_ref

        status = nf90_create('orbit_refcoords_comparison.nc', nf90_netcdf4, ncid)
        call check_nc(status, 'nf90_create')

        status = nf90_def_dim(ncid, 'time', n_points, dimid_time)
        call check_nc(status, 'nf90_def_dim time')

        status = nf90_def_var(ncid, 'time', nf90_double, [dimid_time], varid_time)
        call check_nc(status, 'nf90_def_var time')
        status = nf90_put_att(ncid, varid_time, 'units', 'normalized')
        call check_nc(status, 'nf90_put_att time units')

        status = nf90_def_var(ncid, 's_vmec', nf90_double, [dimid_time], varid_s_vmec)
        call check_nc(status, 'nf90_def_var s_vmec')
        status = nf90_def_var(ncid, 'theta_vmec', nf90_double, [dimid_time], varid_th_vmec)
        call check_nc(status, 'nf90_def_var theta_vmec')
        status = nf90_def_var(ncid, 'phi_vmec', nf90_double, [dimid_time], varid_phi_vmec)
        call check_nc(status, 'nf90_def_var phi_vmec')
        status = nf90_def_var(ncid, 'p_vmec', nf90_double, [dimid_time], varid_p_vmec)
        call check_nc(status, 'nf90_def_var p_vmec')
        status = nf90_def_var(ncid, 'lambda_vmec', nf90_double, [dimid_time], varid_lambda_vmec)
        call check_nc(status, 'nf90_def_var lambda_vmec')
        status = nf90_def_var(ncid, 'mu_vmec', nf90_double, [dimid_time], varid_mu_vmec)
        call check_nc(status, 'nf90_def_var mu_vmec')

        status = nf90_def_var(ncid, 's_refcoords', nf90_double, [dimid_time], varid_s_ref)
        call check_nc(status, 'nf90_def_var s_refcoords')
        status = nf90_def_var(ncid, 'theta_refcoords', nf90_double, [dimid_time], varid_th_ref)
        call check_nc(status, 'nf90_def_var theta_refcoords')
        status = nf90_def_var(ncid, 'phi_refcoords', nf90_double, [dimid_time], varid_phi_ref)
        call check_nc(status, 'nf90_def_var phi_refcoords')
        status = nf90_def_var(ncid, 'p_refcoords', nf90_double, [dimid_time], varid_p_ref)
        call check_nc(status, 'nf90_def_var p_refcoords')
        status = nf90_def_var(ncid, 'lambda_refcoords', nf90_double, [dimid_time], varid_lambda_ref)
        call check_nc(status, 'nf90_def_var lambda_refcoords')
        status = nf90_def_var(ncid, 'mu_refcoords', nf90_double, [dimid_time], varid_mu_ref)
        call check_nc(status, 'nf90_def_var mu_refcoords')

        status = nf90_put_att(ncid, nf90_global, 'description', &
            'RK45 orbit comparison: magfie_vmec vs magfie_refcoords')
        call check_nc(status, 'nf90_put_att description')

        status = nf90_enddef(ncid)
        call check_nc(status, 'nf90_enddef')

        status = nf90_put_var(ncid, varid_time, time_arr)
        call check_nc(status, 'nf90_put_var time')

        status = nf90_put_var(ncid, varid_s_vmec, traj_vmec(1, :))
        call check_nc(status, 'nf90_put_var s_vmec')
        status = nf90_put_var(ncid, varid_th_vmec, traj_vmec(2, :))
        call check_nc(status, 'nf90_put_var theta_vmec')
        status = nf90_put_var(ncid, varid_phi_vmec, traj_vmec(3, :))
        call check_nc(status, 'nf90_put_var phi_vmec')
        status = nf90_put_var(ncid, varid_p_vmec, traj_vmec(4, :))
        call check_nc(status, 'nf90_put_var p_vmec')
        status = nf90_put_var(ncid, varid_lambda_vmec, traj_vmec(5, :))
        call check_nc(status, 'nf90_put_var lambda_vmec')
        status = nf90_put_var(ncid, varid_mu_vmec, mu_vmec)
        call check_nc(status, 'nf90_put_var mu_vmec')

        status = nf90_put_var(ncid, varid_s_ref, traj_refcoords(1, :))
        call check_nc(status, 'nf90_put_var s_refcoords')
        status = nf90_put_var(ncid, varid_th_ref, traj_refcoords(2, :))
        call check_nc(status, 'nf90_put_var theta_refcoords')
        status = nf90_put_var(ncid, varid_phi_ref, traj_refcoords(3, :))
        call check_nc(status, 'nf90_put_var phi_refcoords')
        status = nf90_put_var(ncid, varid_p_ref, traj_refcoords(4, :))
        call check_nc(status, 'nf90_put_var p_refcoords')
        status = nf90_put_var(ncid, varid_lambda_ref, traj_refcoords(5, :))
        call check_nc(status, 'nf90_put_var lambda_refcoords')
        status = nf90_put_var(ncid, varid_mu_ref, mu_refcoords)
        call check_nc(status, 'nf90_put_var mu_refcoords')

        status = nf90_close(ncid)
        call check_nc(status, 'nf90_close')
    end subroutine write_orbits_netcdf

end program test_orbit_refcoords_rk45
