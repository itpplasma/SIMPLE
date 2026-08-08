program test_cuda_native_bridge
    use, intrinsic :: iso_fortran_env, only: dp => real64, sp => real32, int64
    use simple_cuda_native, only: trace_orbits_cuda_native_landreman, &
        CUDA_NATIVE_PROFILE_COUNT
    use field_can_mod, only: field_can_t
    use orbit_symplectic_base, only: symplectic_integrator_t, &
        CASH_KARP, DORMAND_PRINCE
    use boozer_sub, only: boozer_state
    use boozer_rk_tables, only: rk_field_table, rk_profile_table, &
        rk_num_points, rk_x_min, rk_inv_h_step, rk_period, rk_inv_period, &
        rk_tables_ready

    implicit none

    integer, parameter :: npart = 3
    type(symplectic_integrator_t) :: si(npart)
    type(field_can_t) :: f(npart)
    integer :: loss_step(npart), nfev(npart), method, i
    real(dp) :: loss_time(npart), zend(4, npart), observed_duration
    real(dp) :: energy_loss_fraction
    real(dp) :: profile_ms(CUDA_NATIVE_PROFILE_COUNT)
    real(dp), parameter :: duration = 0.02_dp
    integer(int64) :: warp_nfev_slots

    allocate (rk_field_table(2*4*4*4), rk_profile_table(6*4))
    rk_field_table = 0.0_sp
    do i = 1, size(rk_field_table)/2
        rk_field_table(2*i - 1) = 1.0_sp
    end do
    rk_profile_table = 0.0_sp
    do i = 1, 4
        rk_profile_table(6*(i - 1) + 5) = 1.0_sp
    end do
    rk_num_points = [4, 4, 4]
    rk_x_min = 0.0_dp
    rk_inv_h_step = 1.0_dp
    rk_period = [0.0_dp, 3.0_dp, 3.0_dp]
    rk_inv_period = [0.0_dp, 1.0_dp/3.0_dp, 1.0_dp/3.0_dp]
    rk_tables_ready = .true.
    boozer_state%torflux = 1.0_dp

    si(1)%z = [0.2_dp, 0.3_dp, 0.4_dp, 0.5_dp]
    si(2)%z = [0.4_dp, -0.2_dp, 1.1_dp, -0.7_dp]
    si(3)%z = [0.6_dp, 0.8_dp, -0.5_dp, 1.25_dp]
    do i = 1, npart
        si(i)%rtol = 1.0e-10_dp
        f(i)%mu = 0.1_dp*real(i, dp)
        f(i)%ro0 = real(i + 1, dp)
    end do

    do method = CASH_KARP, DORMAND_PRINCE
        if (method /= CASH_KARP .and. method /= DORMAND_PRINCE) cycle
        call check_method(method)
    end do
    print *, 'native CUDA Fortran bridge oracle passed'

contains

    subroutine check_method(method)
        integer, intent(in) :: method
        integer :: j
        real(dp) :: initial_z(4, npart)

        initial_z = reshape([ &
            0.2_dp, 0.3_dp, 0.4_dp, 0.5_dp, &
            0.4_dp, -0.2_dp, 1.1_dp, -0.7_dp, &
            0.6_dp, 0.8_dp, -0.5_dp, 1.25_dp], [4, npart])
        do j = 1, npart
            si(j)%z = initial_z(:, j)
        end do
        call trace_orbits_cuda_native_landreman(si, f, npart, method, &
            duration, 0.005_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0e-14_dp, &
            loss_step, loss_time, zend, nfev, observed_duration, &
            energy_loss_fraction, warp_nfev_slots, profile_ms)

        do j = 1, npart
            if (abs(zend(1, j) - initial_z(1, j)) > 1.0e-7_dp .or. &
                abs(zend(2, j) - initial_z(2, j)) > 1.0e-7_dp .or. &
                abs(zend(3, j) - (initial_z(3, j) + &
                initial_z(4, j)*duration)) > 1.0e-7_dp .or. &
                abs(zend(4, j) - initial_z(4, j)) > 1.0e-7_dp) &
                error stop 'native CUDA Fortran bridge trajectory mismatch'
        end do
        if (any(loss_step /= 0) .or. any(abs(loss_time - duration) > 1.0e-14_dp)) &
            error stop 'native CUDA Fortran bridge loss mismatch'
        if (any(nfev <= 0) .or. warp_nfev_slots <= 0_int64) &
            error stop 'native CUDA Fortran bridge work count mismatch'
        if (abs(observed_duration - duration) > 1.0e-14_dp .or. &
            abs(energy_loss_fraction) > epsilon(1.0_dp)) &
            error stop 'native CUDA Fortran bridge accounting mismatch'
        if (profile_ms(3) <= 0.0_dp .or. profile_ms(5) < profile_ms(3)) &
            error stop 'native CUDA Fortran bridge profile mismatch'
    end subroutine check_method

end program test_cuda_native_bridge
