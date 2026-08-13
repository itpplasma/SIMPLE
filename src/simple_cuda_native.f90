module simple_cuda_native
    use, intrinsic :: iso_c_binding, only: c_double, c_float, c_int, &
        c_int64_t, c_size_t
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use field_can_mod, only: field_can_t
    use orbit_symplectic_base, only: symplectic_integrator_t, CASH_KARP, &
        DORMAND_PRINCE
    use boozer_sub, only: boozer_state
    use boozer_rk_tables, only: rk_field_table, rk_profile_table, &
        rk_num_points, rk_x_min, rk_inv_h_step, rk_period, rk_inv_period, &
        rk_tables_ready

    implicit none
    private

    integer(c_int), parameter :: CUDA_CASH_KARP = 1_c_int
    ! The production DOPRI entry point is the tuned fifth-order/PI controller.
    ! Keep the numeric value stable for callers compiled against the old name.
    integer(c_int), parameter :: CUDA_DORMAND_PRINCE_TUNED = 2_c_int
    integer, parameter, public :: CUDA_NATIVE_PROFILE_COUNT = 5

    public :: trace_orbits_cuda_native_landreman
    public :: initialize_cuda_native_particle, evaluate_cuda_native_particle

    interface
        function cuda_native_rk54(method, particle_count, field_table, &
                field_table_count, profile_table, profile_table_count, &
                point_count, x_min, inv_h_step, period, inv_period, torflux, &
                initial_z, mu, ro0, total_duration, block_duration, tolerance, &
                minimum_timestep, loss_decay_rate, maxloss, observed_duration, &
                final_z, loss_time, status, rhs_evaluations, warp_rhs_slots, &
                profile_ms) bind(c, &
                name='simple_cuda_native_rk54') result(ierr)
            import :: c_double, c_float, c_int, c_int64_t, c_size_t
            integer(c_int), value :: method, particle_count
            real(c_float), intent(in) :: field_table(*), profile_table(*)
            integer(c_size_t), value :: field_table_count, profile_table_count
            integer(c_int), intent(in) :: point_count(3)
            real(c_double), intent(in) :: x_min(3), inv_h_step(3), period(3)
            real(c_double), intent(in) :: inv_period(3)
            real(c_double), value :: torflux
            real(c_double), intent(in) :: initial_z(*), mu(*), ro0(*)
            real(c_double), value :: total_duration, block_duration, tolerance
            real(c_double), value :: minimum_timestep, loss_decay_rate, maxloss
            real(c_double), intent(out) :: observed_duration
            real(c_double), intent(out) :: final_z(*), loss_time(*)
            integer(c_int), intent(out) :: status(*)
            integer(c_int64_t), intent(out) :: rhs_evaluations(*)
            integer(c_int64_t), intent(out) :: warp_rhs_slots
            real(c_double), intent(out) :: profile_ms(5)
            integer(c_int) :: ierr
        end function cuda_native_rk54
    end interface

contains

    subroutine evaluate_cuda_native_particle(f, s, theta, phi)
        use boozer_rk_tables, only: splint_boozer_rk_table_device

        type(field_can_t), intent(inout) :: f
        real(dp), intent(in) :: s, theta, phi
        real(dp) :: aphi, daphi, btheta, dbtheta, bphi, dbphi
        real(dp) :: bmod, dbmod(3)

        call splint_boozer_rk_table_device(s, theta, phi, aphi, daphi, &
            btheta, dbtheta, bphi, dbphi, bmod, dbmod)
        f%Ath = boozer_state%torflux*s
        f%Aph = aphi
        f%dAth = 0.0_dp
        f%dAth(1) = boozer_state%torflux
        f%dAph = 0.0_dp
        f%dAph(1) = daphi
        f%Bmod = bmod
        f%dBmod = dbmod
        f%hth = btheta/bmod
        f%hph = bphi/bmod
        f%dhth = 0.0_dp
        f%dhph = 0.0_dp
        f%dhth(1) = dbtheta/bmod - btheta*dbmod(1)/bmod**2
        f%dhph(1) = dbphi/bmod - bphi*dbmod(1)/bmod**2
        f%dhth(2:3) = -btheta*dbmod(2:3)/bmod**2
        f%dhph(2:3) = -bphi*dbmod(2:3)/bmod**2
    end subroutine evaluate_cuda_native_particle

    subroutine initialize_cuda_native_particle(si, f, z, timestep, tolerance)
        use parmot_mod, only: ro0
        use util, only: sqrt2

        type(symplectic_integrator_t), intent(inout) :: si
        type(field_can_t), intent(inout) :: f
        real(dp), intent(in) :: z(5), timestep, tolerance

        call evaluate_cuda_native_particle(f, z(1), z(2), z(3))
        f%mu = z(4)**2*(1.0_dp - z(5)**2)/f%Bmod
        f%ro0 = ro0/sqrt2
        f%vpar = z(4)*z(5)*sqrt2
        si%z(1:3) = z(1:3)
        si%z(4) = f%vpar*f%hph + f%Aph/f%ro0
        si%pabs = z(4)
        si%dt = timestep/sqrt2
        si%rtol = tolerance
        si%atol = tolerance
        si%ntau = 1
    end subroutine initialize_cuda_native_particle

    subroutine trace_orbits_cuda_native_landreman(si, f, npart, method, &
            total_duration, block_duration, time_scale, loss_tau, maxloss, &
            hmin, loss_step, loss_time, zend, nfev, observed_duration, &
            energy_loss_fraction, warp_nfev_slots, profile_ms)
        type(symplectic_integrator_t), intent(in) :: si(npart)
        type(field_can_t), intent(in) :: f(npart)
        integer, intent(in) :: npart, method
        real(dp), intent(in) :: total_duration, block_duration, time_scale
        real(dp), intent(in) :: loss_tau, maxloss, hmin
        integer, intent(out) :: loss_step(npart), nfev(npart)
        real(dp), intent(out) :: loss_time(npart), zend(4, npart)
        real(dp), intent(out) :: observed_duration, energy_loss_fraction
        integer(int64), intent(out) :: warp_nfev_slots
        real(dp), intent(out) :: profile_ms(CUDA_NATIVE_PROFILE_COUNT)

        real(dp), allocatable :: initial_z(:, :), mu(:), ro0(:)
        integer(c_int), allocatable :: status(:)
        integer(c_int64_t), allocatable :: nfev_c(:)
        integer(c_int) :: cuda_method, ierr, point_count(3)
        integer :: i
        real(dp) :: loss_detection_tolerance

        if (.not. rk_tables_ready .or. .not. allocated(rk_field_table) .or. &
            .not. allocated(rk_profile_table)) &
            error stop 'native CUDA RK requires initialized compact Boozer tables'
        if (method == CASH_KARP) then
            cuda_method = CUDA_CASH_KARP
        else if (method == DORMAND_PRINCE) then
            cuda_method = CUDA_DORMAND_PRINCE_TUNED
        else
            error stop 'native CUDA RK supports only Cash-Karp and Dormand-Prince'
        end if

        allocate (initial_z(4, npart), mu(npart), ro0(npart), status(npart), &
            nfev_c(npart))
        do i = 1, npart
            initial_z(:, i) = si(i)%z
            mu(i) = f(i)%mu
            ro0(i) = f(i)%ro0
        end do
        point_count = int(rk_num_points, c_int)

        ierr = cuda_native_rk54(cuda_method, int(npart, c_int), &
            rk_field_table, size(rk_field_table, kind=c_size_t), &
            rk_profile_table, size(rk_profile_table, kind=c_size_t), &
            point_count, rk_x_min, rk_inv_h_step, rk_period, rk_inv_period, &
            boozer_state%torflux, initial_z, mu, ro0, total_duration, &
            block_duration, si(1)%rtol, hmin, time_scale/loss_tau, maxloss, &
            observed_duration, zend, loss_time, status, nfev_c, &
            warp_nfev_slots, profile_ms)
        if (ierr /= 0_c_int) then
            print '(a,i0)', ' native CUDA RK error = ', ierr
            error stop 'native CUDA RK launch failed'
        end if

        loss_detection_tolerance = &
            16.0_dp*epsilon(1.0_dp)*max(total_duration, 1.0_dp)
        energy_loss_fraction = 0.0_dp
        do i = 1, npart
            if (nfev_c(i) > int(huge(nfev(i)), c_int64_t)) &
                error stop 'native CUDA RK field-evaluation count overflow'
            nfev(i) = int(nfev_c(i))
            if (status(i) /= 0_c_int) then
                loss_step(i) = -int(status(i))
            else if (loss_time(i) < total_duration - loss_detection_tolerance) then
                loss_step(i) = 1
                energy_loss_fraction = energy_loss_fraction + &
                    exp(-loss_time(i)*time_scale/loss_tau)
            else
                loss_step(i) = 0
            end if
        end do
        energy_loss_fraction = energy_loss_fraction/real(npart, dp)
    end subroutine trace_orbits_cuda_native_landreman

end module simple_cuda_native
