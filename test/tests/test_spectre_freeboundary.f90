program test_spectre_freeboundary
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use field_spectre, only: spectre_field_t, create_spectre_field
    use field_can_spectre, only: construct_spectre_coordinates, &
        set_spectre_construction_grid, spectre_mvol, cleanup_spectre, &
        evaluate_spectre_can
    use field_can_base, only: field_can_t
    use interface_crossing, only: apply_crossing, crossing_info_t, &
        CROSSING_LEVEL0, CROSS_CROSSING, CROSS_LOSS
    use magfie_sub, only: set_magfie_spectre_field, init_magfie, SPECTRE, M_TO_CM
    use spectre_fo_hybrid, only: spectre_fo_state_t, &
        spectre_fo_advance_until_exit, SPECTRE_FO_OK, SPECTRE_FO_LOSS

    implicit none

    real(dp), parameter :: RO0_BAR = 0.01_dp
    character(len=1024) :: filename
    type(spectre_field_t) :: field
    integer :: ierr

    if (command_argument_count() /= 1) error stop 'expected free-boundary fixture'
    call get_command_argument(1, filename)
    call create_spectre_field(field, trim(filename), ierr)
    if (ierr /= 0) error stop 'failed to load free-boundary fixture'
    if (field%data%Nvol /= 2 .or. field%data%Mvol /= 3) &
        error stop 'fixture does not contain one external vacuum volume'
    call set_magfie_spectre_field(field)
    call init_magfie(SPECTRE)

    call set_spectre_construction_grid(8, 8, 8, 5, 1, 1.0e-4_dp)
    call construct_spectre_coordinates(field, ierr)
    if (ierr == 0) error stop 'ODE step-limit failure was not catchable'
    call set_spectre_construction_grid(8, 8, 8, 5, 1000000, 1.0e-4_dp)
    call construct_spectre_coordinates(field, ierr)
    if (ierr /= 0) error stop 'canonical construction failed'
    if (spectre_mvol /= 3) error stop 'canonical construction omitted vacuum volume'
    call check_repeated_construction(field)

    call check_field_and_gc(field)
    call check_vacuum_fo(field)
    call check_outer_loss(field)

    print *, 'freeboundary plasma-vacuum GC/FO continuation PASS'

contains

    subroutine check_repeated_construction(field)
        type(spectre_field_t), intent(in) :: field
        type(field_can_t) :: first, current
        integer :: i, failures

        call evaluate_spectre_can(first, 1.5_dp, 0.7_dp, 0.3_dp, 0)
        call cleanup_spectre
        call construct_spectre_coordinates(field, ierr)
        if (ierr /= 0) error stop 'repeated canonical construction failed'
        call evaluate_spectre_can(current, 1.5_dp, 0.7_dp, 0.3_dp, 0)
        if (first%Bmod /= current%Bmod .or. first%Ath /= current%Ath .or. &
            first%Aph /= current%Aph .or. any(first%dAth /= current%dAth) .or. &
            any(first%dAph /= current%dAph)) then
            error stop 'repeated canonical construction changed the chart'
        end if

        failures = 0
        !$omp parallel do private(current) reduction(+:failures)
        do i = 1, 64
            call evaluate_spectre_can(current, 1.5_dp, 0.7_dp, 0.3_dp, 0)
            if (.not. ieee_is_finite(current%Bmod) .or. &
                current%Bmod /= first%Bmod) failures = failures + 1
        end do
        !$omp end parallel do
        if (failures /= 0) error stop 'threaded canonical evaluation changed the field'
        call cleanup_spectre
    end subroutine check_repeated_construction

    subroutine check_field_and_gc(field)
        type(spectre_field_t), intent(in) :: field
        type(crossing_info_t) :: info
        real(dp) :: y(5), y_out(5), acov(3), hcov(3), bmod, sqgb(3)

        call field%evaluate([1.75_dp, 0.7_dp, 0.3_dp], acov, hcov, bmod, sqgb)
        if (.not. ieee_is_finite(bmod) .or. bmod <= 0.0_dp) &
            error stop 'plasma field is not finite'
        call field%evaluate([2.25_dp, 0.7_dp, 0.3_dp], acov, hcov, bmod, sqgb)
        if (.not. ieee_is_finite(bmod) .or. bmod <= 0.0_dp) &
            error stop 'vacuum field is not finite'

        y = [2.0_dp, 0.7_dp, 0.3_dp, 1.0_dp, 1.0_dp]
        call apply_crossing(y, 2, 1, field%data%Mvol, CROSSING_LEVEL0, y_out, info)
        if (info%event_type /= CROSS_CROSSING) &
            error stop 'plasma-to-vacuum transition was treated as loss'
        if (info%vol_from /= 2 .or. info%vol_to /= 3) &
            error stop 'plasma-to-vacuum ownership mismatch'
        if (y_out(4) /= y(4) .or. y_out(5) /= y(5)) &
            error stop 'zero-mu outward transition changed momentum'

        y(5) = -1.0_dp
        call apply_crossing(y, 2, -1, field%data%Mvol, CROSSING_LEVEL0, y_out, info)
        if (info%event_type /= CROSS_CROSSING) &
            error stop 'vacuum-to-plasma transition did not cross'
        if (info%vol_from /= 3 .or. info%vol_to /= 2) &
            error stop 'vacuum-to-plasma ownership mismatch'

        y(5) = 1.0_dp
        call apply_crossing(y, 2, 1, 2, CROSSING_LEVEL0, y_out, info)
        if (info%event_type /= CROSS_LOSS) &
            error stop 'fixed-boundary outer interface stopped being absorbing'
    end subroutine check_field_and_gc

    subroutine check_vacuum_fo(field)
        type(spectre_field_t), intent(in) :: field
        type(spectre_fo_state_t) :: state
        real(dp) :: acov(3), hcov(3), bmod, sqgb(3), e_cov(3, 3)
        real(dp) :: g(3, 3), ginv(3, 3), sqrtg, bvec(3), y_exit(5)
        real(dp) :: dt_used, speed0
        integer :: ierr, owner
        logical :: exited

        state%u = [2.5_dp, 0.7_dp, 0.3_dp]
        state%owner = 3
        state%iface = 2
        state%active = .true.
        call field%coords%evaluate_cart(state%u, state%x)
        state%x = state%x*M_TO_CM
        call field%evaluate(state%u, acov, hcov, bmod, sqgb)
        call field%coords%covariant_basis(state%u, e_cov)
        call field%coords%metric_tensor(state%u, g, ginv, sqrtg)
        bvec = matmul(e_cov, sqgb/sqrtg)
        state%v = bvec/norm2(bvec)
        speed0 = norm2(state%v)
        call spectre_fo_advance_until_exit(state, 1.0e-8_dp, RO0_BAR, dt_used, &
            y_exit, owner, exited, ierr)
        if (ierr /= SPECTRE_FO_OK .or. .not. exited .or. owner /= 3) &
            error stop 'full orbit did not continue in the vacuum volume'
        if (dt_used <= 0.0_dp .or. dt_used > 1.0e-8_dp) &
            error stop 'vacuum full-orbit time accounting'
        if (abs(norm2(state%v) - speed0) > 1.0e-12_dp*speed0) &
            error stop 'vacuum Boris step changed kinetic energy'
    end subroutine check_vacuum_fo

    subroutine check_outer_loss(field)
        type(spectre_field_t), intent(in) :: field
        type(spectre_fo_state_t) :: state
        real(dp) :: e_cov(3, 3), radial(3), y_exit(5), dt_used
        integer :: ierr, owner
        logical :: exited

        state%u = [3.0_dp, 0.7_dp, 0.3_dp]
        state%owner = 3
        state%iface = 2
        state%active = .true.
        call field%coords%evaluate_cart(state%u, state%x)
        state%x = state%x*M_TO_CM
        call field%coords%covariant_basis(state%u, e_cov)
        radial = e_cov(:, 1)/norm2(e_cov(:, 1))
        state%v = radial
        call spectre_fo_advance_until_exit(state, 1.0e-4_dp, RO0_BAR, dt_used, &
            y_exit, owner, exited, ierr)
        if (ierr /= SPECTRE_FO_LOSS .or. state%u(1) /= 3.0_dp) &
            error stop 'vacuum outer boundary did not absorb full orbit'
    end subroutine check_outer_loss
end program test_spectre_freeboundary
