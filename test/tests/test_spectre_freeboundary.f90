program test_spectre_freeboundary
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use field_spectre, only: spectre_field_t, create_spectre_field
    use field_can_spectre, only: construct_spectre_coordinates, &
        set_spectre_construction_grid, spectre_mvol, cleanup_spectre, &
        evaluate_spectre_can
    use field_can_base, only: field_can_t
    use magfie_sub, only: set_magfie_spectre_field, init_magfie, SPECTRE

    implicit none

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

    call check_field(field)

    print *, 'freeboundary plasma-vacuum construction and field evaluation PASS'

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

    subroutine check_field(field)
        type(spectre_field_t), intent(in) :: field
        real(dp) :: acov(3), hcov(3), bmod, sqgb(3)

        call field%evaluate([1.75_dp, 0.7_dp, 0.3_dp], acov, hcov, bmod, sqgb)
        if (.not. ieee_is_finite(bmod) .or. bmod <= 0.0_dp) &
            error stop 'plasma field is not finite'
        call field%evaluate([2.25_dp, 0.7_dp, 0.3_dp], acov, hcov, bmod, sqgb)
        if (.not. ieee_is_finite(bmod) .or. bmod <= 0.0_dp) &
            error stop 'vacuum field is not finite'

    end subroutine check_field

end program test_spectre_freeboundary
