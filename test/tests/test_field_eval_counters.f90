program test_field_eval_counters
    ! Guards the field-evaluation counters that the integrator cost benchmark
    ! reads. Without these, "accuracy per field evaluation" cannot be compared
    ! between the explicit Runge-Kutta paths (first derivatives only) and the
    ! implicit symplectic ones (second derivatives for their Jacobians).
    !
    ! The oracle is a hand-computed count, not a recorded number: a known
    ! sequence of eval_field calls at known mode_secders values has an
    ! arithmetically predictable effect on each counter. A test that merely
    ! recorded whatever the counters produced would pass just as happily with
    ! the split wired up backwards.

    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use field_can_base, only: n_field_evaluations, n_field_evaluations_d1, &
                              n_field_evaluations_d2, count_field_evaluation
    use field_can_mod, only: eval_field => evaluate, field_can_t
    use simple, only: tracer_t
    use simple_main, only: init_field
    use params, only: coord_input, field_input

    implicit none

    integer :: nfail

    nfail = 0
    call check_counter_arithmetic(nfail)
    call check_backend_counts(nfail)

    if (nfail > 0) then
        write (error_unit, "(i0,a)") nfail, " test(s) failed"
        stop 1
    end if
    write (*, "(a)") "test_field_eval_counters: all tests passed"

contains

    ! The split must add up to the total, and must land on the right side. Both
    ! directions matter: a swapped split still sums correctly.
    subroutine check_counter_arithmetic(nfail)
        integer, intent(inout) :: nfail

        integer(8) :: t0, d1_0, d2_0
        integer :: i

        t0 = n_field_evaluations
        d1_0 = n_field_evaluations_d1
        d2_0 = n_field_evaluations_d2

        ! 5 first-derivative and 3 second-derivative evaluations.
        do i = 1, 5
            call count_field_evaluation(0)
        end do
        do i = 1, 3
            call count_field_evaluation(2)
        end do
        ! mode_secders = 1 is also a second-derivative evaluation (radial only).
        call count_field_evaluation(1)

        if (n_field_evaluations - t0 /= 9_8) then
            write (error_unit, "(a,i0)") "  total counted ", n_field_evaluations - t0
            nfail = nfail + 1
        end if
        if (n_field_evaluations_d1 - d1_0 /= 5_8) then
            write (error_unit, "(a,i0)") "  d1 counted ", n_field_evaluations_d1 - d1_0
            nfail = nfail + 1
        end if
        if (n_field_evaluations_d2 - d2_0 /= 4_8) then
            write (error_unit, "(a,i0)") "  d2 counted ", n_field_evaluations_d2 - d2_0
            nfail = nfail + 1
        end if
        if (n_field_evaluations - t0 /= &
            (n_field_evaluations_d1 - d1_0) + (n_field_evaluations_d2 - d2_0)) then
            write (error_unit, "(a)") "  split does not sum to the total"
            nfail = nfail + 1
        end if
        write (*, "(a)") "  counter arithmetic: total and split agree"
    end subroutine check_counter_arithmetic

    ! The real backend must route through the same counter, so a known number of
    ! eval_field calls moves the counters by exactly that number.
    subroutine check_backend_counts(nfail)
        integer, intent(inout) :: nfail

        type(tracer_t) :: norb
        type(field_can_t) :: f
        integer(8) :: t0, d1_0, d2_0
        integer :: i
        character(len=256) :: wout
        logical :: found

        wout = 'wout.nc'
        inquire (file=trim(wout), exist=found)
        if (.not. found) then
            wout = 'test/test_data/wout.nc'
            inquire (file=trim(wout), exist=found)
        end if
        if (.not. found) then
            write (*, "(a)") "  backend count: skipped, no wout.nc in working dir"
            return
        end if

        field_input = trim(wout)
        coord_input = trim(wout)
        call init_field(norb, trim(wout), 5, 5, 3, 2)

        t0 = n_field_evaluations
        d1_0 = n_field_evaluations_d1
        d2_0 = n_field_evaluations_d2

        do i = 1, 7
            call eval_field(f, 0.4_dp, 0.3_dp*real(i, dp), 0.2_dp, 0)
        end do
        do i = 1, 4
            call eval_field(f, 0.4_dp, 0.3_dp*real(i, dp), 0.2_dp, 2)
        end do

        if (n_field_evaluations_d1 - d1_0 /= 7_8) then
            write (error_unit, "(a,i0,a)") "  backend d1 counted ", &
                n_field_evaluations_d1 - d1_0, ", expected 7"
            nfail = nfail + 1
        end if
        if (n_field_evaluations_d2 - d2_0 /= 4_8) then
            write (error_unit, "(a,i0,a)") "  backend d2 counted ", &
                n_field_evaluations_d2 - d2_0, ", expected 4"
            nfail = nfail + 1
        end if
        if (n_field_evaluations - t0 /= 11_8) then
            write (error_unit, "(a,i0,a)") "  backend total counted ", &
                n_field_evaluations - t0, ", expected 11"
            nfail = nfail + 1
        end if
        write (*, "(a)") "  backend counts match the hand-computed call sequence"
    end subroutine check_backend_counts

end program test_field_eval_counters
