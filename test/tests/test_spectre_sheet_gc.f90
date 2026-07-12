program test_spectre_sheet_gc
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use field_spectre, only: spectre_field_t, create_spectre_field
    use magfie_sub, only: set_magfie_spectre_field
    use spectre_sheet_gc, only: sheet_gc_sample_t, sheet_gc_state_t, &
        evaluate_sheet_profile, sheet_gc_rhs, sheet_gc_initialize, &
        sheet_gc_advance, sheet_gc_to_y, SHEET_GC_OK, SHEET_GC_DEGENERATE, &
        SHEET_GC_UNINITIALIZED, SHEET_GC_VPAR_ZERO

    implicit none

    integer, parameter :: NSAMPLE = 9, NPERF = 10000
    real(dp), parameter :: ETAS(NSAMPLE) = &
        [-2.0_dp, -1.5_dp, -1.0_dp, -0.5_dp, 0.0_dp, &
        0.5_dp, 1.0_dp, 1.5_dp, 2.0_dp]
    character(len=1024) :: filename
    type(spectre_field_t) :: field
    type(sheet_gc_sample_t) :: sample
    type(sheet_gc_state_t) :: state
    real(dp) :: rhs(4), y(5), energy_rate, scale, max_residual, t0, t1
    real(dp) :: energy_after, dt_used, p_initial, lambda_initial
    integer :: ierr, i
    logical :: exited

    if (command_argument_count() /= 1) error stop 'expected SPECTRE fixture'
    call sheet_gc_rhs(1, 0.0_dp, 0.7_dp, 0.3_dp, 0.5_dp, 0.02_dp, &
        0.01_dp, rhs, energy_rate, ierr)
    if (ierr /= SHEET_GC_UNINITIALIZED) error stop 'uninitialized field status'
    call get_command_argument(1, filename)
    call create_spectre_field(field, trim(filename), ierr)
    if (ierr /= 0) error stop 'failed to load SPECTRE fixture'
    call set_magfie_spectre_field(field)

    call sheet_gc_rhs(1, 0.0_dp, 0.7_dp, 0.3_dp, 0.0_dp, 0.02_dp, &
        0.01_dp, rhs, energy_rate, ierr)
    if (ierr /= SHEET_GC_VPAR_ZERO) error stop 'zero-vpar degeneracy status'

    call evaluate_sheet_profile(1, -4.0_dp, 0.7_dp, 0.3_dp, sample, ierr)
    if (ierr /= SHEET_GC_OK) error stop 'profile evaluation failed'
    y(1:3) = [1.0_dp, 0.7_dp, 0.3_dp]
    y(4) = sqrt(0.5_dp**2 + 2.0_dp*2.0e-5_dp*sample%bmod)
    y(5) = 0.5_dp/y(4)
    p_initial = y(4)
    lambda_initial = y(5)
    call sheet_gc_initialize(1, 1, y, 2.0e-5_dp, sample%bmod, 0.01_dp, &
        state, ierr)
    if (ierr /= SHEET_GC_OK) error stop 'sheet initialization failed'
    if (state%owner /= 1) error stop 'outward sheet entry owner'
    if (abs(state%p - p_initial) > 5.0_dp*epsilon(1.0_dp)*p_initial) &
        error stop 'sheet p normalization'
    if (abs(state%vpar - p_initial*lambda_initial) > &
        5.0_dp*epsilon(1.0_dp)*abs(state%vpar)) &
        error stop 'sheet vbar normalization'
    call sheet_gc_advance(state, 1.0e-7_dp, 0.01_dp, exited, dt_used, ierr)
    if (ierr /= SHEET_GC_OK) error stop 'sheet advance failed'
    if (dt_used <= 0.0_dp .or. dt_used > 1.0e-7_dp) error stop 'sheet time budget'
    if (exited .neqv. (.not. state%active)) error stop 'sheet exit state'
    call sheet_gc_to_y(state, y)
    if (y(1) /= 1.0_dp) error stop 'sheet coordinate moved off interface'
    call evaluate_sheet_profile(1, state%eta, state%theta, state%zeta, sample, &
        ierr)
    energy_after = 0.5_dp*state%vpar**2 + state%mu*sample%bmod
    if (abs(energy_after - state%energy) > &
        5.0_dp*epsilon(1.0_dp)*state%energy) error stop 'projected sheet energy'

    call evaluate_sheet_profile(1, 4.0_dp, 0.7_dp, 0.3_dp, sample, ierr)
    y(4) = sqrt(0.5_dp**2 + 2.0_dp*2.0e-5_dp*sample%bmod)
    y(5) = -0.5_dp/y(4)
    call sheet_gc_initialize(1, -1, y, 2.0e-5_dp, sample%bmod, 0.01_dp, &
        state, ierr)
    if (ierr /= SHEET_GC_OK) error stop 'inward sheet initialization failed'
    if (state%owner /= 2) error stop 'inward sheet entry owner'

    max_residual = 0.0_dp
    do i = 1, NSAMPLE
        call sheet_gc_rhs(1, ETAS(i), 0.7_dp, 0.3_dp, 0.5_dp, 0.02_dp, &
            0.01_dp, rhs, energy_rate, ierr)
        if (ierr /= SHEET_GC_OK) error stop 'sheet RHS degenerate on fixture'
        if (.not. all(ieee_is_finite(rhs))) error stop 'non-finite sheet RHS'
        scale = max(abs(0.5_dp*rhs(4)), &
            abs(0.02_dp)*sum(abs(rhs(1:3))), tiny(1.0_dp))
        max_residual = max(max_residual, abs(energy_rate)/scale)
    end do
    if (max_residual > 5.0e-13_dp) error stop 'sheet energy residual'

    call cpu_time(t0)
    do i = 1, NPERF
        call sheet_gc_rhs(1, 0.5_dp, 0.7_dp, 0.3_dp, 0.5_dp, 0.02_dp, &
            0.01_dp, rhs, energy_rate, ierr)
    end do
    call cpu_time(t1)
    if (ierr /= SHEET_GC_OK) error stop 'sheet benchmark failed'
    print '(A,ES12.4)', 'sheet_gc: max relative energy residual = ', max_residual
    print '(A,F10.3)', 'sheet_gc: microseconds per RHS = ', &
        1.0e6_dp*(t1 - t0)/real(NPERF, dp)

contains

    elemental logical function ieee_is_finite(x)
        real(dp), intent(in) :: x

        ieee_is_finite = abs(x) <= huge(x)
    end function ieee_is_finite

end program test_spectre_sheet_gc
