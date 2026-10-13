program test_fo_symplectic
    !> Exact symplecticity of full-orbit steps that cross the current sheet.
    !>
    !> Markers are traced until steps cross the sheet; the analytic tangent
    !> map M of each such step must satisfy M^T J M = J to round-off.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use sheet_field, only: field_params_t, region_of
    use fo_sheet, only: sheet_step, tangent_map, step_stats_t, step_record_t, &
        SCHEME_LANDING
    use sheet_setup, only: initial_state
    implicit none

    integer, parameter :: NMARKER = 64, MAX_STEPS = 200000, CHECKS_PER_MARKER = 20
    real(dp), parameter :: H_STEP = 0.1_dp, SPEED = 0.05_dp
    real(dp), parameter :: SYMPLECTIC_TOL = 1.0e-12_dp

    type(field_params_t) :: par
    type(step_stats_t) :: stats
    type(step_record_t) :: rec
    real(dp) :: z(4), z_new(4), m(4, 4), jm(4, 4), defect, worst
    integer :: im, k, region, nchecked, nmulti, ncross

    jm = 0.0_dp
    jm(1, 3) = 1.0_dp
    jm(2, 4) = 1.0_dp
    jm(3, 1) = -1.0_dp
    jm(4, 2) = -1.0_dp
    worst = 0.0_dp
    nchecked = 0
    nmulti = 0
    do im = 1, NMARKER
        call initial_state(par, SPEED, im, z)
        region = region_of(par, z(1:2))
        ncross = 0
        do k = 1, MAX_STEPS
            call sheet_step(par, SCHEME_LANDING, region, z, H_STEP, z_new, &
                stats, rec)
            if (rec%n > 1) then
                call tangent_map(par, rec, m)
                defect = maxval(abs(matmul(transpose(m), matmul(jm, m)) - jm))
                worst = max(worst, defect)
                nchecked = nchecked + 1
                if (rec%n > 2) nmulti = nmulti + 1
                ncross = ncross + 1
            end if
            z = z_new
            if (ncross == CHECKS_PER_MARKER) exit
        end do
    end do

    print '(A,I0,A,I0,A)', 'checked ', nchecked, ' crossing steps (', nmulti, &
        ' with several crossings)'
    print '(A,ES10.2,A,I0)', 'worst max|M^T J M - J| = ', worst, &
        '   solver failures ', stats%failures
    if (nchecked == 0) error stop 'no crossing found'
    if (worst >= SYMPLECTIC_TOL) error stop 'crossing step not symplectic'
    print '(A)', 'PASS'

end program test_fo_symplectic
