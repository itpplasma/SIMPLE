program test_field_line_midpoint
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use field_line_spectre, only: field_line_psi, field_line_step
    use spectre_reader, only: spectre_data_t, load_spectre, free_spectre
    implicit none

    type(spectre_data_t) :: data
    character(len=1024) :: field
    real(dp) :: state(4), initial(4), h
    integer :: ierr

    call get_command_argument(1, field)
    call load_spectre(trim(field), data, ierr)
    if (ierr /= 0) error stop 'failed to load SPECTRE fixture'

    state(1:3) = [0.2_dp, 0.3_dp, 0.1_dp]
    state(4) = field_line_psi(data, 2, state(1), state(2), state(3))
    initial = state
    h = 3.141592653589793238462643383279_dp/(64.0_dp*real(data%Nfp, dp))

    call field_line_step(data, 2, h, state(1), state(2), state(3), state(4))
    call field_line_step(data, 2, -h, state(1), state(2), state(3), state(4))
    if (maxval(abs(state - initial)) > 2.0e-11_dp) then
        error stop 'implicit midpoint step is not reversible'
    end if

    call free_spectre(data)
    print *, 'SPECTRE field-line midpoint reversibility PASS'
end program test_field_line_midpoint
