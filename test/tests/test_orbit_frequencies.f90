program test_orbit_frequencies
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use orbit_frequencies, only: frequency_options_t, frequency_result_t, &
                                 compute_canonical_frequencies, FREQ_SUCCESS, &
                                 FREQ_INVALID_INPUT, ORBIT_CLASS_TRAPPED, &
                                 ORBIT_CLASS_PASSING
    use params, only: coord_input, field_input
    use simple, only: tracer_t, init_params
    use simple_main, only: init_field

    implicit none

    type(tracer_t) :: tracer
    type(frequency_options_t) :: options
    type(frequency_result_t) :: result
    real(dp) :: trapped(5), passing(5)

    field_input = 'wout.nc'
    coord_input = 'wout.nc'
    call init_field(tracer, 'wout.nc', 5, 5, 3, 3)
    call init_params(tracer, 1, 2, 5.0e3_dp, 1024, 1, 1.0e-10_dp)
    tracer%integmode = 3

    options%n_periods = 3
    options%max_steps = 2000000

    trapped = [0.4_dp, 0.7_dp, 0.1_dp, 1.0_dp, 0.1_dp]
    call compute_canonical_frequencies(tracer, trapped, options, result)
    call require(result%status == FREQ_SUCCESS, 'trapped frequency status')
    call require(result%orbit_class == ORBIT_CLASS_TRAPPED, 'trapped classification')
    call require(result%n_periods == options%n_periods, 'trapped period count')
    call require(result%period > 0.0_dp, 'trapped positive period')
    call require(result%omega_b > 0.0_dp, 'trapped positive omega_b')

    passing = [0.4_dp, 0.7_dp, 0.1_dp, 1.0_dp, 0.9_dp]
    call compute_canonical_frequencies(tracer, passing, options, result)
    call require(result%status == FREQ_SUCCESS, 'passing frequency status')
    call require(result%orbit_class == ORBIT_CLASS_PASSING, 'passing classification')
    call require(result%n_periods == options%n_periods, 'passing period count')
    call require(result%parallel_direction /= 0, 'passing direction')
    call require(result%period > 0.0_dp, 'passing positive period')
    call require(result%omega_b > 0.0_dp, 'passing positive omega_b')

    options%n_periods = 0
    call compute_canonical_frequencies(tracer, passing, options, result)
    call require(result%status == FREQ_INVALID_INPUT, 'invalid options status')

contains

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(*), intent(in) :: message

        if (.not. condition) then
            write (*, '(A)') 'FAILED: '//message
            error stop 1
        end if
    end subroutine require

end program test_orbit_frequencies
