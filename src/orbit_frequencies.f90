module orbit_frequencies
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use simple, only: tracer_t, init_sympl, tstep
    use util, only: twopi

    implicit none
    private

    integer, parameter, public :: ORBIT_CLASS_UNKNOWN = 0
    integer, parameter, public :: ORBIT_CLASS_TRAPPED = 1
    integer, parameter, public :: ORBIT_CLASS_PASSING = 2

    integer, parameter, public :: FREQ_SUCCESS = 0
    integer, parameter, public :: FREQ_INVALID_INPUT = 1
    integer, parameter, public :: FREQ_ORBIT_LOST = 2
    integer, parameter, public :: FREQ_INTEGRATOR_ERROR = 3
    integer, parameter, public :: FREQ_MAX_STEPS = 4

    type, public :: frequency_options_t
        integer :: n_periods = 1
        integer :: max_steps = 10000000
    end type frequency_options_t

    type, public :: frequency_result_t
        integer :: status = FREQ_INVALID_INPUT
        integer :: orbit_class = ORBIT_CLASS_UNKNOWN
        integer :: parallel_direction = 0
        integer :: n_periods = 0
        integer :: n_steps = 0
        real(dp) :: period = 0.0_dp
        real(dp) :: period_std = 0.0_dp
        real(dp) :: delta_phi = 0.0_dp
        real(dp) :: delta_phi_std = 0.0_dp
        real(dp) :: omega_b = 0.0_dp
        real(dp) :: omega_phi = 0.0_dp
    end type frequency_result_t

    public :: compute_canonical_frequencies

contains

    subroutine compute_canonical_frequencies(tracer, initial_state, options, result)
        type(tracer_t), intent(in) :: tracer
        real(dp), intent(in) :: initial_state(5)
        type(frequency_options_t), intent(in) :: options
        type(frequency_result_t), intent(out) :: result

        type(tracer_t) :: work
        real(dp) :: z(5), z_previous(5)
        real(dp) :: event_time, event_phi, previous_event_time, previous_event_phi
        real(dp) :: target_theta, fraction, step_time
        real(dp), allocatable :: periods(:), displacements(:)
        integer :: ierr, step, n_events, direction
        logical :: have_previous_event, event_found

        result = frequency_result_t()
        if (options%n_periods < 1 .or. options%max_steps < 1) return
        if (tracer%dtaumin <= 0.0_dp .or. tracer%v0 <= 0.0_dp) return

        allocate (periods(options%n_periods), displacements(options%n_periods))
        periods = 0.0_dp
        displacements = 0.0_dp
        work = tracer
        z = initial_state
        step_time = work%dtaumin/work%v0

        if (work%integmode > 0) then
            call init_sympl(work%si, work%f, z, work%dtaumin, work%dtaumin, &
                            work%relerr, work%integmode)
        else
            work%dtau = work%dtaumin
        end if

        have_previous_event = .false.
        n_events = 0
        direction = 0
        target_theta = 0.0_dp

        do step = 1, options%max_steps
            z_previous = z
            if (work%integmode > 0) then
                call tstep(work%si, work%f, z, ierr)
            else
                call tstep(work, z, ierr)
            end if
            result%n_steps = step

            if (ierr /= 0) then
                if (ierr == 1 .or. ierr == 2) then
                    result%status = FREQ_ORBIT_LOST
                else
                    result%status = FREQ_INTEGRATOR_ERROR
                end if
                return
            end if

            event_found = .false.
            fraction = 0.0_dp

            if (result%orbit_class == ORBIT_CLASS_UNKNOWN) then
                if (z_previous(5) < 0.0_dp .and. z(5) >= 0.0_dp) then
                    result%orbit_class = ORBIT_CLASS_TRAPPED
                    result%parallel_direction = 0
                    fraction = zero_crossing_fraction(z_previous(5), z(5))
                    event_found = .true.
                else if (abs(z(2) - initial_state(2)) >= twopi) then
                    result%orbit_class = ORBIT_CLASS_PASSING
                    direction = merge(1, -1, z(2) > initial_state(2))
                    result%parallel_direction = direction
                    target_theta = initial_state(2) + real(direction, dp)*twopi
                    fraction = zero_crossing_fraction(z_previous(2) - target_theta, &
                                                      z(2) - target_theta)
                    event_found = .true.
                    target_theta = target_theta + real(direction, dp)*twopi
                end if
            else if (result%orbit_class == ORBIT_CLASS_TRAPPED) then
                if (z_previous(5) < 0.0_dp .and. z(5) >= 0.0_dp) then
                    fraction = zero_crossing_fraction(z_previous(5), z(5))
                    event_found = .true.
                end if
            else
                if (direction > 0) then
                    if (z_previous(2) < target_theta .and. z(2) >= target_theta) then
                        fraction = zero_crossing_fraction( &
                            z_previous(2) - target_theta, z(2) - target_theta &
                        )
                        event_found = .true.
                    end if
                else
                    if (z_previous(2) > target_theta .and. z(2) <= target_theta) then
                        fraction = zero_crossing_fraction( &
                            z_previous(2) - target_theta, z(2) - target_theta &
                        )
                        event_found = .true.
                    end if
                end if
                if (event_found) target_theta = target_theta + real(direction, dp)*twopi
            end if

            if (.not. event_found) cycle

            event_time = (real(step - 1, dp) + fraction)*step_time
            event_phi = z_previous(3) + fraction*(z(3) - z_previous(3))
            if (have_previous_event) then
                n_events = n_events + 1
                periods(n_events) = event_time - previous_event_time
                displacements(n_events) = event_phi - previous_event_phi
                if (n_events == options%n_periods) exit
            end if
            previous_event_time = event_time
            previous_event_phi = event_phi
            have_previous_event = .true.
        end do

        if (n_events /= options%n_periods) then
            result%status = FREQ_MAX_STEPS
            result%n_periods = n_events
            return
        end if

        result%n_periods = n_events
        result%period = sum(periods)/real(n_events, dp)
        result%delta_phi = sum(displacements)/real(n_events, dp)
        if (n_events > 1) then
            result%period_std = sample_std(periods, result%period)
            result%delta_phi_std = sample_std(displacements, result%delta_phi)
        end if
        result%omega_b = twopi/result%period
        result%omega_phi = sum(displacements)/sum(periods)
        result%status = FREQ_SUCCESS
    end subroutine compute_canonical_frequencies

    pure real(dp) function zero_crossing_fraction(left, right) result(fraction)
        real(dp), intent(in) :: left, right
        real(dp) :: denominator

        denominator = right - left
        if (abs(denominator) <= tiny(denominator)) then
            fraction = 0.5_dp
        else
            fraction = -left/denominator
        end if
        fraction = max(0.0_dp, min(1.0_dp, fraction))
    end function zero_crossing_fraction

    pure real(dp) function sample_std(values, mean_value) result(std_value)
        real(dp), intent(in) :: values(:), mean_value

        std_value = sqrt(sum((values - mean_value)**2)/real(size(values) - 1, dp))
    end function sample_std

end module orbit_frequencies
