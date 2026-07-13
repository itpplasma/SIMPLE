program test_boundary_event_config
    use, intrinsic :: ieee_arithmetic, only: ieee_positive_inf, ieee_quiet_nan, &
        ieee_value
    use params, only: boundary_event_fraction_tolerance, &
        boundary_event_radial_tolerance, validate_boundary_event_tolerances

    implicit none

    character(16) :: mode

    call get_command_argument(1, mode)
    boundary_event_fraction_tolerance = 1.0d-8
    boundary_event_radial_tolerance = 1.0d-9
    select case (trim(mode))
    case ('valid')
        continue
    case ('fraction_nan')
        boundary_event_fraction_tolerance = ieee_value(0.0d0, ieee_quiet_nan)
    case ('fraction_inf')
        boundary_event_fraction_tolerance = ieee_value(0.0d0, ieee_positive_inf)
    case ('radial_nan')
        boundary_event_radial_tolerance = ieee_value(0.0d0, ieee_quiet_nan)
    case default
        error stop 'unknown boundary event configuration test mode'
    end select
    call validate_boundary_event_tolerances
end program test_boundary_event_config
