program test_radial_potential_config
    ! Only the Boozer canonical field applies Phi(s); other paths must reject a
    ! nonzero radial_electric_potential_slope instead of ignoring it.
    use params, only: radial_electric_potential_slope, integmode, orbit_model, &
        gpu_mode, ORBIT_FULL_ORBIT, validate_radial_electric_potential
    use velo_mod, only: isw_field_type
    use magfie_sub, only: BOOZER, VMEC

    implicit none

    character(16) :: mode

    call get_command_argument(1, mode)
    radial_electric_potential_slope = 0.1d0
    isw_field_type = BOOZER
    integmode = 1
    select case (trim(mode))
    case ('valid')
        continue
    case ('zero_slope_rk')
        radial_electric_potential_slope = 0.0d0
        isw_field_type = VMEC
        integmode = 0
    case ('field')
        isw_field_type = VMEC
    case ('rk')
        integmode = 0
    case ('full_orbit')
        orbit_model = ORBIT_FULL_ORBIT
    case ('gpu')
        gpu_mode = 'production'
    case default
        error stop 'unknown radial potential configuration test mode'
    end select
    call validate_radial_electric_potential
end program test_radial_potential_config
