program test_canonical_chart_config
    !> canonical_chart = 'cartesian' is accepted only with integmode 3-7, the
    !> guiding-centre model, a Boozer or test field and gpu_mode = 'off'.
    use params, only: canonical_chart, integmode, orbit_model, gpu_mode, &
                      ORBIT_GC, ORBIT_FULL_ORBIT, validate_canonical_chart
    use velo_mod, only: isw_field_type
    use magfie_sub, only: TEST, CANFLUX, BOOZER, MEISS, ALBERT, SPECTRE

    implicit none

    character(16) :: mode

    call get_command_argument(1, mode)
    canonical_chart = 'cartesian'
    integmode = 3
    isw_field_type = BOOZER
    orbit_model = ORBIT_GC
    gpu_mode = 'off'
    select case (trim(mode))
    case ('polar_euler')
        canonical_chart = 'polar'
        integmode = 1
        isw_field_type = MEISS
    case ('midpoint')
        continue
    case ('gauss4_test')
        integmode = 7
        isw_field_type = TEST
    case ('euler')
        integmode = 1
    case ('rk')
        integmode = 0
    case ('lobatto')
        integmode = 15
    case ('canflux')
        isw_field_type = CANFLUX
    case ('meiss')
        isw_field_type = MEISS
    case ('albert')
        isw_field_type = ALBERT
    case ('spectre')
        isw_field_type = SPECTRE
    case ('full_orbit')
        orbit_model = ORBIT_FULL_ORBIT
    case ('gpu')
        gpu_mode = 'production'
    case ('unknown')
        canonical_chart = 'darboux'
    case default
        error stop 'unknown canonical chart configuration test mode'
    end select
    call validate_canonical_chart
    print '(a)', 'canonical chart configuration accepted'
end program test_canonical_chart_config
