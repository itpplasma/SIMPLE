program test_api_repeated_sampling
    !> The first simple_sample after simple_init and a repeat on the same surface
    !> must return valid start states. reallocate_arrays resets the starting-line
    !> arrays, which used to be left uninitialised (first call) or zero (repeat).
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use simple_api, only: simple_reset_config, simple_init, simple_sample, &
                          simple_sample_surface
    use params, only: field_input, coord_input, isw_field_type, integmode, &
                      startmode, num_surf, sbeg, deterministic, ntimstep, &
                      ntestpart
    use new_vmec_stuff_mod, only: multharm
    use timing, only: init_timer
    implicit none
    integer, parameter :: n = 8
    real(dp) :: z(5, n)
    character(1024) :: filename

    call get_command_argument(1, filename)
    call init_timer
    call simple_reset_config
    field_input = trim(filename)
    coord_input = trim(filename)
    isw_field_type = 1
    integmode = 0
    startmode = 1
    num_surf = 1
    sbeg(1) = 0.5_dp
    deterministic = .true.
    ntimstep = 2
    ntestpart = n
    ! Sampling/reset semantics do not require a production-resolution field.
    multharm = 1
    call simple_init

    call simple_sample(n, z)
    call check(z, 0.5_dp, 'first sample after init')
    call simple_sample(n, z)
    call check(z, 0.5_dp, 'repeat sample on the same surface')
    call simple_sample_surface(n, 0.7_dp, z)
    call check(z, 0.7_dp, 'new surface')
    call simple_sample_surface(n, 0.7_dp, z)
    call check(z, 0.7_dp, 'repeat on the new surface')

contains

    subroutine check(z, s, label)
        real(dp), intent(in) :: z(:, :), s
        character(*), intent(in) :: label

        if (any(z /= z)) error stop label//': NaN in start states'
        if (any(abs(z(1, :) - s) > 1.0e-12_dp)) error stop label//': wrong s'
        if (any(abs(z(2:3, :)) > 1.0e6_dp)) error stop label//': angle out of range'
        if (all(z(2:3, :) == 0.0_dp)) error stop label//': angles all zero'
        if (any(abs(z(4, :) - 1.0_dp) > 1.0e-12_dp)) error stop label//': wrong v/v0'
        if (any(abs(z(5, :)) > 1.0_dp)) error stop label//': |v_par/v| > 1'
    end subroutine check

end program test_api_repeated_sampling
