program test_spectre_fo_gauss2_radial_refinement
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use field_spectre, only: create_spectre_field
    use magfie_sub, only: spectre_field
    use util, only: c, p_mass, e_charge, ev
    use spectre_fo_gauss2, only: gauss2_decode, gauss2_step
    use spectre_fo_hybrid, only: spectre_fo_state_t, spectre_fo_configure, &
        spectre_fo_advance_until_exit, spectre_fo_enter, FO_GAUSS2, SPECTRE_FO_OK
    implicit none
    type(spectre_fo_state_t) :: state
    real(dp) :: ro0, used, y(5), initial, defect, maximum
    real(dp) :: expected(6), clock, frozen(6), actual(6)
    logical :: ok
    integer :: ierr, owner, k
    logical :: exited
    character(1024) :: filename

    call get_command_argument(1, filename)
    allocate(spectre_field)
    call create_spectre_field(spectre_field, trim(filename), ierr)
    if (ierr /= 0) error stop 'field load'
    call spectre_fo_configure(1.0e-9_dp, FO_GAUSS2)
    ro0 = sqrt(2.0_dp*3500.0_dp*ev/(4.0_dp*p_mass))* &
        4.0_dp*p_mass*c/(2.0_dp*e_charge)/sqrt(2.0_dp)
    ! Marker 20, immediately before the inner-layer energy kicks. Rebase E
    ! on the captured velocity so this is an independent on-shell FO start.
    state%z_gauss = [1.0025024863675955_dp, 3.2758489329512472_dp, &
        -3.9738344980129323_dp, -13.708571643795706_dp, &
        84.369367121529450_dp, -16.358194458166761_dp]
    state%layer_rho = 1.8214223148575103e-3_dp
    state%h_tau = 0.5_dp*3.7552110327543201e-2_dp
    state%refine_length = 5.8285514075440331e-2_dp
    state%g_min = 1.953125e-3_dp
    call gauss2_decode(state%z_gauss, state%layer_rho, ro0, state%x, state%v)
    initial = 0.5_dp*dot_product(state%v, state%v)
    state%energy = initial
    state%radial_scale = (state%refine_length/ &
        (32.0_dp*state%h_tau*sqrt(2.0_dp*initial)))**2
    state%u = state%z_gauss(1:3)
    state%owner = 2
    state%iface = 1
    state%exit_rho = 10.0_dp
    state%active = .true.
    maximum = 0.0_dp
    do k = 1, 20
        call spectre_fo_advance_until_exit(state, 0.005_dp, ro0, used, y, &
            owner, exited, ierr)
        if (ierr /= SPECTRE_FO_OK .or. exited) error stop 'FO advance'
        defect = abs(0.5_dp*dot_product(state%v, state%v)/initial - 1.0_dp)
        maximum = max(maximum, defect)
    end do
    print *, 'relative FO energy envelope:', maximum
    if (maximum > 1.0e-5_dp) error stop 'radial refinement energy kicks'
    ! Entry chooses all parameters once. The first accepted step must equal
    ! an independent call of the fixed-parameter map, including its stage clock.
    y = [0.95_dp, 3.2888548092268115_dp, 0.73631077818510704_dp, &
        1.0_dp, -0.18109296850536349_dp]
    call spectre_fo_enter(state, y, 1, 1, ro0, ierr)
    if (ierr /= SPECTRE_FO_OK) error stop 'Gauss entry'
    state%exit_rho = 10.0_dp
    frozen = [state%h_tau, state%energy, state%layer_rho, &
        state%refine_length, state%g_min, state%radial_scale]
    call gauss2_step(state%z_gauss, state%layer_rho, ro0, state%energy, &
        state%refine_length, state%g_min, state%h_tau, expected, clock, ok, &
        radial_scale=state%radial_scale)
    if (.not. ok) error stop 'entry reference step'
    call spectre_fo_advance_until_exit(state, 0.5_dp*clock, ro0, used, y, &
        owner, exited, ierr)
    if (ierr /= SPECTRE_FO_OK .or. exited) error stop 'entry dispatch'
    if (maxval(abs(state%z_gauss - expected)) > 1.0e-12_dp) &
        error stop 'hybrid shortened the canonical step'
    if (abs(used - clock) > 1.0e-14_dp) error stop 'stage clock quadrature'
    actual = [state%h_tau, state%energy, state%layer_rho, &
        state%refine_length, state%g_min, state%radial_scale]
    if (maxval(abs(actual - frozen)) > 0.0_dp) error stop 'entry parameters changed'
end program test_spectre_fo_gauss2_radial_refinement
