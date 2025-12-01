program test_full_orbit
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neo_biotsavart, only: coils_t, load_coils_from_file, &
                               compute_vector_potential, compute_magnetic_field
    use orbit_full, only: FullOrbitState, init_full_orbit_state, &
        timestep_full_orbit, convert_full_to_gc, compute_energy
    use util, only: pi, twopi

    implicit none

    type(coils_t), target :: coils
    type(FullOrbitState) :: state
    real(dp) :: A_cart(3), B_cart(3), Bmod
    real(dp) :: x(3)
    real(dp) :: s, theta, phi_out, lambda, v
    real(dp) :: energy_init, energy_current
    integer :: i, ierr
    character(256) :: coils_file

    ! Get coils file from command line or use default
    if (command_argument_count() >= 1) then
        call get_command_argument(1, coils_file)
    else
        coils_file = 'coils.5C'
    endif

    print *, 'Loading coils from: ', trim(coils_file)
    call load_coils_from_file(coils_file, coils)

    ! Convert from SI to CGS
    print *, 'Converting units SI -> CGS'
    coils%x = coils%x * 100.0d0
    coils%y = coils%y * 100.0d0
    coils%z = coils%z * 100.0d0
    coils%current = coils%current * 2997924536.8431d0

    ! Test field evaluation at a point
    x = [1200.0d0, 0.0d0, 0.0d0]  ! Cartesian coordinates in cm

    print *, 'Testing field evaluation at x=', x
    A_cart = compute_vector_potential(coils, x)
    B_cart = compute_magnetic_field(coils, x)
    Bmod = sqrt(B_cart(1)**2 + B_cart(2)**2 + B_cart(3)**2)

    print *, 'A_cart = ', A_cart
    print *, 'B_cart = ', B_cart
    print *, 'Bmod = ', Bmod, ' G'

    ! Check if field is reasonable
    if (Bmod < 1d3 .or. Bmod > 1d6) then
        print *, 'ERROR: Bmod out of expected range (1e3 - 1e6 G)'
        stop 1
    endif

    ! Initialize a full orbit state
    ! Start from a simple VMEC position
    s = 0.5d0
    theta = 0.0d0
    phi = 0.0d0
    lambda = 0.5d0  ! pitch angle parameter
    v = 1.3d9       ! velocity ~ 3.5 MeV alpha

    print *, ''
    print *, 'Initializing full orbit state'
    print *, 's=', s, ' theta=', theta, ' phi=', phi
    print *, 'lambda=', lambda, ' v=', v

    ! orbit_model=1 for Pauli particle, mass=4 AMU (alpha), charge=2e
    call init_full_orbit_state(state, s, theta, phi, lambda, v, &
                               1, 4.0d0, 2.0d0, 1d-10, coils)

    print *, ''
    print *, 'Initial Cartesian state:'
    print *, 'x=', state%z(1), ' y=', state%z(2), ' z=', state%z(3)
    print *, 'p_x=', state%z(4), ' p_y=', state%z(5), ' p_z=', state%z(6)
    print *, 'mu=', state%mu
    print *, 'dt=', state%dt

    energy_init = compute_energy(state)
    print *, 'Initial energy=', energy_init

    ! Try a few timesteps
    print *, ''
    print *, 'Attempting timesteps...'

    do i = 1, 10
        call timestep_full_orbit(state, ierr)
        if (ierr /= 0) then
            print *, 'Timestep ', i, ' failed with ierr=', ierr
            exit
        endif

        energy_current = compute_energy(state)
        call convert_full_to_gc(state, s, theta, phi_out, lambda, v)

        print '(A,I3,A,F12.3,A,F12.3,A,F10.6,A,ES12.4)', &
            'Step ', i, ': R=', state%z(1), ' Z=', state%z(3), &
            ' s=', s, ' dE/E=', (energy_current - energy_init) / energy_init
    enddo

    if (ierr == 0) then
        print *, ''
        print *, 'SUCCESS: Full orbit integration working'
    else
        print *, ''
        print *, 'FAILURE: Full orbit integration failed'
        stop 1
    endif

end program test_full_orbit
