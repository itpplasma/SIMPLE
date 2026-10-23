program test_collision_species
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use collis_alp, only: compute_collision_coeffs, loacol_alpha, &
                          init_collision_profiles, efcolf, efcolf_grid
    use simple_main, only: init_collisions
    use params, only: n_d, n_e, v0, particle_energy_eV_effective
    use velo_mod, only: isw_field_type
    use magfie_sub, only: TEST
    implicit none
    real(dp), parameter :: energy = 55000.0_dp, temperature = 1000.0_dp
    real(dp), parameter :: proton_mass = 1.6726e-24_dp, ev = 1.6022e-12_dp
    real(dp) :: alpha(3), coefficients(3), ratios(3), energies(3), speed
    real(dp) :: pitch_rate, drag_rate, pitch_norm, drag_norm
    integer :: mass

    call compute_collision_coeffs(1.0_dp, 12.0_dp, 1.0_dp, 6.0_dp, &
                          1.0e13_dp, 1.0e12_dp, temperature, temperature, temperature, &
                                  energy, speed, alpha, ratios, energies)
    do mass = 1, 4
        call compute_collision_coeffs(1.0_dp, 12.0_dp, 1.0_dp, 6.0_dp, &
                          1.0e13_dp, 1.0e12_dp, temperature, temperature, temperature, &
                  energy, speed, coefficients, ratios, energies, real(mass, dp), 1.0_dp)
        ! Kinetic energy and Maxwellian speeds are independent dimensional oracles.
        if (abs(0.5_dp*mass*proton_mass*speed**2/(energy*ev) - 1) > 1e-13_dp) &
            error stop 'Test-particle kinetic energy is wrong'
        if (maxval(abs(energies - energy/temperature)) > 1e-13_dp) &
            error stop 'Collision Maxwellian energy scale is wrong'
        ! At fixed energy, the electron coefficient scales as Z_test^2/sqrt(A).
        if (abs(coefficients(3)/alpha(3) - 0.5_dp/sqrt(real(mass, dp))) > &
            1e-13_dp) error stop 'Electron collision mass/charge scaling is wrong'
        if (abs(ratios(1) - sqrt(energy/(mass*temperature))) > 1e-13_dp) &
            error stop 'Main-ion relative speed is wrong'
        if (any(coefficients <= 0)) error stop 'Missing main-ion or impurity collisions'
        call loacol_alpha(1.0_dp, 12.0_dp, 1.0_dp, 6.0_dp, &
                          1.0e13_dp, 1.0e12_dp, temperature, temperature, temperature, &
                          energy, speed, pitch_rate, drag_rate, pitch_norm, drag_norm, &
                          real(mass, dp), 1.0_dp)
        if (any(efcolf /= coefficients)) error stop 'Scalar initialization lost species'
        call init_collision_profiles(1.0_dp, 12.0_dp, 1.0_dp, 6.0_dp, &
                                     energy, speed, real(mass, dp), 1.0_dp)
        if (any(efcolf_grid <= 0)) error stop 'Profile initialization lost collisions'
        ! Exercise the production initializer, including its orbit-speed guard.
        n_d = real(mass, dp)
        n_e = 1
        particle_energy_eV_effective = energy
        v0 = speed
        isw_field_type = TEST
        call init_collisions
        if (any(efcolf /= efcolf_grid(:, 1))) &
            error stop 'Scalar and profile species disagree'
    end do
    print *, 'PASS: H, D and variable charge Coulomb collision coefficients'
end program test_collision_species
