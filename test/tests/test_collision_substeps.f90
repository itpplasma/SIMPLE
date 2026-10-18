program test_collision_substeps
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use collis_alp, only: stost, coleff, efcolf, velrat, enrat, &
                          efcolf_grid, velrat_grid, enrat_grid, N_S_GRID
    implicit none

    integer, parameter :: npart = 8192
    real(dp), parameter :: resolutions(3) = [0.2_dp, 0.1_dp, 0.05_dp]
    real(dp), parameter :: pi = acos(-1.0_dp)
    integer :: i, j, k, mode, ierr, nseed, thermal_case
    integer, allocatable :: seed(:), seed_after(:)
    real(dp) :: z(5), initial(5), dpp, dhh, drag, duration, h, remaining
    real(dp) :: average1, average2, target1, target2, target4, exact1, exact2
    real(dp) :: variance1, variance2
    real(dp) :: weak_error(3), angle_cos, enr, sigma, energy, fourth, pitch2
    real(dp) :: u(4), velocity(3), eq_energy, eq_fourth, sampling_energy
    real(dp) :: target_energy, target_fourth, variance_component

    call random_seed(size=nseed)
    allocate (seed(nseed), seed_after(nseed))
    efcolf = [1.0_dp, 0.0_dp, 0.0_dp]
    ! Constant-diffusion small-speed limit of the same coefficients: the radial
    ! operator is the norm of a three-dimensional Ornstein-Uhlenbeck process.
    velrat = 0.0_dp
    enrat = 35.0_dp
    call fill_grid()

    initial = [0.5_dp, 0.2_dp, 0.3_dp, 0.015_dp, 0.8_dp]
    call coleff(initial(4), dpp, dhh, drag)
    duration = 0.3_dp/dhh
    exact1 = initial(5)*exp(-2.0_dp*dhh*duration)
    exact2 = legendre2(initial(5))*exp(-6.0_dp*dhh*duration)
    do k = 1, size(resolutions)
        ! The finite-angle kernel has exact conditional Legendre eigenvalues
        ! P_l(cos(angle)); their product independently predicts this algorithm.
        target1 = initial(5)
        target2 = legendre2(initial(5))
        target4 = legendre4(initial(5))
        remaining = duration
        do
            h = min(remaining, resolutions(k)**2/dhh)
            angle_cos = cos(sqrt(4.0_dp*dhh*h))
            target1 = target1*angle_cos
            target2 = target2*legendre2(angle_cos)
            target4 = target4*legendre4(angle_cos)
            remaining = remaining - h
            if (remaining == 0.0_dp) exit
        end do
        weak_error(k) = max(abs(target1 - exact1), abs(target2 - exact2))
        seed = 541
        call random_seed(put=seed)
        average1 = 0.0_dp
        average2 = 0.0_dp
        do i = 1, npart
            z = initial
            call stost(z, duration, 4, ierr, resolutions(k))
            if (ierr /= 0) error stop 'Low-momentum pitch integration failed'
            if (z(4) /= initial(4) .or. abs(z(5)) > 1.0_dp) &
                error stop 'Pitch update left the velocity sphere'
            average1 = average1 + z(5)/npart
            average2 = average2 + legendre2(z(5))/npart
        end do
        variance1 = (1.0_dp + 2.0_dp*target2)/3.0_dp - target1**2
        variance2 = 0.2_dp + 2.0_dp*target2/7.0_dp + &
                    18.0_dp*target4/35.0_dp - target2**2
        write (*, '(a,f6.3,a,4es13.5)') 'pitch fraction=', resolutions(k), &
            ' measured/target P1,P2: ', average1, target1, average2, target2
        if (abs(average1 - target1) > 5.0_dp*sqrt(variance1/npart) .or. &
            abs(average2 - target2) > 5.0_dp*sqrt(variance2/npart)) &
            error stop 'Pitch weak moments disagree with finite-angle kernel'
    end do
    if (weak_error(3) > 0.1_dp*weak_error(1)) &
        error stop 'Finite-angle moments do not converge to heat kernel'

    enr = enrat(1)
    sigma = 1.0_dp/sqrt(2.0_dp*enr)
    eq_energy = 1.5_dp/enr
    eq_fourth = 3.75_dp/enr**2
    ! One relaxation time; the full operator must retain an isotropic
    ! Maxwellian for each numerical resolution, including its low-p tail.
    duration = 1.0_dp/(4.0_dp*dpp*enr)
    do thermal_case = 0, 1
        do k = 1, size(resolutions)
            target_energy = eq_energy
            target_fourth = eq_fourth
            if (thermal_case == 1) then
                ! Exact finite-time radial OU moments from a cold monoenergetic
                ! start. A no-op energy operator cannot pass this transient oracle.
                variance_component = (1.0_dp - exp(-1.0_dp))/(2.0_dp*enr)
                target_energy = eq_energy + (initial(4)**2 - eq_energy)*exp(-1.0_dp)
                target_fourth = initial(4)**4*exp(-2.0_dp) + &
                               10.0_dp*initial(4)**2*exp(-1.0_dp)*variance_component + &
                                15.0_dp*variance_component**2
            end if
            seed = 363
            call random_seed(put=seed)
            energy = 0.0_dp
            fourth = 0.0_dp
            pitch2 = 0.0_dp
            do i = 1, npart
                call random_number(u)
                velocity(1) = sigma*sqrt(-2.0_dp*log(max(u(1), tiny(1.0_dp))))* &
                              cos(2.0_dp*pi*u(2))
                velocity(2) = sigma*sqrt(-2.0_dp*log(max(u(1), tiny(1.0_dp))))* &
                              sin(2.0_dp*pi*u(2))
                velocity(3) = sigma*sqrt(-2.0_dp*log(max(u(3), tiny(1.0_dp))))* &
                              cos(2.0_dp*pi*u(4))
                z = initial
                z(4) = norm2(velocity)
                z(5) = velocity(3)/z(4)
                if (thermal_case == 1) z = initial
                do j = 1, 16
                    call stost(z, duration/16.0_dp, 1, ierr, resolutions(k))
                    if (ierr /= 0) &
                        error stop 'Thermal full collision integration failed'
                end do
                energy = energy + z(4)**2/npart
                fourth = fourth + z(4)**4/npart
                pitch2 = pitch2 + legendre2(z(5))/npart
            end do
            write (*, '(a,i1,a,f6.3,a,5es13.5)') 'radial case=', thermal_case, &
                ' fraction=', resolutions(k), ' <p2>,target,<p4>,target,<P2>: ', &
                energy, target_energy, fourth, target_fourth, pitch2
            sampling_energy = sqrt((target_fourth - target_energy**2)/npart)
            if (abs(energy - target_energy) > 5.0_dp*sampling_energy + &
                resolutions(k)**2*eq_energy) &
                error stop 'Radial weak energy disagrees with OU oracle'
            if (abs(fourth - target_fourth) > &
                5.0_dp*sqrt(45.0_dp/(enr**4*npart)) + &
                2.0_dp*resolutions(k)**2*eq_fourth) &
                error stop 'Radial fourth moment disagrees with OU oracle'
            if (thermal_case == 0 .and. abs(pitch2) > 5.0_dp/sqrt(5.0_dp*npart)) &
                error stop 'Thermal pitch distribution lost isotropy'
        end do
    end do

    ! These are explicit numerical failures, never clipping or random resets.
    z = initial
    call stost(z, duration, 1, ierr, 1.0e-180_dp)
    if (ierr /= 4 .or. any(z /= initial)) &
        error stop 'Substep underflow did not fail and roll back'
    z = initial
    call coleff(z(4), dpp, dhh, drag)
    call stost(z, 2000.0_dp/dhh, 4, ierr)
    if (ierr /= 5 .or. any(z /= initial)) &
        error stop 'Substep budget exhaustion did not fail and roll back'
    efcolf = 0.0_dp
    call fill_grid()
    do mode = 1, 4
        z = initial
        seed = 541
        call random_seed(put=seed)
        call stost(z, duration, mode, ierr)
        call random_seed(get=seed_after)
        if (ierr /= 0 .or. any(z /= initial) .or. any(seed_after /= seed)) &
            error stop 'Vacuum changed state or consumed random numbers'
    end do
    print *, 'PASS: low-p weak moments, heat kernel, Maxwellian and exact vacuum'
contains
    pure real(dp) function legendre2(x)
        real(dp), intent(in) :: x
        legendre2 = (3.0_dp*x*x - 1.0_dp)/2.0_dp
    end function legendre2

    pure real(dp) function legendre4(x)
        real(dp), intent(in) :: x
        legendre4 = (35.0_dp*x**4 - 30.0_dp*x*x + 3.0_dp)/8.0_dp
    end function legendre4

    subroutine fill_grid()
        integer :: node
        do node = 1, N_S_GRID
            efcolf_grid(:, node) = efcolf
            velrat_grid(:, node) = velrat
            enrat_grid(:, node) = enrat
        end do
    end subroutine fill_grid
end program test_collision_substeps
