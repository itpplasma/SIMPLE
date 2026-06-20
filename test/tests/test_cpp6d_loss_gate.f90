program test_cpp6d_loss_gate
    ! Multi-particle regression gate for the production CPP6D loss path on the real
    ! QA equilibrium test_data/wout.nc. It guards the field-direction bug that the
    ! lambda-less vmec_field_metric had: trapped 6D orbits drifted MONOTONICALLY
    ! outward (s only increased) and every trapped particle was lost, while the
    ! confined fraction collapsed. The unit tests at the time passed because they
    ! traced too few steps with one mild pitch.
    !
    ! The robust, rho*-independent signature of a correct field is that a trapped
    ! orbit BOUNCES: its s dips below the start AND rises above it within a bounce,
    ! instead of climbing straight to the edge. Combined with energy conservation
    ! and a multi-particle confined count, this catches the regression without
    ! demanding exact GC agreement (which the large-rho* QA case does not give for
    ! a full/Pauli orbit).
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use simple, only: init_sympl, init_cpp, init_params, tracer_t, &
                      orbit_timestep_cpp_canonical
    use simple_main, only: init_field
    use orbit_cpp_canonical, only: cpp_canon_energy
    use params, only: field_input, coord_input, integmode, relerr, dtaumin, &
                      orbit_coord
    use velo_mod, only: isw_field_type
    use magfie_sub, only: BOOZER
    use boozer_coordinates_mod, only: use_B_r, use_del_tp_B
    use boozer_sub, only: get_boozer_coordinates
    implicit none

    integer, parameter :: ans_s = 5, ans_tp = 5, amultharm = 5
    type(tracer_t) :: norb
    integer :: nfail, i
    real(dp) :: lams_trap(3) = [0.0_dp, 0.15_dp, 0.30_dp]
    real(dp) :: z0(5)

    nfail = 0
    isw_field_type = BOOZER
    field_input = 'wout.nc'; coord_input = 'wout.nc'
    integmode = 1; relerr = 1.0d-13
    call init_field(norb, 'wout.nc', ans_s, ans_tp, amultharm, integmode)
    use_B_r = .true.
    use_del_tp_B = .true.
    call get_boozer_coordinates
    call init_params(norb, 2, 4, 3.5e6_dp, 256, 1, 1.0d-13)
    orbit_coord = 1
    dtaumin = norb%dtaumin

    ! The deepest-trapped orbit (lambda=0) must BOUNCE inward: the field-direction
    ! bug pinned its s_min at the start and drove it monotonically to the edge.
    ! Shallower pitches oscillate outward-first (the banana tip is near the start),
    ! so they only get the oscillation + energy checks.
    do i = 1, 3
        z0 = [0.5_dp, 0.5_dp, 0.2_dp, 1.0_dp, lams_trap(i)]
        call trapped_bounces(z0, lams_trap(i) == 0.0_dp, nfail)
    end do

    ! Multi-particle: a small pitch spread must NOT all be lost, and at least one
    ! trapped pitch must survive the short trace (catches the all-trapped-ejection).
    call multi_particle_retention(nfail)

    if (nfail == 0) then
        print *, 'ALL CPP6D LOSS-GATE TESTS PASSED'
    else
        print *, 'CPP6D LOSS-GATE TESTS FAILED: ', nfail
        error stop 1
    end if

contains

    subroutine trapped_bounces(z0, require_inward, nfail)
        real(dp), intent(in) :: z0(5)
        logical, intent(in) :: require_inward
        integer, intent(inout) :: nfail
        type(tracer_t) :: cpp
        real(dp) :: z(5), E0, E, dEmax, smin, smax
        integer :: it, ierr, nstep
        logical :: lost

        nstep = 20000
        z = z0
        call init_sympl(cpp%si, cpp%f, z, dtaumin, dtaumin, relerr, integmode)
        use_B_r = .true.
        use_del_tp_B = .true.
        call init_cpp(cpp%cpp, cpp%f, z, dtaumin)
        E0 = cpp_canon_energy(cpp%cpp); dEmax = 0.0_dp
        smin = z(1); smax = z(1); lost = .false.
        do it = 1, nstep
            call orbit_timestep_cpp_canonical(cpp%cpp, cpp%f, z, ierr)
            if (ierr /= 0) then; lost = .true.; exit; end if
            smin = min(smin, z(1)); smax = max(smax, z(1))
            E = cpp_canon_energy(cpp%cpp); dEmax = max(dEmax, abs(E - E0)/abs(E0))
        end do

     print '(A,F5.2,A,2F8.4,A,ES10.2,A,L2)', '  lam=', z0(5), ' s band [', smin, smax, &
            '] dE/E=', dEmax, ' lost=', lost
        ! The field-direction bug gave smin = s0 (monotonic outward); the deepest
        ! trapped orbit must dip at least 0.01 below the start (the bounce signature).
        if (require_inward) then
            call check('deepest trapped orbit bounces inward (s_min < s0 - 0.01)', &
                       smin < z0(1) - 0.01_dp, nfail)
        end if
        call check('trapped orbit makes a radial excursion (s_max > s0 + 0.005)', &
                   smax > z0(1) + 0.005_dp, nfail)
        call check('CPP energy conserved over trace (dE/E < 1e-3)', &
                   dEmax < 1.0d-3, nfail)
    end subroutine trapped_bounces

    subroutine multi_particle_retention(nfail)
        integer, intent(inout) :: nfail
        integer, parameter :: np = 8
        type(tracer_t) :: cpp
        real(dp) :: z(5), lam
        integer :: ip, it, ierr, nstep, nconf, ntrap_conf
        logical :: lost, trapped

        nstep = 8000
        nconf = 0; ntrap_conf = 0
        do ip = 1, np
            lam = -0.9_dp + (ip - 1)*1.8_dp/(np - 1)    ! pitch spread -0.9..0.9
            trapped = abs(lam) < 0.4_dp
            z = [0.5_dp, 0.5_dp, 0.2_dp, 1.0_dp, lam]
            call init_sympl(cpp%si, cpp%f, z, dtaumin, dtaumin, relerr, integmode)
            use_B_r = .true.
            use_del_tp_B = .true.
            call init_cpp(cpp%cpp, cpp%f, z, dtaumin)
            lost = .false.
            do it = 1, nstep
                call orbit_timestep_cpp_canonical(cpp%cpp, cpp%f, z, ierr)
                if (ierr /= 0) then; lost = .true.; exit; end if
            end do
            if (.not. lost) then
                nconf = nconf + 1
                if (trapped) ntrap_conf = ntrap_conf + 1
            end if
        end do
        print '(A,I2,A,I2,A,I2)', '  multi: confined ', nconf, '/', np, &
            ', trapped-confined ', ntrap_conf
        ! Bug signature: every trapped orbit is ejected. The Boozer CPP path may
        ! still lose passing samples over this trace, so the gate checks trapped
        ! retention directly.
        call check('multi-particle: at least one trapped retained', &
                   ntrap_conf >= 1, nfail)
    end subroutine multi_particle_retention

    subroutine check(name, ok, nfail)
        character(*), intent(in) :: name
        logical, intent(in) :: ok
        integer, intent(inout) :: nfail
        if (ok) then
            print '(A,A)', 'PASS  ', name
        else
            print '(A,A)', 'FAIL  ', name
            nfail = nfail + 1
        end if
    end subroutine check

end program test_cpp6d_loss_gate
