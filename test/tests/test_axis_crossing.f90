program test_axis_crossing
    !> Axis-crossing regression test for the symplectic Euler1 integrator
    !> (issue #370): orbits that pass the magnetic axis must continue on the
    !> opposite ray, (r, theta) -> (|r|, theta + pi), instead of being kicked
    !> to r = 0.01 and random-walked out of the device.
    !>
    !> A set of deeply trapped near-axis starts is traced on the Boozer field
    !> of the QA test equilibrium. Assertions, for every start:
    !>   - the orbit is never lost (the old clamp produced exactly such
    !>     losses),
    !>   - r stays inside [0, 1),
    !>   - the energy drift stays bounded.
    !> At least one start must trigger the axis-crossing branch
    !> (EVT_R_NEGATIVE), otherwise the test would silently stop covering it.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use util, only: twopi
    use simple, only: tracer_t, init_sympl, init_params
    use simple_main, only: init_field
    use params, only: field_input, coord_input
    use new_vmec_stuff_mod, only: rmajor
    use orbit_symplectic, only: orbit_timestep_sympl
    use orbit_symplectic_base, only: symplectic_integrator_t
    use field_can_mod, only: field_can_t, eval_field => evaluate, get_val
    use diag_counters, only: diag_counters_init, diag_counters_total, &
                             EVT_R_NEGATIVE
    use velo_mod, only: isw_field_type
    use magfie_sub, only: BOOZER

    implicit none

    integer, parameter :: n_start = 6, nsteps = 40000, kcheck = 20
    ! Bounded oscillation of the symplectic Euler Hamiltonian is allowed
    ! (it grows toward the axis); secular drift is not.
    real(dp), parameter :: h_osc_tol = 1.0e-3_dp, h_sec_tol = 1.0e-6_dp

    type(tracer_t) :: norb
    type(field_can_t) :: fcheck
    real(dp) :: z5(5)
    real(dp) :: r0(n_start), vpar0(n_start), th0(n_start)
    real(dp) :: rbig, dtau, h0, hdrift, rmin, hnow, hsum1, hsum2, hsec
    real(dp) :: osc_tol, sec_tol, dth
    integer :: ierr, istart, kt, nh1, nh2, iphase, ksteps
    integer(8) :: ncross0
    integer(8) :: ncross
    logical :: failed

    ! Mix of passing potato orbits (|vpar| large, small r: the orbit circles
    ! a drift-shifted center and sweeps through the axis) and deeply trapped
    ! near-axis bananas.
    r0 = [0.001_dp, 0.002_dp, 0.001_dp, 0.005_dp, 0.003_dp, 0.002_dp]
    vpar0 = [-0.30_dp, 0.90_dp, 0.95_dp, 0.05_dp, -0.60_dp, -0.90_dp]
    th0 = [0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp]

    isw_field_type = BOOZER
    field_input = 'wout.nc'
    coord_input = 'wout.nc'
    call init_field(norb, 'wout.nc', 5, 5, 5, 1)
    call init_params(norb, 2, 4, 3.5e6_dp, 256, 1, 1.0e-12_dp)
    call diag_counters_init

    rbig = rmajor*1.0e2_dp

    failed = .false.
    do iphase = 1, 2
    ! Phase 1: fine step (512 per torus), tight tolerances; orbits graze the
    ! axis. Phase 2: coarse step (64 per torus), where the Newton iterate
    ! genuinely overshoots r < 0 and the chart switch must engage.
    if (iphase == 1) then
        dtau = twopi*rbig/512.0_dp
        ksteps = nsteps
        osc_tol = h_osc_tol
        sec_tol = h_sec_tol
    else
        dtau = twopi*rbig/64.0_dp
        ksteps = nsteps/8
        osc_tol = 1.0e-2_dp
        sec_tol = 1.0e-4_dp
    end if
    do istart = 1, n_start
        z5(1) = r0(istart)
        z5(2) = th0(istart)
        z5(3) = 0.1_dp
        z5(4) = 1.0_dp
        z5(5) = vpar0(istart)
        call init_sympl(norb%si, norb%f, z5, dtau, dtau, 1.0e-12_dp, 1)

        call value_of_h(norb%si, h0)
        hdrift = 0.0_dp
        rmin = norb%si%z(1)
        hsum1 = 0.0_dp
        hsum2 = 0.0_dp
        nh1 = 0
        nh2 = 0

        do kt = 1, ksteps
            call orbit_timestep_sympl(norb%si, norb%f, ierr)
            if (ierr /= 0) then
                print *, 'FAIL: start', istart, 'lost at step', kt, &
                    ' z = ', norb%si%z
                failed = .true.
                exit
            end if
            if (norb%si%z(1) < 0.0_dp .or. norb%si%z(1) >= 1.0_dp) then
                print *, 'FAIL: start', istart, 'r out of range at step', kt, &
                    ' r = ', norb%si%z(1)
                failed = .true.
                exit
            end if
            rmin = min(rmin, norb%si%z(1))
            if (mod(kt, kcheck) == 0) then
                call value_of_h(norb%si, hnow)
                hdrift = max(hdrift, abs(hnow - h0)/abs(h0))
                if (kt <= ksteps/2) then
                    hsum1 = hsum1 + hnow
                    nh1 = nh1 + 1
                else
                    hsum2 = hsum2 + hnow
                    nh2 = nh2 + 1
                end if
            end if
        end do

        if (.not. failed) then
            if (hdrift > osc_tol) then
                print *, 'FAIL: start', istart, 'energy excursion', hdrift
                failed = .true.
            end if
            hsec = abs(hsum2/max(nh2, 1) - hsum1/max(nh1, 1))/abs(h0)
            if (hsec > sec_tol) then
                print *, 'FAIL: start', istart, 'secular energy drift', hsec
                failed = .true.
            end if
        end if

        print '(a,i1,a,i2,a,es10.2,a,es10.2,a,es10.2,a,i12)', ' phase ', &
            iphase, ' start ', istart, &
            '  rmin = ', rmin, '  max|dH|/H = ', hdrift, '  sec = ', hsec, &
            '  crossings so far = ', diag_counters_total(EVT_R_NEGATIVE)
    end do
    end do

    ! Phase 3: white-box trigger. On a clean field ds/dt vanishes at the
    ! axis, so a negative Newton iterate is rare; the spurious-loss bug fired
    ! on corrupted near-axis field data (#370). Force the implicit solve
    ! through the axis by biasing f%pth (the driver copies it to pthold) and
    ! require the chart switch to engage and the orbit to survive.
    dtau = twopi*rbig/64.0_dp
    ncross0 = diag_counters_total(EVT_R_NEGATIVE)
    trigger: do istart = 1, n_start
        do kt = 1, 9
            z5 = [2.0e-5_dp, th0(istart), 0.1_dp, 1.0_dp, vpar0(istart)]
            call init_sympl(norb%si, norb%f, z5, dtau, dtau, 1.0e-12_dp, 1)
            ! Bias pth by a few times the radial pth variation across the
            ! start radius, so the implicit solve lands at a small negative
            ! radius (a realistic crossing magnitude).
            ! Signed: linearized solution r = r0*(1 - kt) crosses for kt >= 2.
            norb%f%pth = norb%f%pth &
                - real(kt, dp)*z5(1)*norb%f%dpth(1)
            call orbit_timestep_sympl(norb%si, norb%f, ierr)
            ! The bias is unphysical by construction, so a loss flag on the
            ! poisoned orbit is not itself a failure; only the crossing
            ! semantics are asserted once the event fires.
            if (diag_counters_total(EVT_R_NEGATIVE) > ncross0) then
                ! Discriminate the chart switch from the historic clamp:
                ! the clamp teleported to exactly r = 0.01 with theta
                ! unchanged; the switch keeps r near the axis (the committed
                ! magnitude is the small overshoot) and shifts theta by pi.
                if (norb%si%z(1) < 0.0_dp) then
                    print *, 'FAIL: trigger left r negative, r = ', &
                        norb%si%z(1)
                    failed = .true.
                end if
                if (norb%si%z(1) >= 5.0e-3_dp) then
                    print *, 'FAIL: trigger left the axis region, r = ', &
                        norb%si%z(1), ' (clamp-like teleport)'
                    failed = .true.
                end if
                dth = modulo(norb%si%z(2) - th0(istart), twopi)
                if (abs(dth - 0.5_dp*twopi) > 1.0_dp) then
                    print *, 'FAIL: trigger theta shift ', dth, &
                        ' not close to pi (clamp-like continuation)'
                    failed = .true.
                end if
                exit trigger
            end if
        end do
    end do trigger

    ncross = diag_counters_total(EVT_R_NEGATIVE)
    if (ncross < 1) then
        print *, 'FAIL: no start triggered the axis-crossing branch; ', &
            'the test no longer covers it'
        failed = .true.
    end if

    if (failed) then
        print *, 'TEST FAILED'
        error stop 1
    end if
    print *, 'test_axis_crossing PASSED, crossings = ', ncross

contains

    subroutine value_of_h(si, h)
        !> Energy at the current integrator state, from a fresh field
        !> evaluation (the integrator's own f may hold extrapolated values).
        type(symplectic_integrator_t), intent(in) :: si
        real(dp), intent(out) :: h

        fcheck = norb%f
        call eval_field(fcheck, si%z(1), si%z(2), si%z(3), 0)
        call get_val(fcheck, si%z(4))
        h = fcheck%H
    end subroutine value_of_h

end program test_axis_crossing
