program test_explicit_integmodes
    ! The explicit high-order integmodes on the 4D canonical chart -- Gauss-Radau
    ! (20), Bulirsch-Stoer (21), Cash-Karp (22), and two-derivative Runge-Kutta
    ! fixed (24) and adaptive (25) -- all integrate the same canonical
    ! Hamiltonian through the same f_ode as the symplectic schemes.
    !
    ! Two oracles, neither of them a recording of this code:
    !
    !   1. The canonical Hamiltonian H is an exact invariant of the flow in a
    !      static field. Every part of its drift is truncation error, so an
    !      adaptive method's drift must shrink when the tolerance is tightened,
    !      and a fixed-step one's must shrink under refinement. This is what
    !      catches a broken error estimate: without control, tightening rtol
    !      changes nothing.
    !
    !   2. The four adaptive methods are independent -- collocation,
    !      extrapolation, an embedded explicit pair, and a two-derivative pair --
    !      so their agreement at tight tolerance is a real cross-check. Four
    !      unrelated methods do not share a bug.
    !
    ! The analytic test field is used, so the test needs no VMEC file.

    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use orbit_symplectic, only: symplectic_integrator_t, orbit_sympl_init, &
        orbit_timestep_sympl
    use orbit_symplectic_base, only: RADAU15, GBS16, CASHKARP45, TDRK24, TDRK24A
    use field_can_mod, only: eval_field => evaluate, field_can_from_name, &
        field_can_t, field_can_init, get_val

    implicit none

    integer, parameter :: NMODE = 4
    integer, parameter :: MODES(NMODE) = [RADAU15, GBS16, CASHKARP45, TDRK24A]
    character(len=16), parameter :: NAMES(NMODE) = &
        [character(len=16) :: 'Gauss-Radau', 'Bulirsch-Stoer', 'Cash-Karp', 'TDRK adaptive']

    real(dp), parameter :: Z0(4) = [0.1_dp, 0.7_dp, 0.1_dp, 0.0_dp]
    real(dp), parameter :: VPAR0 = 0.8_dp
    real(dp), parameter :: DT = 500.0_dp
    integer,  parameter :: NSTEP = 100

    real(dp) :: zloose(4, NMODE), ztight(4, NMODE)
    real(dp) :: hloose(NMODE), htight(NMODE)
    real(dp) :: err(3), rate, spread
    integer  :: nfail, i, j, k
    integer  :: ns(3)

    nfail = 0

    do i = 1, NMODE
        call trace(MODES(i), 1.0e-6_dp, 1, zloose(:, i), hloose(i), nfail, NAMES(i))
        call trace(MODES(i), 1.0e-11_dp, 1, ztight(:, i), htight(i), nfail, NAMES(i))
        write (*, '(a,a16,a,es10.3,a,es10.3)') '  ', NAMES(i), &
            '  |dH/H| at rtol 1e-6 ', hloose(i), '   at 1e-11 ', htight(i)
        ! Gauss-Radau is 15th order, so on this smooth field it is already at
        ! round-off at rtol 1e-6 and has nothing left to gain -- that is the
        ! method working, not the control failing. Everything above round-off
        ! must respond.
        if (.not. (htight(i) < 0.2_dp*hloose(i) .or. hloose(i) < 1.0e-13_dp)) then
            write (error_unit, '(a,a)') '  tightening the tolerance did not ', &
                'reduce the energy drift for '//trim(NAMES(i))
            nfail = nfail + 1
        end if
    end do

    ! Cross-method agreement at tight tolerance, component by component.
    do k = 1, 4
        spread = maxval(ztight(k, :)) - minval(ztight(k, :))
        if (.not. (spread < 1.0e-6_dp*max(1.0_dp, maxval(abs(ztight(k, :)))))) then
            write (error_unit, '(a,i0,a,es10.3)') &
                '  methods disagree on final z(', k, '), spread ', spread
            nfail = nfail + 1
        end if
    end do
    write (*, '(a,es10.3)') '  cross-method spread at rtol 1e-11, max component ', &
        maxval([(maxval(ztight(k, :)) - minval(ztight(k, :)), k = 1, 4)])

    ! Fixed-step TDRK: no tolerance to tighten, so the invariant has to improve
    ! under refinement instead. trace holds the macro-step DT fixed and splits it
    ! into ntau substeps, so raising ntau refines rather than lengthening.
    ns = [16, 32, 64]
    do j = 1, 3
        call trace(TDRK24, 1.0e-6_dp, ns(j), ztight(:, 1), err(j), nfail, 'TDRK fixed')
    end do
    rate = -1.0_dp
    if (err(3) > 0.0_dp .and. err(1) > 0.0_dp) rate = log(err(1)/err(3))/log(4.0_dp)
    write (*, '(a,es10.3,a,es10.3,a,es10.3,a,f6.2)') &
        '  TDRK fixed  |dH/H| at h, h/2, h/4 = ', err(1), ', ', err(2), ', ', &
        err(3), '   observed order ', rate
    if (.not. (rate > 3.0_dp)) then
        write (error_unit, '(a,f6.2)') &
            '  fixed-step TDRK energy drift decays below order 3: ', rate
        nfail = nfail + 1
    end if

    if (nfail > 0) then
        write (error_unit, '(i0,a)') nfail, ' check(s) failed'
        stop 1
    end if
    write (*, '(a)') 'test_explicit_integmodes: all checks passed'

contains

    ! Trace NSTEP macro-steps and return the final state and the largest
    ! relative excursion of H along the way.
    subroutine trace(mode, rtol, ntau, zend, hdrift, nfail, tag)
        integer, intent(in) :: mode, ntau
        real(dp), intent(in) :: rtol
        real(dp), intent(out) :: zend(4)
        real(dp), intent(out) :: hdrift
        integer, intent(inout) :: nfail
        character(*), intent(in) :: tag

        type(symplectic_integrator_t) :: si
        type(field_can_t) :: f
        real(dp) :: h0
        integer :: it, ierr

        call field_can_init(f, 1.0e-5_dp, 1.0_dp, VPAR0)
        call field_can_from_name('test')
        call eval_field(f, Z0(1), Z0(2), Z0(3), 0)
        call get_val(f, Z0(4))

        ! si%dt is the SUBSTEP and si%ntau the number of them, so the macro-step
        ! is ntau*dt. Dividing here keeps the interval fixed under refinement.
        call orbit_sympl_init(si, f, Z0, DT/real(ntau, dp), ntau, rtol, mode)
        h0 = f%H
        hdrift = 0.0_dp
        do it = 1, NSTEP
            call orbit_timestep_sympl(si, f, ierr)
            if (ierr /= 0) then
                write (error_unit, '(a,a,a,i0)') '  ', trim(tag), &
                    ' step failed with ierr = ', ierr
                nfail = nfail + 1
                exit
            end if
            hdrift = max(hdrift, abs((f%H - h0)/h0))
        end do
        zend = si%z(1:4)
    end subroutine trace

end program test_explicit_integmodes
