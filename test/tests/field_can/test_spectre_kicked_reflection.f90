program test_spectre_kicked_reflection
    !> Kicked Level-1 reflection (CA15) on the current-free two-volume sheet
    !> fixture tok2vol_sheet (circular interfaces, B = F_l/R per volume, a
    !> pure-magnitude jump [[B]] > 0 outward). A forbidden crossing moves along
    !> the interface at frozen (v_par, p, lambda, mu) until either the home
    !> field returns to its entry value B* or the target field drops to B*.
    !> With B proportional to 1/R on both sides:
    !>
    !>   * a start near the |B| minimum (outboard, R(theta0)*(1+[[B]]/B) above
    !>     the largest R) reflects on the same side at the up-down conjugate
    !>     2 pi - theta0;
    !>   * a start further inboard reaches the target level B* first and
    !>     transmits into the outer volume at the equal-B point;
    !>   * both exits keep v_par, p and lambda and land on the frozen field
    !>     value, so H = v_par^2/2 + mu B is exact algebraically;
    !>   * an energetically allowed crossing still transmits (regression on the
    !>     radicand0 >= 0 branch).

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use simple, only: tracer_t
    use simple_main, only: init_field
    use params, only: read_config, params_init, netcdffile, ns_s, ns_tp, &
                      multharm, integmode
    use parmot_mod, only: ro0
    use magfie_sub, only: init_magfie, magfie, SPECTRE
    use interface_crossing, only: apply_crossing, crossing_info_t, &
                                  CROSSING_LEVEL1, CROSS_CROSSING, &
                                  CROSS_REFLECTION

    implicit none

    real(dp), parameter :: TWOPI = 6.283185307179586_dp
    real(dp), parameter :: EQUAL_B_TOL = 1.0d-10
    real(dp), parameter :: CONJ_TOL = 5.0d-3

    character(len=1024) :: h5file
    character(len=256) :: config_file
    type(tracer_t) :: norb
    type(crossing_info_t) :: info
    real(dp) :: y(5), y_out(5)
    real(dp) :: th0, mu, vpar, jump_max, dir_jump_max
    logical :: failed

    if (command_argument_count() < 1) then
        print *, 'usage: test_spectre_kicked_reflection <tok2vol_sheet.h5>'
        error stop 1
    end if
    call get_command_argument(1, h5file)

    call write_input(trim(h5file))
    config_file = 'simple.in'
    call read_config(config_file)
    call init_field(norb, netcdffile, ns_s, ns_tp, multharm, integmode)
    call params_init
    call init_magfie(SPECTRE)

    failed = .false.

    call scan_interface(jump_max, dir_jump_max)
    print '(A,ES11.3,A,ES11.3)', 'interface scan: max [[B]]/B = ', jump_max, &
        '  max |1-cos([[h]])| = ', dir_jump_max

    th0 = outboard_start(0.4_dp)
    call forbidden_state(th0, y, mu, vpar)
    call apply_crossing(y, 1, 1, 2, CROSSING_LEVEL1, y_out, info)
    call check_relocation('reflection', CROSS_REFLECTION, 1, &
        bmod_side(1.0_dp - 1.0d-12, y_out(2), y_out(3)))
    if (abs(mod(y_out(2) + 2.0_dp*TWOPI, TWOPI) &
            - mod(TWOPI - th0, TWOPI)) > CONJ_TOL) then
        print '(A,F10.6)', 'FAIL: exit not at the conjugate -theta0, theta = ', &
            y_out(2)
        failed = .true.
    end if

    th0 = outboard_start(0.8_dp)
    call forbidden_state(th0, y, mu, vpar)
    call apply_crossing(y, 1, 1, 2, CROSSING_LEVEL1, y_out, info)
    call check_relocation('transmission', CROSS_CROSSING, 2, &
        bmod_side(1.0_dp + 1.0d-12, y_out(2), y_out(3)))

    call allowed_state(th0, y, mu, vpar)
    call apply_crossing(y, 1, 1, 2, CROSSING_LEVEL1, y_out, info)
    print '(A,I0,A,I0)', 'allowed: event=', info%event_type, ' vol_to=', &
        info%vol_to
    if (info%event_type /= CROSS_CROSSING .or. info%vol_to /= 2) then
        print *, 'FAIL: allowed crossing did not transmit'
        failed = .true.
    end if

    if (failed) error stop 1
    print *, 'kicked Level-1 reflection PASS'

contains

    function outboard_start(th) result(th0)
        !> theta or 2 pi - theta, whichever side the tangential sheet drift
        !> (sign(ro0 h_zeta) in theta) leaves downhill in the home field.
        real(dp), intent(in) :: th
        real(dp) :: th0

        th0 = th
        if (drift_theta_sign(th0) > 0.0_dp) th0 = TWOPI - th
    end function outboard_start

    subroutine check_relocation(label, event_type, vol_to, bm_exit)
        character(*), intent(in) :: label
        integer, intent(in) :: event_type, vol_to
        real(dp), intent(in) :: bm_exit

        real(dp) :: bm_in

        bm_in = bmod_side(1.0_dp - 1.0d-12, y(2), y(3))
        print '(A,A,I0,A,I0,A,ES11.3,A,ES11.3)', label, ': event=', &
            info%event_type, ' vol_to=', info%vol_to, '  dtheta = ', info%dtheta, &
            '  |B_exit - B*|/B* = ', abs(bm_exit - bm_in)/bm_in
        if (info%event_type /= event_type .or. info%vol_to /= vol_to) then
            print *, 'FAIL: ', label, ' has the wrong event type or volume'
            failed = .true.
        end if
        if (abs(info%dtheta) < 1.0d-6) then
            print *, 'FAIL: ', label, ' produced no tangential displacement'
            failed = .true.
        end if
        if (abs(info%vpar_after - vpar) > 0.0_dp) then
            print *, 'FAIL: ', label, ' changed v_par'
            failed = .true.
        end if
        if (any(abs(y_out(4:5) - y(4:5)) > 0.0_dp)) then
            print *, 'FAIL: ', label, ' changed p or lambda'
            failed = .true.
        end if
        if (abs(bm_exit - bm_in)/bm_in > EQUAL_B_TOL) then
            print *, 'FAIL: ', label, ' exit field differs from the entry value'
            failed = .true.
        end if
    end subroutine check_relocation

    subroutine write_input(h5)
        character(*), intent(in) :: h5
        integer :: unit

        open (newunit=unit, file='simple.in', status='replace', action='write')
        write (unit, '(A)') '&config'
        write (unit, '(A)') "  field_input = '"//h5//"'"
        write (unit, '(A)') '  integ_coords = 6'
        write (unit, '(A)') '  integmode = 3'
        write (unit, '(A)') '  spectre_ncon_phi = 32'
        write (unit, '(A)') '  ntestpart = 1'
        write (unit, '(A)') '  ntimstep = 10'
        write (unit, '(A)') '  npoiper2 = 256'
        write (unit, '(A)') '  relerr = 1d-13'
        write (unit, '(A)') '  facE_al = 500.0d0'
        write (unit, '(A)') '  trace_time = 1.0d-6'
        write (unit, '(A)') '  sbeg = 0.97d0'
        write (unit, '(A)') '/'
        close (unit)
    end subroutine write_input

    function drift_theta_sign(th) result(sgn)
        !> Sign of the theta component of the in-sheet grad-B drift,
        !> sign(ro0) * sign(h_zeta), at the interface point (th, 0).
        real(dp), intent(in) :: th
        real(dp) :: sgn

        real(dp) :: bmod, sqrtg, bder(3), hcov(3), hctr(3), hcurl(3)

        call magfie([1.0_dp - 1.0d-12, th, 0.0_dp], bmod, sqrtg, bder, hcov, &
            hctr, hcurl)
        sgn = sign(1.0_dp, ro0)*sign(1.0_dp, hcov(3))
    end function drift_theta_sign

    function bmod_side(rho, th, ze) result(bmod)
        real(dp), intent(in) :: rho, th, ze
        real(dp) :: bmod

        real(dp) :: sqrtg, bder(3), hcov(3), hctr(3), hcurl(3)

        call magfie([rho, th, ze], bmod, sqrtg, bder, hcov, hctr, hcurl)
    end function bmod_side

    subroutine scan_interface(jump_max, dir_jump_max)
        !> Poloidal scan of the two-sided interface fields: the largest
        !> relative [[B]] and the largest direction jump 1 - h_home.h_target.
        real(dp), intent(out) :: jump_max, dir_jump_max

        integer :: k
        real(dp) :: th, bm, bp, cosj
        real(dp) :: sqrtg, bder(3), hcov_m(3), hctr_m(3), hcurl(3)
        real(dp) :: hcov_p(3), hctr_p(3)

        jump_max = 0.0_dp
        dir_jump_max = 0.0_dp
        do k = 0, 63
            th = TWOPI*real(k, dp)/64.0_dp
            call magfie([1.0_dp - 1.0d-12, th, 0.0_dp], bm, sqrtg, bder, &
                hcov_m, hctr_m, hcurl)
            call magfie([1.0_dp + 1.0d-12, th, 0.0_dp], bp, sqrtg, bder, &
                hcov_p, hctr_p, hcurl)
            cosj = dot_product(hctr_m, hcov_p)
            jump_max = max(jump_max, abs(bp - bm)/bm)
            dir_jump_max = max(dir_jump_max, abs(1.0_dp - cosj))
        end do
    end subroutine scan_interface

    subroutine forbidden_state(th0, y, mu, vpar)
        !> Landing state at theta0 with v_par^2 < 2 mu [[B]](theta0): p = 1,
        !> lambda chosen from the measured jump with a factor 1/2 margin.
        real(dp), intent(in) :: th0
        real(dp), intent(out) :: y(5), mu, vpar

        real(dp) :: bm, bp, lam

        bm = bmod_side(1.0_dp - 1.0d-12, th0, 0.0_dp)
        bp = bmod_side(1.0_dp + 1.0d-12, th0, 0.0_dp)
        if (bp <= bm) then
            print *, 'FAIL: fixture has no positive [[B]] at theta0'
            error stop 1
        end if
        ! vpar^2/(1-lam^2) * lam^2 ... with p=1: mu = (1-lam^2)/bm,
        ! vpar = lam; forbidden iff lam^2 < (1-lam^2) (bp-bm)/bm.
        lam = 0.5_dp*sqrt((bp - bm)/(bp - bm + bm))
        y = [1.0_dp, th0, 0.0_dp, 1.0_dp, lam]
        vpar = y(4)*y(5)
        mu = 0.5_dp*y(4)**2*(1.0_dp - y(5)**2)/bm
    end subroutine forbidden_state

    subroutine allowed_state(th0, y, mu, vpar)
        real(dp), intent(in) :: th0
        real(dp), intent(out) :: y(5), mu, vpar

        real(dp) :: bm, bp, lam

        bm = bmod_side(1.0_dp - 1.0d-12, th0, 0.0_dp)
        bp = bmod_side(1.0_dp + 1.0d-12, th0, 0.0_dp)
        ! lambda^2 well above the forbidden threshold (bp-bm)/bp.
        lam = sqrt(min(0.98_dp, 10.0_dp*max(bp - bm, 0.0_dp)/bp + 0.01_dp))
        y = [1.0_dp, th0, 0.0_dp, 1.0_dp, lam]
        vpar = y(4)*y(5)
        mu = 0.5_dp*y(4)**2*(1.0_dp - y(5)**2)/bm
    end subroutine allowed_state

end program test_spectre_kicked_reflection
