program test_darboux_chart
    !> Cartesian (Darboux) canonical chart on the analytic circular tokamak,
    !> where E, mu and p_phi are exact invariants. The oracles are independent
    !> of the generated leaves: the handwritten polar test field and its
    !> polar guiding-centre equations, finite differences, and the exact
    !> invariants.
    !>
    !> Usage: test_darboux_chart.x leaves | step | axis [nsteps]
    !>   leaves  regular test field vs polar field; q = h/s jets, the vector
    !>           field (vs the polar equations under the chart map) and the
    !>           stage Jacobian blocks vs finite differences
    !>   step    tangent-map symplecticity, order 2s for Gauss1-4, and
    !>           agreement with the polar midpoint away from the axis
    !>   axis    axis-grazing trapped markers over nsteps (default 1e6)
    !>           midpoint steps: bounded energy error, exact p_phi, every
    !>           step converged
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use util, only: pi, twopi
    use boozer_regular, only: regular_field_t
    use field_can_base, only: field_can_t, n_field_evaluations
    use field_can_mod, only: field_can_from_name, eval_field => evaluate, &
        get_val, get_derivatives
    use field_can_test, only: evaluate_test_regular
    use orbit_symplectic, only: orbit_sympl_init, orbit_timestep_sympl
    use orbit_symplectic_base, only: symplectic_integrator_t, &
        SYMPLECTIC_STEP_OK, MIDPOINT, GAUSS1, GAUSS2, GAUSS3, GAUSS4
    use orbit_symplectic_darboux, only: darboux_configure, darboux_advance, &
        darboux_point_eval, darboux_stage_eval, darboux_field_point, &
        darboux_stats, darboux_stats_reset, DARBOUX_FIELD_NONE, &
        DARBOUX_FIELD_TEST

    implicit none

    real(dp), parameter :: ro0 = 1.0e-2_dp
    character(len=16) :: what
    character(len=32) :: arg
    integer(int64) :: nsteps
    integer :: nfail

    nfail = 0
    what = 'leaves'
    nsteps = 1000000_int64
    if (command_argument_count() >= 1) call get_command_argument(1, what)
    if (command_argument_count() >= 2) then
        call get_command_argument(2, arg)
        read (arg, *) nsteps
    end if

    call field_can_from_name('test')
    select case (trim(what))
    case ('leaves')
        call darboux_configure(DARBOUX_FIELD_TEST, 1)
        call check_test_field()
        call check_quotient()
        call check_vector_field()
        call check_stage_jacobian()
    case ('step')
        call check_tangent_map([3.0e-4_dp, 1.0e-5_dp, 0.0_dp], 1)
        call check_tangent_map([3.0e-4_dp, 1.0e-5_dp, 0.0_dp], 2)
        call check_tangent_map([2.0e-2_dp, -1.0e-2_dp, 0.3_dp], 1)
        call check_order()
        call check_polar_agreement()
    case ('axis')
        call check_axis(nsteps)
        call check_kick()
    case default
        error stop 'usage: test_darboux_chart.x leaves|step|axis [nsteps]'
    end select

    if (nfail > 0) then
        print '(a,i0,a)', 'FAILED: ', nfail, ' check(s)'
        error stop 1
    end if
    print '(a)', 'PASSED: Darboux chart '//trim(what)

contains

    subroutine check(condition, label, value)
        logical, intent(in) :: condition
        character(*), intent(in) :: label
        real(dp), intent(in) :: value

        if (condition) then
            print '(a,es10.3)', 'PASS: '//label//' ', value
        else
            print '(a,es10.3)', 'FAIL: '//label//' ', value
            nfail = nfail + 1
        end if
    end subroutine check

    !> Polar field_can_t at (r, theta, phi) with first and second derivatives.
    subroutine polar_field(r, th, pphi, mu, f)
        real(dp), intent(in) :: r, th, pphi, mu
        type(field_can_t), intent(inout) :: f

        f%mu = mu
        f%ro0 = ro0
        call eval_field(f, r, th, 0.0_dp, 2)
        call get_derivatives(f, pphi)
    end subroutine polar_field

    subroutine check_test_field()
        type(field_can_t) :: f
        type(regular_field_t) :: rf
        real(dp) :: r, th, x, y, err, pol(5), reg(5), dpol(5, 2), dreg(5, 2)
        integer :: i, j

        err = 0.0_dp
        do i = 1, 4
            do j = 1, 5
                r = 0.08_dp*i
                th = 1.3_dp*j
                x = r*cos(th)
                y = r*sin(th)
                call polar_field(r, th, 0.0_dp, 0.0_dp, f)
                call evaluate_test_regular(x, y, 1, rf)
                pol = [f%Ath, f%Aph, f%hth, f%hph, f%Bmod]
                reg = [rf%Ath, rf%Aph, rf%hth, rf%hph, rf%Bmod]
                dpol(:, 1) = [f%dAth(1), f%dAph(1), f%dhth(1), f%dhph(1), f%dBmod(1)]
                dpol(:, 2) = [f%dAth(2), f%dAph(2), f%dhth(2), f%dhph(2), f%dBmod(2)]
                ! d/dr = (x d/dx + y d/dy)/r, d/dtheta = x d/dy - y d/dx
                dreg(:, 1) = (x*[rf%dAth(1), rf%dAph(1), rf%dhth(1), rf%dhph(1), &
                    rf%dBmod(1)] + y*[rf%dAth(2), rf%dAph(2), rf%dhth(2), &
                    rf%dhph(2), rf%dBmod(2)])/r
                dreg(:, 2) = x*[rf%dAth(2), rf%dAph(2), rf%dhth(2), rf%dhph(2), &
                    rf%dBmod(2)] - y*[rf%dAth(1), rf%dAph(1), rf%dhth(1), &
                    rf%dhph(1), rf%dBmod(1)]
                err = max(err, maxval(abs(pol - reg)), maxval(abs(dpol - dreg)))
            end do
        end do
        call check(err < 1.0e-13_dp, 'regular test field = polar test field', err)
    end subroutine check_test_field

    !> q = h/s jets against central differences of h/s computed directly.
    subroutine check_quotient()
        type(regular_field_t) :: rf
        real(dp) :: u(2), qa(10), qh(10), qp(10), qm(10), fd, err, d
        integer :: k, rho_case

        err = 0.0_dp
        do rho_case = 1, 3
            u = [0.3_dp, -0.2_dp]*10.0_dp**(1 - rho_case)
            d = 1.0e-3_dp*sqrt(u(1)**2 + u(2)**2)
            call darboux_field_point(u, 0.0_dp, 2, rf, qa, qh)
            do k = 1, 2
                call darboux_field_point(u + d*unit2(k), 0.0_dp, 2, rf, qp, qh)
                call darboux_field_point(u - d*unit2(k), 0.0_dp, 2, rf, qm, qh)
                fd = (direct(u + d*unit2(k)) - direct(u - d*unit2(k)))/(2*d)
                err = max(err, abs(qa(1 + k) - fd)/max(abs(fd), 1.0_dp))
                ! Second jets: d/dx of the first jets.
                fd = (qp(2) - qm(2))/(2*d)
                err = max(err, abs(qa(4 + k) - fd)/max(abs(fd), 1.0_dp))
            end do
            err = max(err, abs(qa(1) - direct(u))/abs(direct(u)))
        end do
        call check(err < 1.0e-6_dp, 'q = A_theta/s jets vs finite differences', err)
    end subroutine check_quotient

    real(dp) function direct(u)
        real(dp), intent(in) :: u(2)
        type(regular_field_t) :: rf

        call evaluate_test_regular(u(1), u(2), 0, rf)
        direct = rf%Ath/(u(1)**2 + u(2)**2)
    end function direct

    pure function unit2(k) result(e)
        integer, intent(in) :: k
        real(dp) :: e(2)

        e = 0.0_dp
        e(k) = 1.0_dp
    end function unit2

    !> The generated vector field against the polar guiding-centre equations
    !> (theta_dot = H_r/p_theta,r etc., as in the polar Gauss residual)
    !> pushed through the chart X = sqrt(2 J) cos(theta),
    !> Y = sigma sqrt(2 J) sin(theta), J = sigma p_theta.
    subroutine check_vector_field()
        type(field_can_t) :: f
        type(regular_field_t) :: rf
        real(dp) :: r, th, pphi, mu, u(2), kappa, pth, ham, vpar, chart(2), &
            rhs(4), gu(2, 2), udot(2), thdot, pthdot, phdot, pphdot, sigma, jj, &
            ref(4), err, errh, err_chain
        integer :: i

        err = 0.0_dp
        errh = 0.0_dp
        err_chain = 0.0_dp
        mu = 0.5_dp*(1.0_dp - 0.3_dp**2)
        do i = 1, 6
            r = 0.04_dp*i
            th = 0.9_dp*i
            call polar_field(r, th, 0.0_dp, mu, f)
            pphi = -0.3_dp*f%hph + f%Aph/ro0 + 0.01_dp*i
            call polar_field(r, th, pphi, mu, f)
            thdot = f%dH(1)/f%dpth(1)
            pthdot = -(f%dH(2) - thdot*f%dpth(2))
            phdot = (f%vpar - thdot*f%hth)/f%hph
            pphdot = -(f%dH(3) - thdot*f%dpth(3))
            sigma = sign(1.0_dp, f%pth)
            jj = sigma*f%pth
            ref(1) = sigma*pthdot*cos(th)/sqrt(2*jj) - sqrt(2*jj)*sin(th)*thdot
            ref(2) = pthdot*sin(th)/sqrt(2*jj) + sigma*sqrt(2*jj)*cos(th)*thdot
            ref(3) = phdot
            ref(4) = pphdot
            u = r*[cos(th), sin(th)]
            call darboux_point_eval(u, 0.0_dp, pphi, sigma, mu, ro0, rf, kappa, &
                                    pth, ham, vpar, chart, rhs, gu, udot)
            err = max(err, maxval(abs(rhs - ref))/maxval(abs(ref)))
            errh = max(errh, abs(ham - f%H)/abs(f%H), abs(pth - f%pth)/abs(f%pth), &
                       abs(chart(1) - sqrt(2*jj)*cos(th))/sqrt(2*jj), &
                       abs(chart(2) - sigma*sqrt(2*jj)*sin(th))/sqrt(2*jj))
            ! Chain rule of the chart map: G du/dt + C_phi,p_phi (dphi, dp_phi)
            ! = (dX, dY)/dt; here C does not depend on phi.
            err_chain = max(err_chain, maxval(abs(matmul(gu, udot) - rhs(1:2) &
                + chart_pphi(u, pphi, sigma, mu)*rhs(4)))/maxval(abs(rhs(1:2))))
        end do
        call check(errh < 1.0e-12_dp, 'H, p_theta and chart map vs polar', errh)
        call check(err < 1.0e-10_dp, 'vector field vs polar equations', err)
        call check(err_chain < 1.0e-6_dp, 'du/dt consistent with the chart map', &
                   err_chain)
    end subroutine check_vector_field

    function chart_pphi(u, pphi, sigma, mu) result(d)
        real(dp), intent(in) :: u(2), pphi, sigma, mu
        real(dp) :: d(2)
        real(dp) :: cp(2), cm(2), dp_
        type(regular_field_t) :: rf
        real(dp) :: kappa, pth, ham, vpar, rhs(4), gu(2, 2), udot(2)

        dp_ = 1.0e-6_dp*max(abs(pphi), 1.0_dp)
        call darboux_point_eval(u, 0.0_dp, pphi + dp_, sigma, mu, ro0, rf, kappa, &
                                pth, ham, vpar, cp, rhs, gu, udot)
        call darboux_point_eval(u, 0.0_dp, pphi - dp_, sigma, mu, ro0, rf, kappa, &
                                pth, ham, vpar, cm, rhs, gu, udot)
        d = (cp - cm)/(2*dp_)
    end function chart_pphi

    !> Stage Jacobian blocks against central differences of the first-jet
    !> vector field and chart map, near the axis and in the bulk.
    subroutine check_stage_jacobian()
        type(regular_field_t) :: rf
        real(dp) :: u(2), phi, pphi, sigma, mu, chart(2), rhs(4), fu(4, 2), &
            fw(4, 2), gu(2, 2), gw(2, 2), cp(2), cm(2), rp(4), rm(4), d, err, scl
        real(dp) :: kappa, pth, ham, vpar, g2(2, 2), udot(2)
        integer :: icase, k

        err = 0.0_dp
        mu = 0.5_dp*(1.0_dp - 0.3_dp**2)
        sigma = 1.0_dp
        phi = 0.4_dp
        do icase = 1, 3
            u = [0.2_dp, 0.1_dp]*10.0_dp**(2 - 2*icase)
            pphi = -0.2_dp
            call darboux_stage_eval(u, phi, pphi, sigma, mu, ro0, chart, rhs, fu, &
                                    fw, gu, gw)
            scl = maxval(abs(rhs))
            do k = 1, 2
                d = 1.0e-6_dp*max(sqrt(u(1)**2 + u(2)**2), 1.0e-3_dp)
                call darboux_point_eval(u + d*unit2(k), phi, pphi, sigma, mu, ro0, &
                                        rf, kappa, pth, ham, vpar, cp, rp, g2, udot)
                call darboux_point_eval(u - d*unit2(k), phi, pphi, sigma, mu, ro0, &
                                        rf, kappa, pth, ham, vpar, cm, rm, g2, udot)
                err = max(err, maxval(abs(fu(:, k) - (rp - rm)/(2*d))) &
                    /max(maxval(abs(fu(:, k))), scl))
                err = max(err, maxval(abs(gu(:, k) - (cp - cm)/(2*d))) &
                    /maxval(abs(gu)))
            end do
            d = 1.0e-6_dp*max(abs(pphi), 1.0_dp)
            call darboux_point_eval(u, phi, pphi + d, sigma, mu, ro0, rf, kappa, &
                                    pth, ham, vpar, cp, rp, g2, udot)
            call darboux_point_eval(u, phi, pphi - d, sigma, mu, ro0, rf, kappa, &
                                    pth, ham, vpar, cm, rm, g2, udot)
            err = max(err, maxval(abs(fw(:, 2) - (rp - rm)/(2*d))) &
                /max(maxval(abs(fw(:, 2))), scl))
            err = max(err, maxval(abs(gw(:, 2) - (cp - cm)/(2*d))) &
                /max(maxval(abs(gw(:, 2))), 1.0e-30_dp))
            ! Axisymmetric field: no phi dependence.
            err = max(err, maxval(abs(fw(:, 1)))/scl, maxval(abs(gw(:, 1))))
        end do
        call check(err < 1.0e-6_dp, 'stage Jacobian blocks vs finite differences', &
                   err)
    end subroutine check_stage_jacobian

    !> Field point u and the first-jet data of the Darboux point w by Newton
    !> on the chart map (test helper, independent of the step).
    subroutine chart_inverse(w, sigma, mu, u, rhs, udot, kappa)
        real(dp), intent(in) :: w(4), sigma, mu
        real(dp), intent(inout) :: u(2)
        real(dp), intent(out) :: rhs(4), udot(2), kappa
        type(regular_field_t) :: rf
        real(dp) :: pth, ham, vpar, chart(2), gu(2, 2), r(2), det
        integer :: it

        do it = 1, 50
            call darboux_point_eval(u, w(3), w(4), sigma, mu, ro0, rf, kappa, pth, &
                                    ham, vpar, chart, rhs, gu, udot)
            r = chart - w(1:2)
            if (all(abs(r) <= 1.0e-15_dp*sqrt(2*kappa))) return
            det = gu(1, 1)*gu(2, 2) - gu(1, 2)*gu(2, 1)
            u = u - [gu(2, 2)*r(1) - gu(1, 2)*r(2), gu(1, 1)*r(2) - gu(2, 1)*r(1)] &
                /det
        end do
        error stop 'chart_inverse did not converge'
    end subroutine chart_inverse

    subroutine advance_from(w0, nst, h, stages, w1, status)
        !> nst Gauss steps from the Darboux point w0.
        real(dp), intent(in) :: w0(4), h
        integer, intent(in) :: nst, stages
        real(dp), intent(out) :: w1(4)
        integer, intent(out) :: status
        type(regular_field_t) :: rf
        real(dp) :: w(4), u(2), rhs(4), udot(2), kappa, wn(4), un(2), ham, vpar
        real(dp) :: mu, sigma
        integer :: k

        mu = 0.5_dp*(1.0_dp - 0.3_dp**2)
        sigma = 1.0_dp
        call darboux_configure(DARBOUX_FIELD_TEST, stages)
        w = w0
        u = [w0(1), sigma*w0(2)]/10.0_dp
        call chart_inverse(w, sigma, mu, u, rhs, udot, kappa)
        do k = 1, nst
            call darboux_advance(w, u, rhs, udot, sigma, mu, ro0, h, 1.0e-14_dp, &
                                 sqrt(2*kappa), wn, un, rhs, udot, rf, kappa, ham, &
                                 vpar, status)
            if (status /= SYMPLECTIC_STEP_OK) return
            w = wn
            u = un
        end do
        w1 = w
    end subroutine advance_from

    !> Reference p_phi of a v_par = -0.3 marker at the axis, so that w(4)
    !> offsets are small numbers.
    real(dp) function pphi_axis()
        pphi_axis = -0.3_dp
    end function pphi_axis

    subroutine check_tangent_map(xyphi, stages)
        real(dp), intent(in) :: xyphi(3)
        integer, intent(in) :: stages
        real(dp) :: w0(4), wp(4), wm(4), w1(4), m(4, 4), jmat(4, 4), delta(4), dev, h
        integer :: i, status
        character(len=64) :: label

        h = twopi/256.0_dp
        w0 = [xyphi, pphi_axis()]
        delta = [1.0e-6_dp, 1.0e-6_dp, 1.0e-5_dp, 1.0e-6_dp]
        do i = 1, 4
            wp = w0
            wp(i) = wp(i) + delta(i)
            call advance_from(wp, 1, h, stages, w1, status)
            if (status /= 0) nfail = nfail + 1
            m(:, i) = w1
            wm = w0
            wm(i) = wm(i) - delta(i)
            call advance_from(wm, 1, h, stages, w1, status)
            if (status /= 0) nfail = nfail + 1
            m(:, i) = (m(:, i) - w1)/(2*delta(i))
        end do
        ! omega = dX^dY + dp_phi^dphi in the ordering (X, Y, phi, p_phi)
        jmat = 0.0_dp
        jmat(1, 2) = -1.0_dp
        jmat(2, 1) = 1.0_dp
        jmat(3, 4) = 1.0_dp
        jmat(4, 3) = -1.0_dp
        dev = maxval(abs(matmul(transpose(m), matmul(jmat, m)) - jmat))
        write (label, '(a,i0,a,es8.1)') 'tangent map Gauss', stages, &
            ' symplectic, |X| ~', abs(xyphi(1))
        call check(dev < 1.0e-7_dp, trim(label), dev)
    end subroutine check_tangent_map

    !> Global error at fixed time T through the near-axis region against a
    !> fine-step reference: order 2s for Gauss-s.
    subroutine check_order()
        real(dp) :: w0(4), wref(4), w1(4), err(3), rate, tend
        integer :: stages, k, nbase, status
        integer, parameter :: nbases(4) = [32, 32, 4, 2]
        character(len=64) :: label

        w0 = [1.0e-3_dp, 2.0e-4_dp, 0.0_dp, pphi_axis()]
        tend = 4.0_dp
        do stages = 1, 4
            ! Coarser steps for the high orders keep the errors above roundoff.
            nbase = nbases(stages)
            call advance_from(w0, 64*nbase, tend/(64*nbase), stages, wref, status)
            do k = 1, 3
                call advance_from(w0, nbase*2**(k - 1), tend/(nbase*2**(k - 1)), &
                                  stages, w1, status)
                err(k) = maxval(abs(w1 - wref))
            end do
            rate = log(err(1)/err(2))/log(2.0_dp)
            write (label, '(a,i0,a,i0,a)') 'Gauss', stages, ' order ', 2*stages, &
                ' (observed)'
            print '(a,i0,a,3es10.2)', '  Gauss', stages, ' errors ', err
            call check(rate > 2*stages - 0.3_dp, trim(label), rate)
        end do
    end subroutine check_order

    !> Away from the axis the Cartesian and polar midpoint maps approximate
    !> the same flow: their difference is O(h**2) and falls by 4 when h is
    !> halved.
    subroutine check_polar_agreement()
        real(dp) :: diff(2), rate
        integer :: k

        do k = 1, 2
            diff(k) = chart_difference(256*2**(k - 1), 64*2**(k - 1))
        end do
        rate = log(diff(1)/diff(2))/log(2.0_dp)
        print '(a,2es10.2)', '  polar - cartesian midpoint ', diff
        call check(diff(1) < 1.0e-3_dp .and. rate > 1.7_dp, &
                   'Cartesian = polar midpoint off axis (order)', rate)
    end subroutine check_polar_agreement

    real(dp) function chart_difference(npoiper, nst)
        integer, intent(in) :: npoiper, nst
        type(symplectic_integrator_t) :: si
        type(field_can_t) :: f
        real(dp) :: z(4), zc(4), zp(4), dt
        integer :: k, ierr, pass

        dt = twopi/real(npoiper, dp)
        do pass = 1, 2
            if (pass == 1) then
                call darboux_configure(DARBOUX_FIELD_NONE, 1)
            else
                call darboux_configure(DARBOUX_FIELD_TEST, 1)
            end if
            call start_marker(0.3_dp, 0.5_dp, -0.3_dp, f, z)
            call orbit_sympl_init(si, f, z, dt, 1, 1.0e-13_dp, MIDPOINT)
            do k = 1, nst
                call orbit_timestep_sympl(si, f, ierr)
                if (ierr /= 0) error stop 'chart_difference: step failed'
            end do
            if (pass == 1) then
                zp = si%z
            else
                zc = si%z
            end if
        end do
        call darboux_configure(DARBOUX_FIELD_NONE, 1)
        chart_difference = maxval(abs(zc - zp))
    end function chart_difference

    subroutine start_marker(r0, th0, vpar0, f, z)
        real(dp), intent(in) :: r0, th0, vpar0
        type(field_can_t), intent(inout) :: f
        real(dp), intent(out) :: z(4)

        f%ro0 = ro0
        f%vpar = vpar0
        call eval_field(f, r0, th0, 0.0_dp, 0)
        f%mu = 0.5_dp*(1.0_dp - vpar0**2)/f%Bmod
        z = [r0, th0, 0.0_dp, vpar0*f%hph + f%Aph/ro0]
        call get_val(f, z(4))
    end subroutine start_marker

    !> Trapped markers started at s = r**2 from 5e-4 down to 1e-12 pass the
    !> axis many times; midpoint in the Darboux chart takes every step.
    subroutine check_axis(nst)
        integer(int64), intent(in) :: nst
        integer, parameter :: n_mark = 6
        real(dp), parameter :: s0(n_mark) = [5.0e-4_dp, 1.0e-4_dp, 1.0e-6_dp, &
            1.0e-8_dp, 1.0e-10_dp, 1.0e-12_dp]
        real(dp), parameter :: th0(n_mark) = [0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, &
            4.0_dp, 5.0_dp]
        real(dp), parameter :: vpar0(n_mark) = [-0.3_dp, -0.3_dp, -0.3_dp, &
            -0.28_dp, -0.32_dp, -0.3_dp]
        type(symplectic_integrator_t) :: si
        type(field_can_t) :: f
        real(dp) :: z(4), h0, pph0, e, emax, pphmax, smin, t, tsum, esum, ttsum, &
            tesum, slope, dt
        integer(int64) :: k, nev0, passes
        integer :: imark, ierr
        logical :: inside

        dt = twopi/256.0_dp
        call darboux_configure(DARBOUX_FIELD_TEST, 1)
        print '(a)', ' mark     s_min   passes  max|dE|/E  slope*T/E  max|dpphi|'// &
            '  failures  evals/step'
        do imark = 1, n_mark
            call start_marker(sqrt(s0(imark)), th0(imark), vpar0(imark), f, z)
            h0 = f%H
            pph0 = z(4)
            call orbit_sympl_init(si, f, z, dt, 1, 1.0e-12_dp, MIDPOINT)
            call darboux_stats_reset()
            nev0 = n_field_evaluations
            emax = 0.0_dp
            pphmax = 0.0_dp
            smin = huge(1.0_dp)
            tsum = 0.0_dp
            esum = 0.0_dp
            ttsum = 0.0_dp
            tesum = 0.0_dp
            passes = 0
            inside = .false.
            ierr = 0
            do k = 1, nst
                call orbit_timestep_sympl(si, f, ierr)
                if (ierr /= SYMPLECTIC_STEP_OK) exit
                e = (f%H - h0)/h0
                t = real(k, dp)
                emax = max(emax, abs(e))
                pphmax = max(pphmax, abs(si%z(4) - pph0))
                smin = min(smin, si%z(1)**2)
                if (si%z(1)**2 < 1.0e-5_dp .and. .not. inside) passes = passes + 1
                inside = si%z(1)**2 < 1.0e-5_dp
                tsum = tsum + t
                esum = esum + e
                ttsum = ttsum + t*t
                tesum = tesum + t*e
            end do
            t = real(min(k, nst), dp)
            slope = (t*tesum - tsum*esum)/max(t*ttsum - tsum*tsum, 1.0_dp)*t
            print '(i5,es10.2,i9,3es11.3,i10,f12.2)', imark, smin, passes, emax, &
                slope, pphmax, darboux_stats%failures, &
                real(n_field_evaluations - nev0, dp)/t
            call check(ierr == 0 .and. darboux_stats%failures == 0, &
                       'every step converged', real(darboux_stats%failures, dp))
            call check(emax < 1.0e-3_dp .and. abs(slope) < 0.1_dp*emax + 1.0e-12_dp, &
                       'energy error bounded and non-secular', emax)
            call check(pphmax <= 0.0_dp, 'p_phi exact', pphmax)
        end do
        call darboux_configure(DARBOUX_FIELD_NONE, 1)
    end subroutine check_axis

    !> A collision kick changes v_par and mu at fixed (r, theta, phi) the
    !> way simple_main update_momentum does. The next step re-encodes the
    !> Darboux state, after which the new energy is conserved.
    subroutine check_kick()
        type(symplectic_integrator_t) :: si
        type(field_can_t) :: f
        real(dp) :: z(4), dt, h1, emax, vpar_new, p_new
        integer(int64) :: encodes
        integer :: k, ierr

        dt = twopi/256.0_dp
        call darboux_configure(DARBOUX_FIELD_TEST, 2)
        call start_marker(1.0e-3_dp, 0.3_dp, -0.3_dp, f, z)
        call orbit_sympl_init(si, f, z, dt, 1, 1.0e-12_dp, GAUSS2)
        call darboux_stats_reset()
        do k = 1, 1000
            call orbit_timestep_sympl(si, f, ierr)
            if (ierr /= 0) error stop 'check_kick: step failed'
        end do
        encodes = darboux_stats%encodes
        ! Kick: |v| = 0.9, pitch 0.5 (update_momentum with the regular f).
        p_new = 0.9_dp
        vpar_new = 0.5_dp*p_new
        f%vpar = vpar_new
        f%mu = 0.5_dp*(p_new**2 - vpar_new**2)/f%Bmod
        si%z(4) = f%vpar*f%hph + f%Aph/f%ro0
        call get_val(f, si%z(4))
        si%pthold = f%pth
        h1 = f%H
        emax = 0.0_dp
        do k = 1, 1000
            call orbit_timestep_sympl(si, f, ierr)
            if (ierr /= 0) error stop 'check_kick: step after kick failed'
            emax = max(emax, abs(f%H - h1)/h1)
        end do
        call check(darboux_stats%encodes == encodes + 1, &
                   'kick re-encodes the Darboux state once', &
                   real(darboux_stats%encodes - encodes, dp))
        call check(emax < 1.0e-6_dp, 'energy after kick conserved', emax)
        call darboux_configure(DARBOUX_FIELD_NONE, 1)
    end subroutine check_kick

end program test_darboux_chart
