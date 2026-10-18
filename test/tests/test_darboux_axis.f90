program test_darboux_axis
    !> Axis-grazing guiding-centre orbits on the analytic circular tokamak in
    !> the s = r**2 chart (field "test_s"), where E, mu and p_phi are exact
    !> invariants. The Darboux-chart midpoint (integmode DARBOUX_MIDPOINT) must
    !> take every step without retry, keep |dE|/E bounded and non-secular and
    !> conserve p_phi to roundoff. The flux-chart Euler1 and midpoint schemes
    !> run on the same markers for comparison and are only reported.
    !>
    !> Arguments (all optional): nsteps, ro0, npoiper2, compare (0/1),
    !> artifact (amplitude of a non-regular eps*cos(theta) Bmod axis term;
    !> diagnostic, the regression test uses 0).
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use util, only: pi, twopi
    use field_can_mod, only: field_can_t, field_can_from_name, &
        eval_field => evaluate, get_val, get_derivatives
    use field_can_base, only: n_field_evaluations
    use field_can_test, only: test_s_axis_artifact
    use orbit_symplectic, only: orbit_sympl_init, orbit_timestep_sympl, &
        advance_symplectic_with_retry
    use orbit_symplectic_base, only: symplectic_integrator_t, &
        SYMPLECTIC_STEP_OK, EXPL_IMPL_EULER, MIDPOINT, DARBOUX_MIDPOINT
    use orbit_symplectic_darboux, only: darboux_state_from_field, &
        darboux_eval, darboux_stats, darboux_stats_reset
    use diag_counters, only: diag_counters_init, diag_counters_total, &
        EVT_R_NEGATIVE, EVT_WARNING_MAXIT_ACCEPT, &
        EVT_WARNING_AXIS_CROSSING_ACCEPT

    implicit none

    integer, parameter :: n_mark = 6
    real(dp), parameter :: s_pass = 1.0e-5_dp
    real(dp) :: s0(n_mark), th0(n_mark), vpar0(n_mark)
    real(dp) :: ro0, dt, emax, eslope_rel, pphmax
    integer(int64) :: nsteps
    integer :: npoiper2, imark, compare, imode, nmode, modes(3)
    logical :: failed
    character(len=64) :: arg

    nsteps = 1000000_int64
    ro0 = 1.0e-2_dp
    npoiper2 = 256
    compare = 0
    if (command_argument_count() >= 1) then
        call get_command_argument(1, arg)
        read (arg, *) nsteps
    end if
    if (command_argument_count() >= 2) then
        call get_command_argument(2, arg)
        read (arg, *) ro0
    end if
    if (command_argument_count() >= 3) then
        call get_command_argument(3, arg)
        read (arg, *) npoiper2
    end if
    if (command_argument_count() >= 4) then
        call get_command_argument(4, arg)
        read (arg, *) compare
    end if
    if (command_argument_count() >= 5) then
        call get_command_argument(5, arg)
        read (arg, *) test_s_axis_artifact
    end if

    s0 = [5.0e-4_dp, 1.0e-4_dp, 1.0e-6_dp, 1.0e-8_dp, 1.0e-10_dp, 1.0e-12_dp]
    th0 = [0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp]
    vpar0 = [-0.3_dp, -0.3_dp, -0.3_dp, -0.28_dp, -0.32_dp, -0.3_dp]
    dt = twopi/real(npoiper2, dp)

    call field_can_from_name("test_s")
    call diag_counters_init

    modes = [DARBOUX_MIDPOINT, EXPL_IMPL_EULER, MIDPOINT]
    nmode = 1
    if (compare /= 0) nmode = 3

    print '(a,i0,a,es9.2,a,i0)', 'nsteps=', nsteps, ' ro0=', ro0, &
        ' npoiper2=', npoiper2
    print '(a)', 'mode mark   min_s     passes  max|dE|/E  slope*T/E  '// &
        'max|dpph|  fails  rneg  wacc  evals/step  maxit  us/step'
    failed = .false.
    do imode = 1, nmode
        do imark = 1, n_mark
            call trace(modes(imode), s0(imark), th0(imark), vpar0(imark), &
                imark, emax, eslope_rel, pphmax, failed)
        end do
    end do

    ! Tangent-map symplecticity M^T J M = J by central finite differences in
    ! (X, Y, phi, p_phi): one step at s ~ 1e-9 next to the axis, one in bulk.
    call tangent_check([3.0e-4_dp, 1.0e-5_dp, 0.0_dp, 0.0_dp], failed)
    call tangent_check([2.0e-2_dp, -1.0e-2_dp, 0.3_dp, 0.0_dp], failed)

    if (failed) error stop 'Darboux-chart axis check failed'
    print '(a)', 'PASS'

contains

    subroutine trace(mode, s_start, th_start, vpar_start, imark, emax, &
            eslope_rel, pphmax, failed)
        integer, intent(in) :: mode, imark
        real(dp), intent(in) :: s_start, th_start, vpar_start
        real(dp), intent(out) :: emax, eslope_rel, pphmax
        logical, intent(inout) :: failed

        type(symplectic_integrator_t) :: si
        type(field_can_t) :: f, fe
        real(dp) :: z(4), h0, pph0, e, smin, tsum, esum, ttsum, tesum, t
        real(dp) :: sigma, w(4), r, th, rhs(4), cpu0, cpu1, cpu_step
        integer(int64) :: k, nfail, passes, nev0, nev1
        integer(int64) :: rneg0, wacc0
        integer :: ierr, st
        logical :: inside

        f%mu = 0.5_dp*(1.0_dp - vpar_start**2)
        f%ro0 = ro0
        f%vpar = vpar_start
        call eval_field(f, s_start, th_start, 0.0_dp, 0)
        f%mu = f%mu/f%Bmod
        z = [s_start, th_start, 0.0_dp, vpar_start*f%hph + f%Aph/ro0]
        call get_val(f, z(4))
        h0 = f%H
        pph0 = z(4)

        call orbit_sympl_init(si, f, z, dt, 1, 1.0e-12_dp, mode)
        call eval_field(f, si%z(1), si%z(2), si%z(3), 0)
        call get_derivatives(f, si%z(4))
        call darboux_stats_reset()
        rneg0 = diag_counters_total(EVT_R_NEGATIVE)
        wacc0 = diag_counters_total(EVT_WARNING_MAXIT_ACCEPT) + &
            diag_counters_total(EVT_WARNING_AXIS_CROSSING_ACCEPT)

        nfail = 0
        passes = 0
        inside = .false.
        smin = huge(1.0_dp)
        emax = 0.0_dp
        pphmax = 0.0_dp
        tsum = 0.0_dp
        esum = 0.0_dp
        ttsum = 0.0_dp
        tesum = 0.0_dp
        nev0 = n_field_evaluations
        nev1 = 0
        cpu_step = 0.0_dp
        do k = 1, nsteps
            call cpu_time(cpu0)
            if (mode == DARBOUX_MIDPOINT) then
                call orbit_timestep_sympl(si, f, ierr)
            else
                call advance_symplectic_with_retry(si, f, &
                    orbit_timestep_sympl, ierr)
            end if
            call cpu_time(cpu1)
            cpu_step = cpu_step + (cpu1 - cpu0)
            if (ierr /= SYMPLECTIC_STEP_OK) then
                nfail = nfail + 1
                exit
            end if
            ! Energy oracle: H at the canonical phase point (theta, p_theta,
            ! phi, p_phi), independent of the scheme's stored quadrature r.
            nev1 = n_field_evaluations
            fe = f
            call darboux_state_from_field(f, si%z, sigma, w)
            r = si%z(1)
            call darboux_eval(w, sigma, si%z(2), r, th, fe, rhs, st)
            n_field_evaluations = nev1
            nev0 = nev0 + 0
            e = (fe%H - h0)/h0
            t = real(k, dp)
            emax = max(emax, abs(e))
            pphmax = max(pphmax, abs(si%z(4) - pph0)/max(abs(pph0), 1.0_dp))
            smin = min(smin, r)
            if (r < s_pass .and. .not. inside) passes = passes + 1
            inside = r < s_pass
            tsum = tsum + t
            esum = esum + e
            ttsum = ttsum + t*t
            tesum = tesum + t*e
        end do
        k = min(k, nsteps)
        ! Least-squares slope of dE/E over the run, times the run length.
        t = real(k, dp)
        eslope_rel = (t*tesum - tsum*esum)/max(t*ttsum - tsum*tsum, 1.0_dp)*t

        print '(i4,i5,es11.3,i9,3es11.3,i7,2i6,f11.1,i7,f9.3)', mode, imark, smin, &
            passes, emax, eslope_rel, pphmax, nfail, &
            diag_counters_total(EVT_R_NEGATIVE) - rneg0, &
            diag_counters_total(EVT_WARNING_MAXIT_ACCEPT) + &
            diag_counters_total(EVT_WARNING_AXIS_CROSSING_ACCEPT) - wacc0, &
            real(n_field_evaluations - nev0, dp)/real(max(k, 1_int64), dp), &
            darboux_stats%max_stage_iterations, &
            1.0e6_dp*cpu_step/real(max(k, 1_int64), dp)

        if (mode /= DARBOUX_MIDPOINT) return
        if (nfail /= 0 .or. darboux_stats%failures /= 0) then
            print *, 'FAIL: non-converged or rejected Darboux step, marker', imark
            failed = .true.
        end if
        if (emax > 1.0e-3_dp .or. abs(eslope_rel) > 0.1_dp*emax + 1.0e-12_dp) then
            print *, 'FAIL: energy not bounded/non-secular, marker', imark
            failed = .true.
        end if
        if (pphmax > 1.0e-10_dp) then
            print *, 'FAIL: p_phi not conserved, marker', imark
            failed = .true.
        end if
    end subroutine trace

    subroutine step_from_canonical(w0, pph_shift, w1, ierr)
        !> One Darboux midpoint step from canonical point w0; w0(4) is given
        !> relative to the p_phi of the reference marker (pph_shift).
        real(dp), intent(in) :: w0(4), pph_shift
        real(dp), intent(out) :: w1(4)
        integer, intent(out) :: ierr

        type(symplectic_integrator_t) :: si
        type(field_can_t) :: f
        real(dp) :: w(4), r, th, rhs(4), sigma, z(4)
        integer :: st

        w = w0
        w(4) = w0(4) + pph_shift
        f%ro0 = ro0
        f%mu = 0.5_dp*(1.0_dp - 0.3_dp**2)
        z = [1.0e-6_dp, 0.0_dp, 0.0_dp, w(4)]
        call orbit_sympl_init(si, f, z, dt, 1, 1.0e-14_dp, DARBOUX_MIDPOINT)
        r = 1.0e-6_dp
        call darboux_eval(w, 1.0_dp, 0.0_dp, r, th, f, rhs, st)
        si%z = [r, th, w(3), w(4)]
        f%pth = 0.5_dp*(w(1)**2 + w(2)**2)
        call orbit_timestep_sympl(si, f, ierr)
        call darboux_state_from_field(f, si%z, sigma, w1)
        w1(4) = w1(4) - pph_shift
    end subroutine step_from_canonical

    subroutine tangent_check(w0, failed)
        real(dp), intent(in) :: w0(4)
        logical, intent(inout) :: failed

        real(dp) :: m(4, 4), jmat(4, 4), wp(4), wm(4), w1(4), delta(4), dev
        real(dp) :: pph_shift
        integer :: i, ierr

        ! Reference p_phi of a v_par = -0.3 marker near the axis.
        pph_shift = -0.3_dp*1.0_dp
        delta = [1.0e-6_dp, 1.0e-6_dp, 1.0e-5_dp, 1.0e-5_dp]
        do i = 1, 4
            wp = w0
            wp(i) = wp(i) + delta(i)
            call step_from_canonical(wp, pph_shift, w1, ierr)
            if (ierr /= 0) failed = .true.
            m(:, i) = w1
            wm = w0
            wm(i) = wm(i) - delta(i)
            call step_from_canonical(wm, pph_shift, w1, ierr)
            if (ierr /= 0) failed = .true.
            m(:, i) = (m(:, i) - w1)/(2.0_dp*delta(i))
        end do
        call step_from_canonical(w0, pph_shift, w1, ierr)
        ! omega = dX^dY + dp_phi^dphi: pairs (q, p) = (Y, X) and (phi, p_phi),
        ! so in the ordering (X, Y) the block is [[0, -1], [1, 0]].
        jmat = 0.0_dp
        jmat(1, 2) = -1.0_dp
        jmat(2, 1) = 1.0_dp
        jmat(3, 4) = 1.0_dp
        jmat(4, 3) = -1.0_dp
        dev = maxval(abs(matmul(transpose(m), matmul(jmat, m)) - jmat))
        print '(a,4es11.3,a,4es11.3,a,es10.3)', 'tangent map from', w0, &
            ' to', w1, '  max|M^T J M - J| =', dev
        if (ierr /= 0 .or. dev > 1.0e-6_dp) then
            print *, 'FAIL: tangent map not symplectic'
            failed = .true.
        end if
    end subroutine tangent_check

end program test_darboux_axis
