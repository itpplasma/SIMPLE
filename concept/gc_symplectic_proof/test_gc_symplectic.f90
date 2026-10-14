program test_gc_symplectic
    !> Executable L5 oracle for the guiding-centre symplectic integrators.
    !>
    !> The map under test is the actual production step `orbit_timestep_sympl`
    !> acting on the canonical state w = (theta, phi, p_theta, p_phi), with
    !> p_theta carried in f%pth exactly as the stepper carries it between steps
    !> (si%z(1) is only the radial Newton seed). Modes:
    !>   map     <file> <field_type> <E_kin_eV>   FD defect ||M^T J M - J||_F
    !>   energy  <file> <field_type> <E_kin_eV> <nsteps>   H envelope and slope
    !> <file> = test selects the analytic tokamak field (smoke only).
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use orbit_symplectic, only: orbit_sympl_init, orbit_timestep_sympl
    use orbit_symplectic_base, only: symplectic_integrator_t, EXPL_IMPL_EULER, &
        IMPL_EXPL_EULER, MIDPOINT, GAUSS1, GAUSS2, GAUSS3, GAUSS4, LOBATTO3
    use field_can_mod, only: field_can_t, eval_field => evaluate, get_val, &
        get_derivatives, field_can_from_name, field_can_init
    use simple, only: tracer_t, init_params, encode_symplectic_state
    use simple_main, only: init_field
    use params, only: coord_input, field_input
    use velo_mod, only: isw_field_type

    implicit none

    integer, parameter :: nvar = 8
    integer, parameter :: modes(nvar) = [EXPL_IMPL_EULER, IMPL_EXPL_EULER, &
        MIDPOINT, GAUSS1, GAUSS2, GAUSS3, GAUSS4, LOBATTO3]
    character(len=16), parameter :: names(nvar) = [character(len=16) :: &
        'EXPL_IMPL_EULER', 'IMPL_EXPL_EULER', 'MIDPOINT', 'GAUSS1', 'GAUSS2', &
        'GAUSS3', 'GAUSS4', 'LOBATTO3']
    ! Pre-committed acceptance bounds (scaled canonical Frobenius norm).
    real(dp), parameter :: tol_one = 1.0e-7_dp, tol_multi = 1.0e-5_dp
    integer, parameter :: nmulti = 100
    real(dp), parameter :: rtol = 1.0e-13_dp

    type(field_can_t) :: f, fe
    type(tracer_t) :: norb
    character(len=256) :: arg, mode_arg, file_arg
    integer :: field_type, nsteps, nfail, istart, nstart
    real(dp) :: ekin, dt, pscale, rseed
    real(dp) :: w0(4), starts(4, 2)
    logical :: smoke

    call get_command_argument(1, mode_arg)
    call get_command_argument(2, file_arg)
    call get_command_argument(3, arg); read(arg, *) field_type
    call get_command_argument(4, arg); read(arg, *) ekin
    nsteps = 0
    if (command_argument_count() >= 5) then
        call get_command_argument(5, arg); read(arg, *) nsteps
    end if
    smoke = trim(file_arg) == 'test'

    ! (s, theta, phi, lambda): a passing and a trapped start.
    starts(:, 1) = [0.35_dp, 0.7_dp, 0.1_dp, 0.8_dp]
    starts(:, 2) = [0.5_dp, 1.2_dp, 0.3_dp, 0.2_dp]
    nstart = 2
    call setup_field
    nfail = 0
    do istart = 1, nstart
        call initial_state(starts(:, istart), w0)
        write(*, '(a,i0,a,4es24.16,a,es10.3)') 'start ', istart, ' w0=', w0, &
            ' pscale=', pscale
        select case (trim(mode_arg))
        case ('map')
            call run_map(w0)
        case ('energy')
            call run_energy(w0)
        case default
            error stop 'mode must be map or energy'
        end select
    end do
    if (nfail > 0) then
        write(*, '(a,i0,a)') 'FAIL: ', nfail, ' check(s) failed'
        error stop 1
    end if
    print *, 'PASS'

contains

    subroutine setup_field
        if (smoke) then
            call field_can_from_name('test')
            call field_can_init(f, 1.0e-5_dp, 1.0_dp, 0.0_dp)
            dt = 7800.0_dp/64.0_dp
            nstart = 1
            return
        end if
        field_input = trim(file_arg)
        coord_input = trim(file_arg)
        isw_field_type = field_type
        call init_field(norb, trim(file_arg), 5, 5, 3, 1)
        call init_params(norb, 2, 4, ekin, 128, 1, rtol)
        dt = norb%dtaumin/sqrt(2.0_dp)
    end subroutine setup_field

    subroutine initial_state(start, w)
        real(dp), intent(in) :: start(4)
        real(dp), intent(out) :: w(4)
        real(dp) :: zgc(5), z(4)

        if (smoke) then
            z = [0.1_dp, 1.5_dp, 0.0_dp, 0.0_dp]
            call eval_field(f, z(1), z(2), z(3), 0)
            z(4) = f%vpar*f%hph + f%Aph/f%ro0
        else
            zgc = [start(1), start(2), start(3), 1.0_dp, start(4)]
            call encode_symplectic_state(f, zgc, z)
        end if
        call eval_field(f, z(1), z(2), z(3), 0)
        call get_derivatives(f, z(4))
        rseed = z(1)
        pscale = abs(f%dpth(1))
        w = [z(2), z(3), f%pth, z(4)]
    end subroutine initial_state

    subroutine radius_of(w, r, ok)
        !> Implicit radius R(q, p): p_theta(R, theta, phi, p_phi) = w(3).
        real(dp), intent(in) :: w(4)
        real(dp), intent(out) :: r
        logical, intent(out) :: ok
        real(dp) :: dr
        integer :: it

        fe = f
        r = rseed
        ok = .false.
        do it = 1, 60
            call eval_field(fe, r, w(1), w(2), 0)
            call get_derivatives(fe, w(4))
            dr = (fe%pth - w(3))/fe%dpth(1)
            r = r - dr
            if (abs(dr) < 1.0e-15_dp*max(1.0_dp, abs(r))) then
                ok = .true.
                exit
            end if
        end do
        call eval_field(fe, r, w(1), w(2), 0)
        call get_derivatives(fe, w(4))
    end subroutine radius_of

    subroutine step_map(mode, h, nst, w, wout, ok)
        !> nst production steps of fixed length h from canonical state w.
        integer, intent(in) :: mode, nst
        real(dp), intent(in) :: h, w(4)
        real(dp), intent(out) :: wout(4)
        logical, intent(out) :: ok
        type(symplectic_integrator_t) :: si
        real(dp) :: r
        integer :: k, ierr

        call radius_of(w, r, ok)
        if (.not. ok) return
        call orbit_sympl_init(si, f, [r, w(1), w(2), w(4)], h, 1, rtol, mode)
        f%pth = w(3)
        si%pthold = w(3)
        do k = 1, nst
            call orbit_timestep_sympl(si, f, ierr)
            if (ierr /= 0) then
                ok = .false.
                return
            end if
        end do
        wout = [si%z(2), si%z(3), f%pth, si%z(4)]
    end subroutine step_map

    subroutine state_dependent_map(w, wout, ok)
        !> Negative control: GAUSS2 with h(w) = dt*(1 + 0.05 sin theta), the
        !> structure of an event-selected (landing) step length.
        real(dp), intent(in) :: w(4)
        real(dp), intent(out) :: wout(4)
        logical, intent(out) :: ok
        call step_map(GAUSS2, dt*(1.0_dp + 0.05_dp*sin(w(1))), 1, w, wout, ok)
    end subroutine state_dependent_map

    subroutine fd_jacobian(mode, nst, w, delta, jac, ok, control, halves)
        integer, intent(in) :: mode, nst
        real(dp), intent(in) :: w(4), delta
        real(dp), intent(out) :: jac(4, 4)
        logical, intent(out) :: ok
        logical, intent(in) :: control, halves
        real(dp) :: scl(4), e(4), wp(4, 8), wo(4)
        integer :: j, k
        ! Eighth-order central derivative, retaining the predeclared spacing
        ! pair and defect bounds. The original fourth-order stencil has a
        ! resolved truncation error above 1e-5 for the 100-step NCSX map.
        real(dp), parameter :: coef(8) = [1.0_dp/280, -4.0_dp/105, &
            1.0_dp/5, -4.0_dp/5, 4.0_dp/5, -1.0_dp/5, &
            4.0_dp/105, -1.0_dp/280]
        real(dp), parameter :: offs(8) = &
            [-4.0_dp, -3.0_dp, -2.0_dp, -1.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
        logical :: okk

        ! Scaled coordinates (theta, phi, p_theta/P, p_phi/P): a conformal
        ! rescaling, so symplecticity is unchanged and the defect is O(1)-scaled.
        scl = [1.0_dp, 1.0_dp, pscale, pscale]
        ok = .true.
        do j = 1, 4
            do k = 1, 8
                e = 0.0_dp
                e(j) = offs(k)*delta*scl(j)
                if (control) then
                    call state_dependent_map(w + e, wo, okk)
                else if (halves) then
                    call step_map(mode, 0.5_dp*dt, 2*nst, w + e, wo, okk)
                else
                    call step_map(mode, dt, nst, w + e, wo, okk)
                end if
                ok = ok .and. okk
                wp(:, k) = wo/scl
            end do
            jac(:, j) = matmul(wp, coef)/delta
        end do
    end subroutine fd_jacobian

    real(dp) function sympl_defect(m) result(d)
        real(dp), intent(in) :: m(4, 4)
        real(dp) :: jm(4, 4), om(4, 4)
        ! Omega for w = (q1, q2, p1, p2): dp^dq ordering, Omega = [0 -I; I 0].
        om = 0.0_dp
        om(1, 3) = -1.0_dp; om(2, 4) = -1.0_dp
        om(3, 1) = 1.0_dp; om(4, 2) = 1.0_dp
        jm = matmul(transpose(m), matmul(om, m)) - om
        d = sqrt(sum(jm**2))
    end function sympl_defect

    subroutine check_map(label, mode, nst, w, tol, control, halves)
        character(*), intent(in) :: label
        integer, intent(in) :: mode, nst
        real(dp), intent(in) :: w(4), tol
        logical, intent(in) :: control, halves
        real(dp) :: m1(4, 4), m2(4, 4), d1, d2, dis
        logical :: ok1, ok2, pass

        call fd_jacobian(mode, nst, w, 1.0e-4_dp, m1, ok1, control, halves)
        call fd_jacobian(mode, nst, w, 5.0e-5_dp, m2, ok2, control, halves)
        if (.not. (ok1 .and. ok2)) then
            write(*, '(a,a)') 'map step refused: ', trim(label)
            nfail = nfail + 1
            return
        end if
        d1 = sympl_defect(m1)
        d2 = sympl_defect(m2)
        dis = sqrt(sum((m1 - m2)**2))
        if (control) then
            pass = min(d1, d2) > 100.0_dp*tol .and. dis < min(d1, d2)
        else
            pass = max(d1, d2) <= tol .and. dis <= tol
        end if
        write(*, '(a,a16,a,i4,a,es10.3,a,es10.3,a,es10.3,a,es10.3,a,es9.2,1x,a)') &
            'MAP ', label, ' steps=', nst, ' defect(1e-4)=', d1, &
            ' defect(5e-5)=', d2, ' spacing_disagreement=', dis, &
            ' |M|=', sqrt(sum(m1**2)), ' tol=', tol, merge('PASS', 'FAIL', pass)
        if (.not. pass) nfail = nfail + 1
    end subroutine check_map

    subroutine run_map(w)
        real(dp), intent(in) :: w(4)
        integer :: i

        do i = 1, nvar
            call check_map(names(i), modes(i), 1, w, tol_one, .false., .false.)
            if (.not. smoke) &
                call check_map(names(i), modes(i), nmulti, w, tol_multi, .false., &
                    .false.)
        end do
        if (smoke) return
        ! Prescribed subdivision (the retry/halving path with a fixed pattern).
        call check_map('GAUSS2_2xhalf', GAUSS2, 1, w, tol_one, .false., .true.)
        ! Negative control: a state-selected step length is not symplectic.
        call check_map('STATE_DEP_H', GAUSS2, 1, w, tol_one, .true., .false.)
        call midpoint_carry_diagnostic(w)
    end subroutine run_map

    subroutine midpoint_carry_diagnostic(w)
        !> Report-only diagnostic of the endpoint p_theta carried by MIDPOINT:
        !> f%pth after the step versus p_theta(si%z) evaluated at the endpoint.
        real(dp), intent(in) :: w(4)
        type(symplectic_integrator_t) :: si
        real(dp) :: r, tols(3)
        integer :: k, ierr
        logical :: ok

        tols = [1.0e-6_dp, 1.0e-9_dp, 1.0e-13_dp]
        do k = 1, 3
            call radius_of(w, r, ok)
            call orbit_sympl_init(si, f, [r, w(1), w(2), w(4)], dt, 1, tols(k), &
                MIDPOINT)
            f%pth = w(3)
            si%pthold = w(3)
            call orbit_timestep_sympl(si, f, ierr)
            fe = f
            call eval_field(fe, si%z(1), si%z(2), si%z(3), 0)
            call get_val(fe, si%z(4))
            write(*, '(a,es8.1,a,es10.3)') 'DIAG MIDPOINT rtol=', tols(k), &
                ' |carried p_theta - p_theta(endpoint)|/P=', abs(f%pth - fe%pth)/pscale
        end do
    end subroutine midpoint_carry_diagnostic

    subroutine run_energy(w)
        real(dp), intent(in) :: w(4)
        type(symplectic_integrator_t) :: si
        real(dp) :: r, h0, y, ymin, ymax, half_max(2), sk, sy, skk, sky, kc
        real(dp) :: slope, wcur(4), width, bound, rseed0
        integer :: i, k, ierr
        logical :: ok, pass

        rseed0 = rseed
        do i = 1, nvar
            rseed = rseed0
            call radius_of(w, r, ok)
            h0 = fe%H
            call orbit_sympl_init(si, f, [r, w(1), w(2), w(4)], dt, 1, rtol, modes(i))
            f%pth = w(3)
            si%pthold = w(3)
            ymin = 0.0_dp; ymax = 0.0_dp; half_max = 0.0_dp
            sk = 0.0_dp; sy = 0.0_dp; skk = 0.0_dp; sky = 0.0_dp
            ok = .true.
            do k = 1, nsteps
                call orbit_timestep_sympl(si, f, ierr)
                if (ierr /= 0) then
                    ok = .false.
                    exit
                end if
                wcur = [si%z(2), si%z(3), f%pth, si%z(4)]
                rseed = si%z(1)
                call radius_of(wcur, r, ok)
                if (.not. ok) exit
                y = fe%H/h0 - 1.0_dp
                ymin = min(ymin, y); ymax = max(ymax, y)
                if (2*k <= nsteps) then
                    half_max(1) = max(half_max(1), abs(y))
                else
                    half_max(2) = max(half_max(2), abs(y))
                end if
                kc = real(k, dp) - 0.5_dp*real(nsteps + 1, dp)
                sk = sk + kc; sy = sy + y; skk = skk + kc*kc; sky = sky + kc*y
            end do
            rseed = rseed0
            if (.not. ok) then
                write(*, '(a,a16,a,i0)') 'ENERGY ', names(i), ' refused at step ', k
                nfail = nfail + 1
                cycle
            end if
            slope = (sky - sk*sy/nsteps)/(skk - sk*sk/nsteps)
            width = ymax - ymin
            ! Bounded-envelope criterion: no secular growth. A linear drift
            ! doubles the second-half maximum and makes |slope|*N ~ width.
            bound = 2.0_dp*half_max(1) + 1.0e-12_dp
            pass = half_max(2) <= bound .and. abs(slope)*nsteps <= 0.5_dp*width &
                + 1.0e-12_dp
            write(*, '(a,a16,a,i0,a,es11.3,a,es11.3,a,es11.3,a,es11.3,a,es10.3,1x,a)') &
                'ENERGY ', names(i), ' N=', nsteps, ' min=', ymin, ' max=', ymax, &
                ' slope/step=', slope, ' slope*N/width=', abs(slope)*nsteps/ &
                max(width, tiny(1.0_dp)), ' max2/max1=', half_max(2)/ &
                max(half_max(1), tiny(1.0_dp)), merge('PASS', 'FAIL', pass)
            if (.not. pass) nfail = nfail + 1
        end do
    end subroutine run_energy

end program test_gc_symplectic
