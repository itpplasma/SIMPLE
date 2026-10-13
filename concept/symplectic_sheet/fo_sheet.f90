module fo_sheet
    !> Full-orbit implicit midpoint steps for the two-volume torus of sheet_field.
    !>
    !> A step of length h that meets the sheet s = r - r0 = 0 is split into
    !> segments. Each segment is a midpoint substep with the Hamiltonian of its
    !> region; all but the last end on the sheet. After each landing the state
    !> is kicked along the exact flow of X_s = J grad s (positions fixed,
    !> p -> p - l grad_q s) and continues in the other region. The kick l makes
    !> the discrete energies H(zbar) of the segments before and after the
    !> landing equal; the segment after it ends at the next landing or at the
    !> end of the step. This is the symplectic collision construction with the
    !> identity reset: H is continuous across the sheet, so l = O(h^2) and the
    !> step is exactly symplectic.
    !>
    !> tangent_map returns the exact Jacobian of a recorded step by implicit
    !> differentiation of the midpoint, landing and energy-matching conditions.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use sheet_field, only: field_params_t, hamiltonian, grad_hamiltonian, &
        hess_hamiltonian, sheet_distance, grad_sheet_distance, region_of, bfield
    implicit none
    private

    public :: midpoint_step, sheet_step, tangent_map, poincare_step, &
        poincare_gauss2_step, boris_step
    public :: step_stats_t, step_record_t, poincare_params_t
    public :: SCHEME_LANDING, SCHEME_NAIVE, SCHEME_MIDPOINT_REGION

    integer, parameter :: SCHEME_LANDING = 1
    !> Naive control: plain midpoint with the region of the step start; the
    !> region switches when the step ends on the other side.
    integer, parameter :: SCHEME_NAIVE = 2
    !> Plain midpoint on the true piecewise H: the branch is the region of the
    !> midpoint (z0 + z1)/2 at every Newton iterate.
    integer, parameter :: SCHEME_MIDPOINT_REGION = 3

    integer, parameter :: MAX_NEWTON = 50, MAX_SEGMENTS = 16
    integer, parameter :: MAX_FIXED_POINT = 60
    real(dp), parameter :: NEWTON_TOL = 1.0e-15_dp
    real(dp), parameter :: LANDING_TOL = 1.0e-14_dp
    real(dp), parameter :: ENERGY_TOL = 1.0e-10_dp
    real(dp), parameter :: KICK_TOL = 1.0e-14_dp
    !> Fraction of the remaining time used to step off the sheet when the next
    !> landing is bracketed.
    real(dp), parameter :: LIFTOFF = 1.0e-9_dp

    type :: step_stats_t
        integer :: crossings = 0
        integer :: multi_crossing_steps = 0
        integer :: failures = 0
        real(dp) :: max_kick = 0.0_dp
    end type step_stats_t

    !> Poincare time transformation dt = g(z) dtau with
    !> g = 1 - (1 - g_min) exp(-((r - r0)/length)^2): fixed steps in the
    !> fictitious time tau are short in physical time near the sheet.
    type :: poincare_params_t
        real(dp) :: g_min = 0.05_dp
        real(dp) :: length = 0.01_dp
    end type poincare_params_t

    type :: step_record_t
        integer :: n = 0
        integer :: region(MAX_SEGMENTS) = 0
        real(dp) :: z_start(4, MAX_SEGMENTS) = 0.0_dp
        real(dp) :: t(MAX_SEGMENTS) = 0.0_dp
        real(dp) :: kick(MAX_SEGMENTS) = 0.0_dp
    end type step_record_t

contains

    subroutine midpoint_step(par, region, z0, h, z1, ok, at_midpoint)
        !> z1 = z0 + h J grad H((z0 + z1)/2) by Newton iteration. With
        !> at_midpoint the branch is the region of the midpoint instead of
        !> region.
        type(field_params_t), intent(in) :: par
        integer, intent(in) :: region
        real(dp), intent(in) :: z0(4), h
        real(dp), intent(out) :: z1(4)
        logical, intent(out) :: ok
        logical, intent(in), optional :: at_midpoint

        real(dp) :: g(4), hs(4, 4), res(4), jac(4, 4), dz(4), zm(4)
        integer :: it, reg
        logical :: geometric

        geometric = .false.
        if (present(at_midpoint)) geometric = at_midpoint
        call grad_hamiltonian(par, region, z0, g)
        z1 = z0 + h*j_times(g)
        dz = 0.0_dp
        ok = .false.
        do it = 1, MAX_NEWTON
            zm = 0.5_dp*(z0 + z1)
            reg = region
            if (geometric) reg = region_of(par, zm(1:2))
            call grad_hamiltonian(par, reg, zm, g)
            call hess_hamiltonian(par, reg, zm, hs)
            res = z1 - z0 - h*j_times(g)
            jac = identity4() - 0.5_dp*h*j_times_mat(hs)
            call solve4(jac, -res, dz)
            z1 = z1 + dz
            if (maxval(abs(dz)/max(abs(z1), 1.0_dp)) <= NEWTON_TOL) then
                ok = .true.
                return
            end if
        end do
        ok = maxval(abs(dz)/max(abs(z1), 1.0_dp)) <= 1.0e3_dp*NEWTON_TOL
    end subroutine midpoint_step

    subroutine poincare_step(par, pp, z0, e0, dtau, z1, dt, ok)
        !> Implicit midpoint step of length dtau for K = g(z) (H(z) - e0). K is
        !> smooth for a smooth field, so the step is exactly symplectic; on
        !> K = 0 its flow is the physical flow with dt/dtau = g. dt is the
        !> midpoint quadrature of that physical time.
        type(field_params_t), intent(in) :: par
        type(poincare_params_t), intent(in) :: pp
        real(dp), intent(in) :: z0(4), e0, dtau
        real(dp), intent(out) :: z1(4), dt
        logical, intent(out) :: ok

        real(dp) :: gk(4), hk(4, 4), res(4), jac(4, 4), dz(4), zm(4), g
        integer :: it

        if (par%width <= 0.0_dp) error stop 'poincare_step needs a smooth layer'
        call grad_hess_k(par, pp, z0, e0, gk, hk, g)
        z1 = z0 + dtau*j_times(gk)
        dz = 0.0_dp
        ok = .false.
        do it = 1, MAX_NEWTON
            zm = 0.5_dp*(z0 + z1)
            call grad_hess_k(par, pp, zm, e0, gk, hk, g)
            res = z1 - z0 - dtau*j_times(gk)
            jac = identity4() - 0.5_dp*dtau*j_times_mat(hk)
            call solve4(jac, -res, dz)
            z1 = z1 + dz
            if (maxval(abs(dz)/max(abs(z1), 1.0_dp)) <= NEWTON_TOL) then
                ok = .true.
                exit
            end if
        end do
        if (.not. ok) ok = maxval(abs(dz)/max(abs(z1), 1.0_dp)) &
                           <= 1.0e3_dp*NEWTON_TOL
        call time_scale(par, pp, 0.5_dp*(z0(1:2) + z1(1:2)), g)
        dt = dtau*g
    end subroutine poincare_step

    subroutine poincare_gauss2_step(par, pp, z0, e0, dtau, z1, dt, ok)
        !> Two-stage Gauss-Legendre step (order 4, symplectic) of length dtau
        !> for K = g(z) (H(z) - e0). Its energy error is O(dtau^4), so the
        !> residual K = k acts far more weakly as the spurious potential k/g.
        !> pp%g_min = 1 gives the plain Gauss-Legendre step for H.
        type(field_params_t), intent(in) :: par
        type(poincare_params_t), intent(in) :: pp
        real(dp), intent(in) :: z0(4), e0, dtau
        real(dp), intent(out) :: z1(4), dt
        logical, intent(out) :: ok

        real(dp), parameter :: SQ3 = 1.7320508075688772_dp
        real(dp), parameter :: A(2, 2) = reshape([0.25_dp, &
            0.25_dp + SQ3/6.0_dp, 0.25_dp - SQ3/6.0_dp, 0.25_dp], [2, 2])
        real(dp) :: y(4, 2), f(4, 2), gk(4), hk(4, 4), jac(8, 8), res(8)
        real(dp) :: dy(8), g(2), jhk(4, 4)
        integer :: it, i, j

        if (par%width <= 0.0_dp) error stop 'poincare_gauss2_step needs a layer'
        call grad_hess_k(par, pp, z0, e0, gk, hk, g(1))
        y(:, 1) = z0 + 0.5_dp*dtau*j_times(gk)*(0.5_dp - SQ3/6.0_dp)*2.0_dp
        y(:, 2) = z0 + 0.5_dp*dtau*j_times(gk)*(0.5_dp + SQ3/6.0_dp)*2.0_dp
        dy = 0.0_dp
        ok = .false.
        ! Fixed-point iteration on the stages: gradients only, contractive for
        ! dtau |J Hess K| < 1. Newton below is the fallback.
        do it = 1, MAX_FIXED_POINT
            do j = 1, 2
                call grad_k(par, pp, y(:, j), e0, gk, g(j))
                f(:, j) = j_times(gk)
            end do
            do i = 1, 2
                dy(4*i - 3:4*i) = z0 + dtau*(A(i, 1)*f(:, 1) + A(i, 2)*f(:, 2)) &
                                  - y(:, i)
            end do
            y(:, 1) = y(:, 1) + dy(1:4)
            y(:, 2) = y(:, 2) + dy(5:8)
            if (maxval(abs(dy)) <= NEWTON_TOL*max(maxval(abs(y)), 1.0_dp)) then
                ok = .true.
                exit
            end if
        end do
        if (ok) then
            do j = 1, 2
                call grad_k(par, pp, y(:, j), e0, gk, g(j))
                f(:, j) = j_times(gk)
            end do
            z1 = z0 + 0.5_dp*dtau*(f(:, 1) + f(:, 2))
            dt = 0.5_dp*dtau*(g(1) + g(2))
            return
        end if
        do it = 1, MAX_NEWTON
            jac = 0.0_dp
            do j = 1, 2
                call grad_hess_k(par, pp, y(:, j), e0, gk, hk, g(j))
                f(:, j) = j_times(gk)
                jhk = j_times_mat(hk)
                do i = 1, 2
                    jac(4*i - 3:4*i, 4*j - 3:4*j) = -dtau*A(i, j)*jhk
                end do
            end do
            do i = 1, 8
                jac(i, i) = jac(i, i) + 1.0_dp
            end do
            do i = 1, 2
                res(4*i - 3:4*i) = y(:, i) - z0 &
                    - dtau*(A(i, 1)*f(:, 1) + A(i, 2)*f(:, 2))
            end do
            call solve_n(8, jac, -res, dy)
            y(:, 1) = y(:, 1) + dy(1:4)
            y(:, 2) = y(:, 2) + dy(5:8)
            if (maxval(abs(dy)) <= NEWTON_TOL*max(maxval(abs(y)), 1.0_dp)) then
                ok = .true.
                exit
            end if
        end do
        do j = 1, 2
            call grad_hess_k(par, pp, y(:, j), e0, gk, hk, g(j))
            f(:, j) = j_times(gk)
        end do
        z1 = z0 + 0.5_dp*dtau*(f(:, 1) + f(:, 2))
        dt = 0.5_dp*dtau*(g(1) + g(2))
    end subroutine poincare_gauss2_step

    subroutine boris_step(par, x, v, h)
        !> Boris push in Cartesian (x, y, z) with the field of the local region:
        !> half-step drift, velocity rotation about B(x), half-step drift. It
        !> preserves phase-space volume and |v| exactly, not p_phi.
        type(field_params_t), intent(in) :: par
        real(dp), intent(inout) :: x(3), v(3)
        real(dp), intent(in) :: h

        real(dp) :: rr, c, s, bc(3), b(3), t(3), sv(3), vp(3)

        x = x + 0.5_dp*h*v
        rr = hypot(x(1), x(2))
        c = x(1)/rr
        s = x(2)/rr
        call bfield(par, [rr, x(3)], bc)
        b = [bc(1)*c - bc(2)*s, bc(1)*s + bc(2)*c, bc(3)]
        t = 0.5_dp*h*b
        sv = 2.0_dp*t/(1.0_dp + dot_product(t, t))
        vp = v + cross3(v, t)
        v = v + cross3(vp, sv)
        x = x + 0.5_dp*h*v
    end subroutine boris_step

    pure function cross3(a, b) result(c)
        real(dp), intent(in) :: a(3), b(3)
        real(dp) :: c(3)

        c = [a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), a(1)*b(2) - a(2)*b(1)]
    end function cross3

    subroutine grad_k(par, pp, z, e0, gk, g)
        !> Gradient of K = g (H - e0).
        type(field_params_t), intent(in) :: par
        type(poincare_params_t), intent(in) :: pp
        real(dp), intent(in) :: z(4), e0
        real(dp), intent(out) :: gk(4), g

        real(dp) :: gh(4), gq(2)

        call grad_hamiltonian(par, 1, z, gh)
        call time_scale(par, pp, z(1:2), g, gq)
        gk = g*gh
        gk(1:2) = gk(1:2) + (hamiltonian(par, 1, z) - e0)*gq
    end subroutine grad_k

    subroutine grad_hess_k(par, pp, z, e0, gk, hk, g)
        !> Gradient and Hessian of K = g (H - e0).
        type(field_params_t), intent(in) :: par
        type(poincare_params_t), intent(in) :: pp
        real(dp), intent(in) :: z(4), e0
        real(dp), intent(out) :: gk(4), hk(4, 4), g

        real(dp) :: gh(4), hh(4, 4), dh, gg(4), hg(4, 4)
        integer :: i, j

        dh = hamiltonian(par, 1, z) - e0
        call grad_hamiltonian(par, 1, z, gh)
        call hess_hamiltonian(par, 1, z, hh)
        call time_scale(par, pp, z(1:2), g, gg(1:2), hg(1:2, 1:2))
        gg(3:4) = 0.0_dp
        hg(3:4, :) = 0.0_dp
        hg(:, 3:4) = 0.0_dp
        gk = g*gh + dh*gg
        do j = 1, 4
            do i = 1, 4
                hk(i, j) = g*hh(i, j) + gg(i)*gh(j) + gh(i)*gg(j) + dh*hg(i, j)
            end do
        end do
    end subroutine grad_hess_k

    pure subroutine time_scale(par, pp, q, g, gq, gqq)
        !> g(q), its gradient and Hessian in (R, Z).
        type(field_params_t), intent(in) :: par
        type(poincare_params_t), intent(in) :: pp
        real(dp), intent(in) :: q(2)
        real(dp), intent(out) :: g
        real(dp), intent(out), optional :: gq(2), gqq(2, 2)

        real(dp) :: s, e, d1, d2, n(2), hs(2, 2)

        s = sheet_distance(par, q)
        e = (1.0_dp - pp%g_min)*exp(-(s/pp%length)**2)
        g = 1.0_dp - e
        if (.not. present(gq)) return
        d1 = 2.0_dp*s/pp%length**2*e
        d2 = (2.0_dp/pp%length**2 - 4.0_dp*s*s/pp%length**4)*e
        call grad_sheet_distance(par, q, n)
        call sheet_hessian(par, q, hs)
        gq = d1*n
        if (.not. present(gqq)) return
        gqq(1, 1) = d2*n(1)*n(1) + d1*hs(1, 1)
        gqq(1, 2) = d2*n(1)*n(2) + d1*hs(1, 2)
        gqq(2, 1) = gqq(1, 2)
        gqq(2, 2) = d2*n(2)*n(2) + d1*hs(2, 2)
    end subroutine time_scale

    subroutine sheet_step(par, scheme, region, z0, h, z1, stats, record)
        !> One step of length h from z0 in the given region. On return region
        !> is the region of z1 and record (if present) holds the segments.
        type(field_params_t), intent(in) :: par
        integer, intent(in) :: scheme
        integer, intent(inout) :: region
        real(dp), intent(in) :: z0(4), h
        real(dp), intent(out) :: z1(4)
        type(step_stats_t), intent(inout) :: stats
        type(step_record_t), intent(out), optional :: record

        type(step_record_t) :: rec
        real(dp) :: z(4), t_rem, t_land, z_land(4), e_prev, kick
        logical :: ok, on_sheet, lands
        integer :: nseg

        rec%n = 0
        if (scheme == SCHEME_MIDPOINT_REGION) then
            call midpoint_step(par, region, z0, h, z1, ok, at_midpoint=.true.)
            if (.not. ok) stats%failures = stats%failures + 1
            if (region_of(par, z1(1:2)) /= region) &
                stats%crossings = stats%crossings + 1
            region = region_of(par, z1(1:2))
            call push(rec, region, z0, h, 0.0_dp)
            if (present(record)) record = rec
            return
        end if
        z = z0
        t_rem = h
        on_sheet = .false.
        nseg = 0
        do
            call midpoint_step(par, region, z, t_rem, z1, ok)
            if (.not. ok) stats%failures = stats%failures + 1
            lands = region_of(par, z1(1:2)) /= region
            if (lands .and. scheme == SCHEME_NAIVE) then
                call push(rec, region, z, t_rem, 0.0_dp)
                region = region_of(par, z1(1:2))
                stats%crossings = stats%crossings + 1
                exit
            end if
            if (lands) then
                call land_on_sheet(par, region, z, t_rem, on_sheet, t_land, &
                    z_land, ok)
                if (.not. ok .or. nseg + 1 >= MAX_SEGMENTS) then
                    stats%failures = stats%failures + 1
                    lands = .false.
                end if
            end if
            if (.not. lands) then
                call push(rec, region, z, t_rem, 0.0_dp)
                region = region_of(par, z1(1:2))
                exit
            end if
            e_prev = hamiltonian(par, region, 0.5_dp*(z + z_land))
            call match_kick(par, 3 - region, z_land, t_rem - t_land, e_prev, &
                kick, ok)
            if (.not. ok) stats%failures = stats%failures + 1
            call push(rec, region, z, t_land, kick)
            nseg = nseg + 1
            stats%crossings = stats%crossings + 1
            stats%max_kick = max(stats%max_kick, abs(kick))
            z = kicked(par, z_land, kick)
            t_rem = t_rem - t_land
            region = 3 - region
            on_sheet = .true.
        end do
        if (nseg > 1) stats%multi_crossing_steps = stats%multi_crossing_steps + 1
        if (present(record)) record = rec
    end subroutine sheet_step

    subroutine push(rec, region, z, t, kick)
        type(step_record_t), intent(inout) :: rec
        integer, intent(in) :: region
        real(dp), intent(in) :: z(4), t, kick

        rec%n = rec%n + 1
        rec%region(rec%n) = region
        rec%z_start(:, rec%n) = z
        rec%t(rec%n) = t
        rec%kick(rec%n) = kick
    end subroutine push

    pure function kicked(par, z, kick) result(zk)
        type(field_params_t), intent(in) :: par
        real(dp), intent(in) :: z(4), kick
        real(dp) :: zk(4)

        real(dp) :: gs(2)

        call grad_sheet_distance(par, z(1:2), gs)
        zk = z
        zk(3:4) = zk(3:4) - kick*gs
    end function kicked

    subroutine land_on_sheet(par, region, z0, t_max, from_sheet, t_land, &
            z_land, ok)
        !> Substep length t_land in (0, t_max) with s(z(t_land)) = 0, by
        !> Illinois false position on the midpoint solve. A start on the sheet
        !> brackets from a lift-off time instead of t = 0.
        type(field_params_t), intent(in) :: par
        integer, intent(in) :: region
        real(dp), intent(in) :: z0(4), t_max
        logical, intent(in) :: from_sheet
        real(dp), intent(out) :: t_land, z_land(4)
        logical, intent(out) :: ok

        real(dp) :: a, b, fa, fb, c, fc
        integer :: it, side
        logical :: step_ok

        ok = .false.
        t_land = t_max
        if (from_sheet) then
            a = LIFTOFF*t_max
            call midpoint_step(par, region, z0, a, z_land, step_ok)
            fa = sheet_distance(par, z_land(1:2))
        else
            a = 0.0_dp
            fa = sheet_distance(par, z0(1:2))
        end if
        b = t_max
        call midpoint_step(par, region, z0, b, z_land, step_ok)
        fb = sheet_distance(par, z_land(1:2))
        if (fa*fb > 0.0_dp) return
        side = 0
        do it = 1, 200
            c = (a*fb - b*fa)/(fb - fa)
            call midpoint_step(par, region, z0, c, z_land, step_ok)
            fc = sheet_distance(par, z_land(1:2))
            t_land = c
            if (abs(fc) <= LANDING_TOL) then
                ok = step_ok
                return
            end if
            if (fc*fb > 0.0_dp) then
                b = c
                fb = fc
                if (side == -1) fa = 0.5_dp*fa
                side = -1
            else
                a = c
                fa = fc
                if (side == 1) fb = 0.5_dp*fb
                side = 1
            end if
        end do
    end subroutine land_on_sheet

    subroutine match_kick(par, region, z_land, t_rem, e_prev, kick, ok)
        !> Kick that gives the segment after the landing the discrete energy
        !> e_prev. That segment runs in region from the kicked state to the
        !> next landing, or to the end of the step if there is none.
        type(field_params_t), intent(in) :: par
        integer, intent(in) :: region
        real(dp), intent(in) :: z_land(4), t_rem, e_prev
        real(dp), intent(out) :: kick
        logical, intent(out) :: ok

        real(dp) :: k0, k1, k2, g0, g1
        integer :: it

        k0 = 0.0_dp
        g0 = next_energy(k0) - e_prev
        k1 = 1.0e-9_dp
        g1 = next_energy(k1) - e_prev
        do it = 1, 60
            if (abs(g1) <= tiny(1.0_dp) .or. abs(g1 - g0) <= tiny(1.0_dp)) exit
            k2 = k1 - g1*(k1 - k0)/(g1 - g0)
            k0 = k1
            g0 = g1
            k1 = k2
            g1 = next_energy(k1) - e_prev
            ! H is a small difference of large momenta, so the residual stalls
            ! at round-off; converge on the kick update instead.
            if (abs(k1 - k0) <= KICK_TOL) exit
        end do
        kick = k1
        ok = abs(g1) <= ENERGY_TOL*abs(e_prev)

    contains

        real(dp) function next_energy(k) result(e)
            real(dp), intent(in) :: k
            real(dp) :: zk(4), z_end(4), t_next
            logical :: step_ok

            zk = kicked(par, z_land, k)
            call midpoint_step(par, region, zk, t_rem, z_end, step_ok)
            if (region_of(par, z_end(1:2)) /= region) then
                call land_on_sheet(par, region, zk, t_rem, .true., t_next, &
                    z_end, step_ok)
            end if
            e = hamiltonian(par, region, 0.5_dp*(zk + z_end))
        end function next_energy

    end subroutine match_kick

    subroutine tangent_map(par, rec, m)
        !> Exact Jacobian dz1/dz0 of a recorded landing-scheme step of fixed
        !> total length. Column 5 of the work arrays carries the derivative with
        !> respect to the pending kick, which the energy match of the next
        !> segment eliminates.
        type(field_params_t), intent(in) :: par
        type(step_record_t), intent(in) :: rec
        real(dp), intent(out) :: m(4, 4)

        real(dp) :: dz(4, 5), dt_rem(5), dz_end(4, 5), dt(5), de(5)
        real(dp) :: de_prev(5), p(4, 4), q(4), z0(4), z_end(4), zbar(4)
        real(dp) :: n4(4), g(4), hs_s(2, 2), dl(4), kick
        integer :: k, reg, c
        logical :: ok

        dz = 0.0_dp
        do c = 1, 4
            dz(c, c) = 1.0_dp
        end do
        dt_rem = 0.0_dp
        de_prev = 0.0_dp
        dz_end = 0.0_dp
        do k = 1, rec%n
            reg = rec%region(k)
            z0 = rec%z_start(:, k)
            call midpoint_step(par, reg, z0, rec%t(k), z_end, ok)
            call midpoint_flow_derivatives(par, reg, z0, z_end, rec%t(k), p, q)
            if (k < rec%n) then
                ! Landing: s(z_end) = 0 fixes the segment time.
                n4 = sheet_normal4(par, z_end)
                do c = 1, 5
                    dt(c) = -dot_product(n4, matmul(p, dz(:, c))) &
                            /dot_product(n4, q)
                end do
            else
                dt = dt_rem
            end if
            do c = 1, 5
                dz_end(:, c) = matmul(p, dz(:, c)) + q*dt(c)
            end do
            dt_rem = dt_rem - dt
            zbar = 0.5_dp*(z0 + z_end)
            call grad_hamiltonian(par, reg, zbar, g)
            do c = 1, 5
                de(c) = 0.5_dp*dot_product(g, dz(:, c) + dz_end(:, c))
            end do
            if (k > 1) then
                ! Energy match with the previous segment fixes the pending kick.
                dl = (de_prev(1:4) - de(1:4))/de(5)
                do c = 1, 4
                    dz_end(:, c) = dz_end(:, c) + dz_end(:, 5)*dl(c)
                    de(c) = de(c) + de(5)*dl(c)
                end do
                dt_rem(1:4) = dt_rem(1:4) + dt_rem(5)*dl
            end if
            de_prev = de
            if (k == rec%n) exit
            ! Kick p -> p - l grad_q s(q) at the landing point.
            kick = rec%kick(k)
            call sheet_hessian(par, z_end(1:2), hs_s)
            n4 = sheet_normal4(par, z_end)
            do c = 1, 4
                dz(:, c) = dz_end(:, c)
                dz(3:4, c) = dz(3:4, c) - kick*matmul(hs_s, dz_end(1:2, c))
            end do
            dz(:, 5) = [0.0_dp, 0.0_dp, -n4(1), -n4(2)]
            dt_rem(5) = 0.0_dp
        end do
        m = dz_end(:, 1:4)
    end subroutine tangent_map

    subroutine midpoint_flow_derivatives(par, region, z0, z1, t, p, q)
        !> For z1 = z0 + t J grad H(zbar): dz1 = p dz0 + q dt with
        !> p = (I - t/2 J Hess)^-1 (I + t/2 J Hess), q = (I - t/2 J Hess)^-1 J grad H.
        type(field_params_t), intent(in) :: par
        integer, intent(in) :: region
        real(dp), intent(in) :: z0(4), z1(4), t
        real(dp), intent(out) :: p(4, 4), q(4)

        real(dp) :: hs(4, 4), g(4), jh(4, 4), a(4, 4), b(4, 4)
        integer :: c

        call hess_hamiltonian(par, region, 0.5_dp*(z0 + z1), hs)
        call grad_hamiltonian(par, region, 0.5_dp*(z0 + z1), g)
        jh = j_times_mat(hs)
        a = identity4() - 0.5_dp*t*jh
        b = identity4() + 0.5_dp*t*jh
        do c = 1, 4
            call solve4(a, b(:, c), p(:, c))
        end do
        call solve4(a, j_times(g), q)
    end subroutine midpoint_flow_derivatives

    pure function sheet_normal4(par, z) result(n4)
        !> grad s as a phase-space covector (dR, dZ, dp_R, dp_Z).
        type(field_params_t), intent(in) :: par
        real(dp), intent(in) :: z(4)
        real(dp) :: n4(4)

        real(dp) :: gs(2)

        call grad_sheet_distance(par, z(1:2), gs)
        n4 = [gs(1), gs(2), 0.0_dp, 0.0_dp]
    end function sheet_normal4

    pure subroutine sheet_hessian(par, q, hs)
        !> Hessian of s = r - r0: (I - n n^T)/r.
        type(field_params_t), intent(in) :: par
        real(dp), intent(in) :: q(2)
        real(dp), intent(out) :: hs(2, 2)

        real(dp) :: n(2), r

        call grad_sheet_distance(par, q, n)
        r = hypot(q(1) - par%r_major, q(2))
        hs(1, 1) = (1.0_dp - n(1)*n(1))/r
        hs(1, 2) = -n(1)*n(2)/r
        hs(2, 1) = hs(1, 2)
        hs(2, 2) = (1.0_dp - n(2)*n(2))/r
    end subroutine sheet_hessian

    pure function j_times(g) result(v)
        !> J g for z = (q, p): (dH/dp, -dH/dq).
        real(dp), intent(in) :: g(4)
        real(dp) :: v(4)

        v = [g(3), g(4), -g(1), -g(2)]
    end function j_times

    pure function j_times_mat(a) result(b)
        real(dp), intent(in) :: a(4, 4)
        real(dp) :: b(4, 4)

        b(1, :) = a(3, :)
        b(2, :) = a(4, :)
        b(3, :) = -a(1, :)
        b(4, :) = -a(2, :)
    end function j_times_mat

    pure function identity4() result(a)
        real(dp) :: a(4, 4)
        integer :: i

        a = 0.0_dp
        do i = 1, 4
            a(i, i) = 1.0_dp
        end do
    end function identity4

    pure subroutine solve_n(n, a, b, x)
        !> Gaussian elimination with partial pivoting.
        integer, intent(in) :: n
        real(dp), intent(in) :: a(n, n), b(n)
        real(dp), intent(out) :: x(n)

        real(dp) :: m(n, n + 1), f
        integer :: k, i, p

        m(:, 1:n) = a
        m(:, n + 1) = b
        do k = 1, n
            p = maxloc(abs(m(k:n, k)), dim=1) + k - 1
            if (p /= k) m([k, p], :) = m([p, k], :)
            do i = k + 1, n
                f = m(i, k)/m(k, k)
                m(i, k:n + 1) = m(i, k:n + 1) - f*m(k, k:n + 1)
            end do
        end do
        x = 0.0_dp
        do k = n, 1, -1
            x(k) = (m(k, n + 1) - dot_product(m(k, k + 1:n), x(k + 1:n)))/m(k, k)
        end do
    end subroutine solve_n

    pure subroutine solve4(a, b, x)
        !> Gaussian elimination with partial pivoting for a 4x4 system.
        real(dp), intent(in) :: a(4, 4), b(4)
        real(dp), intent(out) :: x(4)

        real(dp) :: m(4, 5), f
        integer :: k, i, p

        m(:, 1:4) = a
        m(:, 5) = b
        do k = 1, 4
            p = maxloc(abs(m(k:4, k)), dim=1) + k - 1
            if (p /= k) m([k, p], :) = m([p, k], :)
            do i = k + 1, 4
                f = m(i, k)/m(k, k)
                m(i, k:5) = m(i, k:5) - f*m(k, k:5)
            end do
        end do
        x = 0.0_dp
        do k = 4, 1, -1
            x(k) = (m(k, 5) - dot_product(m(k, k + 1:4), x(k + 1:4)))/m(k, k)
        end do
    end subroutine solve4

end module fo_sheet
