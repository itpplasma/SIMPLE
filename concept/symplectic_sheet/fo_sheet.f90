module fo_sheet
    !> Resolved smooth-layer Hamiltonian controls. Sharp-sheet branch and
    !> energy-kick landing experiments were removed; see the rejection ledger.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use sheet_field, only: field_params_t, hamiltonian, grad_hamiltonian, &
        hess_hamiltonian, sheet_distance, grad_sheet_distance
    implicit none
    private

    public :: poincare_step, poincare_gauss2_step, poincare_params_t
    integer, parameter :: MAX_NEWTON = 50, MAX_FIXED_POINT = 60
    real(dp), parameter :: NEWTON_TOL = 1.0e-15_dp

    !> Poincare time transformation dt = g(z) dtau with
    !> g = 1 - (1 - g_min) exp(-((r - r0)/length)^2): fixed steps in the
    !> fictitious time tau are short in physical time near the sheet.
    type :: poincare_params_t
        real(dp) :: g_min = 0.05_dp
        real(dp) :: length = 0.01_dp
    end type poincare_params_t

contains

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
