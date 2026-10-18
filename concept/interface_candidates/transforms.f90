module transforms
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use candidate_field, only: model_t, local_t, geometry, cross, frame, side_of
    implicit none
    private
    public :: lift, inverse, gc_energy, generators, kernel_available
    public :: register_order2, order2_rhs

    abstract interface
        subroutine high_kernel(pm, y, result, ok)
            import :: model_t, dp
            type(model_t), intent(in) :: pm
            real(dp), intent(in) :: y(6)
            real(dp), intent(out) :: result(6)
            logical, intent(out) :: ok
        end subroutine
        function high_energy(pm, g) result(e)
            import :: model_t, dp
            type(model_t), intent(in) :: pm
            real(dp), intent(in) :: g(6)
            real(dp) :: e
        end function
    end interface
    procedure(high_kernel), pointer :: order2_lift => null()
    procedure(high_kernel), pointer :: order2_rhs => null()
    procedure(high_energy), pointer :: order2_energy => null()

contains

    subroutine register_order2(t2, f2, h2)
        procedure(high_kernel) :: t2, f2
        procedure(high_energy) :: h2
        ! Registration is serial, before the OpenMP region. Both kernels must
        ! use this module's (X[3], u, mu, velocity-phase) and q=m=1 convention.
        order2_lift => t2
        order2_rhs => f2
        order2_energy => h2
    end subroutine

    logical function kernel_available(order) result(available)
        integer, intent(in) :: order
        available = order == 0 .or. order == 1
        if (order == 2) available = associated(order2_lift) .and. &
                                   associated(order2_rhs) .and. &
                                   associated(order2_energy)
    end function

    function gc_energy(pm, g, order) result(e)
        type(model_t), intent(in) :: pm
        real(dp), intent(in) :: g(6)
        integer, intent(in), optional :: order
        real(dp) :: e
        type(local_t) :: lf
        if (present(order)) then
            if (order == 2) then
                if (.not. associated(order2_energy)) error stop 'T2 energy missing'
                e = order2_energy(pm, g)
                return
            end if
        end if
        call geometry(pm, g(1:3), side_of(pm, g(1:3)), lf)
        e = 0.5_dp*g(4)**2 + g(5)*lf%bmag
    end function

    subroutine generators(lf, u, mu, phase, rho, gu, gmu, gz, rho1)
        type(local_t), intent(in) :: lf
        real(dp), intent(in) :: u, mu, phase
        real(dp), intent(out) :: rho(3), gu, gmu, gz, rho1(3)
        real(dp) :: e1(3), e2(3), w(3), rh(3), amp, a1, a2
        real(dp) :: jr(3, 3), d1(3), d2(3), dw(3), tmp(3), g2(3), g1dr(3)
        real(dp) :: hnorm
        integer :: j
        call frame(lf%b, e1, e2)
        w = cos(phase)*e1 + sin(phase)*e2
        rh = cross(lf%b, w)
        amp = sqrt(2.0_dp*mu/lf%bmag)
        rho = amp*rh
        a1 = -0.5_dp*(dot_product(rh, matmul(lf%db, w)) &
                      + dot_product(w, matmul(lf%db, rh)))
        a2 = 0.25_dp*(dot_product(w, matmul(lf%db, w)) &
                      - dot_product(rh, matmul(lf%db, rh)))
        gu = -u*dot_product(rho, lf%kappa) + mu*(lf%tau + a1)
        gmu = dot_product(rho, mu*lf%gradb + u*u*lf%kappa)/lf%bmag &
              - u*mu*(lf%tau + a1)/lf%bmag
        gz = -dot_product(rho, lf%gauge) + u*a2/lf%bmag &
             + amp*dot_product(w, lf%gradb &
               + u*u*lf%kappa/(2.0_dp*mu))/lf%bmag
        hnorm = norm2(cross(lf%b, [0.0_dp, 0.0_dp, 1.0_dp]))
        do j = 1, 3
            tmp = cross(lf%db(:, j), [0.0_dp, 0.0_dp, 1.0_dp])
            d1 = (tmp - e1*dot_product(e1, tmp))/hnorm
            d2 = cross(lf%db(:, j), e1) + cross(lf%b, d1)
            dw = cos(phase)*d1 + sin(phase)*d2
            jr(:, j) = -0.5_dp*rho*lf%gradb(j)/lf%bmag &
                + amp*(cross(lf%db(:, j), w) + cross(lf%b, dw))
        end do
        ! zeta=-phase-pi/2, so d rho/d zeta=amp*w. Littlejohn's
        ! Pi_1=-mu*tau*b/2 has Pi_1 cross b=0.
        g2 = (2.0_dp*u*amp*dot_product(w, lf%kappa)/lf%bmag &
              + mu*a2/lf%bmag)*lf%b &
            + 0.5_dp*(u*u*dot_product(rho, lf%kappa)/lf%bmag &
              + mu*u*(3.0_dp*lf%tau - a1)/lf%bmag)*rho/(2.0_dp*mu) &
            + 0.5_dp*(u*a2/lf%bmag + amp*dot_product(w, &
              lf%gradb/lf%bmag + u*u*lf%kappa/(2.0_dp*lf%bmag*mu)))*amp*w
        g1dr = -matmul(jr, rho) + gmu*rho/(2.0_dp*mu) + gz*amp*w
        rho1 = -g2 - 0.5_dp*g1dr
    end subroutine

    subroutine lift(pm, g, order, y, ok)
        type(model_t), intent(in) :: pm
        real(dp), intent(in) :: g(6)
        integer, intent(in) :: order
        real(dp), intent(out) :: y(6)
        logical, intent(out) :: ok
        type(local_t) :: lf, lp
        real(dp) :: e1(3), e2(3), v(3), rho(3), rho1(3), gu, gmu, gz
        real(dp) :: up, mup, phasep, e
        ok = .false.
        y = 0.0_dp
        if (order == 2) then
            if (.not. kernel_available(order)) return
            call order2_lift(pm, g, y, ok)
            return
        end if
        if (order < 0 .or. order > 1 .or. g(5) <= 0.0_dp) return
        call geometry(pm, g(1:3), side_of(pm, g(1:3)), lf)
        call frame(lf%b, e1, e2)
        v = g(4)*lf%b + sqrt(2.0_dp*g(5)*lf%bmag)* &
            (cos(g(6))*e1 + sin(g(6))*e2)
        rho = cross(lf%b, v)/lf%bmag
        y(1:3) = g(1:3) + rho
        if (order == 0) then
            y(4:6) = v
        else
            call generators(lf, g(4), g(5), g(6), rho, gu, gmu, gz, rho1)
            y(1:3) = g(1:3) + rho + rho1
            up = g(4) - gu
            mup = g(5) - gmu
            phasep = g(6) + gz
            if (mup <= 0.0_dp) return
            call geometry(pm, y(1:3), side_of(pm, g(1:3)), lp)
            call frame(lp%b, e1, e2)
            v = up*lp%b + sqrt(2.0_dp*mup*lp%bmag)* &
                (cos(phasep)*e1 + sin(phasep)*e2)
            e = 0.5_dp*g(4)**2 + g(5)*lf%bmag
            ! Energy representation: an O(epsilon^2) completion, explicitly
            ! included in the one lift whose inverse is solved below.
            y(4:6) = v*sqrt(2.0_dp*e/dot_product(v, v))
        end if
        ok = .true.
    end subroutine

    subroutine inverse(pm, y, order, g, ok, residual)
        type(model_t), intent(in) :: pm
        real(dp), intent(in) :: y(6)
        integer, intent(in) :: order
        real(dp), intent(out) :: g(6)
        logical, intent(out) :: ok
        real(dp), intent(out) :: residual
        type(local_t) :: lf
        real(dp) :: e1(3), e2(3), wp(3), z(6), zp(6), zm(6), gp(6), gm(6)
        real(dp) :: out(6), op(6), om(6), jac(6, 6), res(6), dz(6), scales(6)
        real(dp) :: speed, hh, factor, best, nr
        logical :: good, solved, goodp, goodm
        integer :: it, j, ls
        speed = norm2(y(4:6))
        call geometry(pm, y(1:3), side_of(pm, y(1:3)), lf)
        g(1:3) = y(1:3) + cross(y(4:6), lf%b)/lf%bmag
        call geometry(pm, g(1:3), side_of(pm, g(1:3)), lf)
        g(4) = dot_product(y(4:6), lf%b)
        g(5) = max((speed**2 - g(4)**2)/(2.0_dp*lf%bmag), 1.0e-16_dp)
        wp = y(4:6) - g(4)*lf%b
        call frame(lf%b, e1, e2)
        g(6) = atan2(dot_product(wp, e2), dot_product(wp, e1))
        z = g
        z(5) = log(g(5))
        scales = [speed, speed, speed, speed, speed, speed]
        ok = .false.
        residual = huge(1.0_dp)
        do it = 1, 15
            call from_log(z, g)
            call lift(pm, g, order, out, good)
            if (.not. good) return
            res = (out - y)/scales
            residual = maxval(abs(res))
            if (residual < 2.0e-11_dp) then
                ok = .true.
                return
            end if
            do j = 1, 6
                hh = 2.0e-5_dp
                if (j == 4) hh = hh*speed
                zp = z
                zm = z
                zp(j) = zp(j) + hh
                zm(j) = zm(j) - hh
                call from_log(zp, gp)
                call from_log(zm, gm)
                call lift(pm, gp, order, op, goodp)
                call lift(pm, gm, order, om, goodm)
                if (.not. goodp .or. .not. goodm) return
                jac(:, j) = (op - om)/(2.0_dp*hh*scales)
            end do
            call linear_solve(jac, -res, dz, solved)
            if (.not. solved) return
            factor = 1.0_dp
            best = residual
            do ls = 1, 8
                zp = z + factor*dz
                if (abs(zp(5)) < 100.0_dp) then
                    call from_log(zp, gp)
                    call lift(pm, gp, order, op, goodp)
                    if (goodp) then
                        nr = maxval(abs((op - y)/scales))
                        if (nr < best) exit
                    end if
                end if
                factor = factor*0.5_dp
            end do
            z = zp
        end do
    contains
        subroutine from_log(zz, gg)
            real(dp), intent(in) :: zz(6)
            real(dp), intent(out) :: gg(6)
            gg = zz
            gg(5) = exp(zz(5))
        end subroutine
    end subroutine

    subroutine linear_solve(a0, b0, x, ok)
        real(dp), intent(in) :: a0(6, 6), b0(6)
        real(dp), intent(out) :: x(6)
        logical, intent(out) :: ok
        real(dp) :: a(6, 6), b(6), row(6), tmp, f
        integer :: i, j, k
        a = a0
        b = b0
        ok = .false.
        x = 0.0_dp
        do k = 1, 6
            j = maxloc(abs(a(k:6, k)), dim=1) + k - 1
            if (abs(a(j, k)) < 1.0e-20_dp) return
            row = a(k, :)
            a(k, :) = a(j, :)
            a(j, :) = row
            tmp = b(k)
            b(k) = b(j)
            b(j) = tmp
            do i = k + 1, 6
                f = a(i, k)/a(k, k)
                a(i, k:6) = a(i, k:6) - f*a(k, k:6)
                b(i) = b(i) - f*b(k)
            end do
        end do
        do i = 6, 1, -1
            x(i) = (b(i) - dot_product(a(i, i + 1:6), x(i + 1:6)))/a(i, i)
        end do
        ok = .true.
    end subroutine

end module
