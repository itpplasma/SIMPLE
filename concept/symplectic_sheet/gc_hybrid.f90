module gc_hybrid
    !> Guiding-centre segments away from the sheet and the switches between
    !> guiding centre (GC) and full orbit (FO) for the torus of sheet_field.
    !>
    !> GC equations (q = m = 1): B* = B + v_par curl b, B*_par = b . B*,
    !>     dX/dt = (v_par B* + mu b x grad B)/B*_par,
    !>     dv_par/dt = -mu B* . grad B/B*_par,   dzeta/dt = -|B|,
    !> in the orthonormal cylindrical frame (e_R, e_phi, e_Z); axisymmetry
    !> makes phi ignorable. zeta is the gyrophase of v_perp in the frame
    !> (e1, e2 = b x e1), so the gyroradius vector is rho = -(v x b)/|B|.
    !>
    !> FO -> GC is the first-order transform X = x + (v x B)/B^2 with mu chosen
    !> so that v_par^2/2 + mu B(X) equals the FO energy. GC -> FO lifts with a
    !> gyrophase: either the GC one (generic) or the one closest to it that
    !> keeps the exact invariant p_phi (axisymmetric only).
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use sheet_field, only: field_params_t, region_of, psi_derivatives, &
        psi_value, sheet_distance
    implicit none
    private

    public :: gc_state_t, fo_to_gc, gc_to_fo, gc_advance, fo_velocity
    public :: LIFT_GYROPHASE, LIFT_PPHI

    integer, parameter :: LIFT_GYROPHASE = 1, LIFT_PPHI = 2
    real(dp), parameter :: TWOPI = 6.283185307179586_dp

    type :: gc_state_t
        real(dp) :: r = 0.0_dp
        real(dp) :: z = 0.0_dp
        real(dp) :: vpar = 0.0_dp
        real(dp) :: zeta = 0.0_dp
        real(dp) :: mu = 0.0_dp
    end type gc_state_t

    type :: local_field_t
        real(dp) :: b(3), bmod, h(3), grad_bmod(3), curl_h(3)
    end type local_field_t

contains

    function local_field(par, q) result(lf)
        !> B, |B|, b, grad |B| and curl b at q = (R, Z), (R, phi, Z) components.
        type(field_params_t), intent(in) :: par
        real(dp), intent(in) :: q(2)
        type(local_field_t) :: lf

        real(dp) :: psi, gp(2), hp(2, 2), rr, db_r(3), db_z(3), dh_r(3), dh_z(3)

        call psi_derivatives(par, region_of(par, q), q, psi, gp, hp)
        rr = q(1)
        lf%b = [-gp(2)/rr, par%f_tor/rr, gp(1)/rr]
        db_r = [gp(2)/rr**2 - hp(2, 1)/rr, -par%f_tor/rr**2, &
                hp(1, 1)/rr - gp(1)/rr**2]
        db_z = [-hp(2, 2)/rr, 0.0_dp, hp(1, 2)/rr]
        lf%bmod = norm2(lf%b)
        lf%h = lf%b/lf%bmod
        lf%grad_bmod = [dot_product(lf%b, db_r), 0.0_dp, &
                        dot_product(lf%b, db_z)]/lf%bmod
        dh_r = db_r/lf%bmod - lf%b*lf%grad_bmod(1)/lf%bmod**2
        dh_z = db_z/lf%bmod - lf%b*lf%grad_bmod(3)/lf%bmod**2
        lf%curl_h = [-dh_z(2), dh_z(1) - dh_r(3), lf%h(2)/rr + dh_r(2)]
    end function local_field

    pure function cross(a, b) result(c)
        real(dp), intent(in) :: a(3), b(3)
        real(dp) :: c(3)

        c = [a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), a(1)*b(2) - a(2)*b(1)]
    end function cross

    pure subroutine frame(h, e1, e2)
        real(dp), intent(in) :: h(3)
        real(dp), intent(out) :: e1(3), e2(3)

        e1 = cross(h, [0.0_dp, 0.0_dp, 1.0_dp])
        e1 = e1/norm2(e1)
        e2 = cross(h, e1)
    end subroutine frame

    function fo_velocity(par, z) result(v)
        !> Physical (v_R, v_phi, v_Z) of the reduced FO state.
        type(field_params_t), intent(in) :: par
        real(dp), intent(in) :: z(4)
        real(dp) :: v(3)

        v = [z(3), (par%p_phi - psi_value(par, region_of(par, z(1:2)), &
             z(1:2)))/z(1), z(4) + par%f_tor*log(z(1))]
    end function fo_velocity

    function fo_to_gc(par, z) result(gc)
        type(field_params_t), intent(in) :: par
        real(dp), intent(in) :: z(4)
        type(gc_state_t) :: gc

        type(local_field_t) :: lf
        real(dp) :: v(3), d(3), e1(3), e2(3), vperp(3)

        v = fo_velocity(par, z)
        lf = local_field(par, z(1:2))
        d = cross(v, lf%b)/lf%bmod**2
        gc%r = z(1) + d(1)
        gc%z = z(2) + d(3)
        lf = local_field(par, [gc%r, gc%z])
        gc%vpar = dot_product(v, lf%h)
        gc%mu = (0.5_dp*dot_product(v, v) - 0.5_dp*gc%vpar**2)/lf%bmod
        vperp = v - gc%vpar*lf%h
        call frame(lf%h, e1, e2)
        gc%zeta = atan2(dot_product(vperp, e2), dot_product(vperp, e1))
    end function fo_to_gc

    subroutine gc_to_fo(par, gc, lift, pphi_target, z, pphi, ok)
        !> Particle state of a GC. LIFT_PPHI picks the gyrophase nearest to
        !> gc%zeta at which the FO p_phi equals pphi_target; ok = .false. if
        !> there is none (then the GC gyrophase is used).
        type(field_params_t), intent(in) :: par
        type(gc_state_t), intent(in) :: gc
        integer, intent(in) :: lift
        real(dp), intent(in) :: pphi_target
        real(dp), intent(out) :: z(4), pphi
        logical, intent(out) :: ok

        integer, parameter :: NSCAN = 64
        real(dp) :: zeta, a, b, fa, fb, c, fc, best, d
        integer :: k, it

        ok = .true.
        zeta = gc%zeta
        if (lift == LIFT_PPHI) then
            ok = .false.
            best = huge(1.0_dp)
            do k = 0, NSCAN - 1
                a = gc%zeta + TWOPI*real(k, dp)/NSCAN
                b = a + TWOPI/NSCAN
                fa = lifted_pphi(a) - pphi_target
                fb = lifted_pphi(b) - pphi_target
                if (fa*fb > 0.0_dp) cycle
                do it = 1, 60
                    c = 0.5_dp*(a + b)
                    fc = lifted_pphi(c) - pphi_target
                    if (fa*fc <= 0.0_dp) then
                        b = c
                    else
                        a = c
                        fa = fc
                    end if
                end do
                d = abs(modulo(c - gc%zeta + 0.5_dp*TWOPI, TWOPI) - 0.5_dp*TWOPI)
                if (d < best) then
                    best = d
                    zeta = c
                    ok = .true.
                end if
            end do
        end if
        call lifted_state(zeta, z, pphi)

    contains

        real(dp) function lifted_pphi(zt) result(p)
            real(dp), intent(in) :: zt
            real(dp) :: zz(4)

            call lifted_state(zt, zz, p)
        end function lifted_pphi

        subroutine lifted_state(zt, zz, p)
            real(dp), intent(in) :: zt
            real(dp), intent(out) :: zz(4), p

            type(local_field_t) :: lf
            real(dp) :: e1(3), e2(3), v(3), x(3), vperp

            lf = local_field(par, [gc%r, gc%z])
            call frame(lf%h, e1, e2)
            vperp = sqrt(2.0_dp*gc%mu*lf%bmod)
            v = gc%vpar*lf%h + vperp*(cos(zt)*e1 + sin(zt)*e2)
            x = [gc%r, 0.0_dp, gc%z] - cross(v, lf%b)/lf%bmod**2
            zz = [x(1), x(3), v(1), v(3) - par%f_tor*log(x(1))]
            p = x(1)*v(2) + psi_value(par, region_of(par, zz(1:2)), zz(1:2))
        end subroutine lifted_state

    end subroutine gc_to_fo

    subroutine gc_advance(par, gc, t_max, d_enter, tol, t_used, entered)
        !> Adaptive Dormand-Prince 5(4) on (R, Z, v_par, zeta) until the time
        !> t_max is used up or the GC comes within d_enter of the sheet.
        type(field_params_t), intent(in) :: par
        type(gc_state_t), intent(inout) :: gc
        real(dp), intent(in) :: t_max, d_enter, tol
        real(dp), intent(out) :: t_used
        logical, intent(out) :: entered

        real(dp) :: y(4), y_new(4), h, err, h_next

        y = [gc%r, gc%z, gc%vpar, gc%zeta]
        t_used = 0.0_dp
        h = 0.1_dp
        entered = .false.
        do while (t_used < t_max)
            h = min(h, t_max - t_used)
            call dp_step(y, h, y_new, err)
            err = err/tol
            if (err > 1.0_dp) then
                h = h*max(0.2_dp, 0.9_dp*err**(-0.2_dp))
                cycle
            end if
            h_next = h*min(5.0_dp, 0.9_dp*max(err, 1.0e-10_dp)**(-0.2_dp))
            y = y_new
            t_used = t_used + h
            h = h_next
            if (abs(sheet_distance(par, y(1:2))) < d_enter) then
                entered = .true.
                exit
            end if
        end do
        gc%r = y(1)
        gc%z = y(2)
        gc%vpar = y(3)
        gc%zeta = modulo(y(4), TWOPI)

    contains

        function rhs(yy) result(f)
            real(dp), intent(in) :: yy(4)
            real(dp) :: f(4)

            type(local_field_t) :: lf
            real(dp) :: bstar(3), bstar_par, xdot(3)

            lf = local_field(par, yy(1:2))
            bstar = lf%b + yy(3)*lf%curl_h
            bstar_par = dot_product(lf%h, bstar)
            xdot = (yy(3)*bstar + gc%mu*cross(lf%h, lf%grad_bmod))/bstar_par
            f = [xdot(1), xdot(3), -gc%mu*dot_product(bstar, lf%grad_bmod) &
                 /bstar_par, -lf%bmod]
        end function rhs

        subroutine dp_step(y0, hh, y1, e)
            real(dp), intent(in) :: y0(4), hh
            real(dp), intent(out) :: y1(4), e

            real(dp) :: k1(4), k2(4), k3(4), k4(4), k5(4), k6(4), k7(4), ev(4)

            k1 = rhs(y0)
            k2 = rhs(y0 + hh*k1/5)
            k3 = rhs(y0 + hh*(3*k1 + 9*k2)/40)
            k4 = rhs(y0 + hh*(44*k1/45 - 56*k2/15 + 32*k3/9))
            k5 = rhs(y0 + hh*(19372*k1/6561 - 25360*k2/2187 + 64448*k3/6561 &
                - 212*k4/729))
            k6 = rhs(y0 + hh*(9017*k1/3168 - 355*k2/33 + 46732*k3/5247 &
                + 49*k4/176 - 5103*k5/18656))
            y1 = y0 + hh*(35*k1/384 + 500*k3/1113 + 125*k4/192 &
                - 2187*k5/6784 + 11*k6/84)
            k7 = rhs(y1)
            ev = abs(hh*(71*k1/57600 - 71*k3/16695 + 71*k4/1920 &
                - 17253*k5/339200 + 22*k6/525 - k7/40))
            e = max(maxval(ev(1:3)/max(abs(y0(1:3)), 1.0e-3_dp)), ev(4))
        end subroutine dp_step

    end subroutine gc_advance

end module gc_hybrid
