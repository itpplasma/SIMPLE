module candidate_field
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use sheet_field, only: field_params_t, psi_derivatives
    implicit none
    private
    public :: model_t, local_t, field, geometry, potential, distance, radius
    public :: cross, frame, side_of, normal, canonical_pphi, connection

    type :: model_t
        type(field_params_t) :: torus
        real(dp) :: delta = 0.0_dp
        integer :: helicity = 3
    end type

    type :: local_t
        real(dp) :: bvec(3), bmag, b(3), gradb(3), db(3, 3), curlb(3)
        real(dp) :: kappa(3), tau, gauge(3)
    end type

contains

    pure function cross(a, b) result(c)
        real(dp), intent(in) :: a(3), b(3)
        real(dp) :: c(3)
        c = [a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), &
             a(1)*b(2) - a(2)*b(1)]
    end function

    pure subroutine frame(b, e1, e2)
        real(dp), intent(in) :: b(3)
        real(dp), intent(out) :: e1(3), e2(3)
        e1 = cross(b, [0.0_dp, 0.0_dp, 1.0_dp])
        e1 = e1/norm2(e1)
        e2 = cross(b, e1)
    end subroutine

    pure function radius(pm, x) result(r)
        type(model_t), intent(in) :: pm
        real(dp), intent(in) :: x(3)
        real(dp) :: r
        r = hypot(hypot(x(1), x(2)) - pm%torus%r_major, x(3))
    end function

    pure function distance(pm, x) result(d)
        type(model_t), intent(in) :: pm
        real(dp), intent(in) :: x(3)
        real(dp) :: d
        d = radius(pm, x) - pm%torus%r_sheet
    end function

    pure function side_of(pm, x) result(region)
        type(model_t), intent(in) :: pm
        real(dp), intent(in) :: x(3)
        integer :: region
        region = merge(1, 2, distance(pm, x) < 0.0_dp)
    end function

    pure function normal(pm, x) result(n)
        type(model_t), intent(in) :: pm
        real(dp), intent(in) :: x(3)
        real(dp) :: n(3), rr, r, s
        rr = hypot(x(1), x(2))
        s = rr - pm%torus%r_major
        r = hypot(s, x(3))
        n = [s*x(1)/rr, s*x(2)/rr, x(3)]/r
    end function

    subroutine field(pm, x, region, bvec, jac)
        type(model_t), intent(in) :: pm
        real(dp), intent(in) :: x(3)
        integer, intent(in) :: region
        real(dp), intent(out) :: bvec(3)
        real(dp), intent(out), optional :: jac(3, 3)
        real(dp) :: rr, s, z, cp, sp, c, sn, k, h, d, nn
        real(dp) :: psi, gp(2), hp(2, 2), bc(3), dr(3), dz(3), dt(3)
        real(dp) :: rot(3, 3), deriv(3, 3), coords(3, 3)
        rr = hypot(x(1), x(2))
        cp = x(1)/rr
        sp = x(2)/rr
        s = rr - pm%torus%r_major
        z = x(3)
        nn = real(pm%helicity, dp)
        c = cos(nn*atan2(x(2), x(1)))
        sn = sin(nn*atan2(x(2), x(1)))
        k = -s*sn + z*c
        h = s*c + z*sn
        d = pm%delta
        call psi_derivatives(pm%torus, region, [rr, z], psi, gp, hp)
        bc = [-gp(2)/rr + d*nn*z*k/rr, pm%torus%f_tor/rr - d*k, &
              gp(1)/rr - d*nn*s*k/rr]
        rot = reshape([cp, sp, 0.0_dp, -sp, cp, 0.0_dp, &
                       0.0_dp, 0.0_dp, 1.0_dp], [3, 3])
        bvec = matmul(rot, bc)
        if (.not. present(jac)) return
        dr = [-hp(2, 1)/rr + gp(2)/rr**2 - d*nn*z*sn/rr &
              - d*nn*z*k/rr**2, -pm%torus%f_tor/rr**2 + d*sn, &
              hp(1, 1)/rr - gp(1)/rr**2 - d*nn*(k - s*sn)/rr &
              + d*nn*s*k/rr**2]
        dz = [-hp(2, 2)/rr + d*nn*(k + z*c)/rr, -d*c, &
              hp(1, 2)/rr - d*nn*s*c/rr]
        dt = [-d*nn**2*z*h/rr - bc(2), d*nn*h + bc(1), &
              d*nn**2*s*h/rr]
        deriv(:, 1) = matmul(rot, dr)
        deriv(:, 2) = matmul(rot, dt)
        deriv(:, 3) = matmul(rot, dz)
        coords = reshape([cp, -sp/rr, 0.0_dp, sp, cp/rr, 0.0_dp, &
                          0.0_dp, 0.0_dp, 1.0_dp], [3, 3])
        jac = matmul(deriv, coords)
    end subroutine

    subroutine geometry(pm, x, region, lf)
        type(model_t), intent(in) :: pm
        real(dp), intent(in) :: x(3)
        integer, intent(in) :: region
        type(local_t), intent(out) :: lf
        real(dp) :: jac(3, 3), e1(3), e2(3), jnorm
        integer :: i
        call field(pm, x, region, lf%bvec, jac)
        lf%bmag = norm2(lf%bvec)
        lf%b = lf%bvec/lf%bmag
        lf%gradb = matmul(transpose(jac), lf%b)
        do i = 1, 3
            lf%db(:, i) = (jac(:, i) - lf%b*lf%gradb(i))/lf%bmag
        end do
        lf%curlb = [lf%db(3, 2) - lf%db(2, 3), &
                   lf%db(1, 3) - lf%db(3, 1), lf%db(2, 1) - lf%db(1, 2)]
        lf%kappa = matmul(lf%db, lf%b)
        lf%tau = dot_product(lf%b, lf%curlb)
        call frame(lf%b, e1, e2)
        jnorm = norm2(cross(lf%b, [0.0_dp, 0.0_dp, 1.0_dp]))
        do i = 1, 3
            lf%gauge(i) = dot_product(e2, &
                cross(lf%db(:, i), [0.0_dp, 0.0_dp, 1.0_dp]))/jnorm
        end do
    end subroutine

    subroutine connection(pm, x, region, rstar, curlr)
        type(model_t), intent(in) :: pm
        real(dp), intent(in) :: x(3)
        integer, intent(in) :: region
        real(dp), intent(out) :: rstar(3), curlr(3)
        type(local_t) :: lf, lp, lm
        real(dp) :: xp(3), xm(3), jr(3, 3), gt(3), gc(3), rr, cp, sp, a
        real(dp), parameter :: dx = 2.0e-5_dp
        integer :: j
        call geometry(pm, x, region, lf)
        rstar = lf%gauge + 0.5_dp*lf%tau*lf%b
        if (pm%torus%width <= 0.0_dp) then
            rr = hypot(x(1), x(2))
            cp = x(1)/rr
            sp = x(2)/rr
            a = pm%torus%a_in
            if (region == 2) a = pm%torus%a_out
            call torsion_gradient(rr, x(3), atan2(x(2), x(1)), &
                pm%torus%r_major, pm%torus%f_tor, a, pm%delta, &
                real(pm%helicity, dp), gt)
            gc = [cp*gt(1) - sp*gt(2)/rr, &
                  sp*gt(1) + cp*gt(2)/rr, gt(3)]
            curlr = [-dot_product(lf%b, cross(lf%db(:, 2), lf%db(:, 3))), &
                     -dot_product(lf%b, cross(lf%db(:, 3), lf%db(:, 1))), &
                     -dot_product(lf%b, cross(lf%db(:, 1), lf%db(:, 2)))] &
                    + 0.5_dp*(cross(gc, lf%b) + lf%tau*lf%curlb)
            return
        end if
        do j = 1, 3
            xp = x
            xm = x
            xp(j) = xp(j) + dx
            xm(j) = xm(j) - dx
            call geometry(pm, xp, region, lp)
            call geometry(pm, xm, region, lm)
            jr(:, j) = (lp%gauge + 0.5_dp*lp%tau*lp%b &
                      - lm%gauge - 0.5_dp*lm%tau*lm%b)/(2.0_dp*dx)
        end do
        curlr = [jr(3, 2) - jr(2, 3), jr(1, 3) - jr(3, 1), &
                 jr(2, 1) - jr(1, 2)]
    end subroutine

    include 'torsion_kernel.inc'

    function potential(pm, x) result(avec)
        type(model_t), intent(in) :: pm
        real(dp), intent(in) :: x(3)
        real(dp) :: avec(3), rr, s, cp, sp, h, psi, gp(2), hp(2, 2)
        rr = hypot(x(1), x(2))
        cp = x(1)/rr
        sp = x(2)/rr
        s = rr - pm%torus%r_major
        h = s*cos(pm%helicity*atan2(x(2), x(1))) &
            + x(3)*sin(pm%helicity*atan2(x(2), x(1)))
        call psi_derivatives(pm%torus, side_of(pm, x), [rr, x(3)], &
                             psi, gp, hp)
        avec = [pm%delta*h*s*cp - psi*sp/rr, &
                pm%delta*h*s*sp + psi*cp/rr, &
                -pm%torus%f_tor*log(rr) + pm%delta*h*x(3)]
    end function

    function canonical_pphi(pm, y) result(p)
        type(model_t), intent(in) :: pm
        real(dp), intent(in) :: y(6)
        real(dp) :: p, a(3)
        a = potential(pm, y(1:3))
        p = y(1)*(y(5) + a(2)) - y(2)*(y(4) + a(1))
    end function

end module
