module gc_api
    ! Local-jet API; no dependency on SIMPLE or fortsym at runtime.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use gc_generated, only: map0, map1, rhs0, rhs1, diagnostic
    implicit none
    private
    public :: dp, field_derivs_t, gc_transform, gc_inverse, gc_rhs
    public :: gc_transform_jacobian, gc_diagnostics
    public :: GC_OK, GC_ORDER, GC_FIELD, GC_NEWTON, GC_POLAR
    integer, parameter :: GC_OK = 0, GC_ORDER = 1, GC_FIELD = 2, &
                          GC_NEWTON = 3, GC_POLAR = 4

    type :: field_derivs_t
        real(dp) :: origin(3) = 0.0_dp
        real(dp) :: b(3) = 0.0_dp
        real(dp) :: db(3, 3) = 0.0_dp
        real(dp) :: d2b(3, 3, 3) = 0.0_dp
        real(dp) :: epsilon = 1.0_dp
        integer :: gauge_axis = 1
    end type field_derivs_t
contains
    subroutine pack(n, f, X, u, mu, phi, p, ierr)
        integer, intent(in) :: n
        type(field_derivs_t), intent(in) :: f
        real(dp), intent(in) :: X(3), u, mu, phi
        real(dp), intent(out) :: p(49)
        integer, intent(out) :: ierr
        real(dp) :: b(3), db(3, 3), dx(3), bn, norm_ref
        integer :: j, k
        p = 0.0_dp
        ierr = GC_ORDER
        if (n < 0 .or. n > 1) return
        ierr = GC_FIELD
        if (.not. ieee_is_finite(f%epsilon)) return
        if (f%epsilon <= 0.0_dp) return
        if (f%gauge_axis < 1 .or. f%gauge_axis > 3) return
        if (.not. all(ieee_is_finite(X))) return
        if (.not. all(ieee_is_finite([u, mu, phi]))) return
        dx = X - f%origin
        b = f%b + matmul(f%db, dx)
        do k = 1, 3
            do j = 1, 3
                b = b + 0.5_dp*f%d2b(:, j, k)*dx(j)*dx(k)
            end do
        end do
        bn = norm2(b)
        if (.not. ieee_is_finite(bn)) return
        if (bn <= tiny(1.0_dp)) return
        norm_ref = 1.0_dp - (b(f%gauge_axis)/bn)**2
        if (norm_ref < 1.0e-10_dp) return
        ierr = GC_POLAR
        if (mu <= 0.0_dp) return
        p(1:3) = 0.0_dp
        db = f%db
        do k = 1, 3
            db = db + f%d2b(:, :, k)*dx(k)
        end do
        p(4:7) = [u, mu, phi, f%epsilon]
        p(8:10) = b
        p(11:19) = reshape(db, [9])
        p(20:46) = reshape(f%d2b, [27])
        p(46 + f%gauge_axis) = 1.0_dp
        ierr = GC_OK
    end subroutine pack

    subroutine gc_transform(n, field_derivs, X, u, mu, phi, x_particle, v, ierr)
        integer, intent(in) :: n
        type(field_derivs_t), intent(in) :: field_derivs
        real(dp), intent(in) :: X(3), u, mu, phi
        real(dp), intent(out) :: x_particle(3), v(3)
        integer, intent(out) :: ierr
        real(dp) :: p(49), result(42)
        x_particle = 0.0_dp
        v = 0.0_dp
        call pack(n, field_derivs, X, u, mu, phi, p, ierr)
        if (ierr /= GC_OK) return
        if (n == 0) then
            call map0(p, result)
        else
            call map1(p, result)
        end if
        if (.not. all(ieee_is_finite(result))) then
            ierr = GC_FIELD
            return
        end if
        x_particle = result(1:3) + X
        v = result(4:6)
    end subroutine gc_transform

    subroutine gc_transform_jacobian(n, field_derivs, X, u, mu, phi, jac, ierr)
        integer, intent(in) :: n
        type(field_derivs_t), intent(in) :: field_derivs
        real(dp), intent(in) :: X(3), u, mu, phi
        real(dp), intent(out) :: jac(6, 6)
        integer, intent(out) :: ierr
        real(dp) :: p(49), result(42)
        jac = 0.0_dp
        call pack(n, field_derivs, X, u, mu, phi, p, ierr)
        if (ierr /= GC_OK) return
        if (n == 0) then
            call map0(p, result)
        else
            call map1(p, result)
        end if
        if (.not. all(ieee_is_finite(result))) then
            ierr = GC_FIELD
            return
        end if
        jac = reshape(result(7:42), [6, 6])
    end subroutine gc_transform_jacobian

    subroutine gc_rhs(n, field_derivs, X, u, mu, phi, X_dot, u_dot, mu_dot, &
                      phi_dot, bstar, bstar_parallel, drift, ierr)
        integer, intent(in) :: n
        type(field_derivs_t), intent(in) :: field_derivs
        real(dp), intent(in) :: X(3), u, mu, phi
        real(dp), intent(out) :: X_dot(3), u_dot, mu_dot, phi_dot
        real(dp), intent(out) :: bstar(3), bstar_parallel, drift(3)
        integer, intent(out) :: ierr
        real(dp) :: p(49), result(12)
        X_dot = 0.0_dp
        u_dot = 0.0_dp
        mu_dot = 0.0_dp
        phi_dot = 0.0_dp
        bstar = 0.0_dp
        bstar_parallel = 0.0_dp
        drift = 0.0_dp
        call pack(n, field_derivs, X, u, mu, phi, p, ierr)
        if (ierr /= GC_OK) return
        if (n == 0) then
            call rhs0(p, result)
        else
            call rhs1(p, result)
        end if
        if (.not. all(ieee_is_finite(result))) then
            ierr = GC_FIELD
            return
        end if
        if (abs(result(9)) < 1.0e-10_dp*norm2(field_derivs%b)) then
            ierr = GC_FIELD
            return
        end if
        X_dot = result(1:3)
        u_dot = result(4)
        phi_dot = result(5)
        bstar = result(6:8)
        bstar_parallel = result(9)
        drift = result(10:12)
    end subroutine gc_rhs

    subroutine gc_diagnostics(field_derivs, X, u, mu, phi, coefficients, ierr)
        type(field_derivs_t), intent(in) :: field_derivs
        real(dp), intent(in) :: X(3), u, mu, phi
        real(dp), intent(out) :: coefficients(14)
        integer, intent(out) :: ierr
        real(dp) :: p(49)
        coefficients = 0.0_dp
        call pack(1, field_derivs, X, u, mu, phi, p, ierr)
        if (ierr /= GC_OK) return
        call diagnostic(p, coefficients)
    end subroutine gc_diagnostics

    subroutine gc_inverse(n, field_derivs, x_particle, v, X, u, mu, phi, ierr)
        ! Newton in one fixed local quadratic field jet; principal phase branch.
        integer, intent(in) :: n
        type(field_derivs_t), intent(in) :: field_derivs
        real(dp), intent(in) :: x_particle(3), v(3)
        real(dp), intent(out) :: X(3), u, mu, phi
        integer, intent(out) :: ierr
        real(dp) :: z(6), trial(6), target(6), val(6), testval(6), jac(6, 6)
        real(dp) :: b(3), bh(3), ref(3), e1(3), e2(3), wv(3), delta(6)
        real(dp) :: bn, alpha, residual, new_residual, p(49), result(42)
        integer :: it, j, status
        X = 0.0_dp
        u = 0.0_dp
        mu = 0.0_dp
        phi = 0.0_dp
        call pack(n, field_derivs, x_particle, 0.0_dp, 1.0_dp, 0.0_dp, p, ierr)
        if (ierr /= GC_OK) return
        b = p(8:10)
        bn = norm2(b)
        bh = b/bn
        ref = 0.0_dp
        ref(field_derivs%gauge_axis) = 1.0_dp
        e1 = ref - dot_product(ref, bh)*bh
        e1 = e1/norm2(e1)
        e2 = cross(bh, e1)
        z(4) = dot_product(v, bh)
        wv = v - z(4)*bh
        z(5) = dot_product(wv, wv)/(2*bn)
        z(6) = atan2(dot_product(wv, e2), dot_product(wv, e1))
        z(1:3) = x_particle - field_derivs%epsilon*cross(bh, wv)/bn
        target = [x_particle, v]
        do it = 1, 24
            call pack(n, field_derivs, z(1:3), z(4), z(5), z(6), p, ierr)
            if (ierr /= GC_OK) return
            if (n == 0) then
                call map0(p, result)
            else
                call map1(p, result)
            end if
            val = result(1:6)
            val(1:3) = val(1:3) + z(1:3)
            residual = norm2(val - target)
            if (residual < 2.0e-12_dp*(1.0_dp + norm2(target))) then
                X = z(1:3)
                u = z(4)
                mu = z(5)
                phi = z(6)
                ierr = GC_OK
                return
            end if
            jac = reshape(result(7:42), [6, 6])
            call solve6(jac, target - val, delta, status)
            ierr = GC_NEWTON
            if (status /= 0) return
            alpha = 1.0_dp
            do j = 1, 14
                trial = z + alpha*delta
                call gc_transform(n, field_derivs, trial(1:3), trial(4), trial(5), &
                                  trial(6), testval(1:3), testval(4:6), status)
                new_residual = huge(1.0_dp)
                if (status == GC_OK) new_residual = norm2(testval - target)
                if (new_residual < residual) exit
                alpha = alpha/2
            end do
            if (j > 14) return
            z = trial
        end do
        ierr = GC_NEWTON
    end subroutine gc_inverse

    pure function cross(a, b) result(c)
        real(dp), intent(in) :: a(3), b(3)
        real(dp) :: c(3)
        c = [a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), &
             a(1)*b(2) - a(2)*b(1)]
    end function cross

    subroutine solve6(a, b, x, ierr)
        real(dp), intent(in) :: a(6, 6), b(6)
        real(dp), intent(out) :: x(6)
        integer, intent(out) :: ierr
        real(dp) :: aug(6, 7), row(7), factor
        integer :: k, j, pivot
        aug(:, 1:6) = a
        aug(:, 7) = b
        x = 0.0_dp
        ierr = 1
        do k = 1, 6
            pivot = k - 1 + maxloc(abs(aug(k:6, k)), dim=1)
            if (abs(aug(pivot, k)) < 1.0e-14_dp) return
            row = aug(k, :)
            aug(k, :) = aug(pivot, :)
            aug(pivot, :) = row
            aug(k, k:7) = aug(k, k:7)/aug(k, k)
            do j = 1, 6
                if (j == k) cycle
                factor = aug(j, k)
                aug(j, k:7) = aug(j, k:7) - factor*aug(k, k:7)
            end do
        end do
        x = aug(:, 7)
        ierr = 0
    end subroutine solve6
end module gc_api
