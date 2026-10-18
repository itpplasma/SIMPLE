module gc_fo_switch
    ! First-order canonical GC <-> particle switch (item D).
    ! Forward map T (GC -> particle) and its Newton inverse, both evaluated with
    ! the fortsym-generated selector kernels map0/map1 (gen/gen_switch1.f90).
    ! Units: x_dot = v, v_dot = v x B(x)/eps (unit mass, positive charge).
    ! GC variables: X (Cartesian), upar = v.b(X), mu = vperp^2/(2|B(X)|), phi.
    ! Phase: a = cos(phi) e1 + sin(phi) e2, e1 = normalised projection of the
    ! constant vector axis perpendicular to b(X), e2 = b x e1.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use gc_switch_generated, only: map0, map1
    implicit none
    private
    public :: dp, switch_jets_t, jets_provider, switch_forward, switch_jacobian
    public :: switch_inverse, energy_matching_inverse, switch_pack
    public :: SW_OK, SW_FIELD, SW_NEWTON, SW_POLAR

    integer, parameter :: SW_OK = 0, SW_FIELD = 2, SW_NEWTON = 3, SW_POLAR = 4

    !> Local field jets at one point X: b = B(X), db(i,j) = d_j B_i,
    !> d2b(i,j,k) = d_j d_k B_i (only used for the analytic X-columns of the
    !> Jacobian; zero is allowed and then those columns are approximate).
    type :: switch_jets_t
        real(dp) :: b(3) = 0.0_dp
        real(dp) :: db(3, 3) = 0.0_dp
        real(dp) :: d2b(3, 3, 3) = 0.0_dp
        real(dp) :: axis(3) = [0.0_dp, 0.0_dp, 1.0_dp]
        real(dp) :: epsilon = 1.0_dp
    end type switch_jets_t

    abstract interface
        !> Return the jets of the physical field at Cartesian X.
        subroutine jets_provider(X, jets, ierr)
            import :: dp, switch_jets_t
            real(dp), intent(in) :: X(3)
            type(switch_jets_t), intent(inout) :: jets
            integer, intent(out) :: ierr
        end subroutine jets_provider
    end interface

contains

    subroutine switch_pack(jets, upar, mu, phi, p, ierr)
        type(switch_jets_t), intent(in) :: jets
        real(dp), intent(in) :: upar, mu, phi
        real(dp), intent(out) :: p(130)
        integer, intent(out) :: ierr
        real(dp) :: bn, along
        p = 0.0_dp
        ierr = SW_FIELD
        bn = norm2(jets%b)
        if (.not. ieee_is_finite(bn) .or. bn <= tiny(1.0_dp)) return
        if (jets%epsilon <= 0.0_dp) return
        along = dot_product(jets%axis, jets%b)/(bn*norm2(jets%axis))
        if (1.0_dp - along**2 < 1.0e-10_dp) return
        ierr = SW_POLAR
        if (mu <= 0.0_dp) return
        p(4:7) = [upar, mu, phi, jets%epsilon]
        p(8:10) = jets%b
        p(11:19) = reshape(jets%db, [9])
        p(20:46) = reshape(jets%d2b, [27])
        p(47:49) = jets%axis
        ierr = SW_OK
    end subroutine switch_pack

    !> T: GC (X,upar,mu,phi) -> particle (x,v). order 0 or 1.
    subroutine switch_forward(order, jets, X, upar, mu, phi, x_p, v_p, ierr)
        integer, intent(in) :: order
        type(switch_jets_t), intent(in) :: jets
        real(dp), intent(in) :: X(3), upar, mu, phi
        real(dp), intent(out) :: x_p(3), v_p(3)
        integer, intent(out) :: ierr
        real(dp) :: jac(6, 6)
        call switch_jacobian(order, jets, X, upar, mu, phi, x_p, v_p, jac, ierr)
    end subroutine switch_forward

    subroutine switch_jacobian(order, jets, X, upar, mu, phi, x_p, v_p, jac, ierr)
        integer, intent(in) :: order
        type(switch_jets_t), intent(in) :: jets
        real(dp), intent(in) :: X(3), upar, mu, phi
        real(dp), intent(out) :: x_p(3), v_p(3), jac(6, 6)
        integer, intent(out) :: ierr
        real(dp) :: p(130), res(42)
        x_p = 0.0_dp
        v_p = 0.0_dp
        jac = 0.0_dp
        call switch_pack(jets, upar, mu, phi, p, ierr)
        if (ierr /= SW_OK) return
        if (order == 0) then
            call map0(p, res)
        else
            call map1(p, res)
        end if
        if (.not. all(ieee_is_finite(res))) then
            ierr = SW_FIELD
            return
        end if
        x_p = X + res(1:3)
        v_p = res(4:6)
        jac = reshape(res(7:42), [6, 6])
    end subroutine switch_jacobian

    !> Production-style energy-matching inverse (for contrast): Larmor
    !> displacement with b,|B| at the particle, then upar=v.b(X), mu from the
    !> perpendicular energy at X, phase measured in the frame at X.
    subroutine energy_matching_inverse(provider, jets, x_p, v_p, X, upar, mu, phi, ierr)
        procedure(jets_provider) :: provider
        type(switch_jets_t), intent(inout) :: jets
        real(dp), intent(in) :: x_p(3), v_p(3)
        real(dp), intent(out) :: X(3), upar, mu, phi
        integer, intent(out) :: ierr
        real(dp) :: bn, bh(3), vperp(3)
        call provider(x_p, jets, ierr)
        if (ierr /= SW_OK) return
        bn = norm2(jets%b)
        bh = jets%b/bn
        X = x_p - jets%epsilon*cross(bh, v_p)/bn
        call provider(X, jets, ierr)
        if (ierr /= SW_OK) return
        call local_velocity_split(jets, v_p, upar, mu, phi, vperp)
    end subroutine energy_matching_inverse

    subroutine local_velocity_split(jets, v_p, upar, mu, phi, vperp)
        type(switch_jets_t), intent(in) :: jets
        real(dp), intent(in) :: v_p(3)
        real(dp), intent(out) :: upar, mu, phi, vperp(3)
        real(dp) :: bn, bh(3), e1(3), e2(3)
        bn = norm2(jets%b)
        bh = jets%b/bn
        e1 = jets%axis - dot_product(jets%axis, bh)*bh
        e1 = e1/norm2(e1)
        e2 = cross(bh, e1)
        upar = dot_product(v_p, bh)
        vperp = v_p - upar*bh
        mu = dot_product(vperp, vperp)/(2.0_dp*bn)
        phi = atan2(dot_product(vperp, e2), dot_product(vperp, e1))
    end subroutine local_velocity_split

    !> T^{-1} by Newton: solve T(X,upar,mu,phi) = (x_p,v_p), refreshing the
    !> field jets at every iterate. Starts from the energy-matching guess (or
    !> from z on input when warm=.true.). Converged when the residual stops
    !> decreasing below tol (relative to the particle state scales).
    subroutine switch_inverse(order, provider, jets, x_p, v_p, z, ierr, iterations, &
                              residual, warm)
        integer, intent(in) :: order
        procedure(jets_provider) :: provider
        type(switch_jets_t), intent(inout) :: jets
        real(dp), intent(in) :: x_p(3), v_p(3)
        real(dp), intent(inout) :: z(6)
        integer, intent(out) :: ierr
        integer, intent(out), optional :: iterations
        real(dp), intent(out), optional :: residual
        logical, intent(in), optional :: warm
        integer, parameter :: maxit = 60
        real(dp) :: xt(3), vt(3), jac(6, 6), r(6), dz(6), lscale, vscale
        real(dp) :: res_norm, best, phi0
        integer :: it, ipiv(6), info
        logical :: use_warm
        use_warm = .false.
        if (present(warm)) use_warm = warm
        if (.not. use_warm) then
            call energy_matching_inverse(provider, jets, x_p, v_p, z(1:3), z(4), &
                                         z(5), z(6), ierr)
            if (ierr /= SW_OK) return
        end if
        phi0 = z(6)
        vscale = norm2(v_p)
        lscale = jets%epsilon*vscale/norm2(jets%b)
        best = huge(1.0_dp)
        ierr = SW_NEWTON
        do it = 1, maxit
            call provider(z(1:3), jets, info)
            if (info /= SW_OK) then
                ierr = info
                return
            end if
            call switch_jacobian(order, jets, z(1:3), z(4), z(5), z(6), xt, vt, &
                                 jac, info)
            if (info /= SW_OK) then
                ierr = info
                return
            end if
            r(1:3) = xt - x_p
            r(4:6) = vt - v_p
            res_norm = max(maxval(abs(r(1:3)))/lscale, maxval(abs(r(4:6)))/vscale)
            if (present(iterations)) iterations = it
            if (present(residual)) residual = res_norm
            if (res_norm < 1.0e-15_dp .or. (res_norm >= 0.5_dp*best &
                                             .and. res_norm < 1.0e-12_dp)) then
                ierr = SW_OK
                exit
            end if
            best = min(best, res_norm)
            dz = -r
            call dgesv(6, 1, jac, 6, ipiv, dz, 6, info)
            if (info /= 0) then
                ierr = SW_NEWTON
                return
            end if
            z = z + dz
        end do
        if (abs(z(6) - phi0) > 1.0e3_dp) ierr = SW_NEWTON
    end subroutine switch_inverse

    pure function cross(a, b) result(c)
        real(dp), intent(in) :: a(3), b(3)
        real(dp) :: c(3)
        c = [a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), a(1)*b(2) - a(2)*b(1)]
    end function cross
end module gc_fo_switch
