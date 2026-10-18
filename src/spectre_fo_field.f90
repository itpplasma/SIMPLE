module spectre_fo_field
    !> Per-volume splines for the SPECTRE full orbit in logical coordinates.
    !>
    !> Volume k owns the slab rho in [k-1, k] of the stacked SPECTRE chart. On
    !> that slab one quintic batch spline holds the cylindrical embedding R, Z
    !> (cm) and the covariant vector potential A_theta, A_zeta (G cm^2) in
    !> SPECTRE's A_rho = 0 gauge of that volume; zeta is the cylindrical angle,
    !> periodic over one field period. The basis, the Jacobian and B from dA
    !> follow from the jets, so the splined geometry and field are consistent
    !> with each other. A regular (xi, eta, zeta) axis chart carries R, Z and
    !> all three components of A; its jets join the polar representation before
    !> B is evaluated. Each slab is sampled up to its own edges (the upper edge
    !> just inside it), so the sheet at rho = k has the two one-sided limits of
    !> volumes k and k+1. Outside its slab a spline continues its edge cell's
    !> polynomial smoothly; the full orbit uses this only for trial stages.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use interpolate, only: BatchSplineData3D, construct_batch_splines_3d, &
        destroy_batch_splines_3d, evaluate_batch_splines_3d_der
    use magfie_sub, only: spectre_field, TESLA_TO_GAUSS, M_TO_CM
    implicit none
    private

    integer, parameter :: NQ = 4, AXIS_NQ = 5, SPLINE_ORDER = 5
    real(dp), parameter, public :: FO_AXIS_RADIUS = 0.05_dp
    real(dp), parameter :: AXIS_OUTER = 2*FO_AXIS_RADIUS
    real(dp), parameter :: A_SCALE = TESLA_TO_GAUSS*M_TO_CM**2
    !> Upper-edge sample offset: rho = k itself belongs to volume k+1.
    real(dp), parameter :: EDGE_INSET = 1.0e-12_dp
    real(dp), parameter :: TWOPI = 8.0_dp*atan(1.0_dp)

    type(BatchSplineData3D), allocatable :: vol_spl(:)
    type(BatchSplineData3D) :: axis_spl
    logical :: axis_built = .false.
    integer :: fo_field_mvol = 0
    logical, public, protected :: fo_field_axisymmetric = .false.

    public :: spectre_fo_field_build, spectre_fo_field_clear, &
        spectre_fo_field_local, spectre_fo_field_locate, &
        spectre_fo_local_t

    !> Local geometry and field at one logical point of one volume's spline.
    type :: spectre_fo_local_t
        real(dp) :: x(3) = 0.0_dp         !> Cartesian position (cm)
        real(dp) :: e(3, 3) = 0.0_dp      !> covariant basis e_i = dX/du^i
        real(dp) :: dual(3, 3) = 0.0_dp   !> contravariant basis (1/cm)
        real(dp) :: jac = 0.0_dp          !> signed Jacobian e_1.(e_2 x e_3)
        real(dp) :: acov(3) = 0.0_dp      !> covariant A (G cm^2), A_rho = 0
        real(dp) :: b(3) = 0.0_dp         !> unit vector of the interpolated B
        real(dp) :: bmod = 0.0_dp         !> |B| (G)
    end type spectre_fo_local_t

contains

    subroutine spectre_fo_field_build(n_rho, n_theta, n_zeta)
        !> Spline the field that magfie_sub currently holds, every volume.
        integer, intent(in) :: n_rho, n_theta, n_zeta
        real(dp), allocatable :: y(:, :, :, :)
        real(dp) :: u(3), xcyl(3), acov(3), hcov(3), bmod, period, h(3)
        integer :: k, i, j, l, nfp

        if (.not. allocated(spectre_field)) &
            error stop 'spectre_fo_field_build: no SPECTRE field'
        if (n_rho < SPLINE_ORDER + 1 .or. n_theta < SPLINE_ORDER + 2 .or. &
            n_zeta < SPLINE_ORDER + 2) &
            error stop 'spectre_fo_field_build: grid too coarse for quintic splines'
        call spectre_fo_field_clear
        fo_field_mvol = spectre_field%data%Mvol
        fo_field_axisymmetric = all(spectre_field%data%in == 0)
        nfp = max(spectre_field%data%Nfp, 1)
        period = TWOPI/real(nfp, dp)
        allocate (vol_spl(fo_field_mvol))
        allocate (y(n_rho, n_theta, n_zeta, NQ))
        h = [1.0_dp/real(n_rho - 1, dp), TWOPI/real(n_theta - 1, dp), &
            period/real(n_zeta - 1, dp)]
        do k = 1, fo_field_mvol
            !$omp parallel do collapse(2) private(i, j, l, u, xcyl, acov, hcov, bmod)
            do l = 1, n_zeta
                do j = 1, n_theta
                    do i = 1, n_rho
                        u = [real(k - 1, dp) + real(i - 1, dp)*h(1), &
                            real(j - 1, dp)*h(2), real(l - 1, dp)*h(3)]
                        ! Shared sheet geometry is sampled at the sheet itself;
                        ! only A needs a one-sided sample in the owning gauge.
                        call spectre_field%coords%evaluate_cyl(u, xcyl)
                        if (i == n_rho) u(1) = real(k, dp) - EDGE_INSET
                        if (i == 1 .and. k == 1) u(1) = EDGE_INSET
                        call spectre_field%evaluate(u, acov, hcov, bmod)
                        y(i, j, l, :) = [xcyl(1)*M_TO_CM, xcyl(3)*M_TO_CM, &
                            acov(2)*A_SCALE, acov(3)*A_SCALE]
                    end do
                end do
            end do
            !$omp end parallel do
            call construct_batch_splines_3d([real(k - 1, dp), 0.0_dp, 0.0_dp], &
                [real(k, dp), TWOPI, period], y, &
                [SPLINE_ORDER, SPLINE_ORDER, SPLINE_ORDER], &
                [.false., .true., .true.], vol_spl(k))
        end do
        call build_axis_spline(n_rho, n_zeta, period)
    end subroutine spectre_fo_field_build

    subroutine build_axis_spline(n_rho, n_zeta, period)
        ! Regular first-volume chart xi=rho*cos(theta), eta=rho*sin(theta).
        ! Tabulate all three components of the potential in that chart.
        integer, intent(in) :: n_rho, n_zeta
        real(dp), intent(in) :: period
        real(dp), allocatable :: values(:, :, :, :)
        real(dp) :: xi, eta, rho, u(3), cyl(3), a(3), hc(3), bmod, h
        integer :: n, i, j, k

        n = 2*max(3, n_rho/4) + 1
        h = 2*AXIS_OUTER/real(n - 1, dp)
        allocate (values(n, n, n_zeta, AXIS_NQ))
        !$omp parallel do collapse(2) &
        !$omp private(i, j, k, xi, eta, rho, u, cyl, a, hc, bmod)
        do k = 1, n_zeta
            do j = 1, n
                eta = -AXIS_OUTER + (j - 1)*h
                do i = 1, n
                    xi = -AXIS_OUTER + (i - 1)*h
                    rho = sqrt(xi*xi + eta*eta)
                    u = [rho, atan2(eta, xi), (k - 1)*period/(n_zeta - 1)]
                    call spectre_field%coords%evaluate_cyl(u, cyl)
                    if (rho < EDGE_INSET) u(1:2) = [EDGE_INSET, 0.0_dp]
                    call spectre_field%evaluate(u, a, hc, bmod)
                    values(i, j, k, 1:2) = [cyl(1), cyl(3)]*M_TO_CM
                    if (rho < EDGE_INSET) then
                        values(i, j, k, 3:4) = 0.0_dp
                    else
                        values(i, j, k, 3:4) = [-eta, xi]*a(2)*A_SCALE/rho**2
                    end if
                    values(i, j, k, 5) = a(3)*A_SCALE
                end do
            end do
        end do
        !$omp end parallel do
        call construct_batch_splines_3d([-AXIS_OUTER, -AXIS_OUTER, 0.0_dp], &
            [AXIS_OUTER, AXIS_OUTER, period], values, &
            [SPLINE_ORDER, SPLINE_ORDER, SPLINE_ORDER], &
            [.false., .false., .true.], axis_spl)
        axis_built = .true.
    end subroutine build_axis_spline

    subroutine spectre_fo_field_clear
        integer :: k

        fo_field_axisymmetric = .false.
        fo_field_mvol = 0
        if (axis_built) call destroy_batch_splines_3d(axis_spl)
        axis_built = .false.
        if (.not. allocated(vol_spl)) return
        do k = 1, size(vol_spl)
            call destroy_batch_splines_3d(vol_spl(k))
        end do
        deallocate (vol_spl)
    end subroutine spectre_fo_field_clear

    subroutine spectre_fo_field_local(u, owner, loc)
        real(dp), intent(in) :: u(3)
        integer, intent(in) :: owner
        type(spectre_fo_local_t), intent(out) :: loc
        real(dp) :: y(NQ), dy(3, NQ)

        if (owner == 1 .and. abs(u(1)) < AXIS_OUTER) then
            call axis_local(u, loc)
            return
        end if
        call evaluate_batch_splines_3d_der(vol_spl(owner), u, y, dy)
        call local_from_values(u(3), y(1), dy(:, 1), y(2), dy(:, 2), y(3), &
            dy(:, 3), y(4), dy(:, 4), loc)
    end subroutine spectre_fo_field_local

    pure subroutine local_from_values(zeta, r, dr, z, dz, at, dat, az, daz, loc)
        real(dp), intent(in) :: zeta, r, dr(3), z, dz(3), at, dat(3), az, daz(3)
        type(spectre_fo_local_t), intent(out) :: loc
        real(dp) :: sqgb(3)

        call geometry_from_values(zeta, r, dr, z, dz, loc)
        loc%acov = [0.0_dp, at, az]
        ! Match libneo: |J| B^i = eps^ijk d_j A_k with A_rho = 0.
        sqgb = [daz(2) - dat(3), -daz(1), dat(1)]
        call set_magnetic_field(loc, sqgb)
    end subroutine local_from_values

    pure subroutine geometry_from_values(zeta, r, dr, z, dz, loc)
        real(dp), intent(in) :: zeta, r, dr(3), z, dz(3)
        type(spectre_fo_local_t), intent(out) :: loc
        real(dp) :: c, s, cof(3, 3)
        integer :: i

        c = cos(zeta)
        s = sin(zeta)
        loc%x = [r*c, r*s, z]
        do i = 1, 3
            loc%e(:, i) = [dr(i)*c, dr(i)*s, dz(i)]
        end do
        loc%e(1, 3) = loc%e(1, 3) - r*s
        loc%e(2, 3) = loc%e(2, 3) + r*c
        cof(:, 1) = cross(loc%e(:, 2), loc%e(:, 3))
        cof(:, 2) = cross(loc%e(:, 3), loc%e(:, 1))
        cof(:, 3) = cross(loc%e(:, 1), loc%e(:, 2))
        loc%jac = dot_product(loc%e(:, 1), cof(:, 1))
        loc%dual = cof/loc%jac
    end subroutine geometry_from_values

    pure subroutine set_magnetic_field(loc, sqgb)
        type(spectre_fo_local_t), intent(inout) :: loc
        real(dp), intent(in) :: sqgb(3)
        loc%b = matmul(loc%e, sqgb)/abs(loc%jac)
        loc%bmod = norm2(loc%b)
        loc%b = loc%b/max(loc%bmod, tiny(1.0_dp))
    end subroutine set_magnetic_field

    subroutine axis_values(q, y, dy)
        ! Join representations of smooth geometry and the potential, including
        ! the derivative of the weight. Taking curl afterwards preserves the
        ! field consistency and avoids an artificial chart-boundary sheet.
        real(dp), intent(in) :: q(3)
        real(dp), intent(out) :: y(AXIS_NQ), dy(3, AXIS_NQ)
        real(dp) :: ya(AXIS_NQ), da(3, AXIS_NQ), yb(AXIS_NQ), db(3, AXIS_NQ)
        real(dp) :: p(NQ), dpdu(3, NQ), u(3), rho, c, s, f, df(3), t, w, dw(3)
        integer :: j

        rho = norm2(q(1:2))
        if (rho <= FO_AXIS_RADIUS) then
            call evaluate_batch_splines_3d_der(axis_spl, q, y, dy)
            return
        end if
        c = q(1)/rho
        s = q(2)/rho
        u = [rho, atan2(q(2), q(1)), q(3)]
        call evaluate_batch_splines_3d_der(vol_spl(1), u, p, dpdu)
        yb(1:2) = p(1:2)
        db(1, 1:2) = c*dpdu(1, 1:2) - s/rho*dpdu(2, 1:2)
        db(2, 1:2) = s*dpdu(1, 1:2) + c/rho*dpdu(2, 1:2)
        db(3, 1:2) = dpdu(3, 1:2)
        f = p(3)/rho
        df = [c*(dpdu(1, 3) - f)/rho - s*dpdu(2, 3)/rho**2, &
            s*(dpdu(1, 3) - f)/rho + c*dpdu(2, 3)/rho**2, dpdu(3, 3)/rho]
        yb(3:5) = [-s*f, c*f, p(4)]
        db(:, 3) = [s*c*f/rho, -c*c*f/rho, 0.0_dp] - s*df
        db(:, 4) = [s*s*f/rho, -c*s*f/rho, 0.0_dp] + c*df
        db(:, 5) = [c*dpdu(1, 4) - s/rho*dpdu(2, 4), &
            s*dpdu(1, 4) + c/rho*dpdu(2, 4), dpdu(3, 4)]
        if (rho >= AXIS_OUTER) then
            y = yb
            dy = db
            return
        end if
        call evaluate_batch_splines_3d_der(axis_spl, q, ya, da)
        t = (rho - FO_AXIS_RADIUS)/(AXIS_OUTER - FO_AXIS_RADIUS)
        ! C4 weight matches the continuity of the quintic potential splines.
        w = 1 - t**5*(126 - 420*t + 540*t*t - 315*t**3 + 70*t**4)
        dw = -630*t**4*(1 - t)**4/(AXIS_OUTER - FO_AXIS_RADIUS) &
            *[c, s, 0.0_dp]
        y = w*ya + (1 - w)*yb
        do j = 1, AXIS_NQ
            dy(:, j) = w*da(:, j) + (1 - w)*db(:, j) + dw*(ya(j) - yb(j))
        end do
    end subroutine axis_values

    subroutine axis_local(u, loc)
        real(dp), intent(in) :: u(3)
        type(spectre_fo_local_t), intent(out) :: loc
        real(dp) :: y(AXIS_NQ), dy(3, AXIS_NQ), q(3), c, s
        real(dp) :: basis(3, 2), dual(3, 2), curl_a(3)

        c = cos(u(2))
        s = sin(u(2))
        q = [u(1)*c, u(1)*s, u(3)]
        call axis_values(q, y, dy)
        call geometry_from_values(u(3), y(1), dy(:, 1), y(2), dy(:, 2), loc)
        curl_a = [dy(2, 5) - dy(3, 4), dy(3, 3) - dy(1, 5), &
            dy(1, 4) - dy(2, 3)]
        call set_magnetic_field(loc, curl_a)
        ! Expose the existing polar basis to callers; the magnetic field was
        ! evaluated with the nonsingular Cartesian-chart Jacobian above.
        basis = loc%e(:, 1:2)
        dual = loc%dual(:, 1:2)
        loc%e(:, 1) = c*basis(:, 1) + s*basis(:, 2)
        loc%e(:, 2) = u(1)*(-s*basis(:, 1) + c*basis(:, 2))
        loc%dual(:, 1) = c*dual(:, 1) + s*dual(:, 2)
        loc%dual(:, 2) = 0.0_dp
        if (abs(u(1)) > tiny(1.0_dp)) &
            loc%dual(:, 2) = (-s*dual(:, 1) + c*dual(:, 2))/u(1)
        loc%jac = loc%jac*u(1)
        loc%acov = [c*y(3) + s*y(4), u(1)*(-s*y(3) + c*y(4)), y(5)]
    end subroutine axis_local

    subroutine locate_axis(rt, zt, phi, v, status)
        real(dp), intent(in) :: rt, zt, phi
        real(dp), intent(inout) :: v(3)
        integer, intent(out) :: status
        real(dp) :: q(3), y(AXIS_NQ), dy(3, AXIS_NQ), fr, fz, det
        real(dp) :: delta(2), step, theta
        integer :: it

        q = [v(1)*cos(v(2)), v(1)*sin(v(2)), phi]
        do it = 1, 60
            call axis_values(q, y, dy)
            fr = y(1) - rt
            fz = y(2) - zt
            det = dy(1, 1)*dy(2, 2) - dy(2, 1)*dy(1, 2)
            if (abs(det) <= tiny(1.0_dp)) exit
            delta = [(fr*dy(2, 2) - fz*dy(2, 1))/det, &
                (fz*dy(1, 1) - fr*dy(1, 2))/det]
            step = maxval(abs(delta))
            if (step > AXIS_OUTER) delta = delta*AXIS_OUTER/step
            q(1:2) = q(1:2) - delta
            if (maxval(abs(delta)) <= 1.0e-12_dp .and. &
                abs(fr) + abs(fz) <= 1.0e-12_dp*max(rt, 1.0_dp)) then
                theta = v(2)
                if (norm2(q(1:2)) > tiny(1.0_dp)) theta = atan2(q(2), q(1))
                v = [norm2(q(1:2)), &
                    theta + TWOPI*nint((v(2) - theta)/TWOPI), phi]
                status = 0
                return
            end if
        end do
        status = 2
    end subroutine locate_axis

    subroutine spectre_fo_field_locate(x, u, owner, status)
        !> Logical point of the Cartesian x (cm): zeta is the cylindrical angle,
        !> (rho, theta) a Newton root of (R, Z) at that zeta, starting from u in
        !> volume owner. A root outside the owner's slab moves to the adjacent
        !> volume and is refined there. status = 0 on success, 1 if the point is
        !> beyond the outer boundary, 2 if
        !> Newton fails.
        real(dp), intent(in) :: x(3)
        real(dp), intent(inout) :: u(3)
        integer, intent(inout) :: owner
        integer, intent(out) :: status
        real(dp) :: rt, zt, phi, v(3), y(NQ), dy(3, NQ), fr, fz, det, dr, dt
        real(dp) :: scale
        integer :: it, own, switches

        rt = sqrt(x(1)**2 + x(2)**2)
        zt = x(3)
        phi = atan2(x(2), x(1))
        v = u
        v(3) = phi + TWOPI*nint((u(3) - phi)/TWOPI)
        phi = v(3)
        own = owner
        scale = max(rt, 1.0_dp)
        switches = 0
        status = 2
        do it = 1, 60
            if (own == 1 .and. v(1) < AXIS_OUTER) then
                call locate_axis(rt, zt, phi, v, status)
                if (status /= 0) return
                if (v(1) < AXIS_OUTER) then
                    u = v
                    owner = 1
                    return
                end if
                status = 2
            end if
            call evaluate_batch_splines_3d_der(vol_spl(own), v, y, dy)
            fr = y(1) - rt
            fz = y(2) - zt
            det = dy(1, 1)*dy(2, 2) - dy(2, 1)*dy(1, 2)
            if (abs(det) <= tiny(1.0_dp)) return
            dr = (fr*dy(2, 2) - fz*dy(2, 1))/det
            dt = (fz*dy(1, 1) - fr*dy(1, 2))/det
            dr = max(-0.25_dp, min(0.25_dp, dr))
            v(1) = v(1) - dr
            v(2) = v(2) - dt
            if (own == 1 .and. v(1) < 0.0_dp) then
                ! Polar chart flip: the axis-volume embedding uses rho**m
                ! (rho**2 for m=0), so this is the same physical point.
                v(1) = -v(1)
                v(2) = v(2) + 0.5_dp*TWOPI
            end if
            if (v(1) > real(own, dp) .and. own < fo_field_mvol .and. &
                switches < 8) then
                own = own + 1
                switches = switches + 1
            else if (v(1) < real(own - 1, dp) .and. own > 1 .and. &
                switches < 8) then
                own = own - 1
                switches = switches + 1
            end if
            ! Coordinate and physical residuals must both converge. Requiring
            ! a sub-ulp physical displacement can otherwise oscillate forever
            ! after the cylindrical embedding has reached roundoff.
            if (abs(dr) <= 1.0e-12_dp .and. abs(fr) + abs(fz) <= &
                1.0e-12_dp*scale) exit
        end do
        if (it > 60) return
        if (v(1) > real(fo_field_mvol, dp) .or. &
            (own == 1 .and. v(1) <= 0.0_dp)) then
            ! Preserve the converged exterior root so an entry map can
            ! distinguish an initial particle outside the wall from failure.
            u = v
            owner = own
            status = 1
            return
        end if
        ! Roundoff at a sheet: the root lies on the owner's edge.
        if (v(1) < real(own - 1, dp) - 1.0e-12_dp .or. &
            v(1) > real(own, dp) + 1.0e-12_dp) return
        v(1) = max(real(own - 1, dp), min(real(own, dp), v(1)))
        u = v
        owner = own
        status = 0
    end subroutine spectre_fo_field_locate

    pure function cross(a, b) result(c)
        real(dp), intent(in) :: a(3), b(3)
        real(dp) :: c(3)
        c = [a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), &
            a(1)*b(2) - a(2)*b(1)]
    end function cross
end module spectre_fo_field
