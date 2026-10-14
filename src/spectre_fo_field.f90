module spectre_fo_field
    !> Per-volume splines for the SPECTRE full orbit in logical coordinates.
    !>
    !> Volume k owns the slab rho in [k-1, k] of the stacked SPECTRE chart. On
    !> that slab one quintic batch spline holds the cylindrical embedding R, Z
    !> (cm) and the covariant vector potential A_theta, A_zeta (G cm^2) in
    !> SPECTRE's A_rho = 0 gauge of that volume; zeta is the cylindrical angle,
    !> periodic over one field period. The basis, the Jacobian and B from dA
    !> follow from the jets, so the splined geometry and field are consistent
    !> with each other. Each slab is sampled up to its own edges (the upper edge
    !> just inside it), so the sheet at rho = k has the two one-sided limits of
    !> volumes k and k+1. Outside its slab a spline continues its edge cell's
    !> polynomial smoothly; the full orbit uses this only for trial stages.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use interpolate, only: BatchSplineData3D, construct_batch_splines_3d, &
        destroy_batch_splines_3d, evaluate_batch_splines_3d_der
    use magfie_sub, only: spectre_field, TESLA_TO_GAUSS, M_TO_CM
    implicit none
    private

    integer, parameter :: NQ = 4, SPLINE_ORDER = 5
    real(dp), parameter :: A_SCALE = TESLA_TO_GAUSS*M_TO_CM**2
    !> Upper-edge sample offset: rho = k itself belongs to volume k+1.
    real(dp), parameter :: EDGE_INSET = 1.0e-12_dp
    real(dp), parameter :: TWOPI = 8.0_dp*atan(1.0_dp)

    type(BatchSplineData3D), allocatable :: vol_spl(:)
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
    end subroutine spectre_fo_field_build

    subroutine spectre_fo_field_clear
        integer :: k

        fo_field_axisymmetric = .false.
        fo_field_mvol = 0
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

        call evaluate_batch_splines_3d_der(vol_spl(owner), u, y, dy)
        call local_from_values(u(3), y(1), dy(:, 1), y(2), dy(:, 2), y(3), &
            dy(:, 3), y(4), dy(:, 4), loc)
    end subroutine spectre_fo_field_local

    pure subroutine local_from_values(zeta, r, dr, z, dz, at, dat, az, daz, loc)
        real(dp), intent(in) :: zeta, r, dr(3), z, dz(3), at, dat(3), az, daz(3)
        type(spectre_fo_local_t), intent(out) :: loc
        real(dp) :: c, s, sqgb(3), cof(3, 3)
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
        loc%acov = [0.0_dp, at, az]
        ! Match libneo: |J| B^i = eps^ijk d_j A_k with A_rho = 0.
        sqgb = [daz(2) - dat(3), -daz(1), dat(1)]
        loc%b = matmul(loc%e, sqgb)/abs(loc%jac)
        loc%bmod = norm2(loc%b)
        loc%b = loc%b/max(loc%bmod, tiny(1.0_dp))
    end subroutine local_from_values

    subroutine spectre_fo_field_locate(x, u, owner, status)
        !> Logical point of the Cartesian x (cm): zeta is the cylindrical angle,
        !> (rho, theta) a Newton root of (R, Z) at that zeta, starting from u in
        !> volume owner. A root outside the owner's slab moves to the adjacent
        !> volume and is refined there. status = 0 on success, 1 if the point is
        !> beyond the outer boundary or near the axis chart singularity, 2 if
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
        own = owner
        scale = max(rt, 1.0_dp)
        switches = 0
        status = 2
        do it = 1, 60
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
