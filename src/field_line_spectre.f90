module field_line_spectre
    !> Symplectic field-line tracer for a single SPECTRE/SPEC volume.
    !>
    !> The magnetic field line is the flow of a one-degree-of-freedom canonical
    !> system with the toroidal angle zeta as time, the poloidal angle theta as
    !> coordinate, and the covariant vector-potential component psi = A_theta as
    !> the conjugate momentum. The radial coordinate s is recovered implicitly
    !> from psi at fixed (theta, zeta). The derived flow is
    !>     dtheta/dzeta = -dA_zeta/ds / dA_theta/ds
    !>     dpsi/dzeta   =  dA_zeta/dtheta - dA_theta/dtheta * dA_zeta/ds/dA_theta/ds
    !> (the explicit d/dzeta terms cancel exactly). One step advances the pair
    !> (theta, psi) by an implicit midpoint map. Both update signs
    !> set the field-line orientation, so a sign flip reverses iota; keep them
    !> as written.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use spectre_reader, only: spectre_data_t
    use spectre_basis, only: spectre_vecpot_t, eval_spectre_vector_potential
    use interpolate, only: BatchSplineData3D, construct_batch_splines_3d, &
                           destroy_batch_splines_3d, &
                           evaluate_batch_splines_3d_der

    implicit none
    private

    public :: field_line_psi, field_line_step, spectre_rz
    public :: prepare_field_line_spline, cleanup_field_line_spline

    integer, parameter :: max_newton = 20
    real(dp), parameter :: newton_tol = 1.0e-13_dp
    real(dp), parameter :: fd_scale = epsilon(1.0_dp)**(1.0_dp/3.0_dp)
    real(dp), parameter :: twopi = 6.283185307179586476925286766559_dp
    type(BatchSplineData3D) :: hotloop_spline
    integer :: hotloop_lvol = 0

contains

    !> Canonical momentum psi = A_theta at (lvol, s, theta, zeta).
    pure function field_line_psi(data, lvol, s, theta, zeta) result(psi)
        type(spectre_data_t), intent(in) :: data
        integer, intent(in) :: lvol
        real(dp), intent(in) :: s, theta, zeta
        real(dp) :: psi

        type(spectre_vecpot_t) :: av

        call eval_spectre_vector_potential(data, lvol, s, theta, zeta, av, &
                                           second_derivatives=.false.)
        psi = av%Ath
    end function field_line_psi

    subroutine prepare_field_line_spline(data, lvol, n_s, n_theta, n_zeta, order)
        type(spectre_data_t), intent(in) :: data
        integer, intent(in) :: lvol, n_s, n_theta, n_zeta, order

        type(spectre_vecpot_t) :: av
        real(dp), allocatable :: values(:, :, :, :)
        real(dp) :: xmin(3), xmax(3), x(3)
        integer :: i, j, k

        call cleanup_field_line_spline
        if (min(n_s, n_theta, n_zeta) <= order) then
            error stop 'field-line spline grid must exceed its order'
        end if
        xmin = [-1.0_dp, 0.0_dp, 0.0_dp]
        xmax = [1.0_dp, twopi, twopi/real(data%Nfp, dp)]
        allocate (values(n_s, n_theta, n_zeta, 2))
        do k = 1, n_zeta
            x(3) = grid_coordinate(k, n_zeta, xmin(3), xmax(3))
            do j = 1, n_theta
                x(2) = grid_coordinate(j, n_theta, xmin(2), xmax(2))
                do i = 1, n_s
                    x(1) = grid_coordinate(i, n_s, xmin(1), xmax(1))
                    call eval_spectre_vector_potential(data, lvol, x(1), x(2), &
                                                       x(3), av, &
                                                       second_derivatives=.false.)
                    values(i, j, k, :) = [av%Ath, av%Azt]
                end do
            end do
        end do
        call construct_batch_splines_3d(xmin, xmax, values, [order, order, order], &
                                        [.false., .true., .true.], hotloop_spline)
        hotloop_lvol = lvol
    end subroutine prepare_field_line_spline

    subroutine cleanup_field_line_spline
        if (hotloop_lvol /= 0) call destroy_batch_splines_3d(hotloop_spline)
        hotloop_lvol = 0
    end subroutine cleanup_field_line_spline

    pure function grid_coordinate(i, n, xmin, xmax) result(x)
        integer, intent(in) :: i, n
        real(dp), intent(in) :: xmin, xmax
        real(dp) :: x

        x = xmin + (xmax - xmin)*real(i - 1, dp)/real(n - 1, dp)
    end function grid_coordinate

    !> Advance the canonical pair (theta, psi) by implicit midpoint. The two
    !> midpoint residuals are solved in (s_mid, theta_mid); the endpoint radius
    !> is then recovered from A_theta(s1,theta1,zeta1) = psi1. Numerical
    !> Jacobians use only values and first derivatives of the vector potential.
    subroutine field_line_step(data, lvol, h, s, theta, zeta, psi)
        type(spectre_data_t), intent(in) :: data
        integer, intent(in) :: lvol
        real(dp), intent(in) :: h
        real(dp), intent(inout) :: s, theta, zeta, psi

        type(spectre_vecpot_t) :: av
        real(dp) :: s0, theta0, zeta0, psi0, smid, theta_delta
        logical :: midpoint_ok, radius_ok

        s0 = s
        theta0 = theta
        zeta0 = zeta
        psi0 = psi
        smid = s0
        theta_delta = 0.0_dp
        call solve_midpoint(data, lvol, h, s0, theta0, zeta0, psi0, smid, &
                            theta_delta, av, midpoint_ok)
        if (.not. midpoint_ok) then
            write (*, '(A,I0,A,2(ES13.5,1X))') &
                'field_line_step: midpoint Newton failed in volume ', lvol, &
                ', (s, theta) = ', smid, theta0 + theta_delta
            error stop 'field_line_step: midpoint Newton did not converge'
        end if

        theta = theta0 + 2.0_dp*theta_delta
        psi = 2.0_dp*av%Ath - psi0
        zeta = zeta0 + h
        s = 2.0_dp*smid - s0
        call recover_radius(data, lvol, theta, zeta, psi, s, radius_ok)
        if (.not. radius_ok) error stop 'field_line_step: endpoint radius failed'
    end subroutine field_line_step

    subroutine solve_midpoint(data, lvol, h, s0, theta0, zeta0, psi0, &
                              s, theta_delta, av, converged)
        type(spectre_data_t), intent(in) :: data
        integer, intent(in) :: lvol
        real(dp), intent(in) :: h, s0, theta0, zeta0, psi0
        real(dp), intent(inout) :: s, theta_delta
        type(spectre_vecpot_t), intent(out) :: av
        logical, intent(out) :: converged

        real(dp) :: residual(2), jac(2, 2), det, ds, dtheta
        integer :: iter

        converged = .false.
        do iter = 1, max_newton
            call midpoint_residual(data, lvol, h, s0, theta0, zeta0, psi0, &
                                   s, theta_delta, residual, av)
            if (maxval(abs(residual)) < newton_tol) then
                converged = .true.
                return
            end if
            call midpoint_jacobian(data, lvol, h, s0, theta0, zeta0, psi0, &
                                   s, theta_delta, jac)
            det = jac(1, 1)*jac(2, 2) - jac(1, 2)*jac(2, 1)
            if (abs(det) <= tiny(det)) return
            ds = (-residual(1)*jac(2, 2) + residual(2)*jac(1, 2))/det
            dtheta = (residual(1)*jac(2, 1) - residual(2)*jac(1, 1))/det
            s = s + ds
            theta_delta = theta_delta + dtheta
        end do
    end subroutine solve_midpoint

    subroutine midpoint_residual(data, lvol, h, s0, theta0, zeta0, psi0, &
                                 s, theta_delta, residual, av)
        type(spectre_data_t), intent(in) :: data
        integer, intent(in) :: lvol
        real(dp), intent(in) :: h, s0, theta0, zeta0, psi0, s, theta_delta
        real(dp), intent(out) :: residual(2)
        type(spectre_vecpot_t), intent(out) :: av

        real(dp) :: ratio, flow

        call evaluate_hotloop_av(data, lvol, s, theta0 + theta_delta, &
                                 zeta0 + 0.5_dp*h, av)
        ratio = av%dAzt(1)/av%dAth(1)
        flow = av%dAzt(2) - av%dAth(2)*ratio
        residual(1) = theta_delta + 0.5_dp*h*ratio
        residual(2) = av%Ath - psi0 - 0.5_dp*h*flow
    end subroutine midpoint_residual

    subroutine midpoint_jacobian(data, lvol, h, s0, theta0, zeta0, psi0, &
                                 s, theta_delta, jac)
        type(spectre_data_t), intent(in) :: data
        integer, intent(in) :: lvol
        real(dp), intent(in) :: h, s0, theta0, zeta0, psi0, s, theta_delta
        real(dp), intent(out) :: jac(2, 2)

        type(spectre_vecpot_t) :: av
        real(dp) :: plus(2), minus(2), delta

        delta = fd_scale*max(1.0_dp, abs(s))
        call midpoint_residual(data, lvol, h, s0, theta0, zeta0, psi0, &
                               s + delta, theta_delta, plus, av)
        call midpoint_residual(data, lvol, h, s0, theta0, zeta0, psi0, &
                               s - delta, theta_delta, minus, av)
        jac(:, 1) = (plus - minus)/(2.0_dp*delta)

        delta = fd_scale
        call midpoint_residual(data, lvol, h, s0, theta0, zeta0, psi0, &
                               s, theta_delta + delta, plus, av)
        call midpoint_residual(data, lvol, h, s0, theta0, zeta0, psi0, &
                               s, theta_delta - delta, minus, av)
        jac(:, 2) = (plus - minus)/(2.0_dp*delta)
    end subroutine midpoint_jacobian

    subroutine recover_radius(data, lvol, theta, zeta, psi, s, converged)
        type(spectre_data_t), intent(in) :: data
        integer, intent(in) :: lvol
        real(dp), intent(in) :: theta, zeta, psi
        real(dp), intent(inout) :: s
        logical, intent(out) :: converged

        type(spectre_vecpot_t) :: av
        real(dp) :: residual
        integer :: iter

        converged = .false.
        do iter = 1, max_newton
            call evaluate_hotloop_av(data, lvol, s, theta, zeta, av)
            residual = av%Ath - psi
            if (abs(residual) < newton_tol) then
                converged = .true.
                return
            end if
            s = s - residual/av%dAth(1)
        end do
    end subroutine recover_radius

    subroutine evaluate_hotloop_av(data, lvol, s, theta, zeta, av)
        type(spectre_data_t), intent(in) :: data
        integer, intent(in) :: lvol
        real(dp), intent(in) :: s, theta, zeta
        type(spectre_vecpot_t), intent(out) :: av

        real(dp) :: x(3), values(2), derivatives(3, 2)

        if (hotloop_lvol == lvol) then
            x = [s, modulo(theta, twopi), &
                 modulo(zeta, twopi/real(data%Nfp, dp))]
            call evaluate_batch_splines_3d_der(hotloop_spline, x, values, derivatives)
            av%Ath = values(1)
            av%Azt = values(2)
            av%dAth = derivatives(:, 1)
            av%dAzt = derivatives(:, 2)
        else
            call eval_spectre_vector_potential(data, lvol, s, theta, zeta, av, &
                                               second_derivatives=.false.)
        end if
    end subroutine evaluate_hotloop_av

    !> Cylindrical R, Z of a point in volume lvol from the SPEC interface Fourier
    !> blend, mirroring libneo_coordinates_spectre for output only: interior
    !> volumes blend interfaces lvol-1 and lvol linearly in s; the axis volume
    !> (lvol == 1) uses the sbar**m power law (sbar**2 for m = 0) with
    !> sbar = (s + 1)/2.
    pure subroutine spectre_rz(data, lvol, s, theta, zeta, R, Z)
        type(spectre_data_t), intent(in) :: data
        integer, intent(in) :: lvol
        real(dp), intent(in) :: s, theta, zeta
        real(dp), intent(out) :: R, Z

        integer :: ii, m, n
        real(dp) :: arg, carg, sarg, sbar, fj, alss, blss
        real(dp) :: rc, rs, zc, zs

        R = 0.0_dp
        Z = 0.0_dp
        do ii = 1, data%mn
            m = data%im(ii)
            n = data%in(ii)
            arg = real(m, dp)*theta - real(n, dp)*zeta
            carg = cos(arg)
            sarg = sin(arg)
            if (lvol == 1) then
                sbar = 0.5_dp*(s + 1.0_dp)
                if (m == 0) then
                    fj = sbar*sbar
                else
                    fj = sbar**m
                end if
                rc = data%Rbc(ii, 0) + (data%Rbc(ii, 1) - data%Rbc(ii, 0))*fj
                rs = data%Rbs(ii, 0) + (data%Rbs(ii, 1) - data%Rbs(ii, 0))*fj
                zc = data%Zbc(ii, 0) + (data%Zbc(ii, 1) - data%Zbc(ii, 0))*fj
                zs = data%Zbs(ii, 0) + (data%Zbs(ii, 1) - data%Zbs(ii, 0))*fj
            else
                alss = 0.5_dp*(1.0_dp - s)
                blss = 0.5_dp*(1.0_dp + s)
                rc = alss*data%Rbc(ii, lvol - 1) + blss*data%Rbc(ii, lvol)
                rs = alss*data%Rbs(ii, lvol - 1) + blss*data%Rbs(ii, lvol)
                zc = alss*data%Zbc(ii, lvol - 1) + blss*data%Zbc(ii, lvol)
                zs = alss*data%Zbs(ii, lvol - 1) + blss*data%Zbs(ii, lvol)
            end if
            R = R + rc*carg + rs*sarg
            Z = Z + zc*carg + zs*sarg
        end do
    end subroutine spectre_rz

end module field_line_spectre
