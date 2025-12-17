module field_splined
    !> Splined field decorator for magnetic fields.
    !>
    !> Wraps any magnetic_field_t and provides fast splined evaluation.
    !> Splines are built in scaled coordinates (default: r = sqrt(s) for better
    !> resolution near axis, matching libneo convention).
    !>
    !> The reference coordinate system (ref_coords) is used for:
    !>   - Mapping from reference coordinates to Cartesian (evaluate_cart)
    !>   - Computing covariant basis vectors (covariant_basis)
    !> This allows proper support for VMEC, chartmap, or any coordinate system.
    !>
    !> Usage:
    !>   type(coils_field_t) :: raw_coils
    !>   type(splined_field_t) :: splined
    !>   call create_coils_field('coils.dat', raw_coils)
    !>   call create_splined_field(raw_coils, ref_coords, splined)
    !>   call splined%evaluate(x, Acov, hcov, Bmod)
    !>
    !> The scaling parameter controls how the first coordinate is transformed:
    !>   - sqrt_s_scaling_t (default): grid in r = sqrt(s), better axis resolution
    !>
    !> Limitations:
    !>   - Source field must evaluate in Cartesian coordinates
    !>   - sqgBctr output not supported (error stop if requested)

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use interpolate, only: BatchSplineData3D, construct_batch_splines_3d, &
                           evaluate_batch_splines_3d, evaluate_batch_splines_3d_der, &
                           destroy_batch_splines_3d
    use field_base, only: magnetic_field_t
    use field_vmec, only: vmec_field_t
    use libneo_coordinates, only: coordinate_system_t, vmec_coordinate_system_t, &
                                  chartmap_coordinate_system_t, RHO_TOR, RHO_POL
    use coordinate_scaling, only: coordinate_scaling_t, sqrt_s_scaling_t, &
                                  identity_scaling_t
    use cartesian_coordinates, only: cartesian_coordinate_system_t

    implicit none

    type, extends(magnetic_field_t) :: splined_field_t
        type(BatchSplineData3D) :: splines
        logical :: initialized = .false.
    contains
        procedure :: evaluate => splined_evaluate
        procedure :: evaluate_with_der => splined_evaluate_with_der
        final :: splined_field_cleanup
    end type splined_field_t

contains

    subroutine splined_evaluate(self, x, Acov, hcov, Bmod, sqgBctr)
        !> Evaluate splined field at x in scaled coordinates.
        class(splined_field_t), intent(in) :: self
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: Acov(3), hcov(3), Bmod
        real(dp), intent(out), optional :: sqgBctr(3)

        real(dp) :: y_batch(7)

        call evaluate_batch_splines_3d(self%splines, x, y_batch)

        Acov(1) = y_batch(1)
        Acov(2) = y_batch(2)
        Acov(3) = y_batch(3)

        hcov(1) = y_batch(4)
        hcov(2) = y_batch(5)
        hcov(3) = y_batch(6)

        Bmod = y_batch(7)

        if (present(sqgBctr)) then
            error stop "sqgBctr not implemented for splined_field_t"
        end if
    end subroutine splined_evaluate

    subroutine splined_evaluate_with_der(self, x, Acov, hcov, Bmod, dAcov, dhcov, &
                                         dBmod, sqgBctr)
        !> Evaluate splined field at x in scaled coordinates and return
        !> first derivatives with respect to the same coordinates x.
        class(splined_field_t), intent(in) :: self
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: Acov(3), hcov(3), Bmod
        real(dp), intent(out) :: dAcov(3, 3), dhcov(3, 3), dBmod(3)
        real(dp), intent(out), optional :: sqgBctr(3)

        real(dp) :: y_batch(7)
        real(dp) :: dy_batch(3, 7)

        call evaluate_batch_splines_3d_der(self%splines, x, y_batch, dy_batch)

        Acov(1) = y_batch(1)
        Acov(2) = y_batch(2)
        Acov(3) = y_batch(3)

        hcov(1) = y_batch(4)
        hcov(2) = y_batch(5)
        hcov(3) = y_batch(6)

        Bmod = y_batch(7)

        dAcov(:, 1) = dy_batch(:, 1)
        dAcov(:, 2) = dy_batch(:, 2)
        dAcov(:, 3) = dy_batch(:, 3)

        dhcov(:, 1) = dy_batch(:, 4)
        dhcov(:, 2) = dy_batch(:, 5)
        dhcov(:, 3) = dy_batch(:, 6)

        dBmod = dy_batch(:, 7)

        if (present(sqgBctr)) then
            error stop "sqgBctr not implemented for splined_field_t"
        end if
    end subroutine splined_evaluate_with_der

    subroutine create_splined_field(source, ref_coords, field, scaling, &
                                    n_r, n_th, n_phi, xmin, xmax)
        !> Create splined field by sampling source field onto ref_coords grid.
        !> Grid is built in scaled coordinates (default: r = sqrt(s)).
        class(magnetic_field_t), intent(in) :: source
        class(coordinate_system_t), intent(in), target :: ref_coords
        type(splined_field_t), intent(out) :: field
        class(coordinate_scaling_t), intent(in), optional, target :: scaling
        integer, intent(in), optional :: n_r, n_th, n_phi
        real(dp), intent(in), optional :: xmin(3), xmax(3)

        call build_splines(source, ref_coords, field%splines, scaling, &
                           n_r, n_th, n_phi, xmin, xmax)

        allocate (field%coords, source=ref_coords)
        field%initialized = .true.
    end subroutine create_splined_field

    subroutine build_splines(source, ref_coords, spl, scaling, &
                             n_r_in, n_th_in, n_phi_in, xmin_in, xmax_in)
        !> Build splines by sampling source field on uniform grid.
        !> Grid is in scaled coordinates. Default scaling: r = sqrt(s).
        class(magnetic_field_t), intent(in) :: source
        class(coordinate_system_t), intent(in) :: ref_coords
        type(BatchSplineData3D), intent(out) :: spl
        class(coordinate_scaling_t), intent(in), optional, target :: scaling
        integer, intent(in), optional :: n_r_in, n_th_in, n_phi_in
        real(dp), intent(in), optional :: xmin_in(3), xmax_in(3)

        class(coordinate_scaling_t), pointer :: scaling_ptr
        type(sqrt_s_scaling_t), target :: sqrt_s_scaling
        type(identity_scaling_t), target :: identity_scaling

        integer :: n_r, n_th, n_phi
        real(dp) :: xmin(3), xmax(3)

        integer, parameter :: order(3) = [5, 5, 5]
        logical, parameter :: periodic(3) = [.False., .True., .True.]

        real(dp) :: h_r, h_th, h_phi

        real(dp), dimension(:, :, :), allocatable :: Ar, Ath, Aphi, hr, hth, hphi
        real(dp), dimension(:, :, :), allocatable :: Bmod_arr
        real(dp), dimension(:, :, :, :), allocatable :: y_batch
        integer :: dims(3)

        call configure_spline_grid(ref_coords, n_r_in, n_th_in, n_phi_in, xmin_in, &
                                   xmax_in, n_r, n_th, n_phi, xmin, xmax)

        if (present(scaling)) then
            scaling_ptr => scaling
        else
            scaling_ptr => get_default_scaling(ref_coords, sqrt_s_scaling, &
                                               identity_scaling)
        end if

        h_r = (xmax(1) - xmin(1))/(n_r - 1)
        h_th = (xmax(2) - xmin(2))/(n_th - 1)
        h_phi = (xmax(3) - xmin(3))/(n_phi - 1)

        call sample_spline_arrays(source, ref_coords, scaling_ptr, n_r, n_th, n_phi, &
                                  xmin, h_r, h_th, h_phi, Ar, Ath, Aphi, hr, hth, &
                                  hphi, Bmod_arr)

        dims = shape(Ar)
        allocate (y_batch(dims(1), dims(2), dims(3), 7))

        y_batch(:, :, :, 1) = Ar
        y_batch(:, :, :, 2) = Ath
        y_batch(:, :, :, 3) = Aphi
        y_batch(:, :, :, 4) = hr
        y_batch(:, :, :, 5) = hth
        y_batch(:, :, :, 6) = hphi
        y_batch(:, :, :, 7) = Bmod_arr

        call construct_batch_splines_3d(xmin, xmax, y_batch, order, periodic, spl)

    end subroutine build_splines

    subroutine configure_spline_grid(ref_coords, n_r_in, n_th_in, n_phi_in, xmin_in, &
                                     xmax_in, n_r, n_th, n_phi, xmin, xmax)
        use new_vmec_stuff_mod, only: nper
        use util, only: twopi
        class(coordinate_system_t), intent(in) :: ref_coords
        integer, intent(in), optional :: n_r_in, n_th_in, n_phi_in
        real(dp), intent(in), optional :: xmin_in(3), xmax_in(3)
        integer, intent(out) :: n_r, n_th, n_phi
        real(dp), intent(out) :: xmin(3), xmax(3)

        n_r = 62
        if (present(n_r_in)) n_r = n_r_in
        n_th = 63
        if (present(n_th_in)) n_th = n_th_in
        n_phi = 64
        if (present(n_phi_in)) n_phi = n_phi_in

        select type (cs => ref_coords)
        type is (chartmap_coordinate_system_t)
            if (.not. present(n_r_in)) n_r = cs%nrho
            if (.not. present(n_th_in)) n_th = cs%ntheta
            if (.not. present(n_phi_in)) n_phi = cs%nzeta
        class default
            continue
        end select

        xmin = [1d-12, 0d0, 0d0]
        xmax = [1d0, twopi, twopi/nper]
        if (present(xmin_in)) xmin = xmin_in
        if (present(xmax_in)) xmax = xmax_in

        select type (cs => ref_coords)
        type is (chartmap_coordinate_system_t)
            if (.not. present(xmin_in)) xmin(1) = max(xmin(1), 0.1d0)
            if (.not. present(xmax_in)) xmax(1) = min(xmax(1), 0.9d0)
        class default
            continue
        end select
    end subroutine configure_spline_grid

    subroutine print_progress(i, n)
        integer, intent(in) :: i, n
!$omp critical
        if (mod(i, max(1, n/10)) == 0) then
            write (*, '(a,f6.2,a)', advance='no') char(13), 100.0*i/n, ' %'
            call flush (6)
        end if
!$omp end critical
    end subroutine print_progress

    subroutine sample_spline_arrays(source, ref_coords, scaling_ptr, n_r, n_th, n_phi, &
                                    xmin, h_r, h_th, h_phi, Ar, Ath, Aphi, hr, hth, &
                                    hphi, Bmod_arr)
        class(magnetic_field_t), intent(in) :: source
        class(coordinate_system_t), intent(in) :: ref_coords
        class(coordinate_scaling_t), pointer, intent(in) :: scaling_ptr
        integer, intent(in) :: n_r, n_th, n_phi
        real(dp), intent(in) :: xmin(3), h_r, h_th, h_phi
        real(dp), allocatable, intent(out) :: Ar(:, :, :), Ath(:, :, :), Aphi(:, :, :)
        real(dp), allocatable, intent(out) :: hr(:, :, :), hth(:, :, :), hphi(:, :, :)
        real(dp), allocatable, intent(out) :: Bmod_arr(:, :, :)

        real(dp) :: x_scaled(3), x_ref(3), x_cart(3)
        real(dp) :: e_cov(3, 3), scaling_jac(3)
        real(dp) :: Bmod
        real(dp) :: Acov(3), hcov(3)
        logical :: needs_scaling
        integer :: i_r, i_th, i_phi, i_ctr

        allocate (Ar(n_r, n_th, n_phi), Ath(n_r, n_th, n_phi))
        allocate (Aphi(n_r, n_th, n_phi))
        allocate (hr(n_r, n_th, n_phi), hth(n_r, n_th, n_phi))
        allocate (hphi(n_r, n_th, n_phi))
        allocate (Bmod_arr(n_r, n_th, n_phi))

        i_ctr = 0
!$omp parallel private(i_r, i_th, i_phi, x_scaled, x_ref, x_cart, &
!$omp&                  e_cov, scaling_jac, Bmod, Acov, hcov, needs_scaling)
!$omp do
        do i_phi = 1, n_phi
!$omp atomic
            i_ctr = i_ctr + 1
            call print_progress(i_ctr, n_phi)
            do i_th = 1, n_th
                do i_r = 1, n_r
                    x_scaled(1) = xmin(1) + h_r*dble(i_r - 1)
                    x_scaled(2) = xmin(2) + h_th*dble(i_th - 1)
                    x_scaled(3) = xmin(3) + h_phi*dble(i_phi - 1)

                    call scaling_ptr%inverse(x_scaled, x_ref, scaling_jac)
                    call evaluate_at_ref_coords(source, ref_coords, x_scaled, &
                                                x_ref, x_cart, e_cov, Acov, &
                                                hcov, Bmod, needs_scaling)
                    if (needs_scaling) then
                        Acov(1) = Acov(1)*scaling_jac(1)
                        hcov(1) = hcov(1)*scaling_jac(1)
                    end if

                    Ar(i_r, i_th, i_phi) = Acov(1)
                    Ath(i_r, i_th, i_phi) = Acov(2)
                    Aphi(i_r, i_th, i_phi) = Acov(3)

                    hr(i_r, i_th, i_phi) = hcov(1)
                    hth(i_r, i_th, i_phi) = hcov(2)
                    hphi(i_r, i_th, i_phi) = hcov(3)

                    Bmod_arr(i_r, i_th, i_phi) = Bmod
                end do
            end do
        end do
!$omp end do
!$omp end parallel

        Ar = Ar - Ar(1, 1, 1)
        Ath = Ath - Ath(1, 1, 1)
        Aphi = Aphi - Aphi(1, 1, 1)
    end subroutine sample_spline_arrays

    subroutine evaluate_at_ref_coords(source, ref_coords, x_scaled, x_ref, x_cart, &
                                      e_cov, Acov, hcov, Bmod, needs_scaling)
        !> Evaluate source field for spline sampling.
        !>
        !> Supported source modes:
        !>   1. VMEC source: evaluate at x_ref = (s, theta, phi) and return covariant
        !>      components in reference coordinates. Caller must apply scaling.
        !>   2. Cartesian source: evaluate at x_cart and transform to reference
        !>      coordinates (s, theta, phi). Caller must apply scaling.
        !>
        !> All sources now return components in reference coordinates (s, theta, phi).
        !> Caller applies coordinate scaling to obtain components in scaled coords.
        class(magnetic_field_t), intent(in) :: source
        class(coordinate_system_t), intent(in) :: ref_coords
        real(dp), intent(in) :: x_scaled(3)
        real(dp), intent(in) :: x_ref(3)
        real(dp), intent(out) :: x_cart(3), e_cov(3, 3)
        real(dp), intent(out) :: Acov(3), hcov(3), Bmod
        logical, intent(out) :: needs_scaling

        real(dp) :: A_cart(3), h_cart(3)
        real(dp) :: x_src(3)

        select type (source)
        type is (vmec_field_t)
            select type (ref_coords)
            type is (vmec_coordinate_system_t)
                ! Same coordinate system: evaluate directly in VMEC coordinates.
                call source%evaluate(x_ref, Acov, hcov, Bmod)
                x_cart = 0d0
                e_cov = 0d0
                needs_scaling = .true.
                return
            class default
                ! Different reference coordinates: evaluate VMEC field at the same
                ! physical point and transform covariant components into ref_coords.
                call ref_coords%evaluate_cart(x_ref, x_cart)
                call find_vmec_coords_for_cart_point(source%coords, ref_coords, &
                                                     x_ref, x_cart, x_src)

                call source%evaluate(x_src, Acov, hcov, Bmod)

                call source%coords%cov_to_cart(x_src, Acov, A_cart)
                call source%coords%cov_to_cart(x_src, hcov, h_cart)

                call ref_coords%covariant_basis(x_ref, e_cov)
                Acov = matmul(A_cart, e_cov)
                hcov = matmul(h_cart, e_cov)
            end select
            needs_scaling = .true.
            return
        class default
            continue
        end select

        select type (coords => source%coords)
        type is (cartesian_coordinate_system_t)
            call ref_coords%evaluate_cart(x_ref, x_cart)
            call ref_coords%covariant_basis(x_ref, e_cov)
            call source%evaluate(x_cart, A_cart, h_cart, Bmod)
            Acov = matmul(A_cart, e_cov)
            hcov = matmul(h_cart, e_cov)
            needs_scaling = .true.
        class default
            error stop "evaluate_at_ref_coords: unsupported source coordinate system"
        end select
    end subroutine evaluate_at_ref_coords

    subroutine find_vmec_coords_for_cart_point(vmec_coords, ref_coords, u_ref, &
                                               x_target, u_vmec)
        !> Find VMEC coordinates u_vmec such that
        !> vmec_coords%evaluate_cart(u_vmec) == x_target.
        !> Uses a damped Newton iteration with finite-difference Jacobian.
        class(coordinate_system_t), intent(in) :: vmec_coords
        class(coordinate_system_t), intent(in) :: ref_coords
        real(dp), intent(in) :: u_ref(3)
        real(dp), intent(in) :: x_target(3)
        real(dp), intent(out) :: u_vmec(3)

        integer, parameter :: max_iter = 40
        real(dp), parameter :: tol = 1.0e-10_dp
        real(dp), parameter :: eps = 1.0e-6_dp

        real(dp) :: x_test(3), err(3), err_norm
        real(dp) :: J(3, 3), du(3), rhs(3)
        integer :: iter, info

        call initial_guess_vmec_from_refcoords(ref_coords, u_ref, u_vmec)

        do iter = 1, max_iter
            call vmec_coords%evaluate_cart(u_vmec, x_test)
            err = x_test - x_target
            err_norm = sqrt(sum(err**2))
            if (err_norm < tol) then
                return
            end if

            call compute_jacobian_fd(vmec_coords, u_vmec, eps, J)
            rhs = -err
            call solve_3x3(J, rhs, du, info)
            if (info /= 0) then
                print *, 'find_vmec_coords_for_cart_point: solve_3x3 failed, info=', &
                    info
                print *, '  u_ref    = ', u_ref
                print *, '  u_vmec   = ', u_vmec
                print *, '  x_target = ', x_target
                error stop
            end if

            call damped_newton_update(vmec_coords, u_vmec, du, x_target, err_norm)
        end do

        print *, 'find_vmec_coords_for_cart_point: Newton did not converge'
        print *, '  Final err_norm = ', err_norm
        print *, '  u_ref          = ', u_ref
        print *, '  u_vmec         = ', u_vmec
        error stop
    end subroutine find_vmec_coords_for_cart_point

    subroutine initial_guess_vmec_from_refcoords(ref_coords, u_ref, u_vmec)
        class(coordinate_system_t), intent(in) :: ref_coords
        real(dp), intent(in) :: u_ref(3)
        real(dp), intent(out) :: u_vmec(3)

        select type (ref_coords)
        type is (chartmap_coordinate_system_t)
            if (ref_coords%rho_convention == RHO_TOR .or. &
                ref_coords%rho_convention == RHO_POL) then
                u_vmec = [u_ref(1)**2, u_ref(2), u_ref(3)]
            else
                u_vmec = u_ref
            end if
        class default
            u_vmec = u_ref
        end select
    end subroutine initial_guess_vmec_from_refcoords

    subroutine compute_jacobian_fd(coords, u, eps, J)
        class(coordinate_system_t), intent(in) :: coords
        real(dp), intent(in) :: u(3)
        real(dp), intent(in) :: eps
        real(dp), intent(out) :: J(3, 3)

        real(dp) :: u_plus(3), u_minus(3), x_plus(3), x_minus(3)
        integer :: i

        do i = 1, 3
            u_plus = u
            u_minus = u
            u_plus(i) = u_plus(i) + eps
            u_minus(i) = u_minus(i) - eps

            call coords%evaluate_cart(u_plus, x_plus)
            call coords%evaluate_cart(u_minus, x_minus)
            J(:, i) = (x_plus - x_minus)/(2.0_dp*eps)
        end do
    end subroutine compute_jacobian_fd

    subroutine solve_3x3(A_in, b_in, x, info)
        real(dp), intent(in) :: A_in(3, 3)
        real(dp), intent(in) :: b_in(3)
        real(dp), intent(out) :: x(3)
        integer, intent(out) :: info

        interface
            subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
                integer, intent(in) :: n, nrhs, lda, ldb
                integer, intent(out) :: ipiv(*)
                double precision, intent(inout) :: a(lda, *)
                double precision, intent(inout) :: b(ldb, *)
                integer, intent(out) :: info
            end subroutine dgesv
        end interface

        real(dp) :: A(3, 3), b(3, 1)
        integer :: ipiv(3)

        A = A_in
        b(:, 1) = b_in
        call dgesv(3, 1, A, 3, ipiv, b, 3, info)
        x = b(:, 1)
    end subroutine solve_3x3

    subroutine project_vmec_coords(u)
        use util, only: twopi
        real(dp), intent(inout) :: u(3)
        u(1) = min(1.0_dp, max(0.0_dp, u(1)))
        u(2) = modulo(u(2), twopi)
        u(3) = modulo(u(3), twopi)
    end subroutine project_vmec_coords

    subroutine damped_newton_update(coords, u, du, x_target, err_norm)
        use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
        class(coordinate_system_t), intent(in) :: coords
        real(dp), intent(inout) :: u(3)
        real(dp), intent(in) :: du(3)
        real(dp), intent(in) :: x_target(3)
        real(dp), intent(in) :: err_norm

        real(dp) :: alpha, u_in(3), u_try(3), x_try(3), err_try(3), err_try_norm
        integer :: k

        u_in = u
        alpha = 1.0_dp
        do k = 1, 8
            u_try = u_in + alpha*du
            call project_vmec_coords(u_try)

            call coords%evaluate_cart(u_try, x_try)
            err_try = x_try - x_target
            err_try_norm = sqrt(sum(err_try**2))

            if (.not. ieee_is_nan(err_try_norm)) then
                if (err_try_norm < err_norm) then
                    u = u_try
                    return
                end if
            end if
            alpha = 0.5_dp*alpha
        end do

        u = u_in + 0.1_dp*du
        call project_vmec_coords(u)
    end subroutine damped_newton_update

    subroutine splined_field_cleanup(self)
        type(splined_field_t), intent(inout) :: self

        if (self%initialized) then
            call destroy_batch_splines_3d(self%splines)
            self%initialized = .false.
        end if
    end subroutine splined_field_cleanup

    function get_default_scaling(ref_coords, sqrt_s_scaling, identity_scaling) &
        result(scaling_ptr)
        !> Select default scaling based on coordinate system rho_convention.
        !> RHO_TOR/RHO_POL: identity (coords already in rho = sqrt(s))
        !> PSI_TOR_NORM/other: sqrt_s (coords in s, need to transform)
        class(coordinate_system_t), intent(in) :: ref_coords
        type(sqrt_s_scaling_t), intent(in), target :: sqrt_s_scaling
        type(identity_scaling_t), intent(in), target :: identity_scaling
        class(coordinate_scaling_t), pointer :: scaling_ptr

        select type (cs => ref_coords)
        type is (chartmap_coordinate_system_t)
            if (cs%rho_convention == RHO_TOR .or. &
                cs%rho_convention == RHO_POL) then
                scaling_ptr => identity_scaling
            else
                scaling_ptr => sqrt_s_scaling
            end if
        class default
            scaling_ptr => sqrt_s_scaling
        end select
    end function get_default_scaling

end module field_splined
