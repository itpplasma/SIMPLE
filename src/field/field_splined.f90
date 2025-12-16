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
    use libneo_coordinates, only: coordinate_system_t, chartmap_coordinate_system_t, &
        RHO_TOR, RHO_POL
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

        allocate(field%coords, source=ref_coords)
        field%initialized = .true.
    end subroutine create_splined_field


    subroutine build_splines(source, ref_coords, spl, scaling, &
                            n_r_in, n_th_in, n_phi_in, xmin_in, xmax_in)
        !> Build splines by sampling source field on uniform grid.
        !> Grid is in scaled coordinates. Default scaling: r = sqrt(s).
        use new_vmec_stuff_mod, only: nper
        use util, only: twopi

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
        real(dp) :: x_scaled(3), x_ref(3), x_cart(3)
        real(dp) :: e_cov(3, 3), scaling_jac(3)
        real(dp) :: Bmod
        real(dp) :: Acov(3), hcov(3)
        logical :: needs_scaling

        real(dp), dimension(:,:,:), allocatable :: Ar, Ath, Aphi, hr, hth, hphi
        real(dp), dimension(:,:,:), allocatable :: Bmod_arr
        real(dp), dimension(:,:,:,:), allocatable :: y_batch
        integer :: i_r, i_th, i_phi, i_ctr, dims(3)

        n_r = 62; if (present(n_r_in)) n_r = n_r_in
        n_th = 63; if (present(n_th_in)) n_th = n_th_in
        n_phi = 64; if (present(n_phi_in)) n_phi = n_phi_in

        xmin = [1d-12, 0d0, 0d0]
        xmax = [1d0, twopi, twopi/nper]
        if (present(xmin_in)) xmin = xmin_in
        if (present(xmax_in)) xmax = xmax_in

        if (present(scaling)) then
            scaling_ptr => scaling
        else
            scaling_ptr => get_default_scaling(ref_coords, sqrt_s_scaling, &
                                               identity_scaling)
        end if

        h_r = (xmax(1) - xmin(1)) / (n_r - 1)
        h_th = (xmax(2) - xmin(2)) / (n_th - 1)
        h_phi = (xmax(3) - xmin(3)) / (n_phi - 1)

        allocate(Ar(n_r, n_th, n_phi), Ath(n_r, n_th, n_phi))
        allocate(Aphi(n_r, n_th, n_phi))
        allocate(hr(n_r, n_th, n_phi), hth(n_r, n_th, n_phi))
        allocate(hphi(n_r, n_th, n_phi))
        allocate(Bmod_arr(n_r, n_th, n_phi))

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
                    x_scaled(1) = xmin(1) + h_r * dble(i_r - 1)
                    x_scaled(2) = xmin(2) + h_th * dble(i_th - 1)
                    x_scaled(3) = xmin(3) + h_phi * dble(i_phi - 1)

                    call scaling_ptr%inverse(x_scaled, x_ref, scaling_jac)

                    call evaluate_at_ref_coords(source, ref_coords, x_scaled, &
                                               x_ref, x_cart, e_cov, Acov, hcov, Bmod, &
                                               needs_scaling)

                    if (needs_scaling) then
                        Acov(1) = Acov(1) * scaling_jac(1)
                        hcov(1) = hcov(1) * scaling_jac(1)
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

        dims = shape(Ar)
        allocate(y_batch(dims(1), dims(2), dims(3), 7))

        y_batch(:,:,:,1) = Ar
        y_batch(:,:,:,2) = Ath
        y_batch(:,:,:,3) = Aphi
        y_batch(:,:,:,4) = hr
        y_batch(:,:,:,5) = hth
        y_batch(:,:,:,6) = hphi
        y_batch(:,:,:,7) = Bmod_arr

        call construct_batch_splines_3d(xmin, xmax, y_batch, order, periodic, spl)

    contains

        subroutine print_progress(i, n)
            integer, intent(in) :: i, n
            !$omp critical
            if (mod(i, max(1, n / 10)) == 0) then
                write(*, '(a,f6.2,a)', advance='no') char(13), 100.0*i/n, ' %'
                call flush(6)
            end if
            !$omp end critical
        end subroutine print_progress

    end subroutine build_splines


    subroutine evaluate_at_ref_coords(source, ref_coords, x_scaled, x_ref, x_cart, &
                                      e_cov, Acov, hcov, Bmod, needs_scaling)
        !> Evaluate source field for spline sampling.
        !>
        !> Two supported source modes:
        !>   1. Cartesian source: evaluate at x_ref and return covariant components
        !>      in reference coordinates (s, theta, phi). Caller must apply scaling
        !>      to obtain covariant components in scaled coordinates.
        !>   2. VMEC source: evaluate directly at x_scaled = (rho, theta, phi) and
        !>      return covariant components already in scaled coordinates. Caller
        !>      must not apply scaling.
        class(magnetic_field_t), intent(in) :: source
        class(coordinate_system_t), intent(in) :: ref_coords
        real(dp), intent(in) :: x_scaled(3)
        real(dp), intent(in) :: x_ref(3)
        real(dp), intent(out) :: x_cart(3), e_cov(3, 3)
        real(dp), intent(out) :: Acov(3), hcov(3), Bmod
        logical, intent(out) :: needs_scaling

        real(dp) :: A_cart(3), h_cart(3)

        select type (source)
        type is (vmec_field_t)
            call source%evaluate(x_scaled, Acov, hcov, Bmod)
            x_cart = 0d0
            e_cov = 0d0
            needs_scaling = .false.
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
