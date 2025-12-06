module field_splined
    !> Splined field wrapper.
    !> Single responsibility: wrap any magnetic_field_t and provide splined
    !> evaluation in reference coordinates.

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use interpolate, only: BatchSplineData3D, construct_batch_splines_3d, &
        evaluate_batch_splines_3d, destroy_batch_splines_3d
    use field_base, only: magnetic_field_t
    use libneo_coordinates, only: coordinate_system_t
    use coordinate_scaling, only: coordinate_scaling_t, sqrt_s_scaling_t
    use cartesian_coordinates, only: cartesian_coordinate_system_t
    use simple_coordinates, only: transform_vmec_to_cart

    implicit none

    type, extends(magnetic_field_t) :: splined_field_t
        type(BatchSplineData3D) :: splines
        logical :: initialized = .false.
    contains
        procedure :: evaluate => splined_evaluate
        final :: splined_field_cleanup
    end type splined_field_t

contains

    subroutine splined_evaluate(self, x, Acov, hcov, Bmod, sqgBctr)
        !> Evaluate splined field at x in reference coordinates.
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


    subroutine create_splined_field(source, ref_coords, field, scaling, &
                                    n_r, n_th, n_phi, xmin, xmax)
        !> Create splined field by sampling source field onto ref_coords grid.
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
        !> Build splines by sampling source field on ref_coords grid.
        use new_vmec_stuff_mod, only: nper
        use util, only: twopi

        class(magnetic_field_t), intent(in) :: source
        class(coordinate_system_t), intent(in) :: ref_coords
        type(BatchSplineData3D), intent(out) :: spl
        class(coordinate_scaling_t), intent(in), optional, target :: scaling
        integer, intent(in), optional :: n_r_in, n_th_in, n_phi_in
        real(dp), intent(in), optional :: xmin_in(3), xmax_in(3)

        class(coordinate_scaling_t), pointer :: scaling_ptr
        type(sqrt_s_scaling_t), target :: default_scaling

        integer :: n_r, n_th, n_phi
        real(dp) :: xmin(3), xmax(3)

        integer, parameter :: order(3) = [5, 5, 5]
        logical, parameter :: periodic(3) = [.False., .True., .True.]

        real(dp) :: h_r, h_th, h_phi
        real(dp) :: x_spline(3), x_ref(3), x_cart(3)
        real(dp) :: dxcart_dxvmec(3, 3), scaling_jac(3)
        real(dp) :: A_cart(3), B_cart(3), Bmod
        real(dp) :: Acov(3), hcov(3)

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
            scaling_ptr => default_scaling
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
        !$omp parallel private(i_r, i_th, i_phi, x_spline, x_ref, x_cart, &
        !$omp&                  dxcart_dxvmec, scaling_jac, A_cart, B_cart, Bmod, &
        !$omp&                  Acov, hcov)
        !$omp do
        do i_phi = 1, n_phi
            !$omp atomic
            i_ctr = i_ctr + 1
            call print_progress(i_ctr, n_phi)
            do i_th = 1, n_th
                do i_r = 1, n_r
                    x_spline(1) = xmin(1) + h_r * dble(i_r - 1)
                    x_spline(2) = xmin(2) + h_th * dble(i_th - 1)
                    x_spline(3) = xmin(3) + h_phi * dble(i_phi - 1)

                    call scaling_ptr%inverse(x_spline, x_ref, scaling_jac)

                    call evaluate_at_ref_point(source, x_ref, x_cart, &
                                              dxcart_dxvmec, A_cart, B_cart, Bmod)

                    call transform_to_ref(A_cart, B_cart, dxcart_dxvmec, &
                                         scaling_jac, Bmod, Acov, hcov)

                    if (i_r == 32 .and. i_th == 32 .and. i_phi == 32) then
                        print '(A,3ES22.15)', ' SPLINED_DIAG: x_spline = ', x_spline
                        print '(A,3ES22.15)', ' SPLINED_DIAG: x_ref = ', x_ref
                        print '(A,3ES22.15)', ' SPLINED_DIAG: x_cart = ', x_cart
                        print '(A,3ES22.15)', ' SPLINED_DIAG: A_cart = ', A_cart
                        print '(A,3ES22.15)', ' SPLINED_DIAG: B_cart = ', B_cart
                        print '(A,ES22.15)', ' SPLINED_DIAG: Bmod = ', Bmod
                        print '(A,ES22.15)', ' SPLINED_DIAG: scaling_jac(1) = ', &
                            scaling_jac(1)
                        print '(A,3ES22.15)', ' SPLINED_DIAG: Acov = ', Acov
                        print '(A,3ES22.15)', ' SPLINED_DIAG: hcov = ', hcov
                        print '(A,3ES22.15)', ' SPLINED_DIAG: dxcart_dxvmec(:,1) = ', &
                            dxcart_dxvmec(:,1)
                        print '(A,3ES22.15)', ' SPLINED_DIAG: dxcart_dxvmec(:,2) = ', &
                            dxcart_dxvmec(:,2)
                        print '(A,3ES22.15)', ' SPLINED_DIAG: dxcart_dxvmec(:,3) = ', &
                            dxcart_dxvmec(:,3)
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

        print '(A,ES22.15)', ' SPLINED_DIAG: hth(32,32,32) = ', hth(32,32,32)
        print '(A,ES22.15)', ' SPLINED_DIAG: hr(32,32,32) = ', hr(32,32,32)

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

        ! Diagnostic: print VMEC spline coefficient to verify identical splines
        block
            use new_vmec_stuff_mod, only: sR, sZ
            if (size(sR) > 0) then
                print '(A,ES22.15)', ' VMEC_SPLINE_DIAG: sR(1,1,1,32,32,32) = ', &
                    sR(1,1,1,32,32,32)
                print '(A,ES22.15)', ' VMEC_SPLINE_DIAG: sZ(1,1,1,32,32,32) = ', &
                    sZ(1,1,1,32,32,32)
            end if
        end block

    contains

        subroutine print_progress(i, n)
            integer, intent(in) :: i, n
            !$omp critical
            if (mod(i, n / 10) == 0) then
                write(*, '(a,f6.2,a)', advance='no') char(13), 100.0*i/n, ' %'
                call flush(6)
            else
                write(*, *)
            end if
            !$omp end critical
        end subroutine print_progress

    end subroutine build_splines


    subroutine evaluate_at_ref_point(source, x_ref, x_cart, dxcart_dxvmec, &
                                      A_cart, B_cart, Bmod)
        !> Evaluate source field at a point given in VMEC ref coords (s, theta, phi).
        !> Returns Cartesian field values and the Jacobian for transformation.
        !> Note: Returns raw B_cart (not normalized) to match main branch operation order.
        class(magnetic_field_t), intent(in) :: source
        real(dp), intent(in) :: x_ref(3)
        real(dp), intent(out) :: x_cart(3), dxcart_dxvmec(3, 3)
        real(dp), intent(out) :: A_cart(3), B_cart(3), Bmod

        real(dp) :: h_unused(3)

        select type (coords => source%coords)
        type is (cartesian_coordinate_system_t)
            call transform_vmec_to_cart(x_ref, x_cart, dxcart_dxvmec)
            call source%evaluate(x_cart, A_cart, h_unused, Bmod, B_cart)
        class default
            error stop "evaluate_at_ref_point: non-Cartesian source not yet supported"
        end select
    end subroutine evaluate_at_ref_point


    subroutine transform_to_ref(A_cart, B_cart, dxcart_dxvmec, scaling_jac, &
                               Bmod, Acov, hcov)
        !> Transform covariant vectors from Cartesian to ref coords.
        !> Uses the correct transformation: A_ref_j = sum_i A_cart_i * dxcart_dxvmec(i,j)
        !> For h, we use matmul(B, J)/Bmod to match main branch operation order.
        real(dp), intent(in) :: A_cart(3), B_cart(3)
        real(dp), intent(in) :: dxcart_dxvmec(3, 3)
        real(dp), intent(in) :: scaling_jac(3)
        real(dp), intent(in) :: Bmod
        real(dp), intent(out) :: Acov(3), hcov(3)

        Acov = matmul(A_cart, dxcart_dxvmec)
        Acov(1) = Acov(1) * scaling_jac(1)

        hcov = matmul(B_cart, dxcart_dxvmec) / Bmod
        hcov(1) = hcov(1) * scaling_jac(1)
    end subroutine transform_to_ref


    subroutine splined_field_cleanup(self)
        type(splined_field_t), intent(inout) :: self

        if (self%initialized) then
            call destroy_batch_splines_3d(self%splines)
            self%initialized = .false.
        end if
    end subroutine splined_field_cleanup

end module field_splined
